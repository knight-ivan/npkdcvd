# N2_anova_tukey.R
# Study 3 replicated under the formal ANOVA/Tukey detection rule (Algorithm 2)
# Same design as N2_simulation.R (Z-score version), identical seed structure.
# Purpose: demonstrate that P(R^hat_y = R_y) -> 1 under the formal rule,
#          closing the gap between the Z-score empirical results and Theorem 3.2.
#
# Design: C=5, d=25, r_y=3 (disjoint), sigma_rel=0.5, sigma_irrel=1
# n_y in {30,50,75,100,150,200,300,500}, B=1000 replications

library(parallel)

set.seed(20260510)

C          <- 5
D          <- 25
N_Y_GRID   <- c(30, 50, 75, 100, 150, 200, 300, 500)
B          <- 1000
ALPHA      <- 0.05   # Tukey family-wise error level

R_Y       <- rep(3, C)
rel_start <- cumsum(c(1, R_Y[-C]))          # 1, 4, 7, 10, 13
REL_SETS  <- lapply(1:C, function(y) rel_start[y]:(rel_start[y] + R_Y[y] - 1))


gen_train <- function(n_y) {
  lapply(1:C, function(y) {
    X        <- matrix(rnorm(n_y * D), nrow = n_y, ncol = D)
    rel      <- REL_SETS[[y]]
    X[, rel] <- matrix(rnorm(n_y * R_Y[y], mean = 0.5, sd = 0.5),
                       nrow = n_y, ncol = R_Y[y])
    X
  })
}

bw_silverman <- function(x) {
  s <- sd(x); if (s < 1e-12) return(1e-6)
  1.06 * s * length(x)^(-0.2)
}


# ---- ANOVA/Tukey detection (Algorithm 2) ------------------------------------
# For each class y:
#   1. Compute mean Silverman bandwidth h_bar[j] for each variable j
#   2. One-way ANOVA F-test across variables; if H0 not rejected -> R_hat_y = {}
#   3. Sort means; apply Tukey HSD to find the gap between the smallest group
#      and the rest; variables below the split are declared relevant.
detect_anova_tukey <- function(train_list, alpha = ALPHA) {
  # Compute per-observation Silverman bandwidth for each (class, var)
  # We use the mean over n_y observations as the "group" mean for ANOVA
  bw_sil <- matrix(0, nrow = C, ncol = D)
  for (y in 1:C)
    bw_sil[y, ] <- apply(train_list[[y]], 2L, bw_silverman)

  detected <- matrix(0L, nrow = C, ncol = D)

  for (y in 1:C) {
    h_bar <- bw_sil[y, ]   # length-D vector of mean bandwidths for class y
    n_y   <- nrow(train_list[[y]])

    # Step 1: One-way ANOVA across variables
    # Each variable j contributes n_y pseudo-observations (LOO bandwidths) —
    # but since we only have one mean per variable, we use the means as data
    # and the within-variable variance is approximated by the Silverman rule
    # applied to leave-one-out bandwidth estimates.
    # Simpler and equivalent for detection: treat each h_bar[j] as a group mean
    # with n_y observations; the within-group variance is pooled from
    # the per-observation bandwidths (recompute using leave-one-out).

    # Compute per-observation bandwidths for ANOVA
    # bw_obs[i, j] = Silverman bw for class y variable j leaving out obs i
    # For speed, approximate: bw_obs[i,j] ≈ bw_sil[y,j] (constant within class)
    # This gives within-group MSW = 0 -> use a small regulariser
    # Better: compute LOO bandwidth variance directly

    X_y <- train_list[[y]]
    # Per-observation Silverman bandwidth for each variable
    bw_obs <- matrix(0, nrow = n_y, ncol = D)
    for (j in 1:D) {
      x_j <- X_y[, j]
      # LOO bandwidth estimates
      bw_obs[, j] <- sapply(seq_len(n_y), function(i) bw_silverman(x_j[-i]))
    }

    # ANOVA: groups = variables, observations = per-observation bw estimates
    group_means <- colMeans(bw_obs)    # length D
    grand_mean  <- mean(bw_obs)
    SS_between  <- n_y * sum((group_means - grand_mean)^2)
    SS_within   <- sum((bw_obs - matrix(group_means, n_y, D, byrow = TRUE))^2)
    df_between  <- D - 1
    df_within   <- D * (n_y - 1)
    MS_between  <- SS_between / df_between
    MS_within   <- SS_within  / df_within
    if (MS_within < 1e-15) {
      # Degenerate: all bandwidths identical -> no detection possible
      detected[y, ] <- 0L
      next
    }
    F_stat <- MS_between / MS_within
    p_val  <- pf(F_stat, df_between, df_within, lower.tail = FALSE)

    if (p_val > alpha) {
      detected[y, ] <- 0L
      next
    }

    # Step 2: Tukey HSD — find the split at the MAXIMUM gap in sorted means
    # Sort group means ascending
    ord      <- order(group_means)
    h_sorted <- group_means[ord]

    # Tukey critical value: q / sqrt(2) gives the HSD for pairwise comparisons
    q_crit <- qtukey(1 - alpha, nmeans = D, df = df_within) / sqrt(2)
    SE      <- sqrt(MS_within / n_y)

    # Locate the position of the largest gap
    gaps      <- diff(h_sorted)
    split_idx <- which.max(gaps)   # index o: variables 1..o are "low", o+1..D are "high"

    # Accept the split only if the maximum gap exceeds the Tukey critical value
    if (gaps[split_idx] <= q_crit * SE) {
      detected[y, ] <- 0L
    } else {
      rel_vars <- ord[seq_len(split_idx)]
      detected[y, rel_vars] <- 1L
    }
  }
  detected
}


attr_metrics_full <- function(detected) {
  prec <- numeric(C); rec <- numeric(C); f1 <- numeric(C); exact <- logical(C)
  for (y in 1:C) {
    rel   <- REL_SETS[[y]]
    irrel <- setdiff(seq_len(D), rel)
    tp    <- sum(detected[y, rel]   == 1L)
    fp    <- sum(detected[y, irrel] == 1L)
    fn    <- sum(detected[y, rel]   == 0L)
    prec[y]  <- if (tp + fp > 0) tp / (tp + fp) else 0
    rec[y]   <- if (tp + fn > 0) tp / (tp + fn) else 0
    f1[y]    <- if (prec[y] + rec[y] > 0) 2*prec[y]*rec[y]/(prec[y]+rec[y]) else 0
    exact[y] <- (tp == R_Y[y]) && (fp == 0L)
  }
  list(prec = prec, rec = rec, f1 = f1, exact = as.numeric(exact))
}

one_rep <- function(b, n_y) {
  set.seed(b * 1000 + n_y)
  train <- gen_train(n_y)
  det   <- detect_anova_tukey(train)
  attr_metrics_full(det)
}

n_cores <- max(1L, detectCores() - 1L)
cat(sprintf("Study 3 (ANOVA/Tukey): C=%d, d=%d, r_y=%d, alpha=%.2f\n",
            C, D, R_Y[1], ALPHA))
cat("Using", n_cores, "cores\n\n")

results <- vector("list", length(N_Y_GRID))
names(results) <- as.character(N_Y_GRID)

for (k in seq_along(N_Y_GRID)) {
  n_y <- N_Y_GRID[k]
  cat(sprintf("n_y = %3d ... ", n_y))
  t0   <- proc.time()[["elapsed"]]
  reps <- mclapply(seq_len(B), function(b) one_rep(b, n_y), mc.cores = n_cores)

  prec_mat  <- sapply(reps, `[[`, "prec")
  rec_mat   <- sapply(reps, `[[`, "rec")
  f1_mat    <- sapply(reps, `[[`, "f1")
  exact_mat <- sapply(reps, `[[`, "exact")

  results[[k]] <- list(
    n_y      = n_y,
    prec_mn  = rowMeans(prec_mat),  prec_sd  = apply(prec_mat,  1, sd),
    rec_mn   = rowMeans(rec_mat),   rec_sd   = apply(rec_mat,   1, sd),
    f1_mn    = rowMeans(f1_mat),    f1_sd    = apply(f1_mat,    1, sd),
    exact_mn = rowMeans(exact_mat), exact_sd = apply(exact_mat, 1, sd)
  )
  cat(sprintf("done (%.1f s) | Recall=%.3f  Prec=%.3f  F1=%.3f  P(exact)=%.3f\n",
              proc.time()[["elapsed"]] - t0,
              mean(results[[k]]$rec_mn), mean(results[[k]]$prec_mn),
              mean(results[[k]]$f1_mn),  mean(results[[k]]$exact_mn)))
}

BASEDIR <- file.path("..", "data")
saveRDS(results, file.path(BASEDIR, "N2_anova_results.rds"), compress = FALSE)
cat("\nResults saved to N2_anova_results.rds\n\n")

# ---- LaTeX table ----
cat("========== LATEX TABLE (Supplement Table S.6) ==========\n")
cat("\\begin{tabular}{ccccc}\\hline\n")
cat("$n_y$ & Recall & Precision & F1 & $P(\\hat{R}_y = R_y)$ \\\\\n")
cat("\\hline\\hline\n")
for (k in seq_along(N_Y_GRID)) {
  r <- results[[k]]
  cat(sprintf("%d & %.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f) \\\\\n",
              r$n_y,
              mean(r$rec_mn),   mean(r$rec_sd),
              mean(r$prec_mn),  mean(r$prec_sd),
              mean(r$f1_mn),    mean(r$f1_sd),
              mean(r$exact_mn), mean(r$exact_sd)))
}
cat("\\hline\n\\end{tabular}\n")
