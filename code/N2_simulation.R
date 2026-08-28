# N2_simulation.R
# Phase 2, Study N2: Detection quality vs n_y — empirical validation of Theorem 2
# Shows recall converging to 1 as n_y → ∞ for each class y
#
# Design: C=5, d=25 total vars, r_y=3 per class (disjoint, vars 1-15 relevant,
#         vars 16-25 always irrelevant)
# Relevant vars: N(0.5, 0.5²) — moderate signal, so detection quality improves
#               visibly as n_y grows (Z-score variance decreases toward 0)
# Irrelevant vars: N(0,1)
# Detection: cross-class Silverman Z-score, threshold DETECT_TAU = -1.5
#   For C=5, Z_relevant → -4/√5 ≈ -1.789 < -1.5 (detectable in limit)
#   With σ_rel=0.5, Z variance is large at small n_y → visible convergence curve
# B = 1000 replications per n_y; n_y ∈ {30,50,75,100,150,200,300,500}

FIGURES_DIR <- file.path("..", "manuscript", "figures")
DATA_DIR    <- file.path("..", "data")
library(parallel)

set.seed(20260510)

C          <- 5
D          <- 25
N_Y_GRID   <- c(30, 50, 75, 100, 150, 200, 300, 500)
B          <- 1000
DETECT_TAU <- -1.5

R_Y       <- rep(3, C)
rel_start <- cumsum(c(1, R_Y[-C]))          # 1, 4, 7, 10, 13
REL_SETS  <- lapply(1:C, function(y) rel_start[y]:(rel_start[y] + R_Y[y] - 1))
# vars 16-25 always irrelevant for all classes


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

detect_vars <- function(train_list) {
  bw_sil <- matrix(0, nrow = C, ncol = D)
  for (y in 1:C)
    bw_sil[y, ] <- apply(train_list[[y]], 2L, bw_silverman)
  detected <- matrix(0L, nrow = C, ncol = D)
  for (j in 1:D) {
    h_j <- bw_sil[, j]; s_j <- sd(h_j)
    if (s_j < 1e-10) next
    z_j <- (h_j - mean(h_j)) / s_j
    detected[, j] <- as.integer(z_j < DETECT_TAU)
  }
  detected
}

# Returns named list: prec, rec, f1, exact — each a length-C vector
attr_metrics_full <- function(detected) {
  prec  <- numeric(C); rec <- numeric(C); f1 <- numeric(C); exact <- logical(C)
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
  list(prec = prec, rec = rec, f1 = f1, exact = exact)
}

one_rep <- function(b, n_y) {
  set.seed(b * 1000 + n_y)
  train <- gen_train(n_y)
  det   <- detect_vars(train)
  attr_metrics_full(det)
}

# ---- Main loop ----
n_cores <- max(1L, detectCores() - 1L)
cat(sprintf("C=%d | d=%d | r_y=%d per class | σ_rel=0.5 | DETECT_TAU=%.1f\n",
            C, D, R_Y[1], DETECT_TAU))
cat("Using", n_cores, "cores\n\n")

results <- vector("list", length(N_Y_GRID))
names(results) <- as.character(N_Y_GRID)

for (k in seq_along(N_Y_GRID)) {
  n_y <- N_Y_GRID[k]
  cat(sprintf("n_y = %3d ... ", n_y))
  t0 <- proc.time()[["elapsed"]]
  reps <- mclapply(seq_len(B), function(b) one_rep(b, n_y), mc.cores = n_cores)

  # Correct aggregation: sapply over reps, not array reshape
  prec_mat  <- sapply(reps, `[[`, "prec")   # C x B matrix
  rec_mat   <- sapply(reps, `[[`, "rec")
  f1_mat    <- sapply(reps, `[[`, "f1")
  exact_mat <- sapply(reps, `[[`, "exact")  # C x B logical matrix

  results[[k]] <- list(
    n_y      = n_y,
    prec_mn  = rowMeans(prec_mat),  prec_sd  = apply(prec_mat,  1, sd),
    rec_mn   = rowMeans(rec_mat),   rec_sd   = apply(rec_mat,   1, sd),
    f1_mn    = rowMeans(f1_mat),    f1_sd    = apply(f1_mat,    1, sd),
    exact_mn = rowMeans(exact_mat), exact_sd = apply(exact_mat, 1, sd)
  )

  cat(sprintf("done (%.1f s) | Recall=%.3f  Prec=%.3f  F1=%.3f  P(exact)=%.3f\n",
              proc.time()[["elapsed"]] - t0,
              mean(results[[k]]$rec_mn),
              mean(results[[k]]$prec_mn),
              mean(results[[k]]$f1_mn),
              mean(results[[k]]$exact_mn)))
}

# ---- Save ----
out_path <- file.path("..", "data", "N2_results.rds")
saveRDS(results, out_path, compress = FALSE)
cat("\nResults saved to N2_results.rds\n\n")


# ---- LaTeX table ----
cat("========== LATEX TABLE ==========\n")
cat("\\begin{tabular}{ccccc}\\hline\n")
cat("$n_y$ & Recall & Precision & F1 & $P(\\hat R_y = R_y)$ \\\\\n")
cat("\\hline\\hline\n")
for (k in seq_along(N_Y_GRID)) {
  r <- results[[k]]
  cat(sprintf("%d & %.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f) \\\\\n",
              r$n_y,
              mean(r$rec_mn),  mean(r$rec_sd),
              mean(r$prec_mn), mean(r$prec_sd),
              mean(r$f1_mn),   mean(r$f1_sd),
              mean(r$exact_mn), mean(r$exact_sd)))
}
cat("\\hline\n\\end{tabular}\n\n")


# ---- Figure: recall, precision, F1 vs n_y ----
rec_agg   <- sapply(results, function(r) mean(r$rec_mn))
prec_agg  <- sapply(results, function(r) mean(r$prec_mn))
f1_agg    <- sapply(results, function(r) mean(r$f1_mn))
exact_agg <- sapply(results, function(r) mean(r$exact_mn))

fig_path <- file.path("..", "figures", "N2_detprob.eps")
postscript(fig_path, width = 5.5, height = 4.2,
           horizontal = FALSE, onefile = FALSE, paper = "special")

cols <- c("#1b7837", "#762a83", "#4393c3", "#d6604d")
pchs <- c(16, 17, 15, 4)
ltys <- c(1, 2, 3, 1)
lwds <- c(1.8, 1.8, 1.8, 1.5)

par(mar = c(4, 4.2, 1.2, 1))
plot(N_Y_GRID, rec_agg, type = "b",
     ylim = c(0.5, 1.02),
     xlab = expression(n[y] ~ "(per-class training size)"),
     ylab = "Detection quality (mean over classes and reps)",
     col = cols[1], pch = pchs[1], lty = ltys[1], lwd = lwds[1],
     cex = 1.1, cex.axis = 0.9, cex.lab = 0.95, las = 1)
abline(h = 1, col = "grey70", lty = 3)
lines(N_Y_GRID, prec_agg,  type = "b", col = cols[2], pch = pchs[2],
      lty = ltys[2], lwd = lwds[2], cex = 1.1)
lines(N_Y_GRID, f1_agg,    type = "b", col = cols[3], pch = pchs[3],
      lty = ltys[3], lwd = lwds[3], cex = 1.1)
lines(N_Y_GRID, exact_agg, type = "b", col = cols[4], pch = pchs[4],
      lty = ltys[4], lwd = lwds[4], cex = 1.1)
legend("bottomright", bty = "n",
       legend = c("Recall", "Precision", "F1", expression(P(hat(R)[y]==R[y]))),
       col = cols, pch = pchs, lty = ltys, lwd = lwds, cex = 0.85)
dev.off()
cat("Figure saved to N2_detprob.eps\n")
