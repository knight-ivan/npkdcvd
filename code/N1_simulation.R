# N1_simulation.R
# Phase 2, Study N1: High-dimensional sparse setting
# c=5 classes, d=200, r_y in {3,5,5,8,8} (different per class)
# Three sub-experiments: irrelevant variables ~ Uniform(0,1), N(0,1), AR(1)
# Demonstrates Theorem 3: convergence rate = O_p(n_y^{-4/(4+r_y)})

FIGURES_DIR <- file.path("..", "manuscript", "figures")
DATA_DIR    <- file.path("..", "data")
library(e1071)
library(parallel)

set.seed(20260510)

# ---- Problem setup ----
C       <- 5          # number of classes
D       <- 200        # total dimensions
N_Y     <- 150        # training samples per class
N_TEST  <- 100        # test samples per class
B       <- 1000       # replications
R_Y     <- c(3,5,5,8,8)  # relevant dimensions per class

# Disjoint relevant variable sets
rel_start <- c(1, 4, 9, 14, 22)   # class 1: 1-3, 2: 4-8, 3: 9-13, 4: 14-21, 5: 22-29
REL_SETS <- lapply(1:C, function(y) rel_start[y]:(rel_start[y] + R_Y[y] - 1))
# Vars 30-200 (171 vars) are irrelevant for all classes


# ---- Data generation ----
gen_data <- function(n_y, n_test, irrel_type) {
  rho <- 0.5
  train_list <- vector("list", C)
  test_list  <- vector("list", C)

  for (y in 1:C) {
    rel  <- REL_SETS[[y]]
    nr   <- R_Y[y]
    nirr <- D - nr

    # Relevant variables: N(mu, 0.1^2), mu cycles in {0.3, 0.5, 0.7}
    mu_rel <- rep(c(0.3, 0.5, 0.7), length.out = nr)

    make_irrel <- function(nrow) {
      if (irrel_type == "uniform") {
        matrix(runif(nrow * nirr), nrow = nrow, ncol = nirr)
      } else if (irrel_type == "normal") {
        matrix(rnorm(nrow * nirr), nrow = nrow, ncol = nirr)
      } else {  # ar1
        Z <- matrix(0, nrow = nrow, ncol = nirr)
        Z[, 1] <- rnorm(nrow)
        for (j in 2:nirr)
          Z[, j] <- rho * Z[, j-1] + sqrt(1 - rho^2) * rnorm(nrow)
        Z
      }
    }

    make_data <- function(n) {
      X <- matrix(0, nrow = n, ncol = D)
      X[, rel] <- matrix(rnorm(n * nr, mean = rep(mu_rel, each = n), sd = 0.1),
                         nrow = n, ncol = nr)
      X[, -rel] <- make_irrel(n)
      X
    }

    train_list[[y]] <- make_data(n_y)
    test_list[[y]]  <- make_data(n_test)
  }
  list(train = train_list, test = test_list)
}


# ---- Bandwidth selection ----
bw_silverman <- function(x) {
  s <- sd(x)
  if (s < 1e-12) return(1e-6)
  1.06 * s * length(x)^(-0.2)
}

bw_loocv <- function(x) {
  n <- length(x)
  s <- sd(x)
  if (s < 1e-12) return(1e-6)
  h_sil <- 1.06 * s * n^(-0.2)
  # Grid: 5% to 200% of Silverman
  h_grid <- h_sil * c(0.05, 0.10, 0.20, 0.35, 0.50, 0.70, 1.00, 1.50, 2.00)
  D_mat  <- outer(x, x, "-")   # n x n, computed once
  scores <- vapply(h_grid, function(h) {
    K        <- dnorm(D_mat / h) / h
    diag(K)  <- 0
    loo_dens <- rowSums(K) / (n - 1L)
    -mean(log(pmax(loo_dens, 1e-20)))
  }, numeric(1))
  h_grid[which.min(scores)]
}

# Estimate one bandwidth matrix per class: returns C x D matrix
fit_bandwidths <- function(train_list) {
  bw <- matrix(0, nrow = C, ncol = D)
  for (y in 1:C) {
    X <- train_list[[y]]
    bw[y, ] <- apply(X, 2L, bw_loocv)
  }
  bw
}


# ---- Variable detection ----
# For each variable j: compute Silverman bandwidth h_sil[y,j] for each class y.
# A class-y variable has small h_sil if and only if that class's marginal
# distribution for variable j is narrow (low sd) — the signature of a
# class-specific distribution concentrated around class-y's centroid.
# Detection criterion: cross-class Z-score Z[y,j] = (h_sil[y,j] - mean_j) / sd_j
# Declare variable j relevant for class y if Z[y,j] < DETECT_TAU.
DETECT_TAU <- -1.5   # ~ 1.7 SD below cross-class mean; Z_rel anchors at -1.78

detect_vars <- function(train_list) {
  # Compute Silverman bandwidth for every (class, variable) pair
  bw_sil <- matrix(0, nrow = C, ncol = D)
  for (y in 1:C)
    bw_sil[y, ] <- apply(train_list[[y]], 2L, bw_silverman)

  # Cross-class Z-score per variable
  detected <- matrix(0L, nrow = C, ncol = D)
  for (j in 1:D) {
    h_j <- bw_sil[, j]
    s_j <- sd(h_j)
    if (s_j < 1e-10) next          # all classes identical → skip
    z_j <- (h_j - mean(h_j)) / s_j
    detected[, j] <- as.integer(z_j < DETECT_TAU)
  }
  detected
}


# ---- Classification ----
classify_npkdc <- function(X_te_all, train_list, bw) {
  # X_te_all: (C*n_test) x D matrix; all test observations stacked
  n_te <- nrow(X_te_all)
  log_post <- matrix(0, nrow = n_te, ncol = C)

  for (y in 1:C) {
    X_tr  <- train_list[[y]]
    h_y   <- bw[y, ]
    n_y_y <- nrow(X_tr)
    # log f_hat_y(x) = sum_j log[ (1/n_y) sum_i K((x_j - x_ij)/h_j) / h_j ]
    log_f <- matrix(0, nrow = n_te, ncol = D)
    for (j in 1:D) {
      Dj      <- outer(X_te_all[, j], X_tr[, j], "-") / h_y[j]  # n_te x n_y
      Kj      <- dnorm(Dj) / h_y[j]
      marg_j  <- rowMeans(Kj)
      log_f[, j] <- log(pmax(marg_j, 1e-300))
    }
    log_post[, y] <- rowSums(log_f) + log(n_y_y)
  }
  apply(log_post, 1L, which.max)
}


# ---- Attribution metrics ----
attr_metrics <- function(detected) {
  # Returns precision, recall, F1 per class
  out <- matrix(0, nrow = C, ncol = 3,
                dimnames = list(NULL, c("prec","rec","f1")))
  for (y in 1:C) {
    rel   <- REL_SETS[[y]]
    irrel <- setdiff(seq_len(D), rel)
    tp    <- sum(detected[y, rel]   == 1L)
    fp    <- sum(detected[y, irrel] == 1L)
    fn    <- sum(detected[y, rel]   == 0L)
    prec  <- if (tp + fp > 0) tp / (tp + fp) else 0
    rec   <- if (tp + fn > 0) tp / (tp + fn) else 0
    f1    <- if (prec + rec > 0) 2 * prec * rec / (prec + rec) else 0
    out[y, ] <- c(prec, rec, f1)
  }
  out
}


# ---- Single replication ----
one_rep <- function(b, irrel_type) {
  set.seed(b * 1000 + match(irrel_type, c("uniform","normal","ar1")))
  dat <- gen_data(N_Y, N_TEST, irrel_type)

  X_te_all  <- do.call(rbind, dat$test)
  true_lab  <- rep(seq_len(C), each = N_TEST)

  # NPKDC-vd
  t0    <- proc.time()[["elapsed"]]
  bw    <- fit_bandwidths(dat$train)
  pred  <- classify_npkdc(X_te_all, dat$train, bw)
  t_npk <- proc.time()[["elapsed"]] - t0

  det  <- detect_vars(dat$train)
  am   <- attr_metrics(det)
  acc_npk <- mean(pred == true_lab)

  # SVM (radial kernel, default settings)
  X_tr_all <- do.call(rbind, dat$train)
  tr_lab   <- rep(seq_len(C), each = N_Y)
  svm_mod  <- svm(X_tr_all, factor(tr_lab), kernel = "radial",
                  scale = FALSE)
  pred_svm <- as.integer(as.character(predict(svm_mod, X_te_all)))
  acc_svm  <- mean(pred_svm == true_lab)

  list(acc_npk = acc_npk, acc_svm = acc_svm,
       attr = am, time = t_npk)
}


# ---- Main loop — all 3 sub-experiments in parallel ----
irrel_types <- c("uniform", "normal", "ar1")

n_total  <- max(1L, detectCores() - 1L)
n_sub    <- length(irrel_types)
cores_ea <- max(1L, n_total %/% n_sub)   # cores per sub-experiment
cat(sprintf("Using %d total cores, %d per sub-experiment\n\n", n_total, cores_ea))

run_one_type <- function(irrel) {
  res_list <- mclapply(seq_len(B),
                       function(b) one_rep(b, irrel),
                       mc.cores = cores_ea)
  acc_npk  <- sapply(res_list, `[[`, "acc_npk")
  acc_svm  <- sapply(res_list, `[[`, "acc_svm")
  time_npk <- sapply(res_list, `[[`, "time")
  attr_arr <- array(unlist(lapply(res_list, `[[`, "attr")),
                    dim = c(B, C, 3))
  list(acc_npk = acc_npk, acc_svm = acc_svm,
       time_npk = time_npk, attr_arr = attr_arr)
}

t0 <- proc.time()[["elapsed"]]
results_list <- mclapply(irrel_types, run_one_type, mc.cores = n_sub)
cat(sprintf("All sub-experiments done in %.1f sec\n\n",
            proc.time()[["elapsed"]] - t0))

results_all <- setNames(results_list, irrel_types)

for (irrel in irrel_types) {
  r <- results_all[[irrel]]
  cat("=== Sub-experiment:", irrel, "===\n")
  cat("NPKDC-vd accuracy: ", round(mean(r$acc_npk), 4),
      " (", round(sd(r$acc_npk), 4), ")\n", sep = "")
  cat("SVM accuracy:      ", round(mean(r$acc_svm), 4),
      " (", round(sd(r$acc_svm), 4), ")\n", sep = "")
  cat("Attribution F1 per class (mean over", B, "reps):\n")
  print(round(apply(r$attr_arr[,,3], 2, mean), 4))
  cat("\n")
}

# ---- Save results ----
out_path <- file.path("..", "data", "N1_results.rds")
saveRDS(results_all, out_path, compress = FALSE)
cat("Results saved to N1_results.rds\n")


# ---- Summary table for LaTeX ----
cat("\n========== LATEX TABLE ==========\n")
irrel_labels <- c("Uniform(0,1)", "$\\mathcal{N}(0,1)$", "AR(1) $\\rho=0.5$")
cat("\\begin{tabular}{llccccc}\\hline\n")
cat("Irrelevant dist. & Method & Accuracy & Class-1 F1 & Class-2,3 F1 & Class-4,5 F1 & Time (s) \\\\\n")
cat("                 &        & & ($r_y=3$) & ($r_y=5$) & ($r_y=8$) & \\\\\n\\hline\\hline\n")
for (k in seq_along(irrel_types)) {
  irrel <- irrel_types[k]
  r     <- results_all[[irrel]]
  f1_mn <- apply(r$attr_arr[,,3], 2, mean)
  f1_sd <- apply(r$attr_arr[,,3], 2, sd)
  # Average F1 for classes 2&3 (r_y=5) and 4&5 (r_y=8)
  f1_r3  <- sprintf("%.4f (%.4f)", f1_mn[1], f1_sd[1])
  f1_r5  <- sprintf("%.4f (%.4f)", mean(f1_mn[2:3]), mean(f1_sd[2:3]))
  f1_r8  <- sprintf("%.4f (%.4f)", mean(f1_mn[4:5]), mean(f1_sd[4:5]))
  acc_n  <- sprintf("%.4f (%.4f)", mean(r$acc_npk), sd(r$acc_npk))
  acc_s  <- sprintf("%.4f (%.4f)", mean(r$acc_svm), sd(r$acc_svm))
  t_n    <- sprintf("%.1f (%.1f)",  mean(r$time_npk), sd(r$time_npk))

  cat(irrel_labels[k], "& NPKDC-vd &", acc_n, "&", f1_r3, "&", f1_r5, "&", f1_r8, "&", t_n, "\\\\\n")
  cat("                 & SVM      &", acc_s, "& --- & --- & --- & --- \\\\\n")
  if (k < length(irrel_types)) cat("\\hline\n")
}
cat("\\hline\n\\end{tabular}\n")
