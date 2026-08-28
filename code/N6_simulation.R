# N6_simulation.R
# Synthetic Study 6: Generative vs. Discriminative Attribution
#
# Key experiment: X3 strongly SEPARATES Class 1 from Class 2
#   (Class 2 has X3 ~ N(3, 0.2^2); Class 1 has X3 ~ N(0,1) — flat within Class 1)
# but X3 does NOT belong to Class 1's generative character.
# Ground truth: R_1 = {1,2}, R_2 = {3}, R_3 = {4,5}.
#
# Methods compared:
#   A. NPKDC-vd: bandwidth-based class-specific detection (generative)
#   B. OvR-RF: binary random forest Class y vs. rest, variable importance (discriminative)
#   C. OvR-RF top-r: take top r_y variables by OvR-RF importance per class
#
# The key diagnostic: Class 1 attribution
#   NPKDC-vd correctly identifies {X1, X2} (small bandwidths within Class 1)
#   OvR-RF inflates X3 (X3 discriminates Class 1 from Class 2 powerfully)
#
# Design: C=3, d=15, n_y=200, n_test=100, B=1000 replications
# Uses only: parallel, ranger

FIGURES_DIR <- file.path("..", "manuscript", "figures")
DATA_DIR    <- file.path("..", "data")
library(parallel)

set.seed(20260515)

C      <- 3
D      <- 15
N_Y    <- 200
N_TEST <- 100
B      <- 1000
DETECT_TAU <- -1.5

# True relevant sets and sizes
R_Y       <- c(2L, 1L, 2L)   # |R_1|=2, |R_2|=1, |R_3|=2
REL_SETS  <- list(c(1L, 2L), c(3L), c(4L, 5L))


# ---- Data generation ----------------------------------------------------------
gen_data <- function(n_y, n_test) {
  gen_class <- function(y, n) {
    X <- matrix(rnorm(n * D), nrow = n, ncol = D)  # all ~ N(0,1) by default
    if (y == 1L) {
      X[, 1L] <- rnorm(n, mean = 0.5, sd = 0.2)   # X1: concentrated
      X[, 2L] <- rnorm(n, mean = 0.5, sd = 0.2)   # X2: concentrated
      # X3 stays N(0,1) — FLAT within Class 1 (discriminative but not generative)
    } else if (y == 2L) {
      X[, 3L] <- rnorm(n, mean = 3.0, sd = 0.2)   # X3: concentrated at 3
    } else if (y == 3L) {
      X[, 4L] <- rnorm(n, mean = 0.5, sd = 0.2)   # X4: concentrated
      X[, 5L] <- rnorm(n, mean = 0.5, sd = 0.2)   # X5: concentrated
    }
    X
  }
  list(
    train = lapply(1:C, function(y) gen_class(y, n_y)),
    test  = lapply(1:C, function(y) gen_class(y, n_test))
  )
}


# ---- Silverman bandwidth -------------------------------------------------------
bw_silverman <- function(x) {
  s <- sd(x)
  if (s < 1e-12) return(1e-6)
  1.06 * s * length(x)^(-0.2)
}


# ---- NPKDC-vd detection: cross-class Z-score (ANOVA/Tukey spirit) -------------
detect_npkdc <- function(train_list) {
  bw_sil <- matrix(0, C, D)
  for (y in 1:C)
    bw_sil[y, ] <- apply(train_list[[y]], 2L, bw_silverman)
  detected <- matrix(0L, C, D)
  for (j in 1:D) {
    h_j <- bw_sil[, j]; s_j <- sd(h_j)
    if (s_j < 1e-10) next
    z_j <- (h_j - mean(h_j)) / s_j
    detected[, j] <- as.integer(z_j < DETECT_TAU)
  }
  list(detected = detected, bw = bw_sil)
}


# ---- OvR-RF: binary RF for each class, variable importance --------------------
# Requires: ranger
detect_ovr_rf <- function(train_list, R_Y_vec) {
  if (!requireNamespace("ranger", quietly = TRUE)) return(NULL)
  D_loc <- ncol(train_list[[1L]])
  X_all <- do.call(rbind, train_list)
  ylab  <- rep(1:C, times = vapply(train_list, nrow, integer(1)))

  imp_mat <- matrix(0, C, D_loc)
  for (y in 1:C) {
    ybin <- factor(as.integer(ylab == y), levels = c(0L, 1L))
    df   <- data.frame(y = ybin, X_all)
    mod  <- ranger::ranger(y ~ ., data = df, num.trees = 500,
                           importance = "impurity",
                           probability = FALSE,
                           verbose = FALSE)
    imp_mat[y, ] <- mod$variable.importance
  }

  # Native: top variables by importance — threshold at mean + 1 SD per class
  selected_native <- matrix(0L, C, D_loc)
  for (y in 1:C) {
    thr <- mean(imp_mat[y, ]) + sd(imp_mat[y, ])
    selected_native[y, ] <- as.integer(imp_mat[y, ] >= thr)
  }

  # Top-r: take exactly r_y variables per class
  selected_topr <- matrix(0L, C, D_loc)
  for (y in 1:C) {
    idx <- order(imp_mat[y, ], decreasing = TRUE)[1:R_Y_vec[y]]
    selected_topr[y, idx] <- 1L
  }

  list(score = imp_mat, selected_native = selected_native, selected_topr = selected_topr)
}


# ---- Classification: product-KDE rule -----------------------------------------
classify_npkdc <- function(test_list, train_list, det) {
  # Simple product-KDE Bayes rule
  bw <- det$bw  # C x D
  all_test <- do.call(rbind, test_list)
  true <- rep(1:C, times = vapply(test_list, nrow, integer(1)))
  n_te <- nrow(all_test)
  log_prob <- matrix(0, n_te, C)
  for (y in 1:C) {
    tr <- train_list[[y]]
    n_y <- nrow(tr)
    for (j in 1:D) {
      h <- bw[y, j]
      log_dens_j <- vapply(seq_len(n_te), function(i) {
        log(mean(dnorm((all_test[i, j] - tr[, j]) / h) / h) + 1e-300)
      }, numeric(1))
      log_prob[, y] <- log_prob[, y] + log_dens_j
    }
    log_prob[, y] <- log_prob[, y] + log(1 / C)  # equal priors
  }
  pred <- max.col(log_prob, ties.method = "first")
  list(pred = pred, true = true, accuracy = mean(pred == true))
}


# ---- Attribution metrics -------------------------------------------------------
attr_metrics <- function(detected_row, y) {
  rel   <- REL_SETS[[y]]
  irrel <- setdiff(1:D, rel)
  tp <- sum(detected_row[rel]   == 1L)
  fp <- sum(detected_row[irrel] == 1L)
  fn <- sum(detected_row[rel]   == 0L)
  prec  <- if (tp + fp > 0) tp / (tp + fp) else 0
  rec   <- if (tp + fn > 0) tp / (tp + fn) else 0
  f1    <- if (prec + rec > 0) 2 * prec * rec / (prec + rec) else 0
  exact <- as.numeric((tp == R_Y[y]) && (fp == 0L))
  # Also: selection of X3 for Class 1 (the diagnostic variable)
  x3_selected <- if (y == 1L) as.numeric(detected_row[3L] == 1L) else NA_real_
  c(prec = prec, rec = rec, f1 = f1, exact = exact, x3_for_c1 = x3_selected)
}


# ---- RF accuracy via OvR voting -----------------------------------------------
classify_ovr_rf <- function(train_list, test_list) {
  if (!requireNamespace("ranger", quietly = TRUE)) return(NULL)
  X_all <- do.call(rbind, train_list)
  ylab  <- rep(1:C, times = vapply(train_list, nrow, integer(1)))
  Xte   <- do.call(rbind, test_list)
  true  <- rep(1:C, times = vapply(test_list, nrow, integer(1)))
  votes <- matrix(0, nrow(Xte), C)
  for (y in 1:C) {
    ybin <- factor(as.integer(ylab == y), levels = c(0L, 1L))
    df   <- data.frame(y = ybin, X_all)
    mod  <- ranger::ranger(y ~ ., data = df, num.trees = 500,
                           probability = TRUE, verbose = FALSE)
    pred_prob <- predict(mod, data = data.frame(Xte))$predictions[, "1"]
    votes[, y] <- pred_prob
  }
  pred <- max.col(votes, ties.method = "first")
  mean(pred == true)
}


# ---- Single replication -------------------------------------------------------
one_rep <- function(b) {
  set.seed(b * 1313 + 77)
  dat <- gen_data(N_Y, N_TEST)

  # --- NPKDC-vd ---
  det_np <- detect_npkdc(dat$train)
  cl_np  <- classify_npkdc(dat$test, dat$train, det_np)
  np_acc <- cl_np$accuracy
  np_attr <- t(sapply(1:C, function(y) attr_metrics(det_np$detected[y, ], y)))

  # Top-r fallback for NPKDC-vd (if native detection misses)
  det_topr <- matrix(0L, C, D)
  for (y in 1:C) {
    sc <- -det_np$bw[y, ]  # smaller bw = higher "importance"
    idx <- order(sc, decreasing = TRUE)[1:R_Y[y]]
    det_topr[y, idx] <- 1L
  }
  np_topr_attr <- t(sapply(1:C, function(y) attr_metrics(det_topr[y, ], y)))

  # --- OvR-RF ---
  det_rf <- tryCatch(detect_ovr_rf(dat$train, R_Y),
                     error = function(e) NULL)
  if (is.null(det_rf)) {
    rf_acc   <- NA_real_
    rf_attr  <- matrix(NA_real_, C, 5)
    rf_topr  <- matrix(NA_real_, C, 5)
  } else {
    rf_acc  <- tryCatch(classify_ovr_rf(dat$train, dat$test),
                        error = function(e) NA_real_)
    rf_attr <- t(sapply(1:C, function(y) attr_metrics(det_rf$selected_native[y, ], y)))
    rf_topr <- t(sapply(1:C, function(y) attr_metrics(det_rf$selected_topr[y, ], y)))
  }
  colnames(np_attr) <- colnames(np_topr_attr) <- colnames(rf_attr) <- colnames(rf_topr) <-
    c("prec", "rec", "f1", "exact", "x3_for_c1")

  list(
    np_acc    = np_acc,
    rf_acc    = rf_acc,
    np_native = np_attr,
    np_topr   = np_topr_attr,
    rf_native = rf_attr,
    rf_topr   = rf_topr
  )
}


# ---- Main --------------------------------------------------------------------
n_cores <- max(1L, detectCores() - 1L)
cat(sprintf("Study 6: C=%d, d=%d, n_y=%d, B=%d | R_1={1,2}, R_2={3}, R_3={4,5}\n",
            C, D, N_Y, B))
cat(sprintf("KEY: X3 separates Class 1 from Class 2 (discriminative)\n"))
cat(sprintf("     X3 flat within Class 1 (NOT generatively relevant to R_1)\n\n"))
cat("Using", n_cores, "cores\n\n")

t0 <- proc.time()[["elapsed"]]
reps <- mclapply(seq_len(B), one_rep, mc.cores = n_cores)
cat(sprintf("Done in %.1f sec\n\n", proc.time()[["elapsed"]] - t0))


# ---- Aggregate ----------------------------------------------------------------
np_acc_vec  <- sapply(reps, `[[`, "np_acc")
rf_acc_vec  <- sapply(reps, `[[`, "rf_acc")

# Build C x 5 x B arrays for each method/mode
np_nat_arr <- simplify2array(lapply(reps, `[[`, "np_native"))   # C x 5 x B
np_top_arr <- simplify2array(lapply(reps, `[[`, "np_topr"))
rf_nat_arr <- simplify2array(lapply(reps, `[[`, "rf_native"))
rf_top_arr <- simplify2array(lapply(reps, `[[`, "rf_topr"))

mean_sd <- function(arr, y, metric) {
  vals <- arr[y, metric, ]
  vals <- vals[!is.na(vals)]
  c(mean = mean(vals), sd = sd(vals))
}

cat("=== ACCURACY (Stage 1) ===\n")
cat(sprintf("  NPKDC-vd : %.4f (%.4f)\n", mean(np_acc_vec), sd(np_acc_vec)))
cat(sprintf("  OvR-RF   : %.4f (%.4f)\n",
            mean(rf_acc_vec, na.rm=TRUE), sd(rf_acc_vec, na.rm=TRUE)))

cat("\n=== ATTRIBUTION F1 by class (Stage 2) ===\n")
cat("Method         | Class 1 (R={1,2}) | Class 2 (R={3}) | Class 3 (R={4,5})\n")
cat("               |  F1      exact    |  F1     exact   |  F1      exact\n")

for (mode in c("np_native","np_topr","rf_native","rf_topr")) {
  arr <- switch(mode,
    np_native = np_nat_arr, np_topr = np_top_arr,
    rf_native = rf_nat_arr, rf_topr = rf_top_arr)
  lbl <- switch(mode,
    np_native = "NPKDC-vd(nat)",
    np_topr   = "NPKDC-vd(top)",
    rf_native = "OvR-RF(nat)  ",
    rf_topr   = "OvR-RF(top)  ")
  f1 <- sapply(1:C, function(y) mean_sd(arr, y, "f1"))
  ex <- sapply(1:C, function(y) mean_sd(arr, y, "exact"))
  cat(sprintf("%s | %.3f(%.3f) %.3f | %.3f(%.3f) %.3f | %.3f(%.3f) %.3f\n",
    lbl,
    f1["mean",1], f1["sd",1], ex["mean",1],
    f1["mean",2], f1["sd",2], ex["mean",2],
    f1["mean",3], f1["sd",3], ex["mean",3]))
}

cat("\n=== X3 selection rate for Class 1 (key diagnostic) ===\n")
cat("  NPKDC-vd(nat): P(X3 selected for Class 1) =",
    mean(np_nat_arr[1,"x3_for_c1",], na.rm=TRUE), "\n")
cat("  NPKDC-vd(top): P(X3 in top-2 for Class 1) =",
    mean(np_top_arr[1,"x3_for_c1",], na.rm=TRUE), "\n")
cat("  OvR-RF(nat):   P(X3 selected for Class 1) =",
    mean(rf_nat_arr[1,"x3_for_c1",], na.rm=TRUE), "\n")
cat("  OvR-RF(top):   P(X3 in top-2 for Class 1) =",
    mean(rf_top_arr[1,"x3_for_c1",], na.rm=TRUE), "\n\n")


# ---- Save --------------------------------------------------------------------
saveRDS(list(
  reps = reps, B = B, C = C, D = D, R_Y = R_Y, N_Y = N_Y,
  np_acc_vec = np_acc_vec, rf_acc_vec = rf_acc_vec,
  np_nat_arr = np_nat_arr, np_top_arr = np_top_arr,
  rf_nat_arr = rf_nat_arr, rf_top_arr = rf_top_arr
), file.path("..", "data", "N6_results.rds"),
compress = FALSE)
cat("Results saved to N6_results.rds\n\n")


# ---- LaTeX table (Stage 1: accuracy) -----------------------------------------
cat("=== LATEX: Stage 1 accuracy ===\n")
cat(sprintf("NPKDC-vd & %.4f (%.4f) & $\\checkmark$ \\\\\n",
            mean(np_acc_vec), sd(np_acc_vec)))
cat(sprintf("OvR-RF   & %.4f (%.4f) & $\\checkmark$ \\\\\n",
            mean(rf_acc_vec, na.rm=TRUE), sd(rf_acc_vec, na.rm=TRUE)))


# ---- LaTeX table (Stage 2: F1 by class) --------------------------------------
cat("\n=== LATEX: Stage 2 F1 by class (top-r mode) ===\n")
cat("\\begin{tabular}{llccc}\\hline\n")
cat("Class & Method & Precision & Recall & F1 \\\\\n\\hline\\hline\n")
class_names <- c("Class~1 ($R_1=\\{1,2\\}$)", "Class~2 ($R_2=\\{3\\}$)", "Class~3 ($R_3=\\{4,5\\}$)")
for (y in 1:C) {
  np_p <- mean_sd(np_top_arr, y, "prec")
  np_r <- mean_sd(np_top_arr, y, "rec")
  np_f <- mean_sd(np_top_arr, y, "f1")
  rf_p <- mean_sd(rf_top_arr, y, "prec")
  rf_r <- mean_sd(rf_top_arr, y, "rec")
  rf_f <- mean_sd(rf_top_arr, y, "f1")
  cn <- if (y == 1L) paste0("\\multirow{2}{*}{", class_names[y], "}") else ""
  if (y > 1L) cat("\\hline\n")
  cat(sprintf("%s & NPKDC-vd & %.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f) \\\\\n",
              class_names[y], np_p["mean"], np_p["sd"], np_r["mean"], np_r["sd"],
              np_f["mean"], np_f["sd"]))
  cat(sprintf(" & OvR-RF   & %.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f) \\\\\n",
              rf_p["mean"], rf_p["sd"], rf_r["mean"], rf_r["sd"], rf_f["mean"], rf_f["sd"]))
}
cat("\\hline\n\\end{tabular}\n")

cat("\nP(X3 selected for Class 1):\n")
cat(sprintf("  NPKDC-vd (top-r): %.3f\n", mean(np_top_arr[1,"x3_for_c1",], na.rm=TRUE)))
cat(sprintf("  OvR-RF   (top-r): %.3f\n\n", mean(rf_top_arr[1,"x3_for_c1",], na.rm=TRUE)))
