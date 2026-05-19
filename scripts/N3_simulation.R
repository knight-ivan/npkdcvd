# N3_simulation.R
# Phase 2, Study N3: Class-level vs per-observation attribution stability
#
# Core XAI-era argument: NPKDC-vd provides STABLE class-level attribution
# R̂_y (one set per class), while per-observation explanations vary across
# observations of the same class.
#
# Design: C=5, d=20, r_y=3 per class (disjoint, vars 1-15 relevant, 16-20 irrel)
# Relevant: N(0.5, 0.5²); Irrelevant: N(0,1)
# n_train=150 per class, n_test=100 per class
# B=100 replications
#
# Methods compared:
#  A. NPKDC-vd class-level: R̂_y from Silverman Z-score detection (stable per class)
#  B. Per-observation: Â_{y,i} = top-r_y vars by per-obs log-density contrast φ_j
#     φ_j(x_i, y) = log f̂_{y,j}(x_ij) - mean_{y'} log f̂_{y',j}(x_ij)
#
# Metrics:
#  1. Attribution accuracy: P(exact match with R_y) for both methods
#  2. Within-class Jaccard stability: mean Jaccard over same-class pairs
#     For method A: always 1 (by construction). For method B: < 1.

FIGURES_DIR <- file.path("..", "manuscript", "figures")
DATA_DIR    <- file.path("..", "data")
library(parallel)

set.seed(20260510)

C       <- 5
D       <- 20
N_Y     <- 150
N_TEST  <- 100
B       <- 100
DETECT_TAU <- -1.5
R_Y     <- rep(3, C)
rel_start <- cumsum(c(1, R_Y[-C]))   # 1, 4, 7, 10, 13
REL_SETS  <- lapply(1:C, function(y) rel_start[y]:(rel_start[y] + R_Y[y] - 1))
# vars 16-20 always irrelevant


# ---- Data generation ----
gen_data <- function(n_y, n_test) {
  train <- lapply(1:C, function(y) {
    X <- matrix(rnorm(n_y * D), nrow = n_y, ncol = D)
    X[, REL_SETS[[y]]] <- matrix(rnorm(n_y * R_Y[y], mean = 0.5, sd = 0.5),
                                  nrow = n_y, ncol = R_Y[y])
    X
  })
  test <- lapply(1:C, function(y) {
    X <- matrix(rnorm(n_test * D), nrow = n_test, ncol = D)
    X[, REL_SETS[[y]]] <- matrix(rnorm(n_test * R_Y[y], mean = 0.5, sd = 0.5),
                                  nrow = n_test, ncol = R_Y[y])
    X
  })
  list(train = train, test = test)
}


# ---- Silverman bandwidth and Z-score detection ----
bw_silverman <- function(x) {
  s <- sd(x); if (s < 1e-12) return(1e-6)
  1.06 * s * length(x)^(-0.2)
}

detect_class_level <- function(train_list) {
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
  detected   # C x D binary matrix: R̂_y = which(detected[y,]==1)
}


# ---- 1D KDE density at query points ----
kde1d <- function(query, train_x, h) {
  # returns log-density at each query point
  n <- length(train_x)
  log_dens <- apply(matrix(query), 1, function(q) {
    log(mean(dnorm((q - train_x) / h) / h) + 1e-300)
  })
  log_dens
}


# ---- Per-observation log-density contrast ----
# φ_j(x_i, y) = log f̂_{y,j}(x_ij) - mean_{y'} log f̂_{y',j}(x_ij)
per_obs_contrast <- function(x_test_i, y_true, train_list, bw_sil_mat) {
  phi <- numeric(D)
  for (j in 1:D) {
    log_dens_all <- vapply(1:C, function(yp) {
      kde1d(x_test_i[j], train_list[[yp]][, j], bw_sil_mat[yp, j])
    }, numeric(1))
    phi[j] <- log_dens_all[y_true] - mean(log_dens_all)
  }
  phi
}


# ---- Attribution metrics ----
attr_acc <- function(detected_row, y) {
  rel   <- REL_SETS[[y]]
  irrel <- setdiff(1:D, rel)
  tp <- sum(detected_row[rel] == 1L)
  fp <- sum(detected_row[irrel] == 1L)
  fn <- sum(detected_row[rel] == 0L)
  prec  <- if (tp + fp > 0) tp / (tp + fp) else 0
  rec   <- if (tp + fn > 0) tp / (tp + fn) else 0
  f1    <- if (prec + rec > 0) 2 * prec * rec / (prec + rec) else 0
  exact <- (tp == R_Y[y]) && (fp == 0L)
  c(prec = prec, rec = rec, f1 = f1, exact = as.numeric(exact))
}

jaccard <- function(a, b) {
  # a, b: integer vectors of set membership (length D, 0/1)
  inter <- sum(a & b)
  union <- sum(a | b)
  if (union == 0) return(1)
  inter / union
}


# ---- Single replication ----
one_rep <- function(b) {
  set.seed(b * 999 + 42)
  dat <- gen_data(N_Y, N_TEST)

  # ---- A. Class-level detection ----
  cl_det <- detect_class_level(dat$train)  # C x D binary
  cl_acc <- t(sapply(1:C, function(y) attr_acc(cl_det[y, ], y)))
  # cl_acc: C x 4 matrix (prec, rec, f1, exact)

  # ---- Bandwidths for per-obs KDE ----
  bw_sil <- matrix(0, C, D)
  for (y in 1:C)
    bw_sil[y, ] <- apply(dat$train[[y]], 2L, bw_silverman)

  # ---- B. Per-observation attribution ----
  po_acc_list <- vector("list", C)
  jaccard_list <- vector("list", C)

  for (y in 1:C) {
    X_te <- dat$test[[y]]  # n_test x D

    # φ matrix: n_test x D
    phi_mat <- t(apply(X_te, 1, function(xi)
      per_obs_contrast(xi, y, dat$train, bw_sil)))

    # Â_{y,i} = top-r_y by φ (declare relevant)
    attr_mat <- matrix(0L, N_TEST, D)
    for (i in 1:N_TEST) {
      top_idx <- order(phi_mat[i, ], decreasing = TRUE)[1:R_Y[y]]
      attr_mat[i, top_idx] <- 1L
    }

    # Per-obs accuracy
    acc_i <- t(apply(attr_mat, 1, function(r) attr_acc(r, y)))
    po_acc_list[[y]] <- colMeans(acc_i)  # mean over test obs

    # Within-class Jaccard stability: sample 50 pairs for speed
    n_pairs <- min(N_TEST * (N_TEST - 1) / 2, 200L)
    idx_pairs <- matrix(sample(N_TEST, 2 * n_pairs, replace = TRUE), ncol = 2)
    jac_vals <- apply(idx_pairs, 1, function(ij)
      jaccard(attr_mat[ij[1], ], attr_mat[ij[2], ]))
    jaccard_list[[y]] <- mean(jac_vals)
  }

  po_acc  <- do.call(rbind, po_acc_list)   # C x 4
  jac_per_class <- unlist(jaccard_list)     # length C

  list(
    cl_acc = cl_acc,       # C x 4: class-level prec/rec/f1/exact
    po_acc = po_acc,       # C x 4: per-obs mean prec/rec/f1/exact
    jac_po = jac_per_class # C: per-class within-Jaccard for per-obs
    # class-level Jaccard is always 1 (constant across obs)
  )
}


# ---- Main ----
n_cores <- max(1L, detectCores() - 1L)
cat(sprintf("N3: C=%d, d=%d, r_y=%d, n_y=%d, n_test=%d, B=%d\n",
            C, D, R_Y[1], N_Y, N_TEST, B))
cat("Using", n_cores, "cores\n\n")

t0 <- proc.time()[["elapsed"]]
reps <- mclapply(seq_len(B), one_rep, mc.cores = n_cores)
cat(sprintf("Total time: %.1f sec\n\n", proc.time()[["elapsed"]] - t0))

# ---- Aggregate ----
cl_acc_arr  <- simplify2array(lapply(reps, `[[`, "cl_acc"))   # C x 4 x B
po_acc_arr  <- simplify2array(lapply(reps, `[[`, "po_acc"))   # C x 4 x B
jac_po_mat  <- sapply(reps, `[[`, "jac_po")                   # C x B

# Means over reps and classes
cl_mn <- apply(cl_acc_arr, c(1, 2), mean)    # C x 4
po_mn <- apply(po_acc_arr, c(1, 2), mean)    # C x 4
cl_sd <- apply(cl_acc_arr, c(1, 2), sd)
po_sd <- apply(po_acc_arr, c(1, 2), sd)
jac_po_mn <- rowMeans(jac_po_mat)            # length C
jac_po_sd <- apply(jac_po_mat, 1, sd)

cat("=== Results (mean over B reps and C classes) ===\n\n")
cat("Method A — Class-level NPKDC-vd:\n")
cat(sprintf("  Recall    = %.3f (%.3f)\n", mean(cl_mn[,2]), mean(cl_sd[,2])))
cat(sprintf("  Precision = %.3f (%.3f)\n", mean(cl_mn[,1]), mean(cl_sd[,1])))
cat(sprintf("  F1        = %.3f (%.3f)\n", mean(cl_mn[,3]), mean(cl_sd[,3])))
cat(sprintf("  P(exact)  = %.3f (%.3f)\n", mean(cl_mn[,4]), mean(cl_sd[,4])))
cat(sprintf("  Within-class Jaccard = 1.000 (class-level is constant by construction)\n\n"))

cat("Method B — Per-observation attribution (top-3 by φ_j):\n")
cat(sprintf("  Recall    = %.3f (%.3f)\n", mean(po_mn[,2]), mean(po_sd[,2])))
cat(sprintf("  Precision = %.3f (%.3f)\n", mean(po_mn[,1]), mean(po_sd[,1])))
cat(sprintf("  F1        = %.3f (%.3f)\n", mean(po_mn[,3]), mean(po_sd[,3])))
cat(sprintf("  P(exact)  = %.3f (%.3f)\n", mean(po_mn[,4]), mean(po_sd[,4])))
cat(sprintf("  Within-class Jaccard = %.3f (%.3f)\n\n",
            mean(jac_po_mn), mean(jac_po_sd)))

# ---- Save ----
saveRDS(list(reps=reps, cl_mn=cl_mn, cl_sd=cl_sd,
             po_mn=po_mn, po_sd=po_sd,
             jac_po_mn=jac_po_mn, jac_po_sd=jac_po_sd),
        file.path("..", "data", "N3_results.rds"),
        compress = FALSE)
cat("Results saved to N3_results.rds\n\n")


# ---- LaTeX table ----
cat("========== LATEX TABLE ==========\n")
cat("\\begin{tabular}{lcccc}\\hline\n")
cat("Method & Precision & Recall & F1 & Within-class Jaccard \\\\\n")
cat("\\hline\\hline\n")
cat(sprintf("Class-level (NPKDC-vd) & %.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f) & $1.000$ \\\\\n",
            mean(cl_mn[,1]), mean(cl_sd[,1]),
            mean(cl_mn[,2]), mean(cl_sd[,2]),
            mean(cl_mn[,3]), mean(cl_sd[,3])))
cat(sprintf("Per-observation ($\\hat{A}_{y,i}$) & %.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f) \\\\\n",
            mean(po_mn[,1]), mean(po_sd[,1]),
            mean(po_mn[,2]), mean(po_sd[,2]),
            mean(po_mn[,3]), mean(po_sd[,3]),
            mean(jac_po_mn), mean(jac_po_sd)))
cat("\\hline\n\\end{tabular}\n\n")


# ---- Figure: distribution of per-obs φ values for one representative rep ----
# Show φ_j for all test obs of each class (for one rep), overlaid with class-level threshold
rep1 <- reps[[1]]
set.seed(1 * 999 + 42)
dat1 <- gen_data(N_Y, N_TEST)

bw_sil1 <- matrix(0, C, D)
for (y in 1:C) bw_sil1[y, ] <- apply(dat1$train[[y]], 2L, bw_silverman)

# Compute φ matrix for class 1 test observations (all 100)
phi_c1 <- t(apply(dat1$test[[1]], 1, function(xi)
  per_obs_contrast(xi, 1L, dat1$train, bw_sil1)))
# phi_c1: 100 x 20 matrix

# Class-level detection for class 1
cl_det1 <- detect_class_level(dat1$train)
relevant_c1 <- which(cl_det1[1, ] == 1L)  # detected relevant for class 1

fig_path <- file.path("..", "manuscript", "figures", "N3_contrast.eps")
postscript(fig_path, width = 6, height = 4.2,
           horizontal = FALSE, onefile = FALSE, paper = "special")

# Boxplot of φ values per variable for class 1
col_rel <- "#1b7837"; col_irrel <- "#aaaaaa"; col_det <- "#d6604d"
box_cols <- ifelse(1:D %in% REL_SETS[[1]], col_rel, col_irrel)
# Overwrite detected-but-wrong in orange (unlikely with our design)
box_cols[setdiff(relevant_c1, REL_SETS[[1]])] <- col_det

par(mar = c(4, 4.2, 1.5, 1))
boxplot(phi_c1, at = 1:D, col = box_cols, border = box_cols,
        xlab = "Variable index", ylab = expression(phi[j](x[i], y)),
        main = "",
        whisklty = 1, outpch = 20, outcex = 0.4, outcol = box_cols,
        cex.axis = 0.85, cex.lab = 0.95, las = 1)
abline(h = 0, col = "grey60", lty = 2)
# Mark class-level detected variables
rug(relevant_c1, side = 3, ticksize = 0.08, col = "#d6604d", lwd = 2)
legend("topright", bty = "n", cex = 0.80,
       legend = c(expression("True"~R[y]~"(3 vars)"),
                  "Irrelevant",
                  "Class-level detected"),
       fill = c(col_rel, col_irrel, col_det),
       border = NA)
dev.off()
cat("Figure saved to N3_contrast.eps\n")
