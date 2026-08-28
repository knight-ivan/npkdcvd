# rodeo_study4.R -------------------------------------------------------------
# Study 4 (N3 design): class-level vs per-observation attribution stability,
# under the rebuilt rodeo. C=5, d=20, r_y=3; relevant N(0.5,0.5^2), irrel N(0,1).
#
# Answers report section 4: the class-level within-class Jaccard=1 is MECHANICAL
# (one set per class reused for every obs). We replace it with a MEANINGFUL
# measure -- reproducibility of R_hat_y across independently RESAMPLED training
# sets (mean pairwise Jaccard over B fresh draws). Per-observation attribution
# keeps the log-density-contrast baseline (unstable, for contrast).
# ---------------------------------------------------------------------------
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(here) || !nzchar(here)) here <- "."
source(file.path(here, "npkdcvd_core.R"))
source(file.path(here, "rodeo_core.R"))

C <- 5; D <- 20; N_Y <- 150; N_TEST <- 100
R_Y <- rep(3L, C); rel_start <- cumsum(c(1L, R_Y[-C]))
REL_SETS <- lapply(1:C, function(y) rel_start[y]:(rel_start[y] + R_Y[y] - 1))

gen_one <- function(n) lapply(1:C, function(y) {
  X <- matrix(rnorm(n * D), n, D)
  X[, REL_SETS[[y]]] <- matrix(rnorm(n * R_Y[y], 0.5, 0.5), n, R_Y[y]); X
})
bw_sil <- function(x) { s <- sd(x); if (s < 1e-12) 1e-6 else 1.06 * s * length(x)^(-0.2) }
kde1d_log <- function(q, tr, h) log(mean(dnorm((q - tr) / h) / h) + 1e-300)
jacc <- function(a, b) { u <- sum(a | b); if (u == 0) 1 else sum(a & b) / u }

one_rep <- function(b) {
  set.seed(b * 999 + 42)
  tr <- gen_one(N_Y); te <- gen_one(N_TEST)
  cl_det <- rodeo_detect(tr, keep = D, mc.cores = 1L)          # class-level (rodeo)
  cl_acc <- attr_metrics(cl_det, REL_SETS)                     # C x 4
  bw <- t(sapply(1:C, function(y) apply(tr[[y]], 2L, bw_sil))) # C x D
  # per-observation attribution: top-r_y by log-density contrast phi_j
  po_acc <- matrix(0, C, 4); jac_po <- numeric(C)
  for (y in 1:C) {
    Xte <- te[[y]]
    A <- matrix(0L, N_TEST, D)
    for (i in 1:N_TEST) {
      phi <- vapply(1:D, function(j) {
        ld <- vapply(1:C, function(yp) kde1d_log(Xte[i, j], tr[[yp]][, j], bw[yp, j]), numeric(1))
        ld[y] - mean(ld)
      }, numeric(1))
      A[i, order(phi, decreasing = TRUE)[1:R_Y[y]]] <- 1L
    }
    acc_i <- t(apply(A, 1, function(r) {
      rel <- REL_SETS[[y]]; ir <- setdiff(1:D, rel)
      tp <- sum(r[rel]); fp <- sum(r[ir]); fn <- sum(r[rel] == 0)
      pr <- if (tp + fp) tp / (tp + fp) else 0; rc <- if (tp + fn) tp / (tp + fn) else 0
      c(pr, rc, if (pr + rc) 2 * pr * rc / (pr + rc) else 0, as.numeric(tp == R_Y[y] && fp == 0))
    }))
    po_acc[y, ] <- colMeans(acc_i)
    pr <- matrix(sample(N_TEST, 2 * 200, replace = TRUE), ncol = 2)
    jac_po[y] <- mean(apply(pr, 1, function(ij) jacc(A[ij[1], ], A[ij[2], ])))
  }
  list(cl_acc = cl_acc, cl_det = cl_det, po_acc = po_acc, jac_po = jac_po)
}

B    <- as.integer(Sys.getenv("RODEO_B", "120"))
ncor <- as.integer(Sys.getenv("RODEO_CORES", as.character(detectCores())))
cat(sprintf("=== Study 4: class-level vs per-obs stability (rodeo), B=%d on %d cores ===\n", B, ncor))
out <- mclapply(seq_len(B), one_rep, mc.cores = ncor)

cl <- apply(simplify2array(lapply(out, `[[`, "cl_acc")), c(1, 2), mean)  # C x 4
po <- apply(simplify2array(lapply(out, `[[`, "po_acc")), c(1, 2), mean)
jac_po <- mean(sapply(out, `[[`, "jac_po"))

# cross-training-set reproducibility of R_hat_y (the meaningful stability measure)
rep_jac <- numeric(C)
for (y in 1:C) {
  M <- t(sapply(out, function(o) o$cl_det[y, ]))        # B x D binary
  pr <- matrix(sample(nrow(M), 2 * 500, replace = TRUE), ncol = 2)
  pr <- pr[pr[, 1] != pr[, 2], ]
  rep_jac[y] <- mean(apply(pr, 1, function(ij) jacc(M[ij[1], ], M[ij[2], ])))
}

cat(sprintf("\nClass-level (rodeo):   prec=%.3f rec=%.3f F1=%.3f exact=%.3f\n",
            mean(cl[, 1]), mean(cl[, 2]), mean(cl[, 3]), mean(cl[, 4])))
cat(sprintf("Per-observation:       prec=%.3f rec=%.3f F1=%.3f exact=%.3f  within-class Jaccard=%.3f\n",
            mean(po[, 1]), mean(po[, 2]), mean(po[, 3]), mean(po[, 4]), jac_po))
cat(sprintf("\nMEANINGFUL stability -- R_hat_y reproducibility across resampled training sets:\n"))
cat(sprintf("  mean pairwise Jaccard = %.3f  (per class: %s)\n",
            mean(rep_jac), paste(sprintf("%.3f", rep_jac), collapse = " ")))
