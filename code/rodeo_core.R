# rodeo_core.R  --------------------------------------------------------------
# Correct density-rodeo for NPKDC-vd (Route 2 rebuild, 2026-08-24).
#
# This replaces the never-implemented Algorithm 1 with a genuine rodeo, fixing
# the exact objections in the referee reports:
#
#   * R1.C2  -- uses the DENSITY derivative  Z_j = d f_hat / d h_j
#               (NOT d log f_hat / d h_j), so E[Z_j] admits the direct
#               second-order expansion; no E[d log/dh] != d/dh E[.] gap.
#   * R1.C3  -- CROSS-FITTING: the KDE is built on a fit split; the statistic is
#               evaluated on a disjoint eval split, so the per-eval-point
#               statistics are conditionally independent given the fit split and
#               carry no O(K(0)/h) self-contribution. Variance is the genuine
#               contribution-level  Var(Z_{ji})/n_fit.
#   * R2.C3  -- NON-CANCELLING aggregation: variables are ranked by the mean of
#               the SQUARED standardised derivative over eval points, so
#               sign-changing curvature cannot average the signal to zero.
#   * R2.C4  -- selection is a DETERMINISTIC largest-gap rule on that statistic
#               (no ANOVA/Tukey distributional claim).
#   * R2.min2-- an optional marginal PRE-SCREEN keeps the joint product-KDE
#               numerically stable in high d; the sure-screening set is reported
#               so the theory can condition on {A0 \supseteq union_y R_y}.
#
# Numerics: all kernel weights are formed in log-space and stabilised by the
# per-eval-point row maximum, so the STANDARDISED statistic Z_j / sd(Z_j) is
# exactly scale-invariant and never underflows, even for large p.
#
# Parallelism: pure vectorised R (one outer() per coordinate per iteration);
# the driver parallelises over replications AND classes with fork-based
# parallel::mclapply (macOS-friendly, shared memory). See rodeo_detect() and the
# run_* helpers.
# ---------------------------------------------------------------------------

suppressWarnings(suppressMessages({
  library(parallel)
}))

## ---- product-Gaussian-kernel log weights ----------------------------------
# logW[l,i] = sum_k [ log phi((Xeval[l,k]-Xfit[i,k])/h[k]) - log h[k] ].
kde_logW <- function(Xeval, Xfit, h) {
  m <- nrow(Xeval); p <- ncol(Xfit)
  logC <- -0.5 * log(2 * pi)
  logW <- matrix(0.0, m, nrow(Xfit))
  for (k in seq_len(p)) {
    U <- outer(Xeval[, k], Xfit[, k], "-") / h[k]
    logW <- logW + (logC - 0.5 * U * U - log(h[k]))
  }
  logW
}

## ---- standardised rodeo statistic at each eval point ----------------------
# For every coordinate j returns the standardised derivative
#   t_j[l] = Z_j(x_l) / sqrt( Var_i(Z_{ji}(x_l)) / n_fit ),
# where Z_{ji}(x_l) = W_{li} (u_{lij}^2 - 1)/h_j and W is the FULL product kernel.
# Because W is stabilised per eval-point (row max removed), the common exp(M_l)
# scale cancels between Z_j and its sd, so t_j is scale-free (and underflow-free).
rodeo_tstat <- function(Xeval, Xfit, h) {
  m <- nrow(Xeval); n <- nrow(Xfit); p <- ncol(Xfit)
  logW <- kde_logW(Xeval, Xfit, h)
  Mrow <- apply(logW, 1L, max)
  Ws   <- exp(logW - Mrow)                       # m x n, in (0,1]
  tmat <- matrix(0.0, m, p)
  for (j in seq_len(p)) {
    U       <- outer(Xeval[, j], Xfit[, j], "-") / h[j]
    contrib <- Ws * ((U * U - 1) / h[j])         # m x n  (scaled Z_{ji})
    zbar    <- rowMeans(contrib)                 # scaled Z_j
    zvar    <- rowMeans(contrib * contrib) - zbar * zbar
    se      <- sqrt(pmax(zvar, 0) / n)
    tmat[, j] <- zbar / pmax(se, 1e-300)
  }
  tmat                                            # m x p, scale-free
}

## ---- marginal pre-screen (sure-screening set A0) --------------------------
# Keep, per class, the coordinates most likely to carry class-specific
# structure: rank by a robust marginal-concentration score (small within-class
# sd relative to the pooled sd = candidate relevant), keep the top `keep`.
# Returns an integer vector of retained column indices for class y.
rodeo_prescreen <- function(Xy, X_pool_sd, keep = 40) {
  p <- ncol(Xy)
  if (keep >= p) return(seq_len(p))
  sdy   <- apply(Xy, 2L, sd)
  score <- X_pool_sd / pmax(sdy, 1e-12)          # large => class y is tight here
  sort(order(score, decreasing = TRUE)[seq_len(keep)])
}

## ---- one class: cross-fitted iterative rodeo ------------------------------
# X : n_y x p matrix for one class (already pre-screened columns `cols`).
# Returns a length-p logical `selected`, plus the per-coordinate score and the
# final bandwidths (diagnostics).
rodeo_select_class <- function(X, cols = seq_len(ncol(X)),
                               c0 = 1.5, beta = 0.8, max_iter = 40,
                               lambda_mult = 2.0, gap_factor = 3.0,
                               max_frac = 0.5, cf_folds = 2L) {
  n <- nrow(X); p <- ncol(X)
  sdcol <- apply(X, 2L, sd)
  h0    <- rep(c0 * max(sdcol, 1e-6), p)          # common, oversmoothed start
  # cross-fitting folds: average the standardised-squared statistic across folds
  fold  <- ((sample.int(n) - 1L) %% cf_folds) + 1L
  T_acc <- matrix(0.0, cf_folds, p)
  h_acc <- numeric(p); nproc <- 0L                 # accumulate final bandwidths
  for (f in seq_len(cf_folds)) {
    Xfit <- X[fold != f, , drop = FALSE]
    Xev  <- X[fold == f, , drop = FALSE]
    if (nrow(Xfit) < 5L || nrow(Xev) < 2L) { T_acc[f, ] <- 0; next }
    h      <- h0
    active <- rep(TRUE, p)
    lambda <- lambda_mult * sqrt(2 * log(max(nrow(Xfit), 2L)))
    Tbest  <- rep(0.0, p)
    for (it in seq_len(max_iter)) {
      if (!any(active)) break
      tt <- rodeo_tstat(Xev, Xfit[, active, drop = FALSE], h[active])
      Tj <- colMeans(tt * tt)                     # non-cancelling aggregate
      Tbest[active] <- pmax(Tbest[active], Tj)
      keep_shrink <- Tj > (lambda * lambda)       # still significant => shrink
      idx <- which(active)
      h[idx[keep_shrink]] <- h[idx[keep_shrink]] * beta
      active[idx[!keep_shrink]] <- FALSE          # freeze the rest
    }
    T_acc[f, ] <- Tbest
    h_acc <- h_acc + h; nproc <- nproc + 1L
  }
  score <- colMeans(T_acc)                          # per-coordinate rodeo score
  bw    <- if (nproc > 0L) h_acc / nproc else h0    # fold-averaged final bandwidths
  # deterministic largest-gap selection (Decision D): split the sorted scores at
  # their largest gap, and accept the split only when that gap clearly dominates
  # typical gaps (scale-free: gap_factor x median gap) and the relevant cluster
  # is a sparse minority (<= max_frac of the screened coordinates).
  ord  <- order(score, decreasing = TRUE)
  ss   <- score[ord]
  selected <- rep(FALSE, p)
  if (length(ss) >= 2L) {
    gaps <- -diff(ss)                               # >= 0, gaps down the ranking
    kmax <- which.max(gaps)
    posg <- gaps[gaps > 0]
    med  <- if (length(posg)) stats::median(posg) else 0
    if (med <= 0) med <- 1e-12
    if (gaps[kmax] > gap_factor * med &&
        kmax <= max(1L, floor(p * max_frac))) {
      selected[ord[seq_len(kmax)]] <- TRUE
    }
  }
  list(selected = selected, score = score, cols = cols, bw = bw)
}

## ===========================================================================
## BASELINE-CORRECTED (density-ratio) rodeo  --  Milestone 2
## Target: R_y = { j : the ratio f_y / b depends on x_j }, b = shared background.
##
## Method (background whitening): the rodeo detects departure from a UNIFORM
## baseline, so map every coordinate through the pooled-background CDF B_j (the
## probability integral transform). On this quantile scale the shared background
## becomes Uniform, hence "f_y matches b in coordinate j"  <=>  "transformed
## coordinate is uniform"  <=>  rodeo score at the noise level. A location shift
## (same dispersion, shifted mean) whitens to a CONCENTRATION on the quantile
## scale, so it is caught too -- unlike a raw concentration/sd rule. This reuses
## the validated plain rodeo and avoids any unstable ratio of density estimates.
## Assumes a product background b(x)=prod_j b_j(x_j); state that in the paper.
## ===========================================================================

# Probability integral transform: map each column of X through the empirical CDF
# of the corresponding column of Xref (the pooled background). Background-
# distributed values map to ~Uniform(0,1).
pit_transform <- function(X, Xref) {
  p <- ncol(X); U <- matrix(0.0, nrow(X), p)
  for (j in seq_len(p)) {
    ref <- sort(Xref[, j])
    U[, j] <- findInterval(X[, j], ref) / (length(ref) + 1)
  }
  U
}

# Baseline-corrected detection: whiten every class through the shared pooled
# background, then run the plain rodeo. Returns a C x D 0/1 matrix.
rodeo_detect_bc <- function(train_list, keep = 40, mc.cores = 1L, ...) {
  X_all <- do.call(rbind, train_list)
  tl_u  <- lapply(train_list, function(Xy) pit_transform(Xy, X_all))
  rodeo_detect(tl_u, keep = keep, mc.cores = mc.cores, ...)
}

## ---- all classes -> C x D detection matrix (drop-in for the harness) ------
# train_list : list of C matrices (n_y x D). Parallel over classes.
# Returns a C x D 0/1 matrix, same shape as detect_gap()/detect_zscore().
rodeo_detect <- function(train_list, keep = 40, mc.cores = 1L, ...) {
  C <- length(train_list); D <- ncol(train_list[[1]])
  X_pool_sd <- apply(do.call(rbind, train_list), 2L, sd)
  one <- function(y) {
    Xy   <- train_list[[y]]
    cols <- rodeo_prescreen(Xy, X_pool_sd, keep = keep)
    res  <- rodeo_select_class(Xy[, cols, drop = FALSE], cols = cols, ...)
    out  <- integer(D)
    out[cols[res$selected]] <- 1L
    out
  }
  rows <- if (mc.cores > 1L) mclapply(seq_len(C), one, mc.cores = mc.cores)
          else lapply(seq_len(C), one)
  do.call(rbind, rows)
}
