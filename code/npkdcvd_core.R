# npkdcvd_core.R
# ---------------------------------------------------------------------------
# Shared NPKDC-vd primitives, unified from N1_simulation.R and N2_anova_tukey.R
# so every robustness study (Theme E of the revision) reuses ONE detector.
#
# The detector is the marginal-concentration rule the paper is being realigned
# to (Decision B): the per-(class,variable) bandwidth h_{y,j} = Silverman rule
# is PROPORTIONAL to the within-class sd(X_j^y); a variable is "relevant" for
# class y when that bandwidth is anomalously small (class y is concentrated).
# Two selection rules are provided:
#   detect_zscore : cross-class z-score threshold  (graphical/diagnostic rule)
#   detect_gap    : within-class largest-gap rule   (Decision D; theory-backed,
#                   independent of the other classes)
# ---------------------------------------------------------------------------

## ---- marginal bandwidths ---------------------------------------------------
bw_silverman <- function(x) {
  s <- stats::sd(x)
  if (s < 1e-12) return(1e-6)
  1.06 * s * length(x)^(-0.2)
}

bw_loocv <- function(x,
                     grid_mult = c(0.05, 0.10, 0.20, 0.35, 0.50, 0.70, 1.00, 1.50, 2.00)) {
  n <- length(x); s <- stats::sd(x)
  if (s < 1e-12) return(1e-6)
  h_grid <- (1.06 * s * n^(-0.2)) * grid_mult
  Dm <- outer(x, x, "-")
  scores <- vapply(h_grid, function(h) {
    K <- stats::dnorm(Dm / h) / h; diag(K) <- 0
    -mean(log(pmax(rowSums(K) / (n - 1L), 1e-20)))
  }, numeric(1))
  h_grid[which.min(scores)]
}

# C x D matrix of marginal bandwidths (one row per class).
bw_matrix <- function(train_list, bw_fun = bw_silverman) {
  C <- length(train_list); D <- ncol(train_list[[1]])
  bw <- matrix(0, C, D)
  for (y in seq_len(C)) bw[y, ] <- apply(train_list[[y]], 2L, bw_fun)
  bw
}

## ---- selection rules -------------------------------------------------------
# Cross-class z-score rule (diagnostic): flag var j for class y when its
# bandwidth is > |tau| SD below the cross-class mean bandwidth for var j.
detect_zscore <- function(train_list, tau = -1.5, bw_fun = bw_silverman) {
  bw <- bw_matrix(train_list, bw_fun)
  C <- nrow(bw); D <- ncol(bw)
  det <- matrix(0L, C, D)
  for (j in seq_len(D)) {
    hj <- bw[, j]; sj <- stats::sd(hj)
    if (sj < 1e-10) next
    det[, j] <- as.integer((hj - mean(hj)) / sj < tau)
  }
  det
}

# Within-class deterministic largest-gap rule (Decision D). For each class,
# sort its D bandwidths, split at the largest gap, and accept the split only if
# the gap exceeds the order of the bandwidth estimation noise,
# c_gap * sqrt(log D / n_y) * (median bandwidth as scale). Independent of the
# other classes -> matches "identity independent of all other classes".
detect_gap <- function(train_list, c_gap = 1.0, bw_fun = bw_silverman) {
  C <- length(train_list); D <- ncol(train_list[[1]])
  bw <- bw_matrix(train_list, bw_fun)
  det <- matrix(0L, C, D)
  for (y in seq_len(C)) {
    n_y <- nrow(train_list[[y]])
    h   <- bw[y, ]; ord <- order(h); hs <- h[ord]
    gaps <- diff(hs)
    if (!length(gaps)) next
    thresh <- c_gap * sqrt(log(max(D, 2)) / n_y) * stats::median(hs)
    # Relevant = the concentrated (smallest-bandwidth) cluster, so split at the
    # FIRST gap from below that exceeds the noise scale -- NOT the global max,
    # which can fall between the background and an even-more-diffuse outlier.
    big <- which(gaps > thresh)
    if (length(big)) det[y, ord[seq_len(big[1])]] <- 1L
  }
  det
}

## ---- plug-in product-KDE classifier ---------------------------------------
classify_npkdc <- function(X_test, train_list, bw = NULL,
                           bw_fun = bw_loocv, priors = NULL) {
  C <- length(train_list); D <- ncol(X_test)
  if (is.null(bw))     bw     <- bw_matrix(train_list, bw_fun)
  if (is.null(priors)) priors <- vapply(train_list, nrow, numeric(1))
  n_te <- nrow(X_test)
  logpost <- matrix(0, n_te, C)
  for (y in seq_len(C)) {
    Xtr <- train_list[[y]]; hy <- bw[y, ]
    logf <- matrix(0, n_te, D)
    for (j in seq_len(D)) {
      Kj <- stats::dnorm(outer(X_test[, j], Xtr[, j], "-") / hy[j]) / hy[j]
      logf[, j] <- log(pmax(rowMeans(Kj), 1e-300))
    }
    logpost[, y] <- rowSums(logf) + log(priors[y])
  }
  max.col(logpost, ties.method = "first")
}

## ---- attribution metrics ---------------------------------------------------
attr_metrics <- function(detected, rel_sets) {
  C <- nrow(detected); D <- ncol(detected)
  out <- matrix(0, C, 4, dimnames = list(NULL, c("prec", "rec", "f1", "exact")))
  for (y in seq_len(C)) {
    rel <- rel_sets[[y]]; irrel <- setdiff(seq_len(D), rel)
    tp <- sum(detected[y, rel] == 1L); fp <- sum(detected[y, irrel] == 1L)
    fn <- sum(detected[y, rel] == 0L)
    prec <- if (tp + fp > 0) tp / (tp + fp) else 0
    rec  <- if (tp + fn > 0) tp / (tp + fn) else 0
    f1   <- if (prec + rec > 0) 2 * prec * rec / (prec + rec) else 0
    out[y, ] <- c(prec, rec, f1, as.numeric(tp == length(rel) && fp == 0L))
  }
  out
}
