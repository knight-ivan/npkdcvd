# 01_data_and_core_npkdc.R
# Core data-generation, NPKDC-vd detection, KDE classification, and metrics.

make_rel_sets <- function(R_Y) {
  rel_start <- cumsum(c(1L, R_Y[-length(R_Y)]))
  lapply(seq_along(R_Y), function(y) rel_start[y]:(rel_start[y] + R_Y[y] - 1L))
}

gen_highdim_sparse_data <- function(C = 5L, D = 200L, R_Y = c(3,5,5,8,8),
                                    n_y = 150L, n_test = 100L,
                                    irrel_type = c("uniform", "normal", "ar1"),
                                    rel_sd = 0.1, rho = 0.5) {
  irrel_type <- match.arg(irrel_type)
  REL_SETS <- make_rel_sets(R_Y)

  make_irrel <- function(nrow, ncol) {
    if (ncol <= 0L) return(matrix(numeric(0), nrow = nrow, ncol = 0L))
    if (irrel_type == "uniform") {
      matrix(runif(nrow * ncol), nrow = nrow, ncol = ncol)
    } else if (irrel_type == "normal") {
      matrix(rnorm(nrow * ncol), nrow = nrow, ncol = ncol)
    } else {
      Z <- matrix(0, nrow = nrow, ncol = ncol)
      Z[, 1L] <- rnorm(nrow)
      if (ncol >= 2L) {
        for (j in 2L:ncol) Z[, j] <- rho * Z[, j - 1L] + sqrt(1 - rho^2) * rnorm(nrow)
      }
      Z
    }
  }

  make_one_class <- function(n, y) {
    rel <- REL_SETS[[y]]
    nr <- length(rel)
    X <- make_irrel(n, D)
    mu_rel <- rep(c(0.3, 0.5, 0.7), length.out = nr)
    X[, rel] <- matrix(rnorm(n * nr, mean = rep(mu_rel, each = n), sd = rel_sd), nrow = n, ncol = nr)
    X
  }

  train <- lapply(seq_len(C), function(y) make_one_class(n_y, y))
  test  <- lapply(seq_len(C), function(y) make_one_class(n_test, y))
  list(train = train, test = test, rel_sets = REL_SETS, R_Y = R_Y, C = C, D = D,
       n_y = n_y, n_test = n_test, irrel_type = irrel_type)
}

bw_silverman <- function(x) {
  s <- stats::sd(x)
  if (!is.finite(s) || s < 1e-12) return(1e-6)
  1.06 * s * length(x)^(-0.2)
}

npkdc_zscores <- function(train_list) {
  C <- length(train_list); D <- ncol(train_list[[1L]])
  bw_sil <- matrix(0, nrow = C, ncol = D)
  for (y in seq_len(C)) bw_sil[y, ] <- apply(train_list[[y]], 2L, bw_silverman)
  z <- matrix(0, nrow = C, ncol = D)
  for (j in seq_len(D)) {
    h <- bw_sil[, j]
    sj <- stats::sd(h)
    if (!is.finite(sj) || sj < 1e-10) next
    z[, j] <- (h - mean(h)) / sj
  }
  z
}

detect_npkdc_vd <- function(train_list, tau = -1.5) {
  z <- npkdc_zscores(train_list)
  detected <- (z < tau) * 1L
  list(selected_native = detected, score = -z, zscore = z)
}

select_top_r_by_score <- function(score, R_Y) {
  C <- nrow(score); D <- ncol(score)
  out <- matrix(0L, nrow = C, ncol = D)
  for (y in seq_len(C)) {
    rr <- min(R_Y[y], D)
    idx <- order(score[y, ], decreasing = TRUE)[seq_len(rr)]
    out[y, idx] <- 1L
  }
  out
}

kde1d_logdensity <- function(query, train_x, h) {
  h <- max(h, 1e-6)
  # Vectorized over query; memory safe for the study sizes used here.
  Dmat <- outer(query, train_x, "-") / h
  log(pmax(rowMeans(stats::dnorm(Dmat) / h), 1e-300))
}

classify_npkdc_selected <- function(test_list, train_list, selected_by_class, use_fallback_all = TRUE) {
  C <- length(train_list); D <- ncol(train_list[[1L]])
  X_te <- do.call(rbind, test_list)
  true <- rep(seq_len(C), times = vapply(test_list, nrow, integer(1)))
  n_te <- nrow(X_te)
  log_post <- matrix(0, nrow = n_te, ncol = C)
  bw <- lapply(train_list, function(X) apply(X, 2L, bw_silverman))

  for (y in seq_len(C)) {
    vars <- which(selected_by_class[y, ] == 1L)
    if (length(vars) == 0L && use_fallback_all) vars <- seq_len(D)
    if (length(vars) == 0L) {
      log_post[, y] <- log(nrow(train_list[[y]]))
      next
    }
    lf <- numeric(n_te)
    for (j in vars) lf <- lf + kde1d_logdensity(X_te[, j], train_list[[y]][, j], bw[[y]][j])
    log_post[, y] <- lf + log(nrow(train_list[[y]]))
  }
  pred <- max.col(log_post, ties.method = "first")
  list(pred = pred, true = true, accuracy = mean(pred == true), class_recall = classwise_recall(true, pred, C))
}

feature_metrics_one <- function(selected_vec, true_set, D) {
  selected_set <- which(selected_vec == 1L)
  tp <- length(intersect(selected_set, true_set))
  fp <- length(setdiff(selected_set, true_set))
  fn <- length(setdiff(true_set, selected_set))
  precision <- if ((tp + fp) > 0L) tp / (tp + fp) else 0
  recall <- if ((tp + fn) > 0L) tp / (tp + fn) else 0
  f1 <- if ((precision + recall) > 0) 2 * precision * recall / (precision + recall) else 0
  exact <- as.numeric(length(selected_set) == length(true_set) && all(sort(selected_set) == sort(true_set)))
  c(precision = precision, recall = recall, f1 = f1, exact = exact,
    n_selected = length(selected_set))
}

feature_metrics_by_class <- function(selected_mat, rel_sets) {
  C <- length(rel_sets); D <- ncol(selected_mat)
  out <- t(vapply(seq_len(C), function(y) feature_metrics_one(selected_mat[y, ], rel_sets[[y]], D), numeric(5)))
  rownames(out) <- paste0("class", seq_len(C))
  out
}

classwise_recall <- function(true, pred, C) {
  vapply(seq_len(C), function(y) {
    idx <- which(true == y)
    if (length(idx) == 0L) NA_real_ else mean(pred[idx] == y)
  }, numeric(1))
}

selection_frequency_update <- function(counts, selected_mat, method, version) {
  key <- paste(method, version, sep = "|")
  if (is.null(counts[[key]])) counts[[key]] <- matrix(0L, nrow = nrow(selected_mat), ncol = ncol(selected_mat))
  counts[[key]] <- counts[[key]] + selected_mat
  counts
}
