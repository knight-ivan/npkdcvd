# 02_benchmark_methods.R
# Modern class-specific feature-discovery benchmarks in R.

fit_ovr_lasso <- function(train_list, test_list, R_Y, alpha = 1, nfolds = 5) {
  if (!requireNamespace("glmnet", quietly = TRUE)) return(NULL)
  C <- length(train_list); D <- ncol(train_list[[1L]])
  X <- do.call(rbind, train_list)
  ylab <- rep(seq_len(C), times = vapply(train_list, nrow, integer(1)))
  Xte <- do.call(rbind, test_list)
  true <- rep(seq_len(C), times = vapply(test_list, nrow, integer(1)))
  selected <- matrix(0L, C, D)
  score <- matrix(0, C, D)
  prob <- matrix(0, nrow = nrow(Xte), ncol = C)

  for (y in seq_len(C)) {
    yy <- as.integer(ylab == y)
    fold_n <- min(nfolds, table(yy))
    fit <- glmnet::cv.glmnet(X, yy, family = "binomial", alpha = alpha,
                             nfolds = max(3L, as.integer(fold_n)), standardize = TRUE)
    cf <- as.matrix(stats::coef(fit, s = "lambda.1se"))[-1L, 1L]
    score[y, ] <- abs(cf)
    selected[y, ] <- as.integer(abs(cf) > 1e-10)
    prob[, y] <- as.numeric(stats::predict(fit, newx = Xte, s = "lambda.1se", type = "response"))
  }
  pred <- max.col(prob, ties.method = "first")
  list(method = "OvR-LASSO", selected_native = selected, score = score,
       selected_topr = select_top_r_by_score(score, R_Y),
       pred = pred, true = true, accuracy = mean(pred == true),
       class_recall = classwise_recall(true, pred, C))
}

fit_ove_lasso <- function(train_list, test_list, R_Y, alpha = 1, nfolds = 5) {
  if (!requireNamespace("glmnet", quietly = TRUE)) return(NULL)
  C <- length(train_list); D <- ncol(train_list[[1L]])
  X_all <- do.call(rbind, train_list)
  ylab <- rep(seq_len(C), times = vapply(train_list, nrow, integer(1)))
  Xte <- do.call(rbind, test_list)
  true <- rep(seq_len(C), times = vapply(test_list, nrow, integer(1)))
  score <- matrix(0, C, D)
  votes <- matrix(0, nrow = nrow(Xte), ncol = C)

  for (y in seq_len(C)) {
    for (k in seq_len(C)) {
      if (k == y) next
      idx <- which(ylab %in% c(y, k))
      yy <- as.integer(ylab[idx] == y)
      fold_n <- min(nfolds, table(yy))
      fit <- glmnet::cv.glmnet(X_all[idx, , drop = FALSE], yy, family = "binomial",
                               alpha = alpha, nfolds = max(3L, as.integer(fold_n)),
                               standardize = TRUE)
      cf <- as.matrix(stats::coef(fit, s = "lambda.1se"))[-1L, 1L]
      score[y, ] <- score[y, ] + abs(cf) / (C - 1L)
      p_y <- as.numeric(stats::predict(fit, newx = Xte, s = "lambda.1se", type = "response"))
      pair_pred <- ifelse(p_y >= 0.5, y, k)
      votes[cbind(seq_along(pair_pred), pair_pred)] <- votes[cbind(seq_along(pair_pred), pair_pred)] + 1L
    }
  }
  selected_native <- (score > 1e-10) * 1L
  pred <- max.col(votes, ties.method = "first")
  list(method = "OvE-LASSO", selected_native = selected_native, score = score,
       selected_topr = select_top_r_by_score(score, R_Y),
       pred = pred, true = true, accuracy = mean(pred == true),
       class_recall = classwise_recall(true, pred, C))
}

fit_rbf_svm <- function(train_list, test_list, cost = 1, gamma = NULL) {
  if (!requireNamespace("e1071", quietly = TRUE)) return(NULL)
  C <- length(train_list)
  X <- do.call(rbind, train_list)
  ylab <- factor(rep(seq_len(C), times = vapply(train_list, nrow, integer(1))))
  Xte <- do.call(rbind, test_list)
  true <- rep(seq_len(C), times = vapply(test_list, nrow, integer(1)))
  args <- list(x = X, y = ylab, kernel = "radial", cost = cost, scale = TRUE)
  if (!is.null(gamma)) args$gamma <- gamma
  mod <- do.call(e1071::svm, args)
  pred <- as.integer(as.character(stats::predict(mod, Xte)))
  list(method = "RBF-SVM", pred = pred, true = true, accuracy = mean(pred == true),
       class_recall = classwise_recall(true, pred, C))
}

fit_xgb_treeshap <- function(train_list, test_list, R_Y,
                             nrounds = 60, max_depth = 3, eta = 0.1,
                             subsample = 0.9, colsample_bytree = 0.9) {
  if (!requireNamespace("xgboost", quietly = TRUE)) return(NULL)
  C <- length(train_list); D <- ncol(train_list[[1L]])
  X <- do.call(rbind, train_list)
  ylab <- rep(seq_len(C), times = vapply(train_list, nrow, integer(1))) - 1L
  Xte <- do.call(rbind, test_list)
  true <- rep(seq_len(C), times = vapply(test_list, nrow, integer(1)))
  dtrain <- xgboost::xgb.DMatrix(data = X, label = ylab)
  dtest <- xgboost::xgb.DMatrix(data = Xte)
  params <- list(objective = "multi:softprob", num_class = C, eval_metric = "mlogloss",
                 max_depth = max_depth, eta = eta, subsample = subsample,
                 colsample_bytree = colsample_bytree, nthread = 1)
  mod <- xgboost::xgb.train(params = params, data = dtrain, nrounds = nrounds, verbose = 0)
  prob_vec <- stats::predict(mod, dtest)
  prob <- matrix(prob_vec, ncol = C, byrow = TRUE)
  pred <- max.col(prob, ties.method = "first")

  contrib <- stats::predict(mod, dtest, predcontrib = TRUE)
  score <- matrix(0, C, D)
  # xgboost R may return either an array [n, class, feature+1] or a matrix.
  if (length(dim(contrib)) == 3L) {
    # Most common recent format: n x C x (D+1)
    dims <- dim(contrib)
    if (dims[2L] == C && dims[3L] >= D) {
      for (y in seq_len(C)) {
        idx <- which(true == y)
        score[y, ] <- apply(abs(contrib[idx, y, seq_len(D), drop = FALSE]), 3L, mean)
      }
    } else if (dims[3L] == C && dims[2L] >= D) {
      for (y in seq_len(C)) {
        idx <- which(true == y)
        score[y, ] <- apply(abs(contrib[idx, seq_len(D), y, drop = FALSE]), 2L, mean)
      }
    } else {
      stop("Unexpected xgboost predcontrib array dimensions: ", paste(dims, collapse = " x "))
    }
  } else {
    # Older format often flattens classes: n x ((D+1)*C)
    mat <- as.matrix(contrib)
    expected <- (D + 1L) * C
    if (ncol(mat) < expected) stop("Unexpected xgboost predcontrib matrix columns.")
    for (y in seq_len(C)) {
      idx <- which(true == y)
      class_cols <- ((y - 1L) * (D + 1L) + 1L):((y - 1L) * (D + 1L) + D)
      score[y, ] <- colMeans(abs(mat[idx, class_cols, drop = FALSE]))
    }
  }
  selected_topr <- select_top_r_by_score(score, R_Y)
  # TreeSHAP is a dense ranking; native sparse set is not intrinsic.
  selected_native <- selected_topr
  list(method = "XGBoost-TreeSHAP", selected_native = selected_native, score = score,
       selected_topr = selected_topr, pred = pred, true = true, accuracy = mean(pred == true),
       class_recall = classwise_recall(true, pred, C))
}
