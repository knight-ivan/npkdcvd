# reviewer_metrics.R: audit outputs for ASA Discoveries R2.
# Uncertainty is computed from INDEPENDENT Monte Carlo replications, averaging
# across classes within a replication. Classes within a replicate are dependent.
reviewer_script_dir <- function() {
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[1L]))) else getwd()
}
reviewer_output_dir <- function() {
  # RODEO_OUT redirects all outputs (run_R2.sh smoke mode uses it so that a smoke
  # test can never overwrite the authoritative B=1000 results in results_R2/).
  path <- Sys.getenv("RODEO_OUT", file.path(reviewer_script_dir(), "..", "results_R2"))
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  path
}
reviewer_cores <- function() {
  z <- suppressWarnings(as.integer(Sys.getenv("RODEO_CORES", "1")))
  if (is.na(z) || z < 1L) stop("RODEO_CORES must be a positive integer")
  if (.Platform$OS.type == "windows") return(1L)
  min(z, parallel::detectCores())
}
reviewer_B <- function(default = 1000L) {
  z <- suppressWarnings(as.integer(Sys.getenv("RODEO_B", as.character(default))))
  if (is.na(z) || z < 1L) stop("RODEO_B must be a positive integer")
  z
}
reviewer_mc <- function(ids, fn) {
  parallel::mclapply(ids, fn, mc.cores = reviewer_cores(), mc.set.seed = FALSE)
}
reviewer_metric_rows <- function(det, truth, replication, group = "", n_y = NA_integer_) {
  if (nrow(det) != length(truth)) stop("truth/detector dimension mismatch")
  metrics <- attr_metrics(det, truth)
  C <- nrow(metrics); D <- ncol(det)
  data.frame(replication = rep(replication, C), setting = rep(group, C),
    n_y = rep(n_y, C), class = seq_len(C),
    precision = metrics[, "prec"], recall = metrics[, "rec"],
    f1 = metrics[, "f1"], exact = metrics[, "exact"],
    fp = vapply(seq_len(C), function(y)
      sum(det[y, setdiff(seq_len(D), truth[[y]])] != 0), numeric(1)),
    fn = vapply(seq_len(C), function(y)
      sum(det[y, truth[[y]]] == 0), numeric(1)),
    selected = rowSums(det), stringsAsFactors = FALSE)
}
reviewer_summary <- function(rows, by = c("setting", "n_y"), metrics = c("precision", "recall", "f1", "exact")) {
  # Group into independent replication means; then form t-based MC intervals.
  if (!nrow(rows)) return(data.frame())
  keys <- unique(rows[, by, drop = FALSE])
  ans <- list()
  for (k in seq_len(nrow(keys))) {
    use <- rep(TRUE, nrow(rows))
    for (v in by) use <- use & (rows[[v]] == keys[[v]][k])
    sub <- rows[use, , drop = FALSE]
    for (metric in metrics) {
      repmeans <- tapply(sub[[metric]], sub$replication, mean)
      N <- length(repmeans); mn <- mean(repmeans)
      se <- if (N > 1L) sd(repmeans) / sqrt(N) else NA_real_
      crit <- if (N > 1L) qt(0.975, N - 1L) else NA_real_
      ans[[length(ans)+1L]] <- data.frame(keys[k, , drop = FALSE],
        metric = metric, B = N, estimate = mn, mc_se = se,
        ci_low = if (N > 1L) max(0, mn - crit*se) else NA_real_,
        ci_high = if (N > 1L) min(1, mn + crit*se) else NA_real_)
    }
  }
  do.call(rbind, ans)
}
# Newly specified, transparent accuracy gate. Historical labels in older PDFs
# must NOT be represented as if they had been generated with this rule.
# Accuracy is macro-averaged over test classes to guard against imbalance.
stage1_gate <- function(accuracy, threshold = 0.60) {
  if (!is.finite(accuracy) || accuracy < 0 || accuracy > 1 ||
      !is.finite(threshold) || threshold <= 0 || threshold >= 1)
    stop("accuracy/threshold out of range")
  accuracy >= threshold
}
