# 03_run_N4_modern_comparison_B1000.R
# R counterpart of the N4 modern class-specific feature-discovery comparison.
# Usage:
#   Rscript R/03_run_N4_modern_comparison_B1000.R --B 1000 --out results/N4_R_B1000
# For a quick test:
#   Rscript R/03_run_N4_modern_comparison_B1000.R --B 2 --small TRUE --out results/smoke_N4

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0L && !is.na(a)) a else b
cmd_all <- commandArgs(FALSE)
file_arg <- cmd_all[grep("^--file=", cmd_all)]
this_file <- sub("^--file=", "", file_arg[1] %||% "R/03_run_N4_modern_comparison_B1000.R")
script_dir <- dirname(normalizePath(this_file, mustWork = FALSE))
source(file.path(script_dir, "00_setup.R"))
source(file.path(script_dir, "01_data_and_core_npkdc.R"))
source(file.path(script_dir, "02_benchmark_methods.R"))

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(B = 1000L, out = file.path(getwd(), "N4_R_output"), seed = 20260512L,
              small = FALSE, irrel_type = "normal", tau = -1.5)
  if (length(args) == 0L) return(out)
  for (i in seq(1L, length(args), by = 2L)) {
    key <- sub("^--", "", args[i]); val <- args[i + 1L]
    if (key == "B") out$B <- as.integer(val)
    if (key == "out") out$out <- val
    if (key == "seed") out$seed <- as.integer(val)
    if (key == "small") out$small <- tolower(val) %in% c("true", "1", "yes")
    if (key == "irrel_type") out$irrel_type <- val
    if (key == "tau") out$tau <- as.numeric(val)
  }
  out
}

args <- parse_args()
make_output_dirs(args$out)
set.seed(args$seed)
pkg_available <- check_npkdc_packages(stop_if_missing = FALSE)

if (args$small) {
  C <- 3L; D <- 30L; R_Y <- c(3L, 4L, 5L); n_y <- 35L; n_test <- 12L
  nrounds_xgb <- 10L
} else {
  C <- 5L; D <- 200L; R_Y <- c(3L, 5L, 5L, 8L, 8L); n_y <- 150L; n_test <- 100L
  nrounds_xgb <- 60L
}

rep_rows <- list(); class_rows <- list(); sel_counts <- list()
methods_seen <- character(0)

for (b in seq_len(args$B)) {
  if (b %% max(1L, floor(args$B / 10L)) == 0L) message("rep ", b, " / ", args$B)
  set.seed(args$seed + 1000L * b)
  dat <- gen_highdim_sparse_data(C = C, D = D, R_Y = R_Y, n_y = n_y, n_test = n_test,
                                 irrel_type = args$irrel_type)
  rel_sets <- dat$rel_sets

  # NPKDC-vd
  det <- detect_npkdc_vd(dat$train, tau = args$tau)
  np_topr <- select_top_r_by_score(det$score, R_Y)
  pred_np <- classify_npkdc_selected(dat$test, dat$train, np_topr)
  methods <- list(list(method = "NPKDC-vd", selected_native = det$selected_native,
                       selected_topr = np_topr, pred = pred_np$pred, true = pred_np$true,
                       accuracy = pred_np$accuracy, class_recall = pred_np$class_recall))

  # Optional benchmarks
  ovrl <- fit_ovr_lasso(dat$train, dat$test, R_Y)
  if (!is.null(ovrl)) methods[[length(methods) + 1L]] <- ovrl
  ovel <- fit_ove_lasso(dat$train, dat$test, R_Y)
  if (!is.null(ovel)) methods[[length(methods) + 1L]] <- ovel
  xgb <- fit_xgb_treeshap(dat$train, dat$test, R_Y, nrounds = nrounds_xgb)
  if (!is.null(xgb)) methods[[length(methods) + 1L]] <- xgb
  svm <- fit_rbf_svm(dat$train, dat$test)
  if (!is.null(svm)) methods[[length(methods) + 1L]] <- svm

  for (m in methods) {
    methods_seen <- union(methods_seen, m$method)
    # Feature metrics only for methods with selected sets.
    if (!is.null(m$selected_topr)) {
      for (version in c("topr", "native")) {
        smat <- if (version == "topr") m$selected_topr else m$selected_native
        fmet <- feature_metrics_by_class(smat, rel_sets)
        sel_counts <- selection_frequency_update(sel_counts, smat, m$method, version)
        for (y in seq_len(C)) {
          class_rows[[length(class_rows) + 1L]] <- data.frame(
            rep = b, method = m$method, version = version, class = y,
            accuracy = m$accuracy, class_recall = m$class_recall[y],
            precision = fmet[y, "precision"], recall = fmet[y, "recall"],
            f1 = fmet[y, "f1"], exact = fmet[y, "exact"],
            n_selected = fmet[y, "n_selected"]
          )
        }
      }
    } else {
      for (y in seq_len(C)) {
        class_rows[[length(class_rows) + 1L]] <- data.frame(
          rep = b, method = m$method, version = "prediction_only", class = y,
          accuracy = m$accuracy, class_recall = m$class_recall[y],
          precision = NA_real_, recall = NA_real_, f1 = NA_real_, exact = NA_real_,
          n_selected = NA_real_
        )
      }
    }
    rep_rows[[length(rep_rows) + 1L]] <- data.frame(rep = b, method = m$method, accuracy = m$accuracy)
  }
}

rep_df <- do.call(rbind, rep_rows)
class_df <- do.call(rbind, class_rows)
utils::write.csv(rep_df, file.path(args$out, "results", "n4_R_replication_metrics.csv"), row.names = FALSE)
utils::write.csv(class_df, file.path(args$out, "results", "n4_R_class_metrics_long.csv"), row.names = FALSE)

# Summaries
class_summary <- aggregate(cbind(accuracy, class_recall, precision, recall, f1, exact, n_selected) ~ method + version + class,
                           data = class_df, FUN = safe_mean)
utils::write.csv(class_summary, file.path(args$out, "results", "n4_R_classwise_summary.csv"), row.names = FALSE)
macro_summary <- aggregate(cbind(accuracy, class_recall, precision, recall, f1, exact) ~ method + version,
                           data = class_df, FUN = safe_mean)
utils::write.csv(macro_summary, file.path(args$out, "results", "n4_R_macro_summary.csv"), row.names = FALSE)

# Selection frequencies
sel_rows <- list()
for (key in names(sel_counts)) {
  parts <- strsplit(key, "\\|", fixed = FALSE)[[1L]]
  mat <- sel_counts[[key]] / args$B
  for (y in seq_len(nrow(mat))) {
    sel_rows[[length(sel_rows) + 1L]] <- data.frame(method = parts[1L], version = parts[2L], class = y,
                                                     feature = seq_len(ncol(mat)), freq = mat[y, ])
  }
}
sel_df <- do.call(rbind, sel_rows)
utils::write.csv(sel_df, file.path(args$out, "results", "n4_R_selection_frequency.csv"), row.names = FALSE)

metadata <- list(B = args$B, seed = args$seed, C = C, D = D, R_Y = R_Y, n_y = n_y, n_test = n_test,
                 irrel_type = args$irrel_type, tau = args$tau, methods_seen = methods_seen,
                 package_available = as.list(pkg_available))
saveRDS(metadata, file.path(args$out, "results", "n4_R_metadata.rds"))
message("Saved N4 R results to ", args$out)
