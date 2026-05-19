# 00_setup.R
# Reproducibility setup for the NPKDC-vd simulation studies.
# This file intentionally avoids automatic package installation so that
# manuscript results are reproducible in a controlled R environment.

npkdc_required_packages <- c("glmnet", "e1071", "xgboost")

check_npkdc_packages <- function(stop_if_missing = FALSE) {
  availability <- vapply(npkdc_required_packages, requireNamespace, logical(1), quietly = TRUE)
  missing <- names(availability)[!availability]
  if (length(missing) > 0L) {
    msg <- paste0(
      "Missing optional benchmark package(s): ", paste(missing, collapse = ", "),
      ". NPKDC-vd core routines will still run, but the corresponding benchmark(s) will be skipped. ",
      "Install them with install.packages(c(",
      paste(sprintf('"%s"', missing), collapse = ", "), ")) if you want the full benchmark."
    )
    if (stop_if_missing) stop(msg, call. = FALSE) else message(msg)
  }
  availability
}

make_output_dirs <- function(out_dir) {
  dirs <- c(out_dir, file.path(out_dir, "results"), file.path(out_dir, "figures"), file.path(out_dir, "tables"), file.path(out_dir, "logs"))
  for (d in dirs) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  invisible(dirs)
}

safe_mean <- function(x) if (length(x) == 0L || all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
safe_sd <- function(x) if (length(x) <= 1L || all(is.na(x))) NA_real_ else sd(x, na.rm = TRUE)

# Deterministic color palette for base-R figures.
npkdc_method_cols <- c(
  "NPKDC-vd" = "#1b7837",
  "OvR-LASSO" = "#762a83",
  "OvE-LASSO" = "#d6604d",
  "XGBoost-TreeSHAP" = "#4393c3",
  "RBF-SVM" = "#666666"
)
