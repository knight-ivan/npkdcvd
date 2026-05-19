# 05_smoke_test.R
# Small end-to-end R smoke test. Designed to finish quickly.
# Usage from the package root:
#   Rscript R/05_smoke_test.R

this_file <- commandArgs(FALSE)[grep("--file=", commandArgs(FALSE))[1]]
if (length(this_file) == 0L || is.na(this_file)) this_file <- "R/05_smoke_test.R"
script_dir <- dirname(normalizePath(sub("--file=", "", this_file), mustWork = FALSE))
source(file.path(script_dir, "00_setup.R"))
source(file.path(script_dir, "01_data_and_core_npkdc.R"))
source(file.path(script_dir, "02_benchmark_methods.R"))

out <- file.path(dirname(script_dir), "smoke_test_output")
make_output_dirs(out)
log_file <- file.path(out, "logs", "smoke_test.log")
sink(log_file, split = TRUE)
cat("NPKDC-vd R smoke test\n")
cat("R version:", R.version.string, "\n")
cat("Time:", as.character(Sys.time()), "\n\n")

availability <- check_npkdc_packages(stop_if_missing = FALSE)
print(availability)

set.seed(20260512)
dat <- gen_highdim_sparse_data(C = 3, D = 30, R_Y = c(3,4,5), n_y = 35, n_test = 12, irrel_type = "normal")
rel_sets <- dat$rel_sets
cat("Generated small data: C=3, D=30, n_y=35, n_test=12\n")

det <- detect_npkdc_vd(dat$train, tau = -1.5)
np_topr <- select_top_r_by_score(det$score, dat$R_Y)
pred_np <- classify_npkdc_selected(dat$test, dat$train, np_topr)
print(list(NPKDC_accuracy = pred_np$accuracy, NPKDC_class_recall = pred_np$class_recall))
print(feature_metrics_by_class(np_topr, rel_sets))

ovr <- fit_ovr_lasso(dat$train, dat$test, dat$R_Y)
if (!is.null(ovr)) print(list(OvR_accuracy = ovr$accuracy)) else cat("OvR-LASSO skipped: glmnet unavailable.\n")

ove <- fit_ove_lasso(dat$train, dat$test, dat$R_Y)
if (!is.null(ove)) print(list(OvE_accuracy = ove$accuracy)) else cat("OvE-LASSO skipped: glmnet unavailable.\n")

svm <- fit_rbf_svm(dat$train, dat$test)
if (!is.null(svm)) print(list(SVM_accuracy = svm$accuracy)) else cat("RBF-SVM skipped: e1071 unavailable.\n")

xgb <- fit_xgb_treeshap(dat$train, dat$test, dat$R_Y, nrounds = 5)
if (!is.null(xgb)) print(list(XGB_accuracy = xgb$accuracy)) else cat("XGBoost-TreeSHAP skipped: xgboost unavailable.\n")

cat("\nSmoke test completed.\n")
sink()
cat("Smoke-test log written to ", log_file, "\n", sep = "")
