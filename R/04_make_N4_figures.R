# 04_make_N4_figures.R
# Base-R figures for the R N4 outputs.
# Usage: Rscript R/04_make_N4_figures.R --in results/N4_R_B1000

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0L && !is.na(a)) a else b

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(input = file.path(getwd(), "N4_R_output"))
  if (length(args) == 0L) return(out)
  for (i in seq(1L, length(args), by = 2L)) {
    key <- sub("^--", "", args[i]); val <- args[i + 1L]
    if (key %in% c("in", "input")) out$input <- val
  }
  out
}

args <- parse_args()
cmd_all <- commandArgs(FALSE)
file_arg <- cmd_all[grep("^--file=", cmd_all)]
this_file <- sub("^--file=", "", file_arg[1] %||% "R/04_make_N4_figures.R")
script_dir <- dirname(normalizePath(this_file, mustWork = FALSE))
source(file.path(script_dir, "00_setup.R"))
fig_dir <- file.path(args$input, "figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
sel <- utils::read.csv(file.path(args$input, "results", "n4_R_selection_frequency.csv"))
cls <- utils::read.csv(file.path(args$input, "results", "n4_R_classwise_summary.csv"))

plot_heatmap_one <- function(df, main = "") {
  mat <- xtabs(freq ~ class + feature, data = df)
  image(t(mat[nrow(mat):1, ]), axes = FALSE, col = colorRampPalette(c("white", "#deebf7", "#9ecae1", "#3182bd"))(101), zlim = c(0,1), main = main, xlab = "Feature", ylab = "Class")
  axis(1, at = seq(0, 1, length.out = ncol(mat)), labels = colnames(mat), cex.axis = 0.45)
  axis(2, at = seq(0, 1, length.out = nrow(mat)), labels = rev(rownames(mat)), las = 1, cex.axis = 0.8)
  box()
}

# Heatmaps: top-r selection frequencies.
methods <- unique(sel$method[sel$version == "topr"])
pdf(file.path(fig_dir, "N4_R_selection_frequency_heatmaps_topr.pdf"), width = 10, height = 8)
par(mfrow = c(2, ceiling(length(methods)/2)), mar = c(3, 3, 3, 1))
for (m in methods) plot_heatmap_one(subset(sel, method == m & version == "topr"), main = m)
dev.off()

# Accuracy-F1 class-wise scatter for top-r methods.
topr <- subset(cls, version == "topr")
pdf(file.path(fig_dir, "N4_R_accuracy_vs_featureF1_classwise.pdf"), width = 6.5, height = 5)
plot(NA, xlim = c(0,1), ylim = c(0,1), xlab = "Feature-recovery F1", ylab = "Class-wise prediction recall", main = "Class-wise joint performance")
for (m in unique(topr$method)) {
  dd <- subset(topr, method == m)
  col <- npkdc_method_cols[m]; if (is.na(col)) col <- "black"
  points(dd$f1, dd$class_recall, pch = 19, col = col)
  text(dd$f1, dd$class_recall, labels = dd$class, pos = 3, cex = 0.75, col = col)
}
legend("bottomleft", legend = unique(topr$method), col = npkdc_method_cols[unique(topr$method)], pch = 19, bty = "n", cex = 0.8)
dev.off()

message("Figures saved to ", fig_dir)
