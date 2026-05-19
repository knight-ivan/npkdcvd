# 07_make_phase3_visualizations.R
# R counterparts for Phase-3 visualization figures.
# Usage:
#   Rscript R/07_make_phase3_visualizations.R --freq results_from_previous_runs/N4/n4_selection_frequency_topr.csv --out phase3_R_figures
#
# V1 is generated from class-by-feature selection frequencies.
# V2 can be generated from a user-supplied 28x28 digit class-mean CSV with columns:
#   class,pixel,row,col,value
# If no digit CSV is supplied, a clearly labeled synthetic MNIST-resolution template is produced.

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0L && !is.na(a)) a else b
cmd_all <- commandArgs(FALSE)
file_arg <- cmd_all[grep("^--file=", cmd_all)]
this_file <- sub("^--file=", "", file_arg[1] %||% "R/07_make_phase3_visualizations.R")
script_dir <- dirname(normalizePath(this_file, mustWork = FALSE))
source(file.path(script_dir, "00_setup.R"))

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(freq = file.path(dirname(script_dir), "results_from_previous_runs", "N4", "n4_selection_frequency_topr.csv"),
              out = file.path(dirname(script_dir), "phase3_R_figures"), digit_csv = NA_character_)
  if (length(args) == 0L) return(out)
  for (i in seq(1L, length(args), by = 2L)) {
    key <- sub("^--", "", args[i]); val <- args[i + 1L]
    if (key == "freq") out$freq <- val
    if (key == "out") out$out <- val
    if (key == "digit_csv") out$digit_csv <- val
  }
  out
}

args <- parse_args()
if (!dir.exists(args$out)) dir.create(args$out, recursive = TRUE, showWarnings = FALSE)

# ---------- V1: Jaccard overlap matrix ----------
if (file.exists(args$freq)) {
  freq <- utils::read.csv(args$freq)
  if (!"freq" %in% names(freq)) names(freq)[names(freq) == "selection_frequency"] <- "freq"
  methods <- unique(freq$method)
  pdf(file.path(args$out, "V1_R_jaccard_overlap.pdf"), width = 8, height = 8)
  par(mfrow = c(2, ceiling(length(methods)/2)), mar = c(3, 3, 3, 4))
  for (m in methods) {
    dd <- subset(freq, method == m)
    mat <- xtabs(freq ~ class + feature, data = dd)
    # Stable support: selected if frequency >= 0.5.
    support <- mat >= 0.5
    C <- nrow(support)
    jac <- matrix(0, C, C)
    for (i in seq_len(C)) for (j in seq_len(C)) {
      inter <- sum(support[i, ] & support[j, ])
      union <- sum(support[i, ] | support[j, ])
      jac[i, j] <- if (union == 0) 1 else inter / union
    }
    image(t(jac[C:1, ]), axes = FALSE, zlim = c(0,1),
          col = colorRampPalette(c("white", "#bdd7e7", "#3182bd"))(101),
          main = m, xlab = "Class", ylab = "Class")
    axis(1, at = seq(0, 1, length.out = C), labels = seq_len(C))
    axis(2, at = seq(0, 1, length.out = C), labels = rev(seq_len(C)), las = 1)
    for (i in seq_len(C)) for (j in seq_len(C)) {
      val <- jac[i,j]
      text((j-1)/(C-1), (C-i)/(C-1), sprintf("%.2f", val), col = ifelse(val > 0.55, "white", "black"), cex = 0.85)
    }
    box()
  }
  dev.off()
  message("Saved V1 R Jaccard figure to ", args$out)
} else {
  warning("Frequency CSV not found: ", args$freq)
}

# ---------- V2: MNIST-resolution digit heatmap ----------
make_template_digit <- function(class_id) {
  z <- matrix(0, 28, 28)
  rr <- 5:24; cc <- 6:23
  if (class_id %% 5 == 0) { z[5:7, cc] <- 1; z[22:24, cc] <- 1; z[rr, 6:8] <- 0.8; z[rr, 21:23] <- 0.8 }
  if (class_id %% 5 == 1) { z[4:24, 13:16] <- 1; z[22:24, 9:20] <- 0.8 }
  if (class_id %% 5 == 2) { z[5:7, cc] <- 1; z[14:16, cc] <- 0.8; z[22:24, cc] <- 1; z[8:14, 20:23] <- .8; z[16:22, 6:9] <- .8 }
  if (class_id %% 5 == 3) { z[5:7, cc] <- 1; z[14:16, cc] <- .8; z[22:24, cc] <- 1; z[rr, 20:23] <- .8 }
  if (class_id %% 5 == 4) { z[4:15, 6:8] <- .8; z[4:24, 20:23] <- 1; z[14:16, cc] <- .8 }
  z + matrix(runif(28*28, 0, 0.03), 28, 28)
}

if (!is.na(args$digit_csv) && file.exists(args$digit_csv)) {
  pix <- utils::read.csv(args$digit_csv)
  classes <- sort(unique(pix$class))
  pdf(file.path(args$out, "V2_R_digit_28x28_heatmap.pdf"), width = 10, height = 5)
  par(mfrow = c(2, ceiling(length(classes)/2)), mar = c(1,1,2,1))
  for (cl in classes) {
    dd <- subset(pix, class == cl)
    z <- matrix(dd$value[order(dd$row, dd$col)], 28, 28, byrow = TRUE)
    image(t(z[28:1, ]), axes = FALSE, col = gray.colors(101, start = 1, end = 0), main = paste("Class", cl))
  }
  dev.off()
} else {
  pdf(file.path(args$out, "V2_R_digit_28x28_template_heatmap.pdf"), width = 10, height = 5)
  par(mfrow = c(2,5), mar = c(1,1,2,1))
  for (cl in 0:9) {
    z <- make_template_digit(cl)
    image(t(z[28:1, ]), axes = FALSE, col = gray.colors(101, start = 1, end = 0), main = paste("Template", cl))
  }
  dev.off()
  message("No digit CSV supplied; saved clearly labeled synthetic 28x28 template heatmap.")
}
