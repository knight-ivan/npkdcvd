## reproduce_all_figures.R
## Reproduces all figures in the main paper and Supplementary Material.
##
## Usage (from the repository root):
##   Rscript reproduce_all_figures.R
##
## Figures are written to figures/ as EPS and PDF.
## Run reproduce_all_tables.R first to generate the results/ .rds files,
## or set USE_PRECOMPUTED <- TRUE to use the committed pre-computed results.

USE_PRECOMPUTED <- TRUE   # set FALSE to recompute from raw simulation output

root <- normalizePath(dirname(sys.frame(1)$ofile %||% "."), mustWork = FALSE)
`%||%` <- function(a, b) if (!is.null(a)) a else b

cat("== NPKDC-vd: Reproducing all figures ==\n\n")

## -- Figure 1: N1 bandwidth Z-score heatmap -----------------------------------
cat("Figure N1_zscore (Synthetic Study 2 Z-score) ...\n")
source(file.path(root, "scripts", "N1_figure.R"))

## -- Figure 2: N2 detection probability vs sample size -----------------------
cat("Figure N2_detprob (Synthetic Study 3 detection probability) ...\n")
source(file.path(root, "scripts", "N2_figure.R"))

## -- Figure 3: N3 contrast plot -----------------------------------------------
cat("Figure N3_contrast (Synthetic Study 2 contrast) ...\n")
source(file.path(root, "scripts", "N3_figure.R"))

## -- Figure 4: N6 selection frequency (generative vs discriminative) ----------
cat("Figure N6_selfreq (Synthetic Study 5 selection frequency) ...\n")
source(file.path(root, "scripts", "N6_figure_v2.R"))

## -- Figure N4: Modern comparison figures (Study 4) ---------------------------
cat("Figure N4 (modern comparison) ...\n")
cmd <- sprintf('Rscript "%s" --in "%s"',
               file.path(root, "R", "04_make_N4_figures.R"),
               file.path(root, "results", "N4"))
system(cmd)

## -- Supplementary figures: bandwidth box-plots (frogs, digits, wave) --------
cat("Supplementary visualizations (frogs, digits, waveform) ...\n")
source(file.path(root, "R", "07_make_phase3_visualizations.R"))

cat("\n== All figures done. Written to figures/. ==\n")
