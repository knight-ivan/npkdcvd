## reproduce_all_tables.R
## Reproduces all tables in the main paper and Supplementary Material.
##
## Usage (from the repository root):
##   Rscript reproduce_all_tables.R
##
## Each study writes its tabular output to results/<study>/.
## Full runs (B = 500 or 1000 replications) can take several hours.
## Set B_OVERRIDE to a small number (e.g. 5) for a quick smoke test.

B_OVERRIDE <- NULL   # set e.g. B_OVERRIDE <- 5L for a quick test

root <- normalizePath(dirname(sys.frame(1)$ofile %||% "."), mustWork = FALSE)
`%||%` <- function(a, b) if (!is.null(a)) a else b

cat("== NPKDC-vd: Reproducing all tables ==\n\n")

## -- Study N1: Synthetic Study 1 (overlapping class-specific variable sets) ---
cat("Study N1 ...\n")
source(file.path(root, "scripts", "N1_simulation.R"))

## -- Study N2: Synthetic Study 3 (detection quality vs sample size) -----------
cat("Study N2 (detection probability) ...\n")
source(file.path(root, "scripts", "N2_simulation.R"))

## -- Study N3: Synthetic Study 2 (high-dimensional sparse setting) ------------
cat("Study N3 ...\n")
source(file.path(root, "scripts", "N3_simulation.R"))

## -- Study N4: Synthetic Study 4 (modern comparison, B = 1000 reps) ----------
cat("Study N4 (modern comparison, B=1000) ...\n")
B <- B_OVERRIDE %||% 1000L
cmd <- sprintf('Rscript "%s" --B %d --out "%s"',
               file.path(root, "R", "03_run_N4_modern_comparison_B1000.R"),
               B, file.path(root, "results", "N4"))
system(cmd)

## -- Study N6: Synthetic Study 5 (generative vs discriminative attribution) ---
cat("Study N6 ...\n")
source(file.path(root, "scripts", "N6_simulation.R"))

## -- Real Data V1: Anuran species classification (Jaccard / F1 table) --------
cat("Real Data V1 (Anuran) ...\n")
source(file.path(root, "scripts", "V1_jaccard.R"))

## -- Real Data V2: Handwritten digit attribution ------------------------------
cat("Real Data V2 (Digits) ...\n")
source(file.path(root, "scripts", "V2_digit_simulation.R"))

cat("\n== All tables done. Results written to results/. ==\n")
