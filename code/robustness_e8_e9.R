# robustness_e8_e9.R ----------------------------------------------------------
# E8 (baseline-correction ablation on E7) and E9 (no-rodeo marginal-bandwidth
# baseline on E1-E7), written as saved, reproducible R outputs.
#
#   E9: the Silverman-bandwidth first-gap detector of npkdcvd_core.R::detect_gap
#       ("gap" rule), run on the same E1-E7 generators, seeds and n_y grids as the
#       proposed rodeo (robustness_experiments.R).  It answers the same marginal
#       concentration question without the iterative density-derivative path.
#   E8: E7 (shared non-flat nuisance) with baseline correction OFF (plain rodeo,
#       "rodeo" rule) versus ON (pooled empirical PIT, "rodeo_bc" rule).
#
# The proposed-rodeo E1-E7 outputs must already exist (run_R2.sh runs
# robustness_experiments.R E1..E7 --rule=rodeo first); they are read, not rerun.
# Every run uses the scenario registry's seed rule b*1000 + n0 + 100000*index.
#
#   RODEO_CORES=14 RODEO_B=1000 Rscript robustness_e8_e9.R
# ---------------------------------------------------------------------------
.this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
.here <- if (is.na(.this_file)) "." else dirname(normalizePath(.this_file))
# Sourcing is side-effect free: the CLI block of robustness_experiments.R is
# guarded by sys.nframe() == 0.
source(file.path(.here, "robustness_experiments.R"))

B     <- reviewer_B(1000L)
cores <- reviewer_cores()
out   <- reviewer_output_dir()
scen  <- paste0("E", 1:7)

# Summary row at the LARGEST n_y of a saved run (estimate and replication-level MC SE).
largest <- function(name, rule) {
  f <- file.path(out, sprintf("robustness_%s_%s_summary.csv", name, rule))
  if (!file.exists(f)) stop("missing ", basename(f), "; run robustness_experiments.R first")
  s <- read.csv(f)
  s <- s[s$n_y == max(s$n_y), , drop = FALSE]
  g <- function(m, col) s[[col]][s$metric == m]
  data.frame(setting = name, rule = rule, n_y = max(s$n_y), B = g("f1", "B"),
             f1 = g("f1", "estimate"), f1_mc_se = g("f1", "mc_se"),
             exact = g("exact", "estimate"), exact_mc_se = g("exact", "mc_se"))
}

## E9: no-rodeo marginal-bandwidth baseline, all n_y of every design ----------
for (nm in scen) invisible(run_scenario(nm, B = B, rule = "gap", cores = cores))
e9 <- do.call(rbind, lapply(scen, function(nm)
  rbind(largest(nm, "rodeo"), largest(nm, "gap"))))
write.csv(e9, file.path(out, "E9_rodeo_vs_gap_largest_n.csv"), row.names = FALSE)

## E8: baseline correction OFF vs ON on E7, both n_y --------------------------
invisible(run_scenario("E7", B = B, rule = "rodeo_bc", cores = cores))
e8 <- do.call(rbind, lapply(c("rodeo", "rodeo_bc"), function(rl) {
  s <- read.csv(file.path(out, sprintf("robustness_E7_%s_summary.csv", rl)))
  s <- s[s$metric %in% c("f1", "exact"), c("n_y", "metric", "B", "estimate", "mc_se")]
  cbind(baseline_correction = if (rl == "rodeo") "OFF" else "ON (pooled PIT)", s)
}))
write.csv(e8, file.path(out, "E8_E7_baseline_ablation.csv"), row.names = FALSE)

cat(sprintf("\n=== E9: proposed rodeo vs no-rodeo baseline at the largest n_y (B=%d) ===\n", B))
print(e9, row.names = FALSE, digits = 4)
cat(sprintf("\n=== E8: E7 baseline-correction ablation (B=%d) ===\n", B))
print(e8, row.names = FALSE, digits = 4)
