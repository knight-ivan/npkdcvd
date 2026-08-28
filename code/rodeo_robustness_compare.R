# rodeo_robustness_compare.R -------------------------------------------------
# E8 (ablation) + E9 (density competitor) for the robustness grid.
#   E9 density competitor: the marginal-sd concentration detector (Route-1,
#       npkdcvd_core.R detect_gap) is a genuine density-based variable selector;
#       run it head-to-head against the rodeo on E1-E7.
#   E8 ablation on E7 (shared non-flat nuisance): baseline correction OFF
#       (plain rodeo) vs ON (rodeo_bc) vs KDE-no-rodeo (sd detector).
# Sources robustness_experiments.R for the E1-E7 generators + runner.
# ---------------------------------------------------------------------------
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(here) || !nzchar(here)) here <- "."
invisible(capture.output(source(file.path(here, "robustness_experiments.R"))))

B    <- as.integer(Sys.getenv("RODEO_B", "150"))
ncor <- as.integer(Sys.getenv("RODEO_CORES", as.character(detectCores())))
scen <- c("E1", "E2", "E3", "E4", "E5", "E6", "E7")

grab <- function(nm, rule) {
  invisible(capture.output(r <- run_scenario(nm, B = B, rule = rule, cores = ncor)))
  last <- r[[length(r)]]                          # metrics at the largest n_y
  c(f1 = unname(last["f1"]), ex = unname(last["exact"]))
}

cat(sprintf("=== E9: rodeo vs density competitor (sd detector), B=%d, largest n_y ===\n", B))
cat(sprintf("%-4s | %-18s | %-18s\n", "scen", "rodeo", "sd-detector (E9)"))
cat(strrep("-", 46), "\n")
for (nm in scen) {
  a <- grab(nm, "rodeo"); b <- grab(nm, "gap")
  cat(sprintf("%-4s | F1=%.3f ex=%.3f | F1=%.3f ex=%.3f\n",
              nm, a["f1"], a["ex"], b["f1"], b["ex"]))
}

cat("\n=== E8 ablation on E7 (shared non-flat nuisance) ===\n")
for (rl in c("rodeo", "rodeo_bc", "gap")) {
  m <- grab("E7", rl)
  lbl <- switch(rl, rodeo = "rodeo (baseline OFF)", rodeo_bc = "rodeo + baseline (ON)",
                gap = "KDE-no-rodeo (sd detector)")
  cat(sprintf("  %-26s F1=%.3f  exact=%.3f\n", lbl, m["f1"], m["ex"]))
}
