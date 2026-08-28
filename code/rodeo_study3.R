# rodeo_study3.R -------------------------------------------------------------
# Study 3 (N2 design): selection-consistency convergence (Theorem 2) under the
# rebuilt rodeo. C=5, d=25, r_y=3 disjoint; relevant N(0.5,0.5^2) MODERATE signal,
# irrelevant N(0,1). Show recall / F1 / P(exact) -> 1 as n_y grows.
# ---------------------------------------------------------------------------
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(here) || !nzchar(here)) here <- "."
source(file.path(here, "npkdcvd_core.R"))
source(file.path(here, "rodeo_core.R"))

C <- 5; D <- 25
R_Y <- rep(3L, C); rel_start <- cumsum(c(1L, R_Y[-C]))
REL_SETS <- lapply(1:C, function(y) rel_start[y]:(rel_start[y] + R_Y[y] - 1))
N_Y_GRID <- c(30, 50, 75, 100, 150, 200, 300, 500)

gen_train <- function(n_y) lapply(1:C, function(y) {
  X <- matrix(rnorm(n_y * D), n_y, D)
  X[, REL_SETS[[y]]] <- matrix(rnorm(n_y * R_Y[y], 0.5, 0.5), n_y, R_Y[y])
  X
})

B    <- as.integer(Sys.getenv("RODEO_B", "300"))
ncor <- as.integer(Sys.getenv("RODEO_CORES", as.character(detectCores())))
cat(sprintf("=== Study 3: selection convergence (rodeo, adaptive), B=%d/n_y on %d cores ===\n", B, ncor))
cat(sprintf("%5s | %6s %6s %6s %8s\n", "n_y", "recall", "prec", "F1", "P(exact)"))

res <- data.frame()
for (n_y in N_Y_GRID) {
  out <- mclapply(seq_len(B), function(b) {
    set.seed(b * 1000 + n_y)
    attr_metrics(rodeo_detect(gen_train(n_y), keep = D, mc.cores = 1L), REL_SETS)
  }, mc.cores = ncor)
  arr <- array(unlist(out), dim = c(C, 4, B))          # prec,rec,f1,exact
  m   <- apply(arr, 2, mean)                            # avg over classes & reps
  cat(sprintf("%5d | %6.3f %6.3f %6.3f %8.3f\n", n_y, m[2], m[1], m[3], m[4]))
  res <- rbind(res, data.frame(n_y = n_y, recall = m[2], prec = m[1], f1 = m[3], exact = m[4]))
}
saveRDS(res, file.path(here, "..", "data", "rodeo_N2_study3_results.rds"))
cat("\nsaved rodeo_N2_study3_results.rds\n")
