# rodeo_study1.R -------------------------------------------------------------
# Study 1 (V1 design): overlapping class sets + the naming property, under the
# rebuilt rodeo. C=10, d=30, cyclic R_y={y,...,y+5} mod 30 (r_y=6, OVERLAPPING);
# relevant N(0.5,(0.02k)^2), irrelevant Uniform(0,1). Recover the inter-class
# Jaccard overlap matrix J_hat(R_hat_y,R_hat_y') and compare to true J.
# ---------------------------------------------------------------------------
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(here) || !nzchar(here)) here <- "."
source(file.path(here, "rodeo_core.R"))

C <- 10; D <- 30; N_Y <- 150
REL_SETS <- lapply(1:C, function(y) (((y - 1) + 0:5) %% 30) + 1)

true_jac <- matrix(0, C, C)
for (y in 1:C) for (yp in 1:C) {
  i <- length(intersect(REL_SETS[[y]], REL_SETS[[yp]]))
  u <- length(union(REL_SETS[[y]], REL_SETS[[yp]]))
  true_jac[y, yp] <- if (u) i / u else 1
}

gen_train <- function(n_y) lapply(1:C, function(y) {
  X <- matrix(runif(n_y * D), n_y, D)
  for (k in 1:6) X[, REL_SETS[[y]][k]] <- rnorm(n_y, 0.5, 0.02 * k)
  X
})
pairwise_jac <- function(det) {
  J <- matrix(0, C, C)
  for (y in 1:C) { Ry <- which(det[y, ] == 1L)
    for (yp in 1:C) { Ryp <- which(det[yp, ] == 1L)
      u <- length(union(Ry, Ryp)); J[y, yp] <- if (u) length(intersect(Ry, Ryp)) / u else 0 } }
  J
}

B    <- as.integer(Sys.getenv("RODEO_B", "300"))
ncor <- as.integer(Sys.getenv("RODEO_CORES", as.character(detectCores())))
cat(sprintf("=== Study 1: overlapping sets / Jaccard recovery (rodeo), B=%d on %d cores ===\n", B, ncor))
out <- mclapply(seq_len(B), function(b) {
  set.seed(b * 777 + 13)
  det <- rodeo_detect(gen_train(N_Y), keep = D, mc.cores = 1L)
  exact <- all(vapply(1:C, function(y)
    setequal(which(det[y, ] == 1L), REL_SETS[[y]]), logical(1)))
  list(J = pairwise_jac(det), exact = exact)
}, mc.cores = ncor)

est_jac <- Reduce("+", lapply(out, `[[`, "J")) / length(out)
exact_rate <- mean(vapply(out, `[[`, logical(1), "exact"))
maxdiff <- max(abs(est_jac - true_jac))

cat(sprintf("\nall-classes exact recovery rate = %.4f\n", exact_rate))
cat(sprintf("max |J_hat - J| = %.3e\n", maxdiff))
cat("true J (row 1):     ", paste(sprintf("%.3f", true_jac[1, ]), collapse = " "), "\n")
cat("estimated J (row 1):", paste(sprintf("%.3f", est_jac[1, ]), collapse = " "), "\n")
saveRDS(list(true_jac = true_jac, est_jac = est_jac, exact_rate = exact_rate),
        file.path(here, "..", "data", "rodeo_V1_study1_results.rds"))
