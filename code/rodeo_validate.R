# rodeo_validate.R -----------------------------------------------------------
# Correctness gate for the Route-2 density-rodeo (rodeo_core.R).
#  (1) controlled low-d sanity: 1 class, 3 structured + 17 flat coords.
#  (2) Study-1 "uniform" generator at full scale (C=5, D=200, n_y=150),
#      parallel over replications on all cores. Target: attribution F1 ~ 1,
#      matching the original uniform sub-experiment -- but now via a genuine,
#      cross-fitted, non-cancelling rodeo instead of the sd detector.
# ---------------------------------------------------------------------------
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(here) || !nzchar(here)) here <- "."
source(file.path(here, "npkdcvd_core.R"))   # attr_metrics(detected, rel_sets)
source(file.path(here, "rodeo_core.R"))

## ---- (1) controlled low-d sanity ------------------------------------------
cat("=== (1) controlled sanity: d=20, relevant={1,2,3} ===\n")
set.seed(1)
n <- 150; d <- 20; rel <- 1:3
Xc <- matrix(runif(n * d), n, d)                       # flat U(0,1) background
mu <- c(0.3, 0.5, 0.7)
for (j in seq_along(rel)) Xc[, rel[j]] <- rnorm(n, mu[j], 0.10)  # concentrated
res <- rodeo_select_class(Xc, cols = seq_len(d))
cat("selected coords:", paste(which(res$selected), collapse = " "), "\n")
cat("top-6 rodeo scores:\n"); print(round(sort(res$score, decreasing = TRUE)[1:6], 2))
cat("\n")

## ---- Study-1 uniform generator (matches N1_simulation.R) ------------------
C <- 5; D <- 200; N_Y <- 150; N_TEST <- 100
R_Y <- c(3, 5, 5, 8, 8)
rel_start <- c(1, 4, 9, 14, 22)
REL_SETS  <- lapply(1:C, function(y) rel_start[y]:(rel_start[y] + R_Y[y] - 1))

gen_data <- function(n_y, irrel_type, rho = 0.5) {
  make_irrel <- function(nr) {
    if (irrel_type == "uniform") matrix(runif(nr * D), nr, D)
    else if (irrel_type == "normal") matrix(rnorm(nr * D), nr, D)
    else {                                              # ar1
      Z <- matrix(0, nr, D); Z[, 1] <- rnorm(nr)
      for (j in 2:D) Z[, j] <- rho * Z[, j - 1] + sqrt(1 - rho^2) * rnorm(nr)
      Z
    }
  }
  tl <- vector("list", C)
  for (y in 1:C) {
    rel <- REL_SETS[[y]]; nr <- R_Y[y]
    mu_rel <- rep(c(0.3, 0.5, 0.7), length.out = nr)
    X <- make_irrel(n_y)
    X[, rel] <- matrix(rnorm(n_y * nr, rep(mu_rel, each = n_y), 0.1), n_y, nr)
    tl[[y]] <- X
  }
  tl
}

one_rep <- function(b, irrel_type) {
  set.seed(20260824 + b)
  tl  <- gen_data(N_Y, irrel_type)
  det <- rodeo_detect(tl, keep = 40, mc.cores = 1L)    # reps are the parallel axis
  attr_metrics(det, REL_SETS)                          # C x 4: prec, rec, f1, exact
}

B    <- as.integer(Sys.getenv("RODEO_B", "60"))
ncor <- as.integer(Sys.getenv("RODEO_CORES", as.character(detectCores())))
types <- c("uniform", "normal", "ar1")
cat(sprintf("=== (2) Study-1 all sub-experiments: C=%d D=%d n_y=%d, B=%d/type on %d cores ===\n",
            C, D, N_Y, B, ncor))
for (irrel in types) {
  t0  <- proc.time()[["elapsed"]]
  out <- mclapply(seq_len(B), function(b) one_rep(b, irrel), mc.cores = ncor)
  el  <- proc.time()[["elapsed"]] - t0
  ok  <- vapply(out, function(o) is.matrix(o) && all(is.finite(o)), logical(1))
  out <- out[ok]
  arr <- array(unlist(out), dim = c(C, 4, length(out)))
  mn  <- apply(arr, c(1, 2), mean); colnames(mn) <- c("prec","rec","f1","exact")
  cat(sprintf("\n[%s] %d/%d reps in %.1f s (%.3f s/rep) | macro F1 = %.4f | exact = %.4f\n",
              irrel, length(out), B, el, el / B, mean(mn[,"f1"]), mean(mn[,"exact"])))
  print(round(mn, 4))
}
