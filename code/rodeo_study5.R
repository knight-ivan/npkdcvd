# rodeo_study5.R -------------------------------------------------------------
# Milestone 3, flagship: Study 5 (generative vs discriminative) under the
# rebuilt rodeo, ADAPTIVE gap rule (no oracle top-r).
# Design (from N6_simulation.R): C=3, D=15, n_y=200. Truth R_1={1,2}, R_2={3},
# R_3={4,5}. X3 is N(0,1) within class 1 (diffuse; discriminative-not-generative)
# but N(3,0.2^2) in class 2. The paper's claim: NPKDC-vd EXCLUDES X3 from R_1.
#
# We compare: plain rodeo (uniform baseline = generative question) vs the
# baseline-corrected PIT rodeo (pooled background = deviation-from-pool). The key
# diagnostic is P(X3 selected for class 1): it MUST stay low to keep the story.
# ---------------------------------------------------------------------------
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(here) || !nzchar(here)) here <- "."
source(file.path(here, "npkdcvd_core.R"))
source(file.path(here, "rodeo_core.R"))

C <- 3; D <- 15; N_Y <- 200
REL_SETS <- list(c(1L, 2L), 3L, c(4L, 5L))

gen_data <- function(n_y) {
  gen_class <- function(y, n) {
    X <- matrix(rnorm(n * D), n, D)                 # all N(0,1) by default
    if (y == 1L) { X[, 1] <- rnorm(n, 0.5, 0.2); X[, 2] <- rnorm(n, 0.5, 0.2) }  # X3 stays N(0,1)
    else if (y == 2L) X[, 3] <- rnorm(n, 3.0, 0.2)
    else if (y == 3L) { X[, 4] <- rnorm(n, 0.5, 0.2); X[, 5] <- rnorm(n, 0.5, 0.2) }
    X
  }
  lapply(1:C, function(y) gen_class(y, n_y))
}

one_rep <- function(b) {
  set.seed(b * 1313 + 77)
  tl <- gen_data(N_Y)
  det_plain <- rodeo_detect(tl,    keep = D, mc.cores = 1L)
  det_bc    <- rodeo_detect_bc(tl, keep = D, mc.cores = 1L)
  list(plain = attr_metrics(det_plain, REL_SETS), plain_x3 = det_plain[1, 3],
       bc    = attr_metrics(det_bc,    REL_SETS), bc_x3    = det_bc[1, 3])
}

B    <- as.integer(Sys.getenv("RODEO_B", "200"))
ncor <- as.integer(Sys.getenv("RODEO_CORES", as.character(detectCores())))
cat(sprintf("=== Study 5 under the rodeo (adaptive rule), B=%d on %d cores ===\n", B, ncor))
out <- mclapply(seq_len(B), one_rep, mc.cores = ncor)

summ <- function(key, x3key) {
  arr <- array(unlist(lapply(out, `[[`, key)), dim = c(C, 4, length(out)))
  f1  <- apply(arr[, 3, ], 1, mean); ex <- apply(arr[, 4, ], 1, mean)
  x3  <- mean(vapply(out, `[[`, numeric(1), x3key))
  list(f1 = f1, ex = ex, x3 = x3)
}
for (m in c("plain", "bc")) {
  s <- summ(m, paste0(m, "_x3"))
  cat(sprintf("\n[%s]  F1: C1=%.3f C2=%.3f C3=%.3f | exact: %.3f %.3f %.3f | P(X3->C1)=%.3f\n",
              ifelse(m == "plain", "plain rodeo (generative)", "baseline-corrected (pool)"),
              s$f1[1], s$f1[2], s$f1[3], s$ex[1], s$ex[2], s$ex[3], s$x3))
}
cat("\n(narrative needs: C1/C2/C3 F1 ~ 1 AND P(X3->C1) ~ 0)\n")
