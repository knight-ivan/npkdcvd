# rodeo_validate_bc.R --------------------------------------------------------
# Validation gate for Milestone 2 (baseline-corrected / PIT rodeo).
#  (A) HARD design: sharp SHARED nuisance {4,5,6} (must reject) + location shift
#      {7} (must select) + concentrated relevant {1,2,3}. Truth R_1 = {1,2,3,7}.
#      Compare plain vs baseline-corrected over many seeds.
#  (B) REGRESSION: Study-1 (uniform/normal/AR1) must still give F1 ~ 1 under the
#      baseline-corrected detector.
# ---------------------------------------------------------------------------
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(here) || !nzchar(here)) here <- "."
source(file.path(here, "npkdcvd_core.R"))
source(file.path(here, "rodeo_core.R"))

ncor <- as.integer(Sys.getenv("RODEO_CORES", as.character(detectCores())))

## ---- (A) hard design: plain vs baseline-corrected -------------------------
gen_hard <- function(seed) {
  set.seed(seed); C <- 3; d <- 20; n <- 150
  gen <- function(y) {
    X <- matrix(runif(n * d), n, d)
    X[, 4:6] <- matrix(rnorm(n * 3, 0, 0.1), n, 3)        # sharp SHARED nuisance
    X[, 7]   <- rnorm(n, if (y == 1) 1.5 else 0, 1)       # location shift class1
    if (y == 1) { mu <- c(0.3, 0.5, 0.7); for (j in 1:3) X[, j] <- rnorm(n, mu[j], 0.1) }
    X
  }
  lapply(1:C, gen)
}
truth1 <- c(1, 2, 3, 7)
hard_one <- function(seed) {
  tl <- gen_hard(seed)
  s_plain <- which(rodeo_detect(tl,    keep = 20, mc.cores = 1L)[1, ] == 1)
  s_bc    <- which(rodeo_detect_bc(tl, keep = 20, mc.cores = 1L)[1, ] == 1)
  c(plain = identical(as.numeric(s_plain), as.numeric(truth1)),
    bc    = identical(as.numeric(s_bc),    as.numeric(truth1)))
}
S <- 100L
cat(sprintf("=== (A) hard design (sharp shared + location shift), %d seeds ===\n", S))
hres <- simplify2array(mclapply(seq_len(S), hard_one, mc.cores = ncor))
cat(sprintf("  exact-recovery of {1,2,3,7}:  plain = %.2f   baseline-corrected = %.2f\n\n",
            mean(hres["plain", ]), mean(hres["bc", ])))

## ---- (B) Study-1 regression under baseline correction ---------------------
C <- 5; D <- 200; N_Y <- 150
R_Y <- c(3, 5, 5, 8, 8); rel_start <- c(1, 4, 9, 14, 22)
REL_SETS <- lapply(1:C, function(y) rel_start[y]:(rel_start[y] + R_Y[y] - 1))
gen_data <- function(n_y, irrel_type, rho = 0.5) {
  make_irrel <- function(nr) {
    if (irrel_type == "uniform") matrix(runif(nr * D), nr, D)
    else if (irrel_type == "normal") matrix(rnorm(nr * D), nr, D)
    else { Z <- matrix(0, nr, D); Z[, 1] <- rnorm(nr)
           for (j in 2:D) Z[, j] <- rho * Z[, j - 1] + sqrt(1 - rho^2) * rnorm(nr); Z }
  }
  tl <- vector("list", C)
  for (y in 1:C) {
    rel <- REL_SETS[[y]]; nr <- R_Y[y]; mu_rel <- rep(c(0.3, 0.5, 0.7), length.out = nr)
    X <- make_irrel(n_y)
    X[, rel] <- matrix(rnorm(n_y * nr, rep(mu_rel, each = n_y), 0.1), n_y, nr)
    tl[[y]] <- X
  }
  tl
}
B <- as.integer(Sys.getenv("RODEO_B", "60"))
cat(sprintf("=== (B) Study-1 regression (baseline-corrected), B=%d/type on %d cores ===\n", B, ncor))
for (irrel in c("uniform", "normal", "ar1")) {
  out <- mclapply(seq_len(B), function(b) {
    set.seed(20260824 + b)
    attr_metrics(rodeo_detect_bc(gen_data(N_Y, irrel), keep = 40, mc.cores = 1L), REL_SETS)
  }, mc.cores = ncor)
  arr <- array(unlist(out), dim = c(C, 4, B))
  mn  <- apply(arr, c(1, 2), mean)
  cat(sprintf("  [%s]  macro F1 = %.4f   exact-recovery = %.4f\n",
              irrel, mean(mn[, 3]), mean(mn[, 4])))
}
