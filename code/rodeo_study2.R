# Formal-selector analysis for Study 2. The original N1_simulation.R is retained
# as a LEGACY bandwidth Z-score/classification diagnostic and is not replaced.
# DGP is taken from N1_simulation.R: C=5,d=200,n_y=150,r=(3,5,5,8,8);
# relevant ~ N(mu_j,0.1^2), mu_j cycling (0.3,0.5,0.7).
# Noise: independent Uniform(0,1), N(0,1), or AR(1) rho=0.5.
# Main analysis uses the final adaptive rodeo with documented top-40 screen,
# no oracle set size; both screened and conditional recovery are saved.
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])))
source(file.path(script_dir, "npkdcvd_core.R"))
source(file.path(script_dir, "rodeo_core.R"))
source(file.path(script_dir, "reviewer_metrics.R"))
C <- 5L; D <- 200L; N_Y <- 150L
r_y <- c(3L, 5L, 5L, 8L, 8L)
starts <- c(1L, 4L, 9L, 14L, 22L)
truth <- lapply(seq_len(C), function(y) seq.int(starts[y], length.out = r_y[y]))
KEEP <- as.integer(Sys.getenv("RODEO_KEEP", "40"))
if (is.na(KEEP) || KEEP < 2L || KEEP > D) stop("invalid RODEO_KEEP")
generate <- function(noise) {
  lapply(seq_len(C), function(y) {
    X <- matrix(0, N_Y, D)
    rr <- truth[[y]]; m <- D - length(rr)
    mm <- rep(c(.3, .5, .7), length.out = length(rr))
    X[, rr] <- matrix(rnorm(N_Y*length(rr), mean=rep(mm, each=N_Y), sd=.1), N_Y, length(rr))
    if (noise == "uniform") Z <- matrix(runif(N_Y*m), N_Y, m)
    else if (noise == "normal") Z <- matrix(rnorm(N_Y*m), N_Y, m)
    else if (noise == "ar1") {
      Z <- matrix(0, N_Y, m); Z[, 1L] <- rnorm(N_Y)
      for (j in 2:m) Z[, j] <- .5*Z[, j-1L] + sqrt(1-.5^2)*rnorm(N_Y)
    } else stop("unknown noise")
    X[, -rr] <- Z
    X
  })
}
B <- reviewer_B()
noise_types <- c("uniform", "normal", "ar1")
allrows <- list()
for (noise in noise_types) {
  out <- reviewer_mc(seq_len(B), function(b) {
    set.seed(20261022L + 100003L*b + match(noise, noise_types)*1009L)
    X <- generate(noise)
    z <- rodeo_detect_detail(X, keep=KEEP, mc.cores=1L)
    rows <- reviewer_metric_rows(z$detected, truth, b, noise, N_Y)
    rows$screen_all <- vapply(seq_len(C), function(y)
      as.numeric(all(truth[[y]] %in% z$screen[[y]])), numeric(1))
    rows$conditional_exact <- ifelse(rows$screen_all == 1, rows$exact, NA_real_)
    rows$sparsity <- r_y
    rows
  })
  allrows[[noise]] <- do.call(rbind, out)
  cat(sprintf("Study 2 (%s): F1 %.3f, exact %.3f, screen survival %.3f\n",
      noise, mean(allrows[[noise]]$f1), mean(allrows[[noise]]$exact),
      mean(allrows[[noise]]$screen_all)))
}
raw <- do.call(rbind, allrows); rownames(raw) <- NULL
resdir <- reviewer_output_dir()
write.csv(raw, file.path(resdir, "study2_rodeo_replications.csv"), row.names=FALSE)
sm <- reviewer_summary(raw)
write.csv(sm, file.path(resdir, "study2_rodeo_summary.csv"), row.names=FALSE)
# Stratified metrics by r_y: group by noise x sparsity and independent replication.
for (rr in unique(r_y)) {
  sm_rr <- reviewer_summary(raw[raw$sparsity == rr,,drop=FALSE])
  sm_rr$sparsity <- rr
  if (!exists("strata")) strata <- list()
  strata[[as.character(rr)]] <- sm_rr
}
write.csv(do.call(rbind, strata), file.path(resdir,"study2_rodeo_sparsity.csv"),row.names=FALSE)
saveRDS(list(data=raw, summary=sm, B=B, keep=KEEP, truth=truth, noise=noise_types),
        file.path(resdir,"study2_rodeo_results.rds"))
