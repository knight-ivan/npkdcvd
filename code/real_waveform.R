# real_waveform.R — Breiman waveform-40 validation: rebuilt-rodeo attribution.
# Synthetic (reproducible): 3 classes, vars 1-21 = signal (R_y={1..21}), 22-40 = noise.
#   RODEO_CORES=22 RODEO_B=1000 Rscript real_waveform.R
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (length(here) && nzchar(here)) setwd(here)
source("rodeo_core.R"); suppressMessages(library(parallel))
cores <- max(1L, as.integer(Sys.getenv("RODEO_CORES", as.character(detectCores() - 2L))))
B     <- as.integer(Sys.getenv("RODEO_B", "1000"))
n_y   <- 500L

h  <- function(i) pmax(6 - abs(i - 11), 0)                 # Breiman triangular waveform
H  <- cbind(h(1:21), h((1:21) - 4), h((1:21) + 4))        # 21 x 3 base waveforms
pr <- list(c(1, 2), c(1, 3), c(2, 3))
gen_class <- function(n, cl) {
  u <- runif(n); a <- pr[[cl]][1]; b <- pr[[cl]][2]
  base <- outer(u, H[, a]) + outer(1 - u, H[, b])          # n x 21
  cbind(base + matrix(rnorm(n * 21), n, 21), matrix(rnorm(n * 19), n, 19))   # -> n x 40
}
rel <- 1:21
am  <- function(sel, D = 40) {                             # attribution metrics vs rel set
  tp <- sum(sel[rel]); fp <- sum(sel[-rel]); prec <- if (tp + fp > 0) tp / (tp + fp) else 1
  rec <- tp / length(rel); f1 <- if (prec + rec > 0) 2 * prec * rec / (prec + rec) else 0
  truth <- seq_len(D) %in% rel
  c(prec = prec, rec = rec, f1 = f1, exact = as.numeric(all(truth == (sel == 1))),
    n_noise_sel = fp)
}
one <- function(b) {
  set.seed(b * 104729 + 3)
  sel <- rodeo_detect(lapply(1:3, function(cl) gen_class(n_y, cl)), keep = 40, mc.cores = 1L)  # 3 x 40
  t(apply(sel, 1, am))
}
cat(sprintf("Waveform-40: 3 classes, n_y=%d, B=%d, cores=%d\n", n_y, B, cores))
res <- mclapply(seq_len(B), one, mc.cores = cores)
M   <- do.call(rbind, res)
cat(sprintf("Attribution (mean over classes x reps): Prec=%.4f Rec=%.4f F1=%.4f P(exact)=%.4f  mean #noise selected=%.3f\n",
            mean(M[, "prec"]), mean(M[, "rec"]), mean(M[, "f1"]), mean(M[, "exact"]), mean(M[, "n_noise_sel"])))
saveRDS(list(M = M), file.path("..", "data", "real_waveform_results.rds"))
cat("saved real_waveform_results.rds\n")
