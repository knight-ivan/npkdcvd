# real_waveform.R — Waveform-40 generated benchmark: corrected cross-fitted rodeo
# attribution. Controlled generator: 3 classes, vars 1-21 are discriminative
# waveform coordinates by construction; vars 22-40 are independent N(0,1) noise
# shared by all classes. The manuscript does not treat {1,...,21} as a generative
# naming truth. Each replication runs the selector on the raw scale (first, so its
# random stream equals a raw-only run) and after the pooled empirical PIT.
#   RODEO_CORES=14 RODEO_B=1000 Rscript real_waveform.R
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (length(here) && nzchar(here)) setwd(here)
source("rodeo_core.R"); suppressMessages(library(parallel))
cores <- max(1L, as.integer(Sys.getenv("RODEO_CORES", as.character(max(1L, detectCores() - 2L)))))
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
  trl <- lapply(1:3, function(cl) gen_class(n_y, cl))
  sel    <- rodeo_detect(trl, keep = 40, mc.cores = 1L)      # raw scale, 3 x 40
  sel_bc <- rodeo_detect_bc(trl, keep = 40, mc.cores = 1L)   # pooled empirical PIT
  list(raw = t(apply(sel, 1, am)), pit = t(apply(sel_bc, 1, am)),
       sel = sel, sel_bc = sel_bc)
}
cat(sprintf("Waveform-40: 3 classes, n_y=%d, B=%d, cores=%d\n", n_y, B, cores))
res <- mclapply(seq_len(B), one, mc.cores = cores)
M    <- do.call(rbind, lapply(res, `[[`, "raw"))
M_bc <- do.call(rbind, lapply(res, `[[`, "pit"))
freq    <- Reduce(`+`, lapply(res, `[[`, "sel")) / length(res)      # 3 x 40 selection frequency
freq_bc <- Reduce(`+`, lapply(res, `[[`, "sel_bc")) / length(res)
rep_summary <- function(M, lab) {
  relsel <- M[, "rec"] * length(rel); noisel <- M[, "n_noise_sel"]
  cat(sprintf("%-22s Prec=%.4f Rec=%.4f F1=%.4f P(exact)=%.4f | per class: #waveform sel=%.3f #noise sel=%.3f; P(>=1 noise)=%.3f; P(empty)=%.3f\n",
      lab, mean(M[, "prec"]), mean(M[, "rec"]), mean(M[, "f1"]), mean(M[, "exact"]),
      mean(relsel), mean(noisel), mean(noisel > 0), mean(relsel + noisel == 0)))
}
rep_summary(M, "raw scale:")
rep_summary(M_bc, "pooled empirical PIT:")
outdir <- Sys.getenv("RODEO_OUT", file.path("..", "results_R2")); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
saveRDS(list(M = M, M_bc = M_bc, freq = freq, freq_bc = freq_bc, B = B, n_y = n_y),
        file.path(outdir, "real_waveform_results.rds"))
cat(sprintf("\nsaved %s\n", file.path(outdir, "real_waveform_results.rds")))
