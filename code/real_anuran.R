# real_anuran.R — Anuran Calls (MFCCs): NPKDC-vd classification + corrected
# cross-fitted rodeo attribution. Data: UCI Anuran Calls (MFCCs), Frogs_MFCCs.csv.
# Each replication reports the selector on the raw MFCC scale and, as main-paper Algorithm 3
# prescribes when the shared background is non-flat, after the pooled empirical
# PIT (rodeo_detect_bc).  The raw-scale call is made first, so its random stream
# is identical to a raw-only run.  Classification does not use the rodeo.
#   ANURAN_CSV=/path/Frogs_MFCCs.csv  RODEO_CORES=14  RODEO_B=1000  Rscript real_anuran.R
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (length(here) && nzchar(here)) setwd(here)
source("rodeo_core.R")
suppressMessages({library(parallel); library(e1071)})

CSV   <- Sys.getenv("ANURAN_CSV", "/tmp/anuran/Frogs_MFCCs.csv")
cores <- max(1L, as.integer(Sys.getenv("RODEO_CORES", as.character(max(1L, detectCores() - 2L)))))
B     <- as.integer(Sys.getenv("RODEO_B", "1000"))
N_TR  <- 100L; N_TE <- 50L

d    <- read.csv(CSV, check.names = FALSE)
mfcc <- as.matrix(d[, grep("^MFCCs_", names(d))]); colnames(mfcc) <- paste0("M", seq_len(ncol(mfcc)))
sp   <- d$Species
keep_sp   <- sort(names(which(table(sp) >= (N_TR + N_TE))))   # 7 species with enough samples
idx_by_sp <- lapply(keep_sp, function(s) which(sp == s))
C <- length(keep_sp)
cat(sprintf("Anuran: %d species (>=%d samples), d=%d MFCC, B=%d, cores=%d\n",
            C, N_TR + N_TE, ncol(mfcc), B, cores))

bw_silverman <- function(x) { s <- sd(x); if (s < 1e-12) return(1e-4); 1.06 * s * length(x)^(-0.2) }
log_kde_all <- function(Xte, Xtr, bw) {
  inv <- 1 / bw; A <- sweep(Xte, 2, inv, `*`); Bm <- sweep(Xtr, 2, inv, `*`)
  D2 <- pmax(outer(rowSums(A^2), rowSums(Bm^2), `+`) - 2 * tcrossprod(A, Bm), 0)
  Dd <- length(bw); cst <- -0.5 * Dd * log(2 * pi) - sum(log(bw))
  lK <- -0.5 * D2 + cst; rm <- apply(lK, 1, max); log(rowSums(exp(lK - rm))) + rm - log(nrow(Xtr))
}
metrics <- function(yhat, yte, C) {
  pv <- sv <- numeric(C)
  for (k in seq_len(C)) { tp <- sum(yhat == k & yte == k); fp <- sum(yhat == k & yte != k)
    tn <- sum(yhat != k & yte != k)
    pv[k] <- if (tp + fp > 0) tp / (tp + fp) else 1; sv[k] <- if (tn + fp > 0) tn / (tn + fp) else 1 }
  c(acc = mean(yhat == yte), prec = mean(pv), spec = mean(sv))
}
one_rep <- function(b, X) {
  set.seed(b * 7919 + 11); D <- ncol(X)
  spl <- lapply(idx_by_sp, function(ix) sample(ix, N_TR + N_TE))
  Xtr <- do.call(rbind, lapply(spl, function(s) X[s[1:N_TR], , drop = FALSE]))
  Xte <- do.call(rbind, lapply(spl, function(s) X[s[N_TR + 1:N_TE], , drop = FALSE]))
  ytr <- rep(seq_len(C), each = N_TR); yte <- rep(seq_len(C), each = N_TE)
  ld <- matrix(0, nrow(Xte), C)
  for (k in seq_len(C)) { bw <- apply(Xtr[ytr == k, , drop = FALSE], 2, bw_silverman)
    ld[, k] <- log_kde_all(Xte, Xtr[ytr == k, , drop = FALSE], bw) }
  m_npkdc <- metrics(max.col(ld, "first"), yte, C)
  ys <- tryCatch(as.integer(as.character(predict(
          svm(Xtr, factor(ytr), kernel = "radial"), Xte))), error = function(e) rep(1L, nrow(Xte)))
  m_svm <- metrics(ys, yte, C)
  trl <- lapply(seq_len(C), function(k) Xtr[ytr == k, , drop = FALSE])
  sel    <- rodeo_detect(trl, keep = D, mc.cores = 1L)      # raw scale
  sel_bc <- rodeo_detect_bc(trl, keep = D, mc.cores = 1L)   # pooled empirical PIT
  list(npkdc = m_npkdc, svm = m_svm, sel = sel, sel_bc = sel_bc)
}
run <- function(label, X) {
  cat(sprintf("\n=== Anuran %s: d=%d ===\n", label, ncol(X)))
  res <- mclapply(seq_len(B), function(b) one_rep(b, X), mc.cores = cores)
  npk <- do.call(rbind, lapply(res, `[[`, "npkdc")); svm <- do.call(rbind, lapply(res, `[[`, "svm"))
  selfreq    <- Reduce(`+`, lapply(res, `[[`, "sel")) / length(res)
  selfreq_bc <- Reduce(`+`, lapply(res, `[[`, "sel_bc")) / length(res)
  nsel    <- vapply(res, function(r) mean(rowSums(r$sel)), numeric(1))
  nsel_bc <- vapply(res, function(r) mean(rowSums(r$sel_bc)), numeric(1))
  f <- function(v) sprintf("%.4f (%.4f)", mean(v), sd(v))
  cat("NPKDC-vd  acc", f(npk[, "acc"]), " prec", f(npk[, "prec"]), " spec", f(npk[, "spec"]), "\n")
  cat("SVM       acc", f(svm[, "acc"]), " prec", f(svm[, "prec"]), " spec", f(svm[, "spec"]), "\n")
  list(label = label, npkdc = npk, svm = svm, selfreq = selfreq, selfreq_bc = selfreq_bc,
       nsel = nsel, nsel_bc = nsel_bc, species = keep_sp, B = length(res))
}

resO <- run("Original (22 MFCC)", mfcc)
set.seed(2026); Xext <- cbind(mfcc, matrix(rnorm(nrow(mfcc) * 5), nrow(mfcc), 5))
colnames(Xext)[23:27] <- paste0("Noise", 1:5)
resE <- run("Extended (22 MFCC + 5 noise)", Xext)
for (sc in c("selfreq", "selfreq_bc")) {
  cat(sprintf("\n[%s] Extended noise-probe (cols 23-27) selection frequency across species:\n", sc))
  print(round(colMeans(resE[[sc]][, 23:27, drop = FALSE]), 4))
  cat(sprintf("[%s] Per-species MFCC selection frequency (Original):\n", sc))
  rownames(resO[[sc]]) <- keep_sp; print(round(resO[[sc]], 2))
}
outdir <- Sys.getenv("RODEO_OUT", file.path("..", "results_R2")); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
saveRDS(list(Original = resO, Extended = resE, species = keep_sp,
             selector = "cross-fitted rodeo + largest gap; raw scale and pooled empirical PIT"),
        file.path(outdir, "real_anuran_results.rds"))
cat(sprintf("\nsaved %s\n", file.path(outdir, "real_anuran_results.rds")))
