# real_digit_attribution.R — per-digit attribution diagnostics with the corrected
# cross-fitted rodeo, for UCI optdigits (8x8, d=64) and MNIST (28x28).
# Needs the binaries written by code/fetch_digit_data.py:
#   digits8_X.bin, digits8_y.bin, mnist_X.bin, mnist_y.bin  (in DIGIT_DIR).
#   DIGIT_DIR=/path RODEO_CORES=14 RODEO_B=1000 Rscript real_digit_attribution.R
#
# Per class and replication (100 training images per digit) it records
#   * the continuous adapted-bandwidth diagnostic profile reported in the paper,
#     (max(b) - b) / (max(b) - min(b)) for the fold-averaged final rodeo
#     bandwidths b (1 = most concentrated pixel);
#   * the rodeo score T_j and the number of pixels the binary gap rule selects.
# Both datasets use the active-pixel pre-screen of V2_digit_simulation.R (ANOVA
# F > 2 and pooled mean intensity > 0.01 on 300 images per digit; 458 MNIST pixels),
# which also removes the constant border pixels whose point masses would otherwise
# masquerade as maximal concentration.
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (length(here) && nzchar(here)) setwd(here)
source("rodeo_core.R"); suppressMessages(library(parallel))
cores <- max(1L, as.integer(Sys.getenv("RODEO_CORES", as.character(max(1L, detectCores() - 2L)))))
B     <- as.integer(Sys.getenv("RODEO_B", "1000"))
DDIR  <- Sys.getenv("DIGIT_DIR", "/tmp")
N_TR  <- 100L; C <- 10L

read_digits <- function(stem, n, d, scale) {
  X <- matrix(readBin(file.path(DDIR, paste0(stem, "_X.bin")), "numeric", n * d,
                      size = 4L, endian = "little"), nrow = n, ncol = d, byrow = TRUE) / scale
  y <- readBin(file.path(DDIR, paste0(stem, "_y.bin")), "integer", n, size = 4L, endian = "little")
  list(X = X, idx = lapply(0:9, function(k) which(y == k)))
}

profile_run <- function(label, X, idx, active, k_top) {
  one <- function(b) {
    set.seed(b * 2711 + 5)
    rows <- lapply(idx, function(ii) {
      r  <- rodeo_select_class(X[sample(ii, N_TR), active, drop = FALSE])
      bw <- r$bw
      list(prof = (max(bw) - bw) / (max(bw) - min(bw) + 1e-12),
           score = r$score, nsel = sum(r$selected),
           nshrunk = sum(bw < max(bw) * (1 - 1e-9)),     # pixels shrunk at least once
           maxscore = max(r$score))
    })
    list(prof  = do.call(rbind, lapply(rows, `[[`, "prof")),
         score = do.call(rbind, lapply(rows, `[[`, "score")),
         nsel  = vapply(rows, `[[`, numeric(1), "nsel"),
         nshrunk  = vapply(rows, `[[`, numeric(1), "nshrunk"),
         maxscore = vapply(rows, `[[`, numeric(1), "maxscore"))
  }
  res   <- mclapply(seq_len(B), one, mc.cores = cores)
  prof  <- Reduce(`+`, lapply(res, `[[`, "prof"))  / length(res)   # 10 x |active|
  score <- Reduce(`+`, lapply(res, `[[`, "score")) / length(res)
  nsel  <- do.call(rbind, lapply(res, `[[`, "nsel"))               # B x 10
  nshrunk  <- do.call(rbind, lapply(res, `[[`, "nshrunk"))
  maxscore <- do.call(rbind, lapply(res, `[[`, "maxscore"))
  lam2 <- (2.0 * sqrt(2 * log(N_TR / 2)))^2                        # shrinkage threshold, fold size N_TR/2
  share <- function(M) apply(M, 1, function(s) { s <- sort(s, decreasing = TRUE); sum(s[seq_len(k_top)]) / sum(s) })
  out <- list(label = label, B = length(res), N_TR = N_TR, active = active, k_top = k_top,
              uniform_share = k_top / length(active), profile = prof, score = score,
              nsel = nsel, nshrunk = nshrunk, maxscore = maxscore, lambda2 = lam2,
              top_share_profile = share(prof), top_share_score = share(score))
  cat(sprintf("\n=== %s: %d pixels, B=%d, N_TR=%d ===\n", label, length(active), out$B, N_TR))
  cat(sprintf("top-%d share of the adapted-bandwidth profile (uniform = %.4f):\n", k_top, out$uniform_share))
  print(round(out$top_share_profile, 4))
  cat(sprintf("top-%d share of the rodeo score:\n", k_top)); print(round(out$top_share_score, 4))
  cat("mean number of gap-selected pixels per digit:\n"); print(round(colMeans(nsel), 3))
  cat(sprintf("class fits with any pixel shrunk: %.4f; mean #pixels shrunk: %.3f\n",
              mean(nshrunk > 0), mean(nshrunk)))
  cat(sprintf("largest pixel score: mean %.2f, max %.2f (shrinkage threshold lambda^2 = %.2f)\n",
              mean(maxscore), max(maxscore), lam2))
  out
}

screen_active <- function(X, idx, sn = 300L) {
  set.seed(1L)
  sc <- do.call(rbind, lapply(idx, function(ii) X[sample(ii, sn), ]))
  ys <- rep(0:9, each = sn); Nsc <- nrow(sc); d <- ncol(X)
  gm <- colMeans(sc); cm <- do.call(rbind, lapply(0:9, function(k) colMeans(sc[ys == k, , drop = FALSE])))
  SSB <- colSums((cm - matrix(gm, C, d, byrow = TRUE))^2) * sn
  SSW <- colSums(sc^2) - colSums(cm^2) * sn
  Fst <- (SSB / (C - 1)) / pmax(SSW / (Nsc - C), 1e-12)
  which(Fst > 2 & colMeans(sc) > 0.01)
}

uci <- read_digits("digits8", 5620L, 64L, 16)
res_uci <- profile_run("UCI optdigits 8x8", uci$X, uci$idx, screen_active(uci$X, uci$idx), k_top = 8L)

mn <- read_digits("mnist", 70000L, 784L, 255)
res_mn <- profile_run("MNIST 28x28", mn$X, mn$idx, screen_active(mn$X, mn$idx), k_top = 20L)

outdir <- Sys.getenv("RODEO_OUT", file.path("..", "results_R2")); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
saveRDS(list(uci8x8 = res_uci, mnist = res_mn,
             selector = "corrected cross-fitted rodeo; raw pixel scale"),
        file.path(outdir, "real_digit_attribution_results.rds"))
cat(sprintf("\nsaved %s\n", file.path(outdir, "real_digit_attribution_results.rds")))
