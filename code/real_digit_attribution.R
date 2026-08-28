# real_digit_attribution.R — per-digit MNIST attribution under the rebuilt rodeo.
# Regenerates fig:digit_heatmap (V2_digit_28x28.pdf) with rodeo_detect.
# Needs /tmp/mnist_X.bin,/tmp/mnist_y.bin (run code/fetch_digit_data.py first).
#   RODEO_CORES=22 RODEO_B=200 Rscript real_digit_attribution.R
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (length(here) && nzchar(here)) setwd(here)
source("rodeo_core.R"); suppressMessages(library(parallel))
cores <- max(1L, as.integer(Sys.getenv("RODEO_CORES", as.character(detectCores() - 2L))))
B     <- as.integer(Sys.getenv("RODEO_B", "200"))
N_TR  <- 100L; C <- 10L; dB <- 784L; nB <- 70000L

X <- matrix(readBin("/tmp/mnist_X.bin", "numeric", nB * dB, size = 4L, endian = "little"),
            nrow = nB, ncol = dB, byrow = TRUE) / 255.0
y <- readBin("/tmp/mnist_y.bin", "integer", nB, size = 4L, endian = "little")
idx <- lapply(0:9, function(k) which(y == k))

# Pre-screen to active pixels (same rule as V2_digit_simulation.R)
set.seed(1L); sn <- 300L
sc <- do.call(rbind, lapply(idx, function(ii) X[sample(ii, sn), ]))
ys <- rep(0:9, each = sn); Nsc <- nrow(sc)
gm <- colMeans(sc); cm <- do.call(rbind, lapply(0:9, function(k) colMeans(sc[ys == k, , drop = FALSE])))
SSB <- colSums((cm - matrix(gm, C, dB, byrow = TRUE))^2) * sn
SSW <- colSums(sc^2) - colSums(cm^2) * sn
Fst <- (SSB / (C - 1)) / pmax(SSW / (Nsc - C), 1e-12)
active <- which(Fst > 2 & colMeans(sc) > 0.01)
cat(sprintf("MNIST attribution: active pixels %d/%d, B=%d, cores=%d\n", length(active), dB, B, cores))

# Continuous rodeo SCORE T_j per pixel (graphical profile). The binary gap rule is
# for sparse-recovery in the simulations; real-data maps use the continuous score.
one <- function(b) {
  set.seed(b * 2711 + 5)
  do.call(rbind, lapply(idx, function(ii) {
    Xy <- X[sample(ii, N_TR), active, drop = FALSE]
    b  <- rodeo_select_class(Xy)$bw                              # adapted bandwidths (small = concentrated)
    (max(b) - b) / (max(b) - min(b) + 1e-12)                     # -> [0,1], 1 = most concentrated
  }))                                                            # 10 x |active|
}
res  <- mclapply(seq_len(B), one, mc.cores = cores)
freq <- Reduce(`+`, res) / length(res)                           # 10 x |active| mean score
full <- matrix(0.0, C, dB); full[, active] <- freq               # 10 x 784
saveRDS(list(freq_full = full, active = active, B = B), file.path("..", "data", "real_digit_attr_results.rds"))

# 10-panel heatmap (bright = high rodeo score = class density concentrated here)
pdf(file.path("..", "figures", "V2_digit_28x28_rodeo.pdf"), width = 10, height = 4.5)
op <- par(mfrow = c(2, 5), mar = c(1, 1, 2, 1))
for (k in 1:C) {
  img <- matrix(full[k, ], 28, 28)[, 28:1]                       # orient upright
  image(img, axes = FALSE, col = grey.colors(64, start = 0, end = 1), main = paste("digit", k - 1))
}
par(op); dev.off()
cat("saved data/real_digit_attr_results.rds and figures/V2_digit_28x28_rodeo.pdf\n")
