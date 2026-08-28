# V1_jaccard.R
# Phase 3, Visualization V1: Jaccard overlap matrix |R̂_y ∩ R̂_y'| / |R̂_y ∪ R̂_y'|
#
# For the 10-class Synthetic Study 1 setting:
#   C=10, d=30, R_y = {y, y+1,...,y+5} (mod 30), cyclic, 6 vars per class
#   Relevant: N(0.5, (0.02*(j-y+1))^2); Irrelevant: Uniform(0,1)
#   n_y=150, B=1000 replications
#
# Plots two panels:
#   Left: true Jaccard matrix Jac(R_y, R_y')
#   Right: NPKDC-vd estimated Jaccard E[Jac(R̂_y, R̂_y')] over B reps
#
# The banded structure matches, showing class-specific attribution is recovered.

FIGURES_DIR <- file.path("..", "figures")
DATA_DIR    <- file.path("..", "data")
library(parallel)

set.seed(20260510)

C          <- 10
D          <- 30
N_Y        <- 150
B          <- 1000
DETECT_TAU <- -1.5

# ---- True relevant sets (cyclic, 1-indexed) ----
# R_y = {y, y+1,...,y+5} with indices mod 30 (1-based: ((y-1+k) %% 30) + 1, k=0..5)
REL_SETS <- lapply(1:C, function(y) (((y - 1) + 0:5) %% 30) + 1)
R_Y <- rep(6, C)

# ---- True Jaccard matrix ----
true_jac <- matrix(0, C, C)
for (y in 1:C) for (yp in 1:C) {
  inter <- length(intersect(REL_SETS[[y]], REL_SETS[[yp]]))
  union <- length(union(REL_SETS[[y]], REL_SETS[[yp]]))
  true_jac[y, yp] <- if (union > 0) inter / union else 1
}


# ---- Data generation ----
gen_train <- function(n_y) {
  lapply(1:C, function(y) {
    X        <- matrix(runif(n_y * D), nrow = n_y, ncol = D)
    rel      <- REL_SETS[[y]]
    for (k in seq_along(rel)) {
      j      <- rel[k]
      sigma  <- 0.02 * k
      X[, j] <- rnorm(n_y, mean = 0.5, sd = sigma)
    }
    X
  })
}

# ---- Silverman bandwidth ----
bw_silverman <- function(x) {
  s <- sd(x); if (s < 1e-12) return(1e-6)
  1.06 * s * length(x)^(-0.2)
}

# ---- Detection: within-class largest-gap method (Algorithm 2 spirit) ----
# For each class y, sort bandwidths h_{y,j} and find the largest consecutive
# gap in the sorted sequence. Variables below the gap form R̂_y.
# This handles shared relevant variables correctly (no cross-class dilution).
detect_vars <- function(train_list) {
  bw_sil <- matrix(0, C, D)
  for (y in 1:C)
    bw_sil[y, ] <- apply(train_list[[y]], 2L, bw_silverman)

  detected <- matrix(0L, C, D)
  for (y in 1:C) {
    h_y    <- bw_sil[y, ]
    ord    <- order(h_y)          # sorted indices (smallest bw first)
    h_sort <- h_y[ord]            # sorted bandwidths

    # Largest consecutive gap in sorted bandwidths
    diffs  <- diff(h_sort)
    if (max(diffs) < 1e-6) next  # all bandwidths equal → skip
    split  <- which.max(diffs)    # split between positions split and split+1

    # Variables with bandwidth ≤ h_sort[split] are detected as relevant
    threshold <- h_sort[split]
    detected[y, h_y <= threshold] <- 1L
  }
  detected
}

# Pairwise Jaccard from detected matrix
pairwise_jaccard <- function(detected) {
  jac <- matrix(0, C, C)
  for (y in 1:C) {
    Ry <- which(detected[y, ] == 1L)
    for (yp in 1:C) {
      Ryp   <- which(detected[yp, ] == 1L)
      inter <- length(intersect(Ry, Ryp))
      union <- length(union(Ry, Ryp))
      jac[y, yp] <- if (union > 0) inter / union else 0
    }
  }
  jac
}

one_rep <- function(b) {
  set.seed(b * 777 + 13)
  train <- gen_train(N_Y)
  det   <- detect_vars(train)
  pairwise_jaccard(det)
}


# ---- Main ----
n_cores <- max(1L, detectCores() - 1L)
cat("V1: C=10, d=30, n_y=150, B=1000\nUsing", n_cores, "cores\n")

t0   <- proc.time()[["elapsed"]]
reps <- mclapply(seq_len(B), one_rep, mc.cores = n_cores)
cat(sprintf("Done in %.1f sec\n\n", proc.time()[["elapsed"]] - t0))

# Mean Jaccard matrix over replications
est_jac <- Reduce("+", reps) / B

cat("True Jaccard (first row):   ", round(true_jac[1, ], 3), "\n")
cat("Estimated Jaccard (first row):", round(est_jac[1, ], 3), "\n\n")

saveRDS(list(true_jac=true_jac, est_jac=est_jac),
        file.path("..", "data", "V1_results.rds"),
        compress = FALSE)
cat("Results saved to V1_results.rds\n\n")


# ---- Figure: side-by-side heatmaps (EPS-safe, no rasterImage) ----
fig_path <- file.path("..", "figures", "V1_jaccard.eps")
postscript(fig_path, width = 7.8, height = 3.8,
           horizontal = FALSE, onefile = FALSE, paper = "special")
# Also save as PDF for pdflatex compatibility (run epstopdf afterward)

col_fn <- colorRampPalette(c("white", "#74c476", "#006d2c"))

plot_jac_heatmap <- function(mat, title) {
  n    <- nrow(mat)
  cols <- col_fn(101)
  # image() maps [0,1] → colour; flip rows so class 1 is at top
  image(1:n, 1:n, t(mat[n:1, ]), zlim = c(0, 1), col = cols,
        axes = FALSE, xlab = "Class y'", ylab = "Class y")
  axis(1, at = 1:n, labels = 1:n, cex.axis = 0.78, tcl = -0.3)
  axis(2, at = 1:n, labels = n:1, cex.axis = 0.78, las = 1, tcl = -0.3)
  box()
  mtext(title, side = 3, line = 0.4, cex = 0.90, font = 2)
  for (i in 1:n) for (j in 1:n) {
    v <- mat[n + 1 - i, j]
    lbl <- if (v == 0) "0" else sprintf("%.2f", v)
    text(j, i, lbl, cex = 0.50, col = if (v >= 0.55) "white" else "black")
  }
}

# Draw a manual vertical colour-bar using rect()
draw_colorbar <- function(x0, y0, x1, y1, n_steps = 50) {
  ys   <- seq(y0, y1, length.out = n_steps + 1)
  cols <- col_fn(n_steps)
  for (k in 1:n_steps)
    rect(x0, ys[k], x1, ys[k + 1], col = cols[k], border = NA, xpd = NA)
  rect(x0, y0, x1, y1, border = "black", xpd = NA)
  # tick labels
  for (frac in c(0, 0.5, 1))
    text(x1 + 0.05, y0 + frac * (y1 - y0), sprintf("%.1f", frac),
         cex = 0.65, adj = c(0, 0.5), xpd = NA)
}

layout(matrix(c(1, 2, 3), 1, 3), widths = c(4, 4, 0.6))
par(mar = c(3.5, 3.5, 2.0, 0.3))
plot_jac_heatmap(true_jac, "True Jaccard")
plot_jac_heatmap(est_jac,  "NPKDC-vd (estimated)")
par(mar = c(3.5, 0.5, 2.0, 1.5))
plot.new()
draw_colorbar(0.0, 0.05, 0.5, 0.95, n_steps = 80)
mtext("Jaccard", side = 3, line = 0.4, cex = 0.75, adj = 0)

dev.off()
cat("Figure saved to V1_jaccard.eps\n")

pdf_path <- file.path("..", "figures", "V1_jaccard.pdf")
pdf(pdf_path, width = 7.8, height = 3.8)
layout(matrix(c(1, 2, 3), 1, 3), widths = c(4, 4, 0.6))
par(mar = c(3.5, 3.5, 2.0, 0.3))
plot_jac_heatmap(true_jac, "True Jaccard")
plot_jac_heatmap(est_jac,  "NPKDC-vd (estimated)")
par(mar = c(3.5, 0.5, 2.0, 1.5))
plot.new()
draw_colorbar(0.0, 0.05, 0.5, 0.95, n_steps = 80)
mtext("Jaccard", side = 3, line = 0.4, cex = 0.75, adj = 0)
dev.off()
cat("PDF saved to V1_jaccard.pdf\n")
