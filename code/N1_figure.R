# N1_figure.R
# Improved Z-score figure for Study 2 (d=200, C=5)
# Shows variables 1-40 only (all relevant vars are in 1-29; vars 30+ are all Z≈0)

FIGURES_DIR <- file.path("..", "figures")
DATA_DIR    <- file.path("..", "data")
set.seed(42)
C   <- 5L
D   <- 200L
N_Y <- 150L
R_Y <- c(3L, 5L, 5L, 8L, 8L)
rel_start <- c(1L, 4L, 9L, 14L, 22L)
REL_SETS  <- lapply(1:C, function(y) rel_start[y]:(rel_start[y] + R_Y[y] - 1L))
DETECT_TAU <- -1.5

bw_silverman <- function(x) {
  s <- sd(x); if (s < 1e-12) return(1e-6)
  1.06 * s * length(x)^(-0.2)
}

gen_data <- function(irrel_type = "normal") {
  rho <- 0.5
  lapply(1:C, function(y) {
    rel <- REL_SETS[[y]]; nr <- R_Y[y]; nirr <- D - nr
    mu_rel <- rep(c(0.3, 0.5, 0.7), length.out = nr)
    X <- matrix(0, nrow = N_Y, ncol = D)
    X[, rel] <- matrix(rnorm(N_Y * nr, mean = rep(mu_rel, each = N_Y), sd = 0.1),
                       nrow = N_Y, ncol = nr)
    X[, -rel] <- matrix(rnorm(N_Y * nirr), nrow = N_Y, ncol = nirr)
    X
  })
}

train <- gen_data()

# Compute Silverman bandwidths and Z-scores
bw_sil <- matrix(0, C, D)
for (y in 1:C) bw_sil[y,] <- apply(train[[y]], 2L, bw_silverman)

zscore <- matrix(0, C, D)
for (j in 1:D) {
  hj <- bw_sil[, j]; sj <- sd(hj)
  if (sj < 1e-10) next
  zscore[, j] <- (hj - mean(hj)) / sj
}

# ---- Figure ----
SHOW_MAX <- 40L   # show only variables 1-40; vars 30-200 are all Z≈0 (redundant)
col_rel   <- "#d6604d"
col_irrel <- "grey68"
YLIM      <- c(-4.2, 2.8)

do_panel <- function(y) {
  rel <- REL_SETS[[y]]
  irrel_idx <- setdiff(1:SHOW_MAX, rel)
  rel_idx   <- intersect(rel, 1:SHOW_MAX)
  z <- zscore[y, 1:SHOW_MAX]

  plot(irrel_idx, z[irrel_idx], pch = 20, cex = 0.28, col = col_irrel,
       xlim = c(0.5, SHOW_MAX + 0.5), ylim = YLIM, xaxt = "n", yaxt = "n",
       xlab = "", ylab = "", main = "")
  mtext(sprintf("Class %d  (r = %d)", y, R_Y[y]),
        side = 3, line = 0.25, cex = 0.82, font = 2)
  abline(h = 0,          col = "grey80", lty = 3, lwd = 0.8)
  abline(h = DETECT_TAU, col = "#4393c3", lty = 2, lwd = 1.4)
  points(rel_idx, z[rel_idx], pch = 16, cex = 1.35, col = col_rel)
  axis(1, at = c(1, 10, 20, 29, 40), cex.axis = 0.72, tcl = -0.3)
  axis(2, at = c(-4, -3, -2, -1, 0, 1, 2), cex.axis = 0.72, las = 1, tcl = -0.3)
  box()
}

fig_path <- file.path(FIGURES_DIR, "N1_zscore.eps")
postscript(fig_path, width = 9, height = 5,
           horizontal = FALSE, onefile = FALSE, paper = "special")
layout(matrix(c(1,2,3,4,5,6), nrow = 2, ncol = 3, byrow = TRUE),
       widths = c(1,1,1), heights = c(1,1))
par(mar = c(3.2, 3.4, 1.8, 0.6), oma = c(2, 0.5, 0.5, 0.5))
for (y in 1:5) do_panel(y)
# 6th panel: shared legend + note
par(mar = c(3.2, 0.5, 1.8, 0.5))
plot.new()
legend("center", bty = "n", cex = 1.05,
       legend = c("Relevant variable", "Irrelevant variable",
                  "Detection threshold  tau = -1.5"),
       pch    = c(16, 20, NA),
       col    = c(col_rel, col_irrel, "#4393c3"),
       lty    = c(NA, NA, 2), lwd = c(NA, NA, 1.4),
       pt.cex = c(1.35, 0.7, NA))
text(0.5, 0.12, "Note: variables 30-200\n(171 irrelevant) omitted;\nall have Z = 0.",
     cex = 0.90, adj = c(0.5, 0.5), col = "grey40")
mtext("Variable index (1-40 shown)", side = 1, outer = TRUE, line = 0.7, cex = 0.82)
dev.off()
cat("EPS saved:", fig_path, "\n")

pdf_path <- file.path(FIGURES_DIR, "N1_zscore-eps-converted-to.pdf")
pdf(pdf_path, width = 9, height = 5)
layout(matrix(c(1,2,3,4,5,6), nrow = 2, ncol = 3, byrow = TRUE),
       widths = c(1,1,1), heights = c(1,1))
par(mar = c(3.2, 3.4, 1.8, 0.6), oma = c(2, 0.5, 0.5, 0.5))
for (y in 1:5) do_panel(y)
par(mar = c(3.2, 0.5, 1.8, 0.5))
plot.new()
legend("center", bty = "n", cex = 1.05,
       legend = c("Relevant variable", "Irrelevant variable",
                  "Detection threshold  tau = -1.5"),
       pch    = c(16, 20, NA),
       col    = c(col_rel, col_irrel, "#4393c3"),
       lty    = c(NA, NA, 2), lwd = c(NA, NA, 1.4),
       pt.cex = c(1.35, 0.7, NA))
text(0.5, 0.12, "Note: variables 30-200\n(171 irrelevant) omitted;\nall have Z = 0.",
     cex = 0.90, adj = c(0.5, 0.5), col = "grey40")
mtext("Variable index (1-40 shown)", side = 1, outer = TRUE, line = 0.7, cex = 0.82)
dev.off()
cat("PDF saved:", pdf_path, "\n")
