# study1_jaccard_figure.R ------------------------------------------------------
# Main-paper Figure (Study 1): true vs estimated inter-class Jaccard matrices for
# the proposed cross-fitted rodeo + largest-gap selector.  Reads the B=1000 output
# of rodeo_study1.R (results_R2/study1_rodeo_results.rds) and writes
# figures/V1_jaccard_final.pdf.  Base R graphics only.
#   Rscript study1_jaccard_figure.R
# ---------------------------------------------------------------------------
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (length(here) && nzchar(here)) setwd(here)
res <- readRDS(file.path("..", "results_R2", "study1_rodeo_results.rds"))
J0 <- res$true_jac; Jhat <- res$est_jac; C <- nrow(J0)
pal <- colorRampPalette(c("#FFFFFF", "#C7E9C0", "#41AB5D", "#00441B"))(101)

panel <- function(M, main) {
  image(1:C, 1:C, t(M[C:1, , drop = FALSE]), zlim = c(0, 1), col = pal,
        axes = FALSE, xlab = expression("Class " * y * "'"), ylab = "", main = main,
        cex.main = 0.95, font.main = 1)
  axis(1, at = 1:C, cex.axis = 0.8, tick = FALSE, line = -0.6)
  axis(2, at = 1:C, labels = C:1, las = 1, cex.axis = 0.8, tick = FALSE, line = -0.2)
  mtext(expression("Class " * y), side = 2, line = 1.9, cex = 0.8)
  v <- t(M[C:1, , drop = FALSE])
  lab <- ifelse(abs(v) < 5e-3, "0", formatC(v, format = "f", digits = 2))
  text(rep(1:C, C), rep(1:C, each = C), lab, cex = 0.52,
       col = ifelse(v > 0.6, "white", "black"))
  box(col = "grey60")
}

pdf(file.path("..", "figures", "V1_jaccard_final.pdf"), width = 6.6, height = 3.15)
layout(matrix(1:3, 1), widths = c(1, 1, 0.16))
par(mar = c(3.2, 3.4, 2.2, 0.6), mgp = c(1.8, 0.5, 0))
panel(J0, "True Jaccard")
panel(Jhat, sprintf("NPKDC-vd estimate (mean of B = %d)", res$B))
par(mar = c(3.2, 0.4, 2.2, 2.8))
image(1, seq(0, 1, length.out = 101), t(seq(0, 1, length.out = 101)), col = pal,
      axes = FALSE, xlab = "", ylab = "")
axis(4, at = seq(0, 1, 0.25), las = 1, cex.axis = 0.7)
mtext("Jaccard", side = 3, line = 0.5, cex = 0.75)
box(col = "grey60")
invisible(dev.off())
cat(sprintf("wrote figures/V1_jaccard_final.pdf (B=%d; max|Jhat-J|=%.3g)\n",
            res$B, max(abs(Jhat - J0))))
