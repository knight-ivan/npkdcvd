# N2_figure.R  — improved detection-quality vs n_y figure

FIGURES_DIR <- file.path("..", "figures")
DATA_DIR    <- file.path("..", "data")
res <- readRDS(file.path(DATA_DIR, "N2_results.rds"))

N_Y_GRID  <- c(30, 50, 75, 100, 150, 200, 300, 500)
rec_agg   <- sapply(res, function(r) mean(r$rec_mn))
prec_agg  <- sapply(res, function(r) mean(r$prec_mn))
f1_agg    <- sapply(res, function(r) mean(r$f1_mn))
exact_agg <- sapply(res, function(r) mean(r$exact_mn))

cols <- c("#1b7837", "#762a83", "#4393c3", "#d6604d")
pchs <- c(16, 17, 15, 4)
ltys <- c(1, 2, 3, 1)
lwds <- c(2.0, 2.0, 2.0, 1.6)

do_plot <- function() {
  par(mar = c(4.2, 4.4, 1.2, 1))
  plot(N_Y_GRID, rec_agg, type = "b",
       ylim = c(0.30, 1.04),
       xlab = expression(n[y]~"(per-class training size)"),
       ylab = "Detection quality (mean over classes and reps)",
       col = cols[1], pch = pchs[1], lty = ltys[1], lwd = lwds[1],
       cex = 1.2, cex.axis = 0.92, cex.lab = 0.96, las = 1)
  abline(h = 1,   col = "grey75", lty = 3, lwd = 0.8)
  abline(v = 150, col = "grey50", lty = 4, lwd = 1.2)
  mtext(expression(n[y]==150), side = 3, at = 150, cex = 0.72, line = 0.1, col = "grey40")
  lines(N_Y_GRID, prec_agg,  type = "b", col = cols[2], pch = pchs[2],
        lty = ltys[2], lwd = lwds[2], cex = 1.2)
  lines(N_Y_GRID, f1_agg,    type = "b", col = cols[3], pch = pchs[3],
        lty = ltys[3], lwd = lwds[3], cex = 1.2)
  lines(N_Y_GRID, exact_agg, type = "b", col = cols[4], pch = pchs[4],
        lty = ltys[4], lwd = lwds[4], cex = 1.2)
  legend("bottomright", bty = "n",
         legend = c("Recall", "Precision", "F1",
                    expression(P(hat(R)[y]==R[y]))),
         col = cols, pch = pchs, lty = ltys, lwd = lwds, cex = 0.88)
}

fig_path <- file.path(FIGURES_DIR, "N2_detprob.eps")
postscript(fig_path, width = 6.5, height = 4.2,
           horizontal = FALSE, onefile = FALSE, paper = "special")
do_plot()
dev.off()
cat("EPS saved:", fig_path, "\n")

pdf_path <- file.path(FIGURES_DIR, "N2_detprob-eps-converted-to.pdf")
pdf(pdf_path, width = 6.5, height = 4.2)
do_plot()
dev.off()
cat("PDF saved:", pdf_path, "\n")
