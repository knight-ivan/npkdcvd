# N6_figure_v2.R — improved heatmap figure for Study 6

FIGURES_DIR <- file.path("..", "manuscript", "figures")
DATA_DIR    <- file.path("..", "data")
library(parallel)

set.seed(20260515)
C <- 3L; D <- 15L; N_Y <- 200L; B <- 500L
R_Y      <- c(2L, 1L, 2L)
REL_SETS <- list(c(1L, 2L), c(3L), c(4L, 5L))

gen_data <- function(n_y) {
  lapply(1:C, function(y) {
    X <- matrix(rnorm(n_y * D), nrow = n_y, ncol = D)
    if (y == 1L) { X[,1] <- rnorm(n_y,0.5,0.2); X[,2] <- rnorm(n_y,0.5,0.2) }
    if (y == 2L) { X[,3] <- rnorm(n_y,3.0,0.2) }
    if (y == 3L) { X[,4] <- rnorm(n_y,0.5,0.2); X[,5] <- rnorm(n_y,0.5,0.2) }
    X
  })
}
bw_sil <- function(x) { s <- sd(x); if(s<1e-12) return(1e-6); 1.06*s*length(x)^(-0.2) }

one_rep_sel <- function(b) {
  set.seed(b * 1313 + 77)
  train   <- gen_data(N_Y)
  bw      <- matrix(0, C, D)
  for (y in 1:C) bw[y,] <- apply(train[[y]], 2, bw_sil)
  np_sel  <- matrix(0L, C, D)
  for (y in 1:C) { idx <- order(bw[y,])[1:R_Y[y]]; np_sel[y, idx] <- 1L }
  rf_sel  <- matrix(0L, C, D)
  if (requireNamespace("ranger", quietly=TRUE)) {
    X_all <- do.call(rbind, train)
    ylab  <- rep(1:C, times=sapply(train, nrow))
    for (y in 1:C) {
      ybin <- factor(as.integer(ylab==y), levels=c(0L,1L))
      df   <- data.frame(y=ybin, X_all)
      mod  <- ranger::ranger(y~., data=df, num.trees=200, importance="impurity", verbose=FALSE)
      imp  <- mod$variable.importance
      idx  <- order(imp, decreasing=TRUE)[1:R_Y[y]]
      rf_sel[y, idx] <- 1L
    }
  }
  list(np=np_sel, rf=rf_sel)
}

n_cores <- max(1L, detectCores()-1L)
cat("Computing selection frequencies (B=500)...\n")
sels    <- mclapply(1:B, one_rep_sel, mc.cores=n_cores)
np_freq <- Reduce("+", lapply(sels, `[[`, "np")) / B
rf_freq <- Reduce("+", lapply(sels, `[[`, "rf")) / B
cat("Done.\n")

col_fn <- colorRampPalette(c("white","#fc8d59","#d73027"))

plot_freq <- function(mat, title, rel_cols_list) {
  n <- nrow(mat); d <- ncol(mat)
  cols <- col_fn(101)
  image(1:d, 1:n, t(mat[n:1, ]), zlim=c(0,1), col=cols,
        axes=FALSE, xlab="Variable index", ylab="Class")
  axis(1, at=1:d, labels=1:d, cex.axis=0.78, tcl=-0.3)
  axis(2, at=1:n, labels=n:1, cex.axis=0.90, las=1, tcl=-0.3)
  box()
  mtext(title, side=3, line=0.45, cex=0.90, font=2)
  for (i in 1:n) for (j in 1:d) {
    v   <- mat[n+1-i, j]
    lbl <- if(v==0) "0" else sprintf("%.2f",v)
    text(j, i, lbl, cex=0.60, col=if(v>=0.55)"white" else "black",
         font=if(v>0 && v<1) 1 else 1)
  }
  for (y in 1:n) {
    rel <- rel_cols_list[[n+1-y]]
    for (j in rel) rect(j-0.5, y-0.5, j+0.5, y+0.5, border="#1b7837", lwd=1.8)
  }
}

draw_colorbar <- function(x0, y0, x1, y1, n_steps=60) {
  ys   <- seq(y0, y1, length.out=n_steps+1)
  cols <- col_fn(n_steps)
  for (k in 1:n_steps) rect(x0, ys[k], x1, ys[k+1], col=cols[k], border=NA, xpd=NA)
  rect(x0, y0, x1, y1, border="black", xpd=NA)
  for (frac in c(0,0.5,1))
    text(x1+0.08, y0+frac*(y1-y0), sprintf("%.1f",frac), cex=0.68, adj=c(0,0.5), xpd=NA)
}

do_fig <- function() {
  layout(matrix(c(1,2,3), 1, 3), widths=c(4,4,0.65))
  par(mar=c(3.5,3.4,2.2,0.3))
  plot_freq(np_freq, "NPKDC-vd (top-r)", REL_SETS)
  plot_freq(rf_freq, "OvR-RF (top-r)",   REL_SETS)
  par(mar=c(3.5,0.5,2.2,1.6))
  plot.new()
  draw_colorbar(0.0, 0.05, 0.5, 0.95)
  mtext("Freq.", side=3, line=0.45, cex=0.78, adj=0)
}

fig_path <- file.path(FIGURES_DIR, "N6_selfreq.eps")
postscript(fig_path, width=8.5, height=4.5,
           horizontal=FALSE, onefile=FALSE, paper="special")
do_fig()
dev.off()
cat("EPS saved:", fig_path, "\n")

pdf_path <- file.path(FIGURES_DIR, "N6_selfreq-eps-converted-to.pdf")
pdf(pdf_path, width=8.5, height=4.5)
do_fig()
dev.off()
cat("PDF saved:", pdf_path, "\n")

cat(sprintf("\nKey numbers:\n  NPKDC-vd P(X3 for Class1)=%.3f\n  OvR-RF   P(X3 for Class1)=%.3f\n",
            np_freq[1,3], rf_freq[1,3]))
