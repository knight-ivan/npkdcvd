# N3_figure.R — improved boxplot for Study 5

FIGURES_DIR <- file.path("..", "manuscript", "figures")
DATA_DIR    <- file.path("..", "data")
set.seed(1 * 999 + 42)
C  <- 5L; D <- 20L; N_Y <- 150L; N_TEST <- 100L
DETECT_TAU <- -1.5
R_Y <- rep(3L, C)
rel_start <- cumsum(c(1L, R_Y[-C]))
REL_SETS  <- lapply(1:C, function(y) rel_start[y]:(rel_start[y] + R_Y[y] - 1L))

bw_silverman <- function(x) {
  s <- sd(x); if (s < 1e-12) return(1e-6)
  1.06 * s * length(x)^(-0.2)
}

gen_data <- function(n_y, n_test) {
  train <- lapply(1:C, function(y) {
    X <- matrix(rnorm(n_y * D), nrow = n_y, ncol = D)
    X[, REL_SETS[[y]]] <- matrix(rnorm(n_y * R_Y[y], mean = 0.5, sd = 0.5),
                                  nrow = n_y, ncol = R_Y[y])
    X
  })
  test <- lapply(1:C, function(y) {
    X <- matrix(rnorm(n_test * D), nrow = n_test, ncol = D)
    X[, REL_SETS[[y]]] <- matrix(rnorm(n_test * R_Y[y], mean = 0.5, sd = 0.5),
                                  nrow = n_test, ncol = R_Y[y])
    X
  })
  list(train = train, test = test)
}

per_obs_contrast <- function(xi, y, train_list, bw_sil) {
  sapply(1:D, function(j) {
    log_vals <- sapply(1:C, function(yp) {
      x_tr <- train_list[[yp]][, j]
      h    <- bw_sil[yp, j]
      mean(dnorm((xi[j] - x_tr) / h) / h)
    })
    log_vals_safe <- log(pmax(log_vals, 1e-300))
    log_vals_safe[y] - mean(log_vals_safe[-y])
  })
}

detect_class_level <- function(train_list) {
  bw_sil <- matrix(0, C, D)
  for (y in 1:C) bw_sil[y,] <- apply(train_list[[y]], 2L, bw_silverman)
  detected <- matrix(0L, C, D)
  for (j in 1:D) {
    hj <- bw_sil[, j]; sj <- sd(hj)
    if (sj < 1e-10) next
    zj <- (hj - mean(hj)) / sj
    detected[, j] <- as.integer(zj < DETECT_TAU)
  }
  detected
}

dat1 <- gen_data(N_Y, N_TEST)
bw_sil1 <- matrix(0, C, D)
for (y in 1:C) bw_sil1[y,] <- apply(dat1$train[[y]], 2L, bw_silverman)
phi_c1 <- t(apply(dat1$test[[1]], 1, function(xi)
  per_obs_contrast(xi, 1L, dat1$train, bw_sil1)))
cl_det1    <- detect_class_level(dat1$train)
relevant_c1 <- which(cl_det1[1,] == 1L)

col_rel   <- "#1b7837"
col_irrel <- "#aaaaaa"
box_cols  <- ifelse(1:D %in% REL_SETS[[1]], col_rel, col_irrel)

do_plot <- function() {
  par(mar = c(4, 4.4, 2.5, 1))
  boxplot(phi_c1, at = 1:D, col = box_cols, border = box_cols,
          xlab = "Variable index", ylab = expression(phi[j](x[i], y)),
          main = "",
          whisklty = 1, outpch = 20, outcex = 0.35, outcol = box_cols,
          cex.axis = 0.88, cex.lab = 0.96, las = 1, boxwex = 0.62)
  abline(h = 0, col = "grey60", lty = 2, lwd = 0.9)
  abline(v = 3.5, col = "grey40", lty = 3, lwd = 1.0)
  # group labels at top
  mtext("True R_y  (3 vars)", side = 3, at = 2, cex = 0.76, line = 0.9,
        col = col_rel, font = 2)
  mtext("Irrelevant  (17 vars)", side = 3, at = 11.5, cex = 0.76, line = 0.9,
        col = "grey40")
  rug(relevant_c1, side = 3, ticksize = 0.06, col = "#d6604d", lwd = 2)
  legend("topright", bty = "n", cex = 0.82,
         legend = c(expression("True "~R[y]~" (3 vars)"),
                    "Irrelevant",
                    "Class-level detected"),
         fill   = c(col_rel, col_irrel, "#d6604d"),
         border = NA)
}

fig_path <- file.path(FIGURES_DIR, "N3_contrast.eps")
postscript(fig_path, width = 8.5, height = 4.2,
           horizontal = FALSE, onefile = FALSE, paper = "special")
do_plot()
dev.off()
cat("EPS saved:", fig_path, "\n")

pdf_path <- file.path(FIGURES_DIR, "N3_contrast-eps-converted-to.pdf")
pdf(pdf_path, width = 8.5, height = 4.2)
do_plot()
dev.off()
cat("PDF saved:", pdf_path, "\n")
