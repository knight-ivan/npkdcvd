# Does the MNIST attribution map stay blank as N_TR grows, or does structure
# emerge? Blank at N_TR=100 may be a power artifact of the standardized statistic.
# Run from code/ so source("rodeo_core.R") resolves. Reads /tmp/mnist_*.bin.
source("rodeo_core.R"); suppressMessages(library(parallel))
cores <- as.integer(Sys.getenv("RODEO_CORES","12"))
dB <- 784L; nB <- 70000L; C <- 10L
X <- matrix(readBin("/tmp/mnist_X.bin","numeric",nB*dB,size=4,endian="little"),
            nrow=nB, ncol=dB, byrow=TRUE)/255
y <- readBin("/tmp/mnist_y.bin","integer",nB,size=4,endian="little")
idx <- lapply(0:9, function(k) which(y==k))

# same pre-screen as real_digit_attribution.R
set.seed(1); sn <- 300L
sc <- do.call(rbind, lapply(idx, function(ii) X[sample(ii,sn),]))
ys <- rep(0:9, each=sn); Nsc <- nrow(sc)
gm <- colMeans(sc); cm <- do.call(rbind, lapply(0:9, function(k) colMeans(sc[ys==k,,drop=FALSE])))
SSB <- colSums((cm - matrix(gm,C,dB,byrow=TRUE))^2)*sn
SSW <- colSums(sc^2) - colSums(cm^2)*sn
Fst <- (SSB/(C-1))/pmax(SSW/(Nsc-C),1e-12)
active <- which(Fst>2 & colMeans(sc)>0.01)
cat("active pixels:", length(active), "  cores:", cores, "\n\n")

Bd <- as.integer(Sys.getenv("BDIAG","8"))
for (N_TR in c(100L,250L,500L,800L)) {
  t0 <- Sys.time()
  one <- function(b){
    set.seed(b*13+7)
    t(sapply(idx, function(ii){
      Xy <- X[sample(ii, N_TR), active, drop=FALSE]
      r  <- rodeo_select_class(Xy)
      c(nsel = sum(r$selected),
        spread = max(r$bw) - min(r$bw),
        nshr = sum(r$bw < max(r$bw) - 1e-9*max(r$bw)))
    }))
  }
  res <- mclapply(seq_len(Bd), one, mc.cores=cores)
  A <- Reduce(`+`, res)/length(res)                # 10 x 3 mean over splits
  cat(sprintf("N_TR=%4d (n_fit=%3d): gap-selected/digit=%.2f   bw_spread=%.4g   pixels_shrunk/digit=%.1f   [%.0fs]\n",
      N_TR, N_TR%/%2L, mean(A[,1]), mean(A[,2]), mean(A[,3]),
      as.numeric(Sys.time()-t0, units="secs")))
}
cat("\ndone\n")
