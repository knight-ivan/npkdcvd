# Render the MNIST per-digit map from the rodeo SCORE T_j (peak squared
# standardized derivative) -- the honest continuous profile of the CORRECTED
# statistic -- at larger N_TR. Is it stroke-like, or just a few scattered pixels?
source("rodeo_core.R"); suppressMessages(library(parallel))
cores <- as.integer(Sys.getenv("RODEO_CORES","12"))
dB <- 784L; nB <- 70000L; C <- 10L
N_TR <- as.integer(Sys.getenv("NTR","500"))
B    <- as.integer(Sys.getenv("RODEO_B","120"))
SCR  <- "."

X <- matrix(readBin("/tmp/mnist_X.bin","numeric",nB*dB,size=4,endian="little"),
            nrow=nB, ncol=dB, byrow=TRUE)/255
y <- readBin("/tmp/mnist_y.bin","integer",nB,size=4,endian="little")
idx <- lapply(0:9, function(k) which(y==k))
set.seed(1); sn <- 300L
sc <- do.call(rbind, lapply(idx, function(ii) X[sample(ii,sn),]))
ys <- rep(0:9, each=sn); Nsc <- nrow(sc)
gm <- colMeans(sc); cm <- do.call(rbind, lapply(0:9, function(k) colMeans(sc[ys==k,,drop=FALSE])))
SSB <- colSums((cm - matrix(gm,C,dB,byrow=TRUE))^2)*sn
SSW <- colSums(sc^2) - colSums(cm^2)*sn
Fst <- (SSB/(C-1))/pmax(SSW/(Nsc-C),1e-12)
active <- which(Fst>2 & colMeans(sc)>0.01)
cat(sprintf("MNIST score-map: N_TR=%d, B=%d, active=%d, cores=%d\n", N_TR, B, length(active), cores))

one <- function(b){
  set.seed(b*2711+5)
  do.call(rbind, lapply(idx, function(ii){
    Xy <- X[sample(ii, N_TR), active, drop=FALSE]
    rodeo_select_class(Xy)$score            # per-active-pixel score T_j
  }))                                        # 10 x |active|
}
res   <- mclapply(seq_len(B), one, mc.cores=cores)
score <- Reduce(`+`, res)/length(res)        # 10 x |active| mean score
full  <- matrix(0.0, C, dB); full[, active] <- score
saveRDS(list(score_full=full, active=active, N_TR=N_TR, B=B), file.path(SCR,"mnist_scoremap.rds"))

pdf(file.path(SCR,"mnist_scoremap.pdf"), width=10, height=4.5)
op <- par(mfrow=c(2,5), mar=c(1,1,2,1))
for (k in 1:C) {
  v <- full[k,]; mx <- max(v[active]); if (mx<=0) mx <- 1
  img <- matrix(v/mx, 28, 28)[, 28:1]
  image(img, axes=FALSE, col=grey.colors(64, start=0, end=1), main=paste("digit", k-1))
}
par(op); dev.off()

cat("score range:", sprintf("%.3f .. %.3f", min(score), max(score)), "\n")
cat("mean active-pixel score per digit:\n"); print(round(rowMeans(score),3))
# concentration: what fraction of each digit's total score sits in its top-20 pixels?
topfrac <- sapply(1:C, function(k){ s<-sort(full[k,active],decreasing=TRUE); sum(s[1:20])/sum(s) })
cat("top-20-pixel score share per digit (1=very concentrated):\n"); print(round(topfrac,3))
cat("\ndone\n")
