# Does the corrected rodeo find concentrated per-pixel structure on UCI 8x8
# digits (d=64), or is it flat like MNIST? Decides the fate of supp figs S.5-S.6.
source("rodeo_core.R"); suppressMessages(library(parallel))
cores <- as.integer(Sys.getenv("RODEO_CORES","12"))
nA <- 5620L; dA <- 64L; C <- 10L
N_TR <- as.integer(Sys.getenv("NTR","400"))
B    <- as.integer(Sys.getenv("RODEO_B","150"))
SCR  <- "."

X <- matrix(readBin("/tmp/digits8_X.bin","numeric",nA*dA,size=4,endian="little"),
            nrow=nA, ncol=dA, byrow=TRUE)/16
y <- readBin("/tmp/digits8_y.bin","integer",nA,size=4,endian="little")
idx <- lapply(0:9, function(k) which(y==k))
cat(sprintf("UCI 8x8: N_TR=%d, B=%d, d=%d, class sizes %s\n",
            N_TR, B, dA, paste(sapply(idx,length),collapse=" ")))

one <- function(b){
  set.seed(b*2711+5)
  do.call(rbind, lapply(idx, function(ii){
    Xy <- X[sample(ii, N_TR), , drop=FALSE]
    r  <- rodeo_select_class(Xy)
    list(score=r$score, nsel=sum(r$selected))
  }))
}
# collect
sc_acc <- matrix(0.0, C, dA); nsel_acc <- numeric(C); nb <- 0L
res <- mclapply(seq_len(B), function(b){
  set.seed(b*2711+5)
  s <- matrix(0.0, C, dA); ns <- numeric(C)
  for (k in 1:C){ r <- rodeo_select_class(X[sample(idx[[k]],N_TR),,drop=FALSE]); s[k,]<-r$score; ns[k]<-sum(r$selected) }
  list(s=s, ns=ns)
}, mc.cores=cores)
for (rr in res){ sc_acc <- sc_acc + rr$s; nsel_acc <- nsel_acc + rr$ns; nb <- nb+1L }
score <- sc_acc/nb; nsel <- nsel_acc/nb
saveRDS(list(score=score, nsel=nsel, N_TR=N_TR, B=B), file.path(SCR,"uci8x8_attr.rds"))

pdf(file.path(SCR,"uci8x8_attr.pdf"), width=10, height=4.5)
op <- par(mfrow=c(2,5), mar=c(1,1,2,1))
for (k in 1:C){ v<-score[k,]; mx<-max(v); if(mx<=0)mx<-1
  image(matrix(v/mx,8,8,byrow=TRUE)[,8:1], axes=FALSE, col=grey.colors(64,0,1), main=paste("digit",k-1)) }
par(op); dev.off()

cat("score range:", sprintf("%.3f .. %.3f", min(score), max(score)), "\n")
cat("gap-selected per digit:\n"); print(round(nsel,2))
topfrac <- sapply(1:C, function(k){ s<-sort(score[k,],decreasing=TRUE); sum(s[1:8])/sum(s) })
cat("top-8-of-64 score share per digit (uniform baseline = 0.125):\n"); print(round(topfrac,3))
cat("\ndone\n")
