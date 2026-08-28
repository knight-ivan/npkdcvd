# Honest Anuran attribution figure (color): per-species selection frequency of the 22
# MFCC coefficients (Original, B=1000, corrected rodeo). Blue heatmap + grid + colorbar
# so the sparse structure (identity chiefly in MFCC-1) is clearly visible.
a  <- readRDS("../data/real_anuran_results.rds")
sf <- a$Original$selfreq                       # 7 x 22 (matches main text)
sp <- a$species
nsp <- nrow(sf); nco <- ncol(sf)
short <- c("Adenomera andrei","Adenomera hylaed.","Ameerega trivittata",
           "Hyla minuta","Hypsiboas cineras.","Hypsiboas cordobae","Leptodactylus fuscus")
if (length(short) != nsp) short <- sp
pal <- colorRampPalette(c("#f7fbff","#c6dbef","#6baed6","#2171b5","#08306b"))(100)

pdf("../figures/frog_selfreq.pdf", width=8.6, height=3.4)
layout(matrix(c(1,2), 1, 2), widths = c(9, 1.15))
## main heatmap
par(mar=c(4.0, 11, 1.1, 0.6))
image(1:nco, 1:nsp, t(sf[nsp:1, , drop=FALSE]), col=pal, zlim=c(0,1),
      axes=FALSE, xlab="", ylab="")
segments(0.5+0:nco, 0.5, 0.5+0:nco, nsp+0.5, col="grey85", lwd=0.5)   # vertical grid
segments(0.5, 0.5+0:nsp, nco+0.5, 0.5+0:nsp, col="grey85", lwd=0.5)   # horizontal grid
box(col="grey40")
mtext("MFCC coefficient index", side=1, line=2.5, cex=0.95)
axis(1, at=c(1,5,10,15,20), cex.axis=0.85)
axis(2, at=1:nsp, labels=rev(short), las=1, cex.axis=0.82)
for (i in 1:nsp) for (j in 1:nco) { v <- sf[i,j]
  if (v >= 0.05) text(j, nsp-i+1, sprintf("%.2f", v),
                      col=ifelse(v>0.45,"white","black"), cex=0.62) }
## colorbar
par(mar=c(4.0, 0.4, 1.1, 3.2))
zz <- seq(0, 1, length.out=100)
image(1, zz, matrix(zz, 1, 100), col=pal, axes=FALSE, xlab="", ylab="")
axis(4, at=seq(0,1,0.25), las=1, cex.axis=0.75)
mtext("selection frequency", side=4, line=2.1, cex=0.8)
box(col="grey40")
dev.off()
cat("wrote figures/frog_selfreq.pdf (color); MFCC-1:", round(sf[,1],2), "\n")
