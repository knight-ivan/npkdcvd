# Alternative Anuran figure: horizontal bar chart of the first-MFCC selection
# frequency per species (the coefficient that carries each species' identity).
# Much more legible than a mostly-empty heatmap. Writes figures/frog_bar.pdf.
a  <- readRDS("../data/real_anuran_results.rds")
sf <- a$Original$selfreq                       # 7 x 22
m1 <- sf[,1]                                    # MFCC-1 selection frequency
other_max <- apply(sf[,-1,drop=FALSE], 1, max)  # largest among coefficients 2..22
short <- c("Adenomera andrei","Adenomera hylaedactylus","Ameerega trivittata",
           "Hyla minuta","Hypsiboas cinerascens","Hypsiboas cordobae","Leptodactylus fuscus")
ord <- order(m1)                                # ascending so largest on top
m1o <- m1[ord]; nmo <- short[ord]; omo <- other_max[ord]

pdf("../figures/frog_bar.pdf", width=7.6, height=3.6)
par(mar=c(4.0, 12.5, 1.0, 1.2))
bp <- barplot(m1o, horiz=TRUE, names.arg=nmo, las=1, col="#2171b5", border=NA,
              xlim=c(0,0.8), xlab="", cex.names=0.85, cex.axis=0.85)
mtext("selection frequency of the first MFCC coefficient", side=1, line=2.5, cex=0.95)
# value labels
text(m1o + 0.015, bp, sprintf("%.2f", m1o), adj=0, cex=0.78, col="black")
# faint markers for the largest OTHER coefficient per species (shows others are ~0)
points(omo, bp, pch=4, col="grey55", cex=0.9)
legend("bottomright", pch=c(15,4), col=c("#2171b5","grey55"), bty="n", cex=0.8,
       legend=c("first MFCC", "largest other coefficient"))
dev.off()
cat("wrote figures/frog_bar.pdf\n")
cat("MFCC-1:", round(m1,2), "\n")
cat("max other coef per species:", round(other_max,3), "\n")
