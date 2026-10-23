# frog_bar_fig.R — Supplement figure for the Anuran analysis: per-species selection
# frequency of the first MFCC coefficient under the corrected rodeo-and-gap
# selector, after the pooled empirical PIT (primary; main-paper Algorithm 3) and on the raw
# MFCC scale (sensitivity).  Grey crosses mark the largest selection frequency among
# all other coefficients (primary scale).  Reads results_R2/real_anuran_results.rds
# written by real_anuran.R; writes figures/frog_bar.pdf.
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (length(here) && nzchar(here)) setwd(here)
a   <- readRDS(file.path("..", "results_R2", "real_anuran_results.rds"))
pit <- a$Original$selfreq_bc; raw <- a$Original$selfreq          # 7 x 22
short <- c("Adenomera andrei", "Adenomera hylaedactylus", "Ameerega trivittata",
           "Hyla minuta", "Hypsiboas cinerascens", "Hypsiboas cordobae",
           "Leptodactylus fuscus")
stopifnot(length(a$species) == length(short))
m_pit <- pit[, 1]; m_raw <- raw[, 1]
other_pit <- apply(pit[, -1, drop = FALSE], 1, max)
ord <- order(m_pit)                                   # ascending: largest on top
M <- rbind(raw = m_raw[ord], pit = m_pit[ord])        # barplot draws rows bottom-up per group

pdf(file.path("..", "figures", "frog_bar.pdf"), width = 7.6, height = 4.5)
par(mar = c(5.6, 12.5, 1.0, 1.2))
bp <- barplot(M, beside = TRUE, horiz = TRUE, names.arg = short[ord], las = 1,
              col = c("#9ecae1", "#08519c"), border = NA, xlim = c(0, 1.08),
              cex.names = 0.85, cex.axis = 0.85, space = c(0.1, 0.7))
mtext("selection frequency of the first MFCC coefficient (B = 1000)",
      side = 1, line = 2.5, cex = 0.95)
text(M + 0.012, bp, sprintf("%.2f", M), adj = 0, cex = 0.7)
points(other_pit[ord], bp[2, ], pch = 4, col = "#fd8d3c", cex = 0.95, lwd = 1.6)
legend("bottom", inset = c(0, -0.36), xpd = TRUE, ncol = 2,
       pch = c(15, 15, 4), col = c("#08519c", "#9ecae1", "#fd8d3c"),
       pt.lwd = c(1, 1, 1.6), bty = "n", cex = 0.74,
       legend = c("pooled-PIT selector (primary)", "raw scale (sensitivity)",
                  "largest other coefficient (primary)"))
invisible(dev.off())
cat("wrote figures/frog_bar.pdf\n")
cat("MFCC-1 pooled PIT:", round(m_pit, 2), "\nMFCC-1 raw:       ", round(m_raw, 2), "\n")
cat("largest other coefficient (PIT):", round(other_pit, 3), "\n")
