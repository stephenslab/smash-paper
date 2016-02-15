#This script produces the plots for the Poisson simulations in the paper

library(vioplot)
source("code_paper/sim/vioplot_col.R")
load("res_paper/res_pois.RData")


mise.hf.ti.r.bur.1 = colMeans(matrix(c(mise.hf.ti.r.4.bur.1, mise.hf.ti.r.5.bur.1, mise.hf.ti.r.6.bur.1, mise.hf.ti.r.7.bur.1), byrow = TRUE, nr = 4))
mise.hf.ti.r.bur.8 = colMeans(matrix(c(mise.hf.ti.r.4.bur.8, mise.hf.ti.r.5.bur.8, mise.hf.ti.r.6.bur.8, mise.hf.ti.r.7.bur.8), byrow = TRUE, nr = 4))

# #vertical
# pdf("paper/violin_pois_1.pdf", height = 8, width = 12)
# par(cex.axis = 2, cex.lab = 2, cex.sub = 2, mar = c(6.1, 6.1, 4.1, 2.1), mgp = c(3, 1.5, 0))
# vioplot.col(mise.ash.bur.1, mise.BMSM.bur.1, mise.hf.ti.r.bur.1, names = c("SMASH", "BMSM", "Haar-Fisz"), col = c("magenta", "red", "red"))
# title(xlab = "Method", ylab = "MISE", line = 4)
# abline(h = median(mise.ash.bur.1), lty = 3, col = 3)
# dev.off()
# 
# pdf("paper/violin_pois_8.pdf", height = 8, width = 12)
# par(cex.axis = 2, cex.lab = 2, cex.sub = 2, mar = c(6.1, 6.1, 4.1, 2.1), mgp = c(3, 1.5, 0))
# vioplot.col(mise.ash.bur.8, mise.BMSM.bur.8, mise.hf.ti.r.bur.8, names = c("SMASH", "BMSM", "Haar-Fisz"), col = c("magenta", "red", "red"))
# title(xlab = "Method", ylab = "MISE", line = 4)
# abline(h = median(mise.ash.bur.8), lty = 3, col = 3)
# dev.off()

#horizontal
pdf("paper/violin_pois_1.pdf", height = 8, width = 8)
par(cex.axis = 2, cex.lab = 2, cex.sub = 2, mar = c(6.1, 12.1, 4.1, 2.1), mgp = c(3, 1.5, 0))
vioplot.col(mise.ash.bur.1, mise.BMSM.bur.1, mise.hf.ti.r.bur.1, horizontal = TRUE, names = c("SMASH", "BMSM", "Haar-Fisz"), col = c("magenta", "red", "red"))
title(xlab = "MISE", line = 4)
abline(v = median(mise.ash.bur.1), lty = 3, col = 3)
dev.off()

pdf("paper/violin_pois_8.pdf", height = 8, width = 8)
par(cex.axis = 2, cex.lab = 2, cex.sub = 2, mar = c(6.1, 12.1, 4.1, 2.1), mgp = c(3, 1.5, 0))
vioplot.col(mise.ash.bur.8, mise.BMSM.bur.8, mise.hf.ti.r.bur.8, horizontal = TRUE, names = c("SMASH", "BMSM", "Haar-Fisz"), col = c("magenta", "red", "red"))
title(xlab = "MISE", line = 4)
abline(v = median(mise.ash.bur.8), lty = 3, col = 3)
dev.off()