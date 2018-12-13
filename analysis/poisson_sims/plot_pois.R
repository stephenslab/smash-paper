#This script produces the plots for the Poisson simulations in the paper

library(vioplot)
source("code_paper/sim/vioplot_col.R")
load("res_paper/res_pois.RData")

mise.hf.ti.r.bur.1 = colMeans(matrix(c(mise.hf.ti.r.4.bur.1, mise.hf.ti.r.5.bur.1, mise.hf.ti.r.6.bur.1, mise.hf.ti.r.7.bur.1), byrow = TRUE, nr = 4))
mise.hf.ti.r.bur.8 = colMeans(matrix(c(mise.hf.ti.r.4.bur.8, mise.hf.ti.r.5.bur.8, mise.hf.ti.r.6.bur.8, mise.hf.ti.r.7.bur.8), byrow = TRUE, nr = 4))

#horizontal
par(cex.axis = 2, cex.lab = 2, cex.sub = 2, mar = c(6.1, 12.1, 4.1, 2.1), mgp = c(3, 1.5, 0))
vioplot.col(mise.ash.bur.1, mise.BMSM.bur.1, mise.hf.ti.r.bur.1, horizontal = TRUE, names = c("SMASH", "BMSM", "Haar-Fisz"), col = c("magenta", "red", "red"))
title(xlab = "MISE", line = 4)
abline(v = median(mise.ash.bur.1), lty = 3, col = 3)

par(cex.axis = 2, cex.lab = 2, cex.sub = 2, mar = c(6.1, 12.1, 4.1, 2.1), mgp = c(3, 1.5, 0))
vioplot.col(mise.ash.bur.8, mise.BMSM.bur.8, mise.hf.ti.r.bur.8, horizontal = TRUE, names = c("SMASH", "BMSM", "Haar-Fisz"), col = c("magenta", "red", "red"))
title(xlab = "MISE", line = 4)
abline(v = median(mise.ash.bur.8), lty = 3, col = 3)

#plot the bursts test function
n=1024
t=1:n/n

I_1 = exp(-(abs(t-0.2)/0.01)^1.2)*(t<=0.2) + exp(-(abs(t-0.2)/0.03)^1.2)*(t>0.2);
I_2 = exp(-(abs(t-0.3)/0.01)^1.2)*(t<=0.3) + exp(-(abs(t-0.3)/0.03)^1.2)*(t>0.3);
I_3 = exp(-(abs(t-0.4)/0.01)^1.2)*(t<=0.4) + exp(-(abs(t-0.4)/0.03)^1.2)*(t>0.4);
mu.bur = 2.99/4.51804*(4*I_1+3*I_2+4.5*I_3);

plot(mu.bur, type = 'l', xlab = "position", ylab = "intensity")


#########################################
#To obtain violin plots for all test functions/intensities (most figures not shown in paper)

mise.hf.ti.r.s.1 = colMeans(matrix(c(mise.hf.ti.r.4.s.1, mise.hf.ti.r.5.s.1, mise.hf.ti.r.6.s.1, mise.hf.ti.r.7.s.1), byrow = TRUE, nr = 4))
mise.hf.ti.r.s.8 = colMeans(matrix(c(mise.hf.ti.r.4.s.8, mise.hf.ti.r.5.s.8, mise.hf.ti.r.6.s.8, mise.hf.ti.r.7.s.8), byrow = TRUE, nr = 4))
mise.hf.ti.r.s.128 = colMeans(matrix(c(mise.hf.ti.r.4.s.128, mise.hf.ti.r.5.s.128, mise.hf.ti.r.6.s.128, mise.hf.ti.r.7.s.128), byrow = TRUE, nr = 4))


par(cex.axis = 2, cex.lab = 2, cex.sub = 2, mar = c(6.1, 12.1, 4.1, 2.1), mgp = c(3, 1.5, 0))
vioplot.col(mise.ash.s.1, mise.BMSM.s.1, mise.hf.ti.r.s.1, horizontal = TRUE, names = c("SMASH", "BMSM", "Haar-Fisz"), col = c("magenta", "red", "red"))
title(xlab = "MISE", line = 4)
abline(v = median(mise.ash.s.1), lty = 3, col = 3)

par(cex.axis = 2, cex.lab = 2, cex.sub = 2, mar = c(6.1, 12.1, 4.1, 2.1), mgp = c(3, 1.5, 0))
vioplot.col(mise.ash.s.8, mise.BMSM.s.8, mise.hf.ti.r.s.8, horizontal = TRUE, names = c("SMASH", "BMSM", "Haar-Fisz"), col = c("magenta", "red", "red"))
title(xlab = "MISE", line = 4)

par(cex.axis = 2, cex.lab = 2, cex.sub = 2, mar = c(6.1, 12.1, 4.1, 2.1), mgp = c(3, 1.5, 0))
vioplot.col(mise.ash.s.128, mise.BMSM.s.128, mise.hf.ti.r.s.128, horizontal = TRUE, names = c("SMASH", "BMSM", "Haar-Fisz"), col = c("magenta", "red", "red"))
title(xlab = "MISE", line = 4)


mise.hf.ti.r.ang.1 = colMeans(matrix(c(mise.hf.ti.r.4.ang.1, mise.hf.ti.r.5.ang.1, mise.hf.ti.r.6.ang.1, mise.hf.ti.r.7.ang.1), byrow = TRUE, nr = 4))
mise.hf.ti.r.ang.8 = colMeans(matrix(c(mise.hf.ti.r.4.ang.8, mise.hf.ti.r.5.ang.8, mise.hf.ti.r.6.ang.8, mise.hf.ti.r.7.ang.8), byrow = TRUE, nr = 4))
mise.hf.ti.r.ang.128 = colMeans(matrix(c(mise.hf.ti.r.4.ang.128, mise.hf.ti.r.5.ang.128, mise.hf.ti.r.6.ang.128, mise.hf.ti.r.7.ang.128), byrow = TRUE, nr = 4))


par(cex.axis = 2, cex.lab = 2, cex.sub = 2, mar = c(6.1, 12.1, 4.1, 2.1), mgp = c(3, 1.5, 0))
vioplot.col(mise.ash.ang.1, mise.BMSM.ang.1, mise.hf.ti.r.ang.1, horizontal = TRUE, names = c("SMASH", "BMSM", "Haar-Fisz"), col = c("magenta", "red", "red"))
title(xlab = "MISE", line = 4)
abline(v = median(mise.ash.ang.1), lty = 3, col = 3)

par(cex.axis = 2, cex.lab = 2, cex.sub = 2, mar = c(6.1, 12.1, 4.1, 2.1), mgp = c(3, 1.5, 0))
vioplot.col(mise.ash.ang.8, mise.BMSM.ang.8, mise.hf.ti.r.ang.8, horizontal = TRUE, names = c("SMASH", "BMSM", "Haar-Fisz"), col = c("magenta", "red", "red"))
title(xlab = "MISE", line = 4)

par(cex.axis = 2, cex.lab = 2, cex.sub = 2, mar = c(6.1, 12.1, 4.1, 2.1), mgp = c(3, 1.5, 0))
vioplot.col(mise.ash.ang.128, mise.BMSM.ang.128, mise.hf.ti.r.ang.128, horizontal = TRUE, names = c("SMASH", "BMSM", "Haar-Fisz"), col = c("magenta", "red", "red"))
title(xlab = "MISE", line = 4)


mise.hf.ti.r.bur.1 = colMeans(matrix(c(mise.hf.ti.r.4.bur.1, mise.hf.ti.r.5.bur.1, mise.hf.ti.r.6.bur.1, mise.hf.ti.r.7.bur.1), byrow = TRUE, nr = 4))
mise.hf.ti.r.bur.8 = colMeans(matrix(c(mise.hf.ti.r.4.bur.8, mise.hf.ti.r.5.bur.8, mise.hf.ti.r.6.bur.8, mise.hf.ti.r.7.bur.8), byrow = TRUE, nr = 4))
mise.hf.ti.r.bur.128 = colMeans(matrix(c(mise.hf.ti.r.4.bur.128, mise.hf.ti.r.5.bur.128, mise.hf.ti.r.6.bur.128, mise.hf.ti.r.7.bur.128), byrow = TRUE, nr = 4))


par(cex.axis = 2, cex.lab = 2, cex.sub = 2, mar = c(6.1, 12.1, 4.1, 2.1), mgp = c(3, 1.5, 0))
vioplot.col(mise.ash.bur.1, mise.BMSM.bur.1, mise.hf.ti.r.bur.1, horizontal = TRUE, names = c("SMASH", "BMSM", "Haar-Fisz"), col = c("magenta", "red", "red"))
title(xlab = "MISE", line = 4)
abline(v = median(mise.ash.bur.1), lty = 3, col = 3)

par(cex.axis = 2, cex.lab = 2, cex.sub = 2, mar = c(6.1, 12.1, 4.1, 2.1), mgp = c(3, 1.5, 0))
vioplot.col(mise.ash.bur.8, mise.BMSM.bur.8, mise.hf.ti.r.bur.8, horizontal = TRUE, names = c("SMASH", "BMSM", "Haar-Fisz"), col = c("magenta", "red", "red"))
title(xlab = "MISE", line = 4)

par(cex.axis = 2, cex.lab = 2, cex.sub = 2, mar = c(6.1, 12.1, 4.1, 2.1), mgp = c(3, 1.5, 0))
vioplot.col(mise.ash.bur.128, mise.BMSM.bur.128, mise.hf.ti.r.bur.128, horizontal = TRUE, names = c("SMASH", "BMSM", "Haar-Fisz"), col = c("magenta", "red", "red"))
title(xlab = "MISE", line = 4)


mise.hf.ti.r.cb.1 = colMeans(matrix(c(mise.hf.ti.r.4.cb.1, mise.hf.ti.r.5.cb.1, mise.hf.ti.r.6.cb.1, mise.hf.ti.r.7.cb.1), byrow = TRUE, nr = 4))
mise.hf.ti.r.cb.8 = colMeans(matrix(c(mise.hf.ti.r.4.cb.8, mise.hf.ti.r.5.cb.8, mise.hf.ti.r.6.cb.8, mise.hf.ti.r.7.cb.8), byrow = TRUE, nr = 4))
mise.hf.ti.r.cb.128 = colMeans(matrix(c(mise.hf.ti.r.4.cb.128, mise.hf.ti.r.5.cb.128, mise.hf.ti.r.6.cb.128, mise.hf.ti.r.7.cb.128), byrow = TRUE, nr = 4))


par(cex.axis = 2, cex.lab = 2, cex.sub = 2, mar = c(6.1, 12.1, 4.1, 2.1), mgp = c(3, 1.5, 0))
vioplot.col(mise.ash.cb.1, mise.BMSM.cb.1, mise.hf.ti.r.cb.1, horizontal = TRUE, names = c("SMASH", "BMSM", "Haar-Fisz"), col = c("magenta", "red", "red"))
title(xlab = "MISE", line = 4)
abline(v = median(mise.ash.cb.1), lty = 3, col = 3)

par(cex.axis = 2, cex.lab = 2, cex.sub = 2, mar = c(6.1, 12.1, 4.1, 2.1), mgp = c(3, 1.5, 0))
vioplot.col(mise.ash.cb.8, mise.BMSM.cb.8, mise.hf.ti.r.cb.8, horizontal = TRUE, names = c("SMASH", "BMSM", "Haar-Fisz"), col = c("magenta", "red", "red"))
title(xlab = "MISE", line = 4)

par(cex.axis = 2, cex.lab = 2, cex.sub = 2, mar = c(6.1, 12.1, 4.1, 2.1), mgp = c(3, 1.5, 0))
vioplot.col(mise.ash.cb.128, mise.BMSM.cb.128, mise.hf.ti.r.cb.128, horizontal = TRUE, names = c("SMASH", "BMSM", "Haar-Fisz"), col = c("magenta", "red", "red"))
title(xlab = "MISE", line = 4)


mise.hf.ti.r.b.1 = colMeans(matrix(c(mise.hf.ti.r.4.b.1, mise.hf.ti.r.5.b.1, mise.hf.ti.r.6.b.1, mise.hf.ti.r.7.b.1), byrow = TRUE, nr = 4))
mise.hf.ti.r.b.8 = colMeans(matrix(c(mise.hf.ti.r.4.b.8, mise.hf.ti.r.5.b.8, mise.hf.ti.r.6.b.8, mise.hf.ti.r.7.b.8), byrow = TRUE, nr = 4))
mise.hf.ti.r.b.128 = colMeans(matrix(c(mise.hf.ti.r.4.b.128, mise.hf.ti.r.5.b.128, mise.hf.ti.r.6.b.128, mise.hf.ti.r.7.b.128), byrow = TRUE, nr = 4))


par(cex.axis = 2, cex.lab = 2, cex.sub = 2, mar = c(6.1, 12.1, 4.1, 2.1), mgp = c(3, 1.5, 0))
vioplot.col(mise.ash.b.1, mise.BMSM.b.1, mise.hf.ti.r.b.1, horizontal = TRUE, names = c("SMASH", "BMSM", "Haar-Fisz"), col = c("magenta", "red", "red"))
title(xlab = "MISE", line = 4)
abline(v = median(mise.ash.b.1), lty = 3, col = 3)

par(cex.axis = 2, cex.lab = 2, cex.sub = 2, mar = c(6.1, 12.1, 4.1, 2.1), mgp = c(3, 1.5, 0))
vioplot.col(mise.ash.b.8, mise.BMSM.b.8, mise.hf.ti.r.b.8, horizontal = TRUE, names = c("SMASH", "BMSM", "Haar-Fisz"), col = c("magenta", "red", "red"))
title(xlab = "MISE", line = 4)

par(cex.axis = 2, cex.lab = 2, cex.sub = 2, mar = c(6.1, 12.1, 4.1, 2.1), mgp = c(3, 1.5, 0))
vioplot.col(mise.ash.b.128, mise.BMSM.b.128, mise.hf.ti.r.b.128, horizontal = TRUE, names = c("SMASH", "BMSM", "Haar-Fisz"), col = c("magenta", "red", "red"))
title(xlab = "MISE", line = 4)


mise.hf.ti.r.hs.1 = colMeans(matrix(c(mise.hf.ti.r.4.hs.1, mise.hf.ti.r.5.hs.1, mise.hf.ti.r.6.hs.1, mise.hf.ti.r.7.hs.1), byrow = TRUE, nr = 4))
mise.hf.ti.r.hs.8 = colMeans(matrix(c(mise.hf.ti.r.4.hs.8, mise.hf.ti.r.5.hs.8, mise.hf.ti.r.6.hs.8, mise.hf.ti.r.7.hs.8), byrow = TRUE, nr = 4))
mise.hf.ti.r.hs.128 = colMeans(matrix(c(mise.hf.ti.r.4.hs.128, mise.hf.ti.r.5.hs.128, mise.hf.ti.r.6.hs.128, mise.hf.ti.r.7.hs.128), byrow = TRUE, nr = 4))


par(cex.axis = 2, cex.lab = 2, cex.sub = 2, mar = c(6.1, 12.1, 4.1, 2.1), mgp = c(3, 1.5, 0))
vioplot.col(mise.ash.hs.1, mise.BMSM.hs.1, mise.hf.ti.r.hs.1, horizontal = TRUE, names = c("SMASH", "BMSM", "Haar-Fisz"), col = c("magenta", "red", "red"))
title(xlab = "MISE", line = 4)
abline(v = median(mise.ash.hs.1), lty = 3, col = 3)

par(cex.axis = 2, cex.lab = 2, cex.sub = 2, mar = c(6.1, 12.1, 4.1, 2.1), mgp = c(3, 1.5, 0))
vioplot.col(mise.ash.hs.8, mise.BMSM.hs.8, mise.hf.ti.r.hs.8, horizontal = TRUE, names = c("SMASH", "BMSM", "Haar-Fisz"), col = c("magenta", "red", "red"))
title(xlab = "MISE", line = 4)

par(cex.axis = 2, cex.lab = 2, cex.sub = 2, mar = c(6.1, 12.1, 4.1, 2.1), mgp = c(3, 1.5, 0))
vioplot.col(mise.ash.hs.128, mise.BMSM.hs.128, mise.hf.ti.r.hs.128, horizontal = TRUE, names = c("SMASH", "BMSM", "Haar-Fisz"), col = c("magenta", "red", "red"))
title(xlab = "MISE", line = 4)


