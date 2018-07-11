#This script produces the plots for the Gaussian simulations in the paper.


library(dplyr)
library(reshape2)
library(ggplot2)
library(vioplot)


load("res_paper/res_gaus_dscr.RData")


#########################################3
#violin plots
source("code_paper/sim/vioplot_col.R")

homo.data.smash = res[res$.id == "sp.3.v1" & res$method == "smash.s8" , ]
homo.data.smash.homo = res[res$.id == "sp.3.v1" & res$method == "smash.homo.s8" , ]
homo.data.tithresh = res[res$.id == "sp.3.v1" & res$method == "tithresh.homo.s8" , ]
homo.data.ebayes = res[res$.id == "sp.3.v1" & res$method == "ebayesthresh" , ]
homo.data.smash.true = res[res$.id == "sp.3.v1" & res$method == "smash.true.s8" , ]

homo.data = res[res$.id == "sp.3.v1" & (res$method == "smash.s8" | res$method == "ebayesthresh" | res$method == "tithresh.homo.s8"), ]


#horizontal
pdf("paper/violin_gaus_homo.pdf", height = 8, width = 8)
par(cex.axis = 1.5, cex.lab = 1.5, cex.sub = 1.5, mar = c(6.1, 12.1, 4.1, 2.1), mgp = c(3, 1.5, 0))
vioplot.col(homo.data.smash$mise, homo.data.smash.homo$mise, homo.data.tithresh$mise, homo.data.ebayes$mise, homo.data.smash.true$mise, horizontal = TRUE, ylim = c(5, 18), col = c("magenta", "gold", "gold", "gold", "purple"))
axis(2, 1:5, labels = c("SMASH", "SMASH\nhomo", "TI-thresholding\nhomo", "EbayesThresh\nhomo", "SMASH\ntrue variance"), las = 2)
title(xlab = "MISE", line = 3.5)
abline(v = median(homo.data.smash$mise), lty = 3, col = 3)
abline(h = 4.5, lty = 3)
dev.off()


cblocks.fn = function(t, type) {
  pos = c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 0.76, 0.78, 0.81)
  hgt = 2.88/5 * c(4, (-5), 3, (-4), 5, (-4.2), 2.1, 4.3, (-3.1), 2.1, (-4.2))
  fn = rep(0, length(t))
  for (j in 1:length(pos)) {
    fn = fn + (1 + sign(t - pos[j])) * (hgt[j]/2)
  }
  fn[fn < 0] = 0
  if (type == "mean") {
    return(NULL)
  } else if (type == "var") {
    return(0.01 + 1 * (fn - min(fn))/max(fn))
  }
}

doppler.fn = function(t, type) {
  dop.f = function(x) sqrt(x * (1 - x)) * sin((2 * pi * 1.05)/(x + 0.05))
  fn = dop.f(t)
  if (type == "mean") {
    fn = 3/(max(fn) - min(fn)) * (fn - min(fn))
    return((1 + fn)/5)
  } else if (type == "var") {
    fn = 10 * fn
    fn = fn - min(fn)
    return(1e-05 + 2 * fn)
  }
}

cor.fn = function(t, type) {
  fn = 623.87 * t^3 * (1 - 2 * t) * (t >= 0 & t <= 0.5) + 187.161 * (0.125 - t^3) * t^4 * (t > 0.5 & t <= 0.8) + 3708.470441 * (t - 1)^3 * (t > 0.8 & t <= 1)
  fn = (0.6/(max(fn) - min(fn))) * fn
  if (type == "mean") {
    return(fn - min(fn) + 0.2)
  } else if (type == "var") {
    return(NULL)
  }
}

spikes.fn = function(t, type) {
  spike.f = function(x) (0.75 * exp(-500 * (x - 0.23)^2) + 1.5 * exp(-2000 * (x - 0.33)^2) + 3 * exp(-8000 * (x - 0.47)^2) + 2.25 * exp(-16000 * (x - 0.69)^2) + 0.5 * exp(-32000 * (x - 0.83)^2))
  fn = spike.f(t)
  if (type == "mean") {
    return((1 + fn)/5)
  } else if (type == "var") {
    return(NULL)
  }
}

t = (1:1024)/1024

hetero.data.smash = res[res$.id == "sp.3.v5" & res$method == "smash.s8" , ]
hetero.data.smash.homo = res[res$.id == "sp.3.v5" & res$method == "smash.homo.s8" , ]
hetero.data.tithresh.homo = res[res$.id == "sp.3.v5" & res$method == "tithresh.homo.s8" , ]
hetero.data.tithresh.rmad = res[res$.id == "sp.3.v5" & res$method == "tithresh.rmad.s8" , ]
hetero.data.tithresh.smash = res[res$.id == "sp.3.v5" & res$method == "tithresh.smash.s8" , ]
hetero.data.tithresh.true = res[res$.id == "sp.3.v5" & res$method == "tithresh.true.s8" , ]
hetero.data.ebayes = res[res$.id == "sp.3.v5" & res$method == "ebayesthresh" , ]
hetero.data.smash.true = res[res$.id == "sp.3.v5" & res$method == "smash.true.s8" , ]

vioplot(hetero.data.smash$mise, hetero.data.tithresh.rmad$mise, hetero.data.tithresh.smash$mise, hetero.data.smash.homo$mise, hetero.data.tithresh.homo$mise,
        hetero.data.ebayes$mise, hetero.data.smash.true$mise, hetero.data.tithresh.true$mise)

# pdf("paper/violin_gaus_hetero_ti.pdf", height = 8, width = 12)
# par(cex.axis = 1.7, cex.lab = 2, cex.sub = 2, mar = c(8.1, 6.1, 6.1, 2.1), mgp = c(3, 1.5, 0))
# vioplot.col(hetero.data.smash$mise, hetero.data.tithresh.rmad$mise, hetero.data.tithresh.smash$mise, hetero.data.tithresh.true$mise, ylim = c(10, 80), names=NULL, col = c("magenta", rep("red", 3)))
# axis(side = 1, at = 1:4, labels = c("SMASH \n", "TI-thresh \n RMAD variance", "TI-thresh \n SMASH variance", "TI-thresh \n true variance"), mgp = c(3, 3, 0))
# title(ylab = "MISE", line = 4)
# title(xlab = "Method", line = 6)
# abline(h = median(hetero.data.smash$mise), lty = 3, col = 3)
# dev.off()
# 
# pdf("paper/violin_gaus_hetero_homo.pdf", height = 8, width = 12)
# par(cex.axis = 2, cex.lab = 2, cex.sub = 2, mar = c(8.1, 6.1, 6.1, 2.1), mgp = c(3, 1.5, 0))
# #vioplot.col(hetero.data.smash$mise, hetero.data.smash.homo$mise, hetero.data.tithresh.homo$mise,hetero.data.ebayes$mise, ylim = c(10, 80), names = c("SMASH", "SMASH, homo", "TI-thresh, homo", "EbayesThresh, homo"), col = c("magenta", rep("gold", 3)))
# vioplot.col(hetero.data.smash$mise, hetero.data.smash.homo$mise,hetero.data.ebayes$mise, ylim = c(10, 80), col = c("magenta", rep("gold", 2)))
# axis(side = 1, at = 1:3, labels = c("SMASH \n", "SMASH \n homo", "EbayesThresh \n homo"), mgp = c(3, 3, 0))
# title(ylab = "MISE", line = 4)
# title(xlab = "Method", line = 6)
# abline(h = median(hetero.data.smash$mise), lty = 3, col = 3)
# dev.off()
# 
# pdf("paper/violin_gaus_hetero_smash.pdf", height = 8, width = 12)
# par(cex.axis = 2, cex.lab = 2, cex.sub = 2, mar = c(5.1, 5.1, 6.1, 2.1), mgp = c(3, 1.3, 0))
# vioplot.col(hetero.data.smash$mise, hetero.data.smash.true$mise, ylim = c(10, 80), names = c("SMASH", "SMASH, true variance"), col = c("magenta", "purple"))
# title(ylab = "MISE")
# abline(h = median(hetero.data.smash$mise), lty = 3, col = 3)
# dev.off()


mu = spikes.fn(t, "mean")
sigma.ini = sqrt(cblocks.fn(t, "var"))
sd.fn = sigma.ini/mean(sigma.ini) * sd(mu)/3

# #vertical
# pdf("paper/violin_gaus_hetero_1.pdf", height = 8, width = 12)
# par(cex.axis = 0.9, cex.lab = 1.5, cex.sub = 1.5, mar = c(8.1, 6.1, 6.1, 2.1), mgp = c(3, 1.5, 0))
# vioplot.col(hetero.data.smash$mise, hetero.data.tithresh.rmad$mise, hetero.data.tithresh.smash$mise, hetero.data.tithresh.true$mise, hetero.data.smash.homo$mise, hetero.data.ebayes$mise, hetero.data.smash.true$mise, ylim = c(10, 80), names=NULL, col = c("magenta", rep("red", 3), rep("gold", 2), "purple"))
# axis(side = 1, at = 1:7, labels = c("SMASH \n", "TI-thresh \n RMAD variance", "TI-thresh \n SMASH variance", "TI-thresh \n true variance", "SMASH \n homo", "EbayesThresh \n homo", "SMASH \n true variance"), mgp = c(3, 3, 0))
# title(ylab = "MISE", line = 4)
# title(xlab = "Method", line = 6)
# abline(h = median(hetero.data.smash$mise), lty = 3, col = 3)
# abline(v = 4.5, lty = 3)
# abline(v = 6.5, lty = 3)
# dev.off()

#horizontal
pdf("paper/violin_gaus_hetero_1.pdf", height = 8, width = 8)
par(cex.axis = 1.5, cex.lab = 1.5, cex.sub = 1.5, mar = c(6.1, 11.1, 4.1, 0.5), mgp = c(3, 1.5, 0))
vioplot.col(hetero.data.smash$mise, hetero.data.tithresh.rmad$mise, hetero.data.tithresh.smash$mise, hetero.data.tithresh.true$mise, hetero.data.smash.homo$mise, hetero.data.ebayes$mise, hetero.data.smash.true$mise, ylim = c(10, 80), horizontal = TRUE, names=NULL, col = c("magenta", rep("red", 3), rep("gold", 2), "purple"))
par(lheight = 0.8)
axis(side = 2, at = 1:7, labels = c("SMASH", "TI-thresh\nRMAD variance", "TI-thresh\nSMASH variance", "TI-thresh\ntrue variance", "SMASH\nhomo", "EbayesThresh\nhomo", "SMASH\ntrue variance"), las = 2)
title(xlab = "MISE", line = 4)
abline(v = median(hetero.data.smash$mise), lty = 3, col = 3)
abline(h = 4.5, lty = 3)
abline(h = 6.5, lty = 3)
dev.off()


pdf("paper/violin_gaus_hetero_sd_1.pdf", height = 8, width = 8)
par(cex.axis = 1.5, cex.lab = 1.5, cex.sub = 1.5, mar = c(6.1, 3.1, 4.1, 2.1))
plot(mu, type = 'l', ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
lines(mu + 2* sd.fn, col = 2, lty = 5)
lines(mu - 2* sd.fn, col = 2, lty = 5)
axis(1, labels = FALSE, tick = FALSE)
axis(2)
dev.off()


hetero.data.smash.2 = res[res$.id == "cor.3.v3" & res$method == "smash.s8" , ]
hetero.data.smash.homo.2 = res[res$.id == "cor.3.v3" & res$method == "smash.homo.s8" , ]
hetero.data.tithresh.homo.2 = res[res$.id == "cor.3.v3" & res$method == "tithresh.homo.s8" , ]
hetero.data.tithresh.rmad.2 = res[res$.id == "cor.3.v3" & res$method == "tithresh.rmad.s8" , ]
hetero.data.tithresh.smash.2 = res[res$.id == "cor.3.v3" & res$method == "tithresh.smash.s8" , ]
hetero.data.tithresh.true.2 = res[res$.id == "cor.3.v3" & res$method == "tithresh.true.s8" , ]
hetero.data.ebayes.2 = res[res$.id == "cor.3.v3" & res$method == "ebayesthresh" , ]
hetero.data.smash.true.2 = res[res$.id == "cor.3.v3" & res$method == "smash.true.s8" , ]

# hetero.data.smash.2 = res[res$.id == "blk.3.v4" & res$method == "smash.s8" , ]
# hetero.data.smash.homo.2 = res[res$.id == "blk.3.v4" & res$method == "smash.homo.s8" , ]
# hetero.data.tithresh.homo.2 = res[res$.id == "blk.3.v4" & res$method == "tithresh.homo.s8" , ]
# hetero.data.tithresh.rmad.2 = res[res$.id == "blk.3.v4" & res$method == "tithresh.rmad.s8" , ]
# hetero.data.tithresh.smash.2 = res[res$.id == "blk.3.v4" & res$method == "tithresh.smash.s8" , ]
# hetero.data.tithresh.true.2 = res[res$.id == "blk.3.v4" & res$method == "tithresh.true.s8" , ]
# hetero.data.ebayes.2 = res[res$.id == "blk.3.v4" & res$method == "ebayesthresh" , ]
# hetero.data.smash.true.2 = res[res$.id == "blk.3.v4" & res$method == "smash.true.s8" , ]

# pdf("paper/violin_gaus_hetero_ti_2.pdf", height = 8, width = 12)
# par(cex.axis = 1.7, cex.lab = 2, cex.sub = 2, mar = c(8.1, 6.1, 4.1, 2.1), mgp = c(3, 1.5, 0))
# vioplot.col(hetero.data.smash.2$mise, hetero.data.tithresh.rmad.2$mise, hetero.data.tithresh.smash.2$mise, hetero.data.tithresh.true.2$mise, ylim = c(0.5,5), col = c("magenta", rep("red", 3)))
# #vioplot.col(hetero.data.smash.2$mise, hetero.data.tithresh.rmad.2$mise, hetero.data.tithresh.smash.2$mise, hetero.data.tithresh.true.2$mise, ylim = c(50, 200), col = c("magenta", rep("red", 3)))
# axis(side = 1, at = 1:4, labels = c("SMASH \n", "TI-thresh \n RMAD variance", "TI-thresh \n SMASH variance", "TI-thresh \n true variance"), mgp = c(3, 3, 0))
# title(ylab = "MISE", line = 4)
# title(xlab = "Method", line = 6)
# abline(h = median(hetero.data.smash.2$mise), lty = 3, col = 3)
# dev.off()
# 
# pdf("paper/violin_gaus_hetero_homo_2.pdf", height = 8, width = 12)
# par(cex.axis = 2, cex.lab = 2, cex.sub = 2, mar = c(8.1, 6.1, 4.1, 2.1), mgp = c(3, 1.5, 0))
# #vioplot.col(hetero.data.smash$mise, hetero.data.smash.homo$mise, hetero.data.tithresh.homo$mise,hetero.data.ebayes$mise, ylim = c(10, 80), names = c("SMASH", "SMASH, homo", "TI-thresh, homo", "EbayesThresh, homo"), col = c("magenta", rep("gold", 3)))
# vioplot.col(hetero.data.smash.2$mise, hetero.data.smash.homo.2$mise,hetero.data.ebayes.2$mise, ylim = c(0.5, 5), col = c("magenta", rep("gold", 2)))
# #vioplot.col(hetero.data.smash.2$mise, hetero.data.smash.homo.2$mise,hetero.data.ebayes.2$mise, ylim = c(50, 200), col = c("magenta", rep("gold", 2)))
# axis(side = 1, at = 1:3, labels = c("SMASH \n", "SMASH \n homo", "EbayesThresh \n homo"), mgp = c(3, 3, 0))
# title(ylab = "MISE", line = 4)
# title(xlab = "Method", line = 6)
# abline(h = median(hetero.data.smash.2$mise), lty = 3, col = 3)
# dev.off()
# 
# pdf("paper/violin_gaus_hetero_smash_2.pdf", height = 8, width = 12)
# par(cex.axis = 2, cex.lab = 2, cex.sub = 2, mar = c(5.1, 5.1, 4.1, 2.1), mgp = c(3, 1.3, 0))
# vioplot.col(hetero.data.smash.2$mise, hetero.data.smash.true.2$mise, ylim = c(50, 200), names = c("SMASH", "SMASH, true variance"), col = c("magenta", "purple"))
# title(ylab = "MISE")
# abline(h = median(hetero.data.smash.2$mise), lty = 3, col = 3)
# dev.off()

mu = cor.fn(t, "mean")
sigma.ini = sqrt(doppler.fn(t, "var"))
sd.fn.2 = sigma.ini/mean(sigma.ini) * sd(mu)/3

# #vertical
# pdf("paper/violin_gaus_hetero_2.pdf", height = 8, width = 12)
# par(cex.axis = 0.9, cex.lab = 1.5, cex.sub = 1.5, mar = c(8.1, 6.1, 4.1, 2.1), mgp = c(3, 1.5, 0))
# vioplot.col(hetero.data.smash.2$mise, hetero.data.tithresh.rmad.2$mise, hetero.data.tithresh.smash.2$mise, hetero.data.tithresh.true.2$mise, hetero.data.smash.homo.2$mise, hetero.data.ebayes.2$mise, hetero.data.smash.true.2$mise, ylim = c(0.5,5), col = c("magenta", rep("red", 3), rep("gold", 2), "purple"))
# #vioplot.col(hetero.data.smash.2$mise, hetero.data.tithresh.rmad.2$mise, hetero.data.tithresh.smash.2$mise, hetero.data.tithresh.true.2$mise, ylim = c(50, 200), col = c("magenta", rep("red", 3)))
# axis(side = 1, at = 1:7, labels = c("SMASH \n", "TI-thresh \n RMAD variance", "TI-thresh \n SMASH variance", "TI-thresh \n true variance", "SMASH \n homo", "EbayesThresh \n homo", "SMASH \n true variance"), mgp = c(3, 3, 0))
# title(ylab = "MISE", line = 4)
# title(xlab = "Method", line = 6)
# abline(h = median(hetero.data.smash.2$mise), lty = 3, col = 3)
# abline(v = 4.5, lty = 3)
# abline(v = 6.5, lty = 3)
# dev.off()

#horizontal
pdf("paper/violin_gaus_hetero_2.pdf", height = 8, width = 8)
par(cex.axis = 1.5, cex.lab = 1.5, cex.sub = 1.5, mar = c(6.1, 11.1, 4.1, 0.5), mgp = c(3, 1.5, 0))
vioplot.col(hetero.data.smash.2$mise, hetero.data.tithresh.rmad.2$mise, hetero.data.tithresh.smash.2$mise, hetero.data.tithresh.true.2$mise, hetero.data.smash.homo.2$mise, hetero.data.ebayes.2$mise, hetero.data.smash.true.2$mise, horizontal = TRUE, ylim = c(0.5,5), col = c("magenta", rep("red", 3), rep("gold", 2), "purple"))
par(lheight = 0.8)
axis(side = 2, at = 1:7, labels = c("SMASH", "TI-thresh\nRMAD variance", "TI-thresh\nSMASH variance", "TI-thresh\ntrue variance", "SMASH\nhomo", "EbayesThresh\nhomo", "SMASH\ntrue variance"), las = 2)
title(xlab = "MISE", line = 4)
abline(v = median(hetero.data.smash.2$mise), lty = 3, col = 3)
abline(h = 4.5, lty = 3)
abline(h = 6.5, lty = 3)
dev.off()


pdf("paper/violin_gaus_hetero_sd_2.pdf", height = 8, width = 8)
par(cex.axis = 1.5, cex.lab = 1.5, cex.sub = 1.5, mar = c(6.1, 3.1, 4.1, 2.1))
plot(mu, type = 'l', ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
lines(mu + 2* sd.fn.2, col = 2, lty = 5)
lines(mu - 2* sd.fn.2, col = 2, lty = 5)
axis(1, labels = FALSE, tick = FALSE)
axis(2)
dev.off()
