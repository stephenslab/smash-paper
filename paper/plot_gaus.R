library(dplyr)
library(reshape2)
library(ggplot2)
library(vioplot)


# 
# ##barplots in main body of paper
# mean.table = melt(table.1.v1)
# names(mean.table) = c("method", "testfunction", "mse")
# 
# pdf("paper/smash_gaus_1_v1.pdf",width = 10, height = 5)
# ggplot(mean.table %>% filter(method != "SURE, homo, Symm8"
#                              & method != "SMASH, truth, Symm8"
#                              & method != "TI-thresh, truth, Symm8"
#                              & method != "TI-thresh, SMASH, Symm8"
#                              & method != "SMASH, joint, Haar"), aes(x = testfunction, y = mse, fill = method)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   xlab("Test Function") + 
#   ylab("Relative MSE") + 
#   ggtitle("MSE of various methods relative to \"SMASH, joint, Haar\", for i.i.d. Gaussian noise, SNR = 1") +
#   geom_abline(intercept = 1, slope = 0, linetype = 2, size = 0.3)
# dev.off()
# 
# 
# mean.table = melt(table.3.v1)
# names(mean.table) = c("method", "testfunction", "mse")
# 
# pdf("paper/smash_gaus_3_v1.pdf",width = 10, height = 5)
# ggplot(mean.table %>% filter(method != "SURE, homo, Symm8"
#                              & method != "SMASH, truth, Symm8"
#                              & method != "TI-thresh, truth, Symm8"
#                              & method != "TI-thresh, SMASH, Symm8"
#                              & method != "SMASH, joint, Haar"), aes(x = testfunction, y = mse, fill = method)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   xlab("Test Function") + 
#   ylab("Relative MSE") + 
#   ggtitle("MSE of various methods relative to \"SMASH, joint, Haar\", for i.i.d. Gaussian noise, SNR = 3") +
#   geom_abline(intercept = 1, slope = 0, linetype = 2, size = 0.3)
# dev.off()
# 
# 
# 
# # qplot(x = testfunction, y = mse, fill = method, xlab = "Test Function", ylab = "Relative MSE", main = "MSE of various methods relative to \"SMASH, joint, Haar\", for i.i.d. Gaussian noise, SNR = 1",
# #                        data = mean.table %>% filter(method != "SURE, homo, Symm8"
# #                                                     & method != "SMASH, truth, Symm8"
# #                                                     & method != "TI-thresh, truth, Symm8"
# #                                                     & method != "TI-thresh, SMASH, Symm8"
# #                                                     & method != "SMASH, joint, Haar"
# #                        )
# #                        , geom = "bar", stat = "identity",
# #                        position = "dodge")
# 
# 
# 
# mean.table = melt(table.1.v5)
# names(mean.table) = c("method", "testfunction", "mse")
# 
# 
# pdf("paper/smash_gaus_1_v5.pdf",width = 10, height = 5)
# ggplot(mean.table %>% filter(method != "SURE, homo, Symm8"
#                              & method != "SMASH, truth, Symm8"
#                              & method != "TI-thresh, truth, Symm8"
#                              & method != "TI-thresh, SMASH, Symm8"
#                              & method != "SMASH, joint, Haar"
#                              & method != "Ebayes, homo, Symm8"
#                              & method != "PostMean, homo, Symm8"), aes(x = testfunction, y = mse, fill = method)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   xlab("Test Function") + 
#   ylab("Relative MSE") + 
#   ggtitle("MSE of various methods relative to \"SMASH, joint, Haar\", for Clipped Blocks variance function, SNR = 1") +
#   theme(plot.title = element_text(size = rel(0.9))) +
#   geom_abline(intercept = 1, slope = 0, linetype = 2, size = 0.3)
# dev.off()
# 
# 
# 
# mean.table = melt(table.3.v5)
# names(mean.table) = c("method", "testfunction", "mse")
# 
# 
# pdf("paper/smash_gaus_3_v5.pdf",width = 10, height = 5)
# ggplot(mean.table %>% filter(method != "SURE, homo, Symm8"
#                              & method != "SMASH, truth, Symm8"
#                              & method != "TI-thresh, truth, Symm8"
#                              & method != "TI-thresh, SMASH, Symm8"
#                              & method != "SMASH, joint, Haar"
#                              & method != "Ebayes, homo, Symm8"
#                              & method != "PostMean, homo, Symm8"), aes(x = testfunction, y = mse, fill = method)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   xlab("Test Function") + 
#   ylab("Relative MSE") + 
#   ggtitle("MSE of various methods relative to \"SMASH, joint, Haar\", for Clipped Blocks variance function, SNR = 3") +
#   theme(plot.title = element_text(size = rel(0.9))) +
#   geom_abline(intercept = 1, slope = 0, linetype = 2, size = 0.3)
# dev.off()
# 
# 
# 
# ##error bars?
# 
# 
# 
# ##boxplot?
# 
# 
# ##lineplot
# ggplot(mean.table, aes(x = testfunction, y = mse, color = method, group = method)) + 
#        geom_line(linetype = 2) + 
#        geom_point(shape = 19) + 
#        xlab("Test Function") + 
#        ylab("Relative MSE")
# 



#########################################3
#violin plots
source("paper/vioplot_col.R")

homo.data.smash = res[res$.id == "sp.3.v1" & res$method == "smash.s8" , ]
homo.data.tithresh = res[res$.id == "sp.3.v1" & res$method == "tithresh.homo.s8" , ]
homo.data.ebayes = res[res$.id == "sp.3.v1" & res$method == "ebayesthresh" , ]
homo.data.smash.true = res[res$.id == "sp.3.v1" & res$method == "smash.true.s8" , ]

homo.data = res[res$.id == "sp.3.v1" & (res$method == "smash.s8" | res$method == "ebayesthresh" | res$method == "tithresh.homo.s8"), ]


# pdf("paper/violin_gaus_homo_comp.pdf", height = 8, width = 12)
# par(cex.axis = 2, cex.lab = 2, cex.sub = 2, mar = c(6.1, 6.1, 4.1, 2.1), mgp = c(3, 1.5, 0))
# vioplot.col(homo.data.smash$mise, homo.data.tithresh$mise, homo.data.ebayes$mise, ylim = c(5, 18), names = c("SMASH", "TI-thresholding", "EbayesThresh"), col = c("magenta", "gold", "gold"))
# title(xlab = "Method", ylab = "MISE", line = 4)
# abline(h = median(homo.data.smash$mise), lty = 3, col = 3)
# dev.off()
# 
# pdf("paper/violin_gaus_homo_smash.pdf", height = 8, width = 12)
# par(cex.axis = 2, cex.lab = 2, cex.sub = 2, mar = c(6.1, 6.1, 4.1, 2.1), mgp = c(3, 1.5, 0))
# vioplot.col(homo.data.smash$mise, homo.data.smash.true$mise, ylim = c(5, 18), names = c("SMASH", "SMASH, true variance"), col = c("magenta", "purple"))
# title(xlab = "Method", ylab = "MISE", line = 4)
# abline(h = median(homo.data.smash$mise), lty = 3, col = 3)
# dev.off()

# #vertical
# pdf("paper/violin_gaus_homo.pdf", height = 8, width = 12)
# par(cex.axis = 1.5, cex.lab = 1.5, cex.sub = 1.5, mar = c(6.1, 6.1, 4.1, 2.1), mgp = c(3, 1.5, 0))
# vioplot.col(homo.data.smash$mise, homo.data.tithresh$mise, homo.data.ebayes$mise, homo.data.smash.true$mise, ylim = c(5, 18), names = c("SMASH", "TI-thresholding", "EbayesThresh", "SMASH, true variance"), col = c("magenta", "gold", "gold", "purple"))
# title(xlab = "Method", ylab = "MISE", line = 4)
# abline(h = median(homo.data.smash$mise), lty = 3, col = 3)
# abline(v = 3.5, lty = 3)
# dev.off()

#horizontal
pdf("paper/violin_gaus_homo.pdf", height = 8, width = 8)
par(cex.axis = 1.5, cex.lab = 1.5, cex.sub = 1.5, mar = c(6.1, 12.1, 4.1, 2.1), mgp = c(3, 1.5, 0))
vioplot.col(homo.data.smash$mise, homo.data.tithresh$mise, homo.data.ebayes$mise, homo.data.smash.true$mise, horizontal = TRUE, ylim = c(5, 18), col = c("magenta", "gold", "gold", "purple"))
axis(2, 1:4, labels = c("SMASH", "TI-thresholding", "EbayesThresh", "SMASH\ntrue variance"), las = 2)
title(xlab = "MISE", line = 3.5)
abline(v = median(homo.data.smash$mise), lty = 3, col = 3)
abline(h = 3.5, lty = 3)
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
    return(1e-05 + 1 * (fn - min(fn))/max(fn))
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
plot(sd.fn, type = 'l', axes = FALSE, xlab = "", ylab = "")
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
plot(sd.fn.2, type = 'l', axes = FALSE, xlab = "", ylab = "")
axis(1, labels = FALSE, tick = FALSE)
axis(2)
dev.off()

####################
#smash + ash illustration
library(smash)
library(scales)
source("paper/code_plot.R")
source("paper/image.scale.R")

n = 1024
t = 1:n/n
spike.f = function(x) (0.75 * exp(-500 * (x - 0.23)^2) + 1.5 * exp(-2000 * (x - 0.33)^2) + 3 * exp(-8000 * (x - 0.47)^2) + 2.25 * exp(-16000 * (x - 0.69)^2) + 0.5 * exp(-32000 * (x - 0.83)^2))
mu.sp = spike.f(t)
mu.sp = (1 + mu.sp)/5

pos = c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 0.76, 0.78, 0.81)
hgt = 2.88/5 * c(4, (-5), 3, (-4), 5, (-4.2), 2.1, 4.3, (-3.1), 2.1, (-4.2))
sig.cb = rep(0, length(t))
for (j in 1:length(pos)) {
  sig.cb = sig.cb + (1 + sign(t - pos[j])) * (hgt[j]/2)
}
sig.cb[sig.cb < 0] = 0
sig.cb = 0.1 + (sig.cb - min(sig.cb))/max(sig.cb)
rsnr = sqrt(3)
sig.cb = sig.cb/mean(sig.cb) * sd(mu.sp)/rsnr^2

set.seed(70915)
x.sim = rnorm(n, mu.sp, sig.cb)

mu.smash = ashsmooth.gaus(x.sim, family = "DaubLeAsymm", filter.number = 8)
mu.ti = ti.thresh(x.sim, method = "rmad", family = "DaubLeAsymm", filter.number = 8)
sig.est = sqrt(2/(3 * (n - 2)) * sum((1/2 * x.sim[1:(n - 2)] - x.sim[2:(n - 1)] + 1/2 * x.sim[3:n])^2))
mu.ti.homo = ti.thresh(x.sim, sig = sig.est, family = "DaubLeAsymm", filter.number = 8)


##get wavelet coefficients and their variances
wc.sim = titable(x.sim)$difftable
#wc.var.sim = titable(sig.cb^2)$sumtable
wc.var.sim = titable(sig.cb^2)$sumtable

wc.true = titable(mu.sp)$difftable

##get shrunken estimates
wc.sim.shrunk = list()
wc.pres = list()
for(j in 0:(log2(n) - 1)){
  wc.sim.shrunk[[j+1]] = ash(wc.sim[j+2,],sqrt(wc.var.sim[j+2,]),prior="nullbiased",multiseqoutput=TRUE,pointmass=TRUE,nullcheck=TRUE,VB=FALSE,mixsd=NULL,mixcompdist="normal",gridmult=2,lambda1=1,lambda2=0,df=NULL,control=(list(trace=FALSE)))
  wc.pres[[j+1]] = 1/sqrt(wc.var.sim[j+2,])
}

rbPal = colorRampPalette(c('#191970', '#4169E1','#87CEEB'))

#plot the wc and their shrunken estimates for different resolutions
col.3 <- rev(rbPal(100))[as.numeric(cut(wc.pres[[3]],breaks = 100))]
wc.sig.3 = 1/wc.pres[[3]]
col.bw.3 = wc.sig.3*0.7/(max(wc.sig.3) - min(wc.sig.3)) - (0.7/(max(wc.sig.3) - min(wc.sig.3)))*min(wc.sig.3)
#plot(wc.sim[4, ], wc.sim.shrunk[[3]]$PosteriorMean, xlim = c(-1, 1), pch = 20, ylim = c(-1, 1), col = col.3)
plot(wc.sim[4, ], wc.sim.shrunk[[3]]$PosteriorMean, xlim = c(-1, 1), pch = 20, ylim = c(-1, 1), col = grey(col.bw.3))

col.9 <- rev(rbPal(100))[as.numeric(cut(wc.pres[[9]],breaks = 100))]
wc.sig.9 = 1/wc.pres[[9]]
col.bw.9 = wc.sig.9*0.7/(max(wc.sig.9) - min(wc.sig.9)) - (0.7/(max(wc.sig.9) - min(wc.sig.9)))*min(wc.sig.9)
#plot(wc.sim[8, ], wc.sim.shrunk[[7]]$PosteriorMean, xlim = c(-14, 14), pch = 20, ylim = c(-14, 14), col = col.7)
plot(wc.sim[10, ], wc.sim.shrunk[[9]]$PosteriorMean, xlim = c(-14, 14), pch = 20, ylim = c(-14, 14), col = grey(col.bw.9))


breaks = seq(min(col.bw.3), max(col.bw.3), length.out = 49)
breaks.ori = seq(min(wc.sig.3), max(wc.sig.3), length.out = 50)


pdf("paper/simple_eg_1.pdf", height = 6, width = 8)
par(cex.axis = 1.8, cex.sub = 1.8, cex.lab = 1.8, mar = c(5.1, 5.1, 2.1, 2.1), mgp = c(3, 1.3, 0))
plot(mu.sp, type = 'l', ylim = c(-0.05, 0.85), xlab = "position", ylab = "", lwd =1.7)
lines(mu.sp + 2*sig.cb, col = 3, lty = 2, lwd = 1.7)
lines(mu.sp - 2*sig.cb, col = 3, lty = 2, lwd = 1.7)
points(x.sim, cex = 0.7, pch = 16, col = alpha("black", 0.5))
dev.off()

pdf("paper/simple_eg_2.pdf", height = 6, width = 8)
par(cex.axis = 1.8, cex.sub = 1.8, cex.lab = 1.8, mar = c(5.1, 5.1, 2.1, 2.1), mgp = c(3, 1.3, 0))
hist(wc.true[4, ], breaks = 2, xlab = "wavelet coefficients", xlim = c(-25, 25), ylim = c(0, 800), col = "red", main = "")
hist(wc.true[10, ], breaks = 40, add = TRUE, col = rgb(0, 1, 0, 0.5))
dev.off()

pdf("paper/simple_eg_3.pdf", height = 6, width = 8)
par(cex.axis = 1.8, cex.sub = 1.8, cex.lab = 1.8, mar = c(5.1, 5.1, 2.1, 2.1), mgp = c(3, 1.3, 0))
plot(wc.sim[10, ], wc.sim.shrunk[[9]]$PosteriorMean, xlab = "empirical wavelet coefficients", ylab = "estimated wavelet coefficients", xlim = c(-4, 4), pch = 20, cex = 0.8, ylim = c(-4, 4), col = '#228B22')
points(wc.sim[4, ], wc.sim.shrunk[[3]]$PosteriorMean, pch = 20, cex = 0.8, col = grey(0.1))
abline(0, 1, lty = 3, col = grey(0.3))
legend('bottomright', legend = c("wavelet coefficients, scale 1", "wavelet coefficients, scale 7"), cex = 1.5, pch = c(20, 20), col = c(grey(0.1), '#228B22'))
dev.off()

pdf("paper/simple_eg_4.pdf", height = 5.8, width = 8)
layout(matrix(c(1, 1, 1, 3, 2, 0), nrow = 3), widths = c(8, 2), heights = c(2, 4, 2))
par(mar = c(7.1, 7.1, 2.1, 1.1), cex.axis = 2.5, cex.sub = 2.5, cex.lab = 2.5, mgp = c(3, 1.5, 0))
plot(wc.sim[4, ], wc.sim.shrunk[[3]]$PosteriorMean, xlim = c(-1, 1), pch = 20, cex = 0.8, xlab = "", ylab = "", ylim = c(-1, 1), col = grey(col.bw.3))
title(xlab = "empirical wavelet coefficients", ylab = "estimated wavelet coefficients", line = 4)
par(mar = c(1, 1, 0.1, 5.1), cex.axis = 2, cex.sub = 2, cex.lab = 2, mgp = c(3, 1.3, 0))
image.scale(volcano, col = grey(breaks), breaks = breaks.ori, horiz = FALSE, yaxt = "n")
axis(4, at = seq(0, 0.25, 0.05), las = 2)
par(mar = c(0.1, 0.1, 0.1 ,0.1))
plot(1, 1, xlim = c(0, 1), ylim = c(0, 1), type = "n", axes = FALSE)
legend(x = -0.25, y = 0.17, legend = "standard deviation", cex = 1.5, bty = "n")
dev.off()

pdf("paper/simple_eg_5.pdf", height = 6, width = 8)
par(cex.axis = 1.8, cex.sub = 1.8, cex.lab = 1.8, mar = c(6.1, 6.1, 2.1, 2.1), mgp = c(3, 1.3, 0))
plot(mu.sp, type = "l", lty = 2, ylim = c(-0.05, 0.85), xlab = "", ylab = "")
title(xlab = "position", ylab = "reconstructed mean", line = 4)
lines(mu.ti, col = "blue")
lines(mu.smash, col = "red")
legend(x = 550, y = 0.9, legend = c("SMASH", "TI-thresh \n RMAD variance"), col = c("red", "blue"), lty = c(1, 1), cex = 1.6, pt.cex = 0.5, bty = "n")
dev.off()


# ############
# #simple illustration
# library(smash)
# 
# n = 1024
# t = 1:n/n
# spike.f = function(x) (0.75 * exp(-500 * (x - 0.23)^2) + 1.5 * exp(-2000 * (x - 0.33)^2) + 3 * exp(-8000 * (x - 0.47)^2) + 2.25 * exp(-16000 * (x - 0.69)^2) + 0.5 * exp(-32000 * (x - 0.83)^2))
# mu.sp = spike.f(t)
# mu.sp = (1 + mu.sp)/5
# 
# pos = c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 0.76, 0.78, 0.81)
# hgt = 2.88/5 * c(4, (-5), 3, (-4), 5, (-4.2), 2.1, 4.3, (-3.1), 2.1, (-4.2))
# sig.cb = rep(0, length(t))
# for (j in 1:length(pos)) {
#   sig.cb = sig.cb + (1 + sign(t - pos[j])) * (hgt[j]/2)
# }
# sig.cb[sig.cb < 0] = 0
# sig.cb = 0.1 + (sig.cb - min(sig.cb))/max(sig.cb)
# rsnr = sqrt(3)
# sig.cb = sig.cb/mean(sig.cb) * sd(mu.sp)/rsnr^2
# 
# set.seed(70915)
# x.sim = rnorm(n, mu.sp, sig.cb)
# 
# mu.smash = ashsmooth.gaus(x.sim, family = "DaubLeAsymm", filter.number = 8)
# mu.ti = ti.thresh(x.sim, method = "rmad", family = "DaubLeAsymm", filter.number = 8)
# sig.est = sqrt(2/(3 * (n - 2)) * sum((1/2 * x.sim[1:(n - 2)] - x.sim[2:(n - 1)] + 1/2 * x.sim[3:n])^2))
# mu.ti.homo = ti.thresh(x.sim, sig = sig.est, family = "DaubLeAsymm", filter.number = 8)
# 
# pdf("paper/simple_eg.pdf", height = 16, width = 16)
# par(mfrow = c(2, 2), cex.axis = 2, cex.sub = 2, cex.lab = 2, mar = c(7.1, 7.1, 4.1, 3.1), mgp = c(3, 1.3, 0))
# plot(mu.sp, type = 'l', ylim = c(-0.05, 0.85), xlab = "", ylab = "")
# title(xlab = "position", ylab = "true mean", line = 4)
# plot(sig.cb, type = 'l', ylim = c(0, 0.1), xlab = "", ylab = "")
# title(xlab = "position", ylab = "standard deviation", line = 4)
# plot(x.sim, ylim = c(-0.05, 0.85), xlab = "", ylab = "")
# title(xlab = "position", ylab = "observed value", line = 4)
# plot(mu.sp, type = "l", lty = 2, ylim = c(-0.05, 0.85), xlab = "", ylab = "")
# title(xlab = "position", ylab = "reconstructed mean", line = 4)
# lines(mu.ti, col = "blue")
# lines(mu.smash, col = "red")
# legend(x = 600, y = 0.83, legend = c("SMASH", "TI-thresh \n RMAD variance"), col = c("red", "blue"), lty = c(1, 1), cex = 1.6, pt.cex = 0.5, bty = "n")
# dev.off()


# #######################
# ##ash illustration
# source("paper/code_plot.R")
# 
# n = 1024
# t = 1:n/n
# spike.f = function(x) (0.75 * exp(-500 * (x - 0.23)^2) + 1.5 * exp(-2000 * (x - 0.33)^2) + 3 * exp(-8000 * (x - 0.47)^2) + 2.25 * exp(-16000 * (x - 0.69)^2) + 0.5 * exp(-32000 * (x - 0.83)^2))
# mu.sp = spike.f(t)
# mu.sp = (1 + mu.sp)/5
# 
# 
# mu.wave=0.5+0.2*cos(4*pi*t)+0.1*cos(24*pi*t)
# 
# mu.sine = 0.01 * sin(2*pi*t)
# 
# mu.cor=623.87*t^3*(1-2*t)*(t>=0&t<=0.5)+187.161*(0.125-t^3)*t^4*(t>0.5&t<=0.8)+3708.470441*(t-1)^3*(t>0.8&t<=1)
# mu.cor=(0.6/(max(mu.cor)-min(mu.cor)))*mu.cor
# mu.cor=mu.cor-min(mu.cor)+0.2
# 
# 
# pos = c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 0.76, 0.78, 0.81)
# hgt = 2.88/5 * c(4, (-5), 3, (-4), 5, (-4.2), 2.1, 4.3, (-3.1), 2.1, (-4.2))
# sig.cb = rep(0, length(t))
# for (j in 1:length(pos)) {
#   sig.cb = sig.cb + (1 + sign(t - pos[j])) * (hgt[j]/2)
# }
# sig.cb[sig.cb < 0] = 0
# sig.cb = 0.1 + (sig.cb - min(sig.cb))/max(sig.cb)
# rsnr = sqrt(3)
# sig.cb = sig.cb/mean(sig.cb) * sd(mu.sp)/rsnr^2
# 
# sig.texp = 1e-04 + 4 * (exp(-550 * (t - 0.2)^2) + exp(-200 * (t - 0.5)^2) + exp(-950 * (t - 0.8)^2))
# 
# set.seed(70915)
# #x.sim = rnorm(n, mu.sp, sig.cb)
# x.sim = rnorm(n, mu.wave, sig.cb)
# 
# ##get wavelet coefficients and their variances
# wc.sim = titable(x.sim)$difftable
# #wc.var.sim = titable(sig.cb^2)$sumtable
# wc.var.sim = titable(sig.cb^2)$sumtable
# 
# wc.true = titable(mu.wave)$difftable
# 
# ##get shrunken estimates
# wc.sim.shrunk = list()
# wc.pres = list()
# for(j in 0:(log2(n) - 1)){
#   wc.sim.shrunk[[j+1]] = ash(wc.sim[j+2,],sqrt(wc.var.sim[j+2,]),prior="nullbiased",multiseqoutput=TRUE,pointmass=TRUE,nullcheck=TRUE,VB=FALSE,mixsd=NULL,mixcompdist="normal",gridmult=2,lambda1=1,lambda2=0,df=NULL,control=(list(trace=FALSE)))
#   wc.pres[[j+1]] = 1/sqrt(wc.var.sim[j+2,])
# }
# 
# rbPal = colorRampPalette(c('red', 'white','blue'))
# 
# #plot the wc and their shrunken estimates for different resolutions
# col.3 <- rbPal(10)[as.numeric(cut(wc.pres[[3]],breaks = 10))]
# wc.sig.3 = 1/wc.pres[[3]]
# col.bw.3 = wc.sig.3*0.7/(max(wc.sig.3) - min(wc.sig.3)) - (0.7/(max(wc.sig.3) - min(wc.sig.3)))*min(wc.sig.3)
# plot(wc.sim[4, ], wc.sim.shrunk[[3]]$PosteriorMean, xlim = c(-1, 1), pch = 20, ylim = c(-1, 1), col = col.3)
# plot(wc.sim[4, ], wc.sim.shrunk[[3]]$PosteriorMean, xlim = c(-1, 1), pch = 20, ylim = c(-1, 1), col = grey(col.bw.3))
# 
# col.7 <- rbPal(10)[as.numeric(cut(wc.pres[[7]],breaks = 10))]
# wc.sig.7 = 1/wc.pres[[7]]
# col.bw.7 = wc.sig.7*0.7/(max(wc.sig.7) - min(wc.sig.7)) - (0.7/(max(wc.sig.7) - min(wc.sig.7)))*min(wc.sig.7)
# plot(wc.sim[8, ], wc.sim.shrunk[[7]]$PosteriorMean, xlim = c(-14, 14), pch = 20, ylim = c(-14, 14), col = col.7)
# plot(wc.sim[8, ], wc.sim.shrunk[[7]]$PosteriorMean, xlim = c(-14, 14), pch = 20, ylim = c(-14, 14), col = grey(col.bw.7))
# 
# 
# 
# set.seed(70915)
# #x.sim = rnorm(n, mu.sp, sig.cb)
# x.sim.2 = rnorm(n, mu.sine, sig.cb)
# 
# ##get wavelet coefficients and their variances
# wc.sim.2 = titable(x.sim.2)$difftable
# #wc.var.sim = titable(sig.cb^2)$sumtable
# wc.var.sim.2 = titable(sig.cb^2)$sumtable
# 
# wc.true.2 = titable(mu.sine)$difftable
# 
# ##get shrunken estimates
# wc.sim.shrunk.2 = list()
# wc.pres.2 = list()
# for(j in 0:(log2(n) - 1)){
#   wc.sim.shrunk.2[[j+1]] = ash(wc.sim.2[j+2,],sqrt(wc.var.sim.2[j+2,]),prior="nullbiased",multiseqoutput=TRUE,pointmass=TRUE,nullcheck=TRUE,VB=FALSE,mixsd=NULL,mixcompdist="normal",gridmult=2,lambda1=1,lambda2=0,df=NULL,control=(list(trace=FALSE)))
#   wc.pres.2[[j+1]] = 1/sqrt(wc.var.sim.2[j+2,])
# }
# 
# rbPal = colorRampPalette(c('red', 'white','blue'))
# 
# #plot the wc and their shrunken estimates for different resolutions
# col.3 <- rbPal(10)[as.numeric(cut(wc.pres[[3]],breaks = 10))]
# wc.sig.3 = 1/wc.pres.2[[3]]
# col.bw.3 = wc.sig.3*0.7/(max(wc.sig.3) - min(wc.sig.3)) - (0.7/(max(wc.sig.3) - min(wc.sig.3)))*min(wc.sig.3)
# plot(wc.sim.2[4, ], wc.sim.shrunk.2[[3]]$PosteriorMean, xlim = c(-0.6, 0.6), pch = 20, ylim = c(-0.6, 0.6), col = col.3)
# plot(wc.sim.2[4, ], wc.sim.shrunk.2[[3]]$PosteriorMean, xlim = c(-0.6, 0.6), pch = 20, ylim = c(-0.6, 0.6), col = grey(col.bw.3))
# 
# col.7 <- rbPal(10)[as.numeric(cut(wc.pres[[7]],breaks = 10))]
# wc.sig.7 = 1/wc.pres.2[[7]]
# col.bw.7 = wc.sig.7*0.7/(max(wc.sig.7) - min(wc.sig.7)) - (0.7/(max(wc.sig.7) - min(wc.sig.7)))*min(wc.sig.7)
# plot(wc.sim.2[8, ], wc.sim.shrunk.2[[7]]$PosteriorMean, xlim = c(-14, 14), pch = 20, ylim = c(-14, 14), col = col.7)
# plot(wc.sim.2[8, ], wc.sim.shrunk.2[[7]]$PosteriorMean, xlim = c(-14, 14), pch = 20, ylim = c(-14, 14), col = grey(col.bw.7))
# abline(0,1)
# 
# 
# hist(wc.true[4, ], breaks = 50, xlim = c(-0.2, 0.2), ylim = c(0, 500), col = "red")
# hist(wc.true.2[4, ], breaks = c(-0.005, 0, 0.005), add = TRUE, col = rgb(0, 1, 0, 0.5))
# 
# hist(wc.true[8, ], breaks = 50, xlim = c(-14, 14), ylim = c(0, 500), col = "red")
# hist(wc.true.2[8, ], breaks = c(-0.5, 0, 0.5), add = TRUE, col = rgb(0, 1, 0, 0.5))
# 
# 
# par(mfrow = c(2, 2))
# col.3 <- rbPal(10)[as.numeric(cut(wc.pres[[3]],breaks = 10))]
# wc.sig.3 = 1/wc.pres[[3]]
# col.bw.3 = wc.sig.3*0.7/(max(wc.sig.3) - min(wc.sig.3)) - (0.7/(max(wc.sig.3) - min(wc.sig.3)))*min(wc.sig.3)
# #plot(wc.sim[4, ], wc.sim.shrunk[[3]]$PosteriorMean, xlim = c(-1, 1), pch = 20, ylim = c(-1, 1), col = col.3)
# plot(wc.sim[4, ], wc.sim.shrunk[[3]]$PosteriorMean, xlim = c(-1, 1), pch = 20, ylim = c(-1, 1), col = grey(col.bw.3), xlab = "Observed wavelet coefficients", ylab = "Shrunken wavelet coefficients", main = "Scale 7, not-so-peaky prior")
# abline(0, 1)
# 
# col.7 <- rbPal(10)[as.numeric(cut(wc.pres[[7]],breaks = 10))]
# wc.sig.7 = 1/wc.pres[[7]]
# col.bw.7 = wc.sig.7*0.7/(max(wc.sig.7) - min(wc.sig.7)) - (0.7/(max(wc.sig.7) - min(wc.sig.7)))*min(wc.sig.7)
# #plot(wc.sim[8, ], wc.sim.shrunk[[7]]$PosteriorMean, xlim = c(-14, 14), pch = 20, ylim = c(-14, 14), col = col.8)
# plot(wc.sim[8, ], wc.sim.shrunk[[7]]$PosteriorMean, xlim = c(-14, 14), pch = 20, ylim = c(-14, 14), col = grey(col.bw.7), xlab = "Observed wavelet coefficients", ylab = "Shrunken wavelet coefficients", main = "Scale 3, not-so-peaky prior")
# abline(0, 1)
# 
# col.3 <- rbPal(10)[as.numeric(cut(wc.pres[[3]],breaks = 10))]
# wc.sig.3 = 1/wc.pres.2[[3]]
# col.bw.3 = wc.sig.3*0.7/(max(wc.sig.3) - min(wc.sig.3)) - (0.7/(max(wc.sig.3) - min(wc.sig.3)))*min(wc.sig.3)
# #plot(wc.sim.2[4, ], wc.sim.shrunk.2[[3]]$PosteriorMean, xlim = c(-0.6, 0.6), pch = 20, ylim = c(-0.6, 0.6), col = col.3)
# plot(wc.sim.2[4, ], wc.sim.shrunk.2[[3]]$PosteriorMean, xlim = c(-0.6, 0.6), pch = 20, ylim = c(-0.6, 0.6), col = grey(col.bw.3), xlab = "Observed wavelet coefficients", ylab = "Shrunken wavelet coefficients", main = "Scale 7, peaky prior")
# abline(0, 1)
# 
# col.7 <- rbPal(10)[as.numeric(cut(wc.pres[[7]],breaks = 10))]
# wc.sig.7 = 1/wc.pres.2[[7]]
# col.bw.7 = wc.sig.7*0.7/(max(wc.sig.7) - min(wc.sig.7)) - (0.7/(max(wc.sig.7) - min(wc.sig.7)))*min(wc.sig.7)
# #plot(wc.sim.2[8, ], wc.sim.shrunk.2[[7]]$PosteriorMean, xlim = c(-14, 14), pch = 20, ylim = c(-14, 14), col = col.8)
# plot(wc.sim.2[8, ], wc.sim.shrunk.2[[7]]$PosteriorMean, xlim = c(-3, 3), pch = 20, ylim = c(-3, 3), col = grey(col.bw.7), xlab = "Observed wavelet coefficients", ylab = "Shrunken wavelet coefficients", main = "Scale 7, peaky prior")
# abline(0, 1)


