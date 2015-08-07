library(dplyr)
library(reshape2)
library(ggplot2)
library(vioplot)


##barplots in main body of paper
mean.table = melt(table.1.v1)
names(mean.table) = c("method", "testfunction", "mse")

pdf("paper/smash_gaus_1_v1.pdf",width = 10, height = 5)
ggplot(mean.table %>% filter(method != "SURE, homo, Symm8"
                             & method != "SMASH, truth, Symm8"
                             & method != "TI-thresh, truth, Symm8"
                             & method != "TI-thresh, SMASH, Symm8"
                             & method != "SMASH, joint, Haar"), aes(x = testfunction, y = mse, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  xlab("Test Function") + 
  ylab("Relative MSE") + 
  ggtitle("MSE of various methods relative to \"SMASH, joint, Haar\", for i.i.d. Gaussian noise, SNR = 1") +
  geom_abline(intercept = 1, slope = 0, linetype = 2, size = 0.3)
dev.off()


mean.table = melt(table.3.v1)
names(mean.table) = c("method", "testfunction", "mse")

pdf("paper/smash_gaus_3_v1.pdf",width = 10, height = 5)
ggplot(mean.table %>% filter(method != "SURE, homo, Symm8"
                             & method != "SMASH, truth, Symm8"
                             & method != "TI-thresh, truth, Symm8"
                             & method != "TI-thresh, SMASH, Symm8"
                             & method != "SMASH, joint, Haar"), aes(x = testfunction, y = mse, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  xlab("Test Function") + 
  ylab("Relative MSE") + 
  ggtitle("MSE of various methods relative to \"SMASH, joint, Haar\", for i.i.d. Gaussian noise, SNR = 3") +
  geom_abline(intercept = 1, slope = 0, linetype = 2, size = 0.3)
dev.off()



# qplot(x = testfunction, y = mse, fill = method, xlab = "Test Function", ylab = "Relative MSE", main = "MSE of various methods relative to \"SMASH, joint, Haar\", for i.i.d. Gaussian noise, SNR = 1",
#                        data = mean.table %>% filter(method != "SURE, homo, Symm8"
#                                                     & method != "SMASH, truth, Symm8"
#                                                     & method != "TI-thresh, truth, Symm8"
#                                                     & method != "TI-thresh, SMASH, Symm8"
#                                                     & method != "SMASH, joint, Haar"
#                        )
#                        , geom = "bar", stat = "identity",
#                        position = "dodge")



mean.table = melt(table.1.v5)
names(mean.table) = c("method", "testfunction", "mse")


pdf("paper/smash_gaus_1_v5.pdf",width = 10, height = 5)
ggplot(mean.table %>% filter(method != "SURE, homo, Symm8"
                             & method != "SMASH, truth, Symm8"
                             & method != "TI-thresh, truth, Symm8"
                             & method != "TI-thresh, SMASH, Symm8"
                             & method != "SMASH, joint, Haar"
                             & method != "Ebayes, homo, Symm8"
                             & method != "PostMean, homo, Symm8"), aes(x = testfunction, y = mse, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  xlab("Test Function") + 
  ylab("Relative MSE") + 
  ggtitle("MSE of various methods relative to \"SMASH, joint, Haar\", for Clipped Blocks variance function, SNR = 1") +
  theme(plot.title = element_text(size = rel(0.9))) +
  geom_abline(intercept = 1, slope = 0, linetype = 2, size = 0.3)
dev.off()



mean.table = melt(table.3.v5)
names(mean.table) = c("method", "testfunction", "mse")


pdf("paper/smash_gaus_3_v5.pdf",width = 10, height = 5)
ggplot(mean.table %>% filter(method != "SURE, homo, Symm8"
                             & method != "SMASH, truth, Symm8"
                             & method != "TI-thresh, truth, Symm8"
                             & method != "TI-thresh, SMASH, Symm8"
                             & method != "SMASH, joint, Haar"
                             & method != "Ebayes, homo, Symm8"
                             & method != "PostMean, homo, Symm8"), aes(x = testfunction, y = mse, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  xlab("Test Function") + 
  ylab("Relative MSE") + 
  ggtitle("MSE of various methods relative to \"SMASH, joint, Haar\", for Clipped Blocks variance function, SNR = 3") +
  theme(plot.title = element_text(size = rel(0.9))) +
  geom_abline(intercept = 1, slope = 0, linetype = 2, size = 0.3)
dev.off()



##error bars?



##boxplot?


##lineplot
ggplot(mean.table, aes(x = testfunction, y = mse, color = method, group = method)) + 
       geom_line(linetype = 2) + 
       geom_point(shape = 19) + 
       xlab("Test Function") + 
       ylab("Relative MSE")




#########################################3
#violin plots
homo.data.smash = res[res$.id == "sp.3.v1" & res$method == "smash.s8" , ]
homo.data.tithresh = res[res$.id == "sp.3.v1" & res$method == "tithresh.homo.s8" , ]
homo.data.ebayes = res[res$.id == "sp.3.v1" & res$method == "ebayesthresh" , ]
homo.data.smash.true = res[res$.id == "sp.3.v1" & res$method == "smash.true.s8" , ]

homo.data = res[res$.id == "sp.3.v1" & (res$method == "smash.s8" | res$method == "ebayesthresh" | res$method == "tithresh.homo.s8"), ]


pdf("paper/violin_gaus_homo_comp.pdf", height = 8, width = 12)
vioplot(homo.data.smash$mise, homo.data.tithresh$mise, homo.data.ebayes$mise, names = c("SMASH", "TI-thresholding", "EbayesThresh"))
title(ylab = "MISE", xlab = "Method")
abline(h = median(homo.data.smash$mise), lty = 3, col = 3)
dev.off()

pdf("paper/violin_gaus_homo_smash.pdf", height = 8, width = 12)
vioplot(homo.data.smash$mise, homo.data.smash.true$mise, names = c("SMASH", "SMASH, true variance"))
title(ylab = "MISE", xlab = "Method")
abline(h = median(homo.data.smash$mise), lty = 3, col = 3)
dev.off()


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

pdf("paper/violin_gaus_hetero_ti.pdf", height = 8, width = 12)
vioplot(hetero.data.smash$mise, hetero.data.tithresh.rmad$mise, hetero.data.tithresh.smash$mise,  hetero.data.tithresh.true$mise, names = c("SMASH", "TI-thresh, RMAD", "TI-thresh, SMASH", "TI-thresh, true variance"))
abline(h = median(hetero.data.smash$mise), lty = 3, col = 3)
dev.off()

pdf("paper/violin_gaus_hetero_homo.pdf", height = 8, width = 12)
vioplot(hetero.data.smash$mise, hetero.data.smash.homo$mise, hetero.data.tithresh.homo$mise,hetero.data.ebayes$mise, names = c("SMASH", "SMASH, homo", "TI-thresh, homo", "EbayesThresh, homo"))
abline(h = median(hetero.data.smash$mise), lty = 3, col = 3)
dev.off()

pdf("paper/violin_gaus_hetero_smash.pdf", height = 8, width = 12)
vioplot(hetero.data.smash$mise, hetero.data.smash.true$mise, names = c("SMASH", "SMASH, true variance"))
abline(h = median(hetero.data.smash$mise), lty = 3, col = 3)
dev.off()


############
#simple illustration
library(smash)

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

pdf("paper/simple_eg.pdf", height = 16, width = 16)
par(mfrow = c(2, 2))
plot(mu.sp, type = 'l', ylim = c(0.15, 0.85), xlab = "position", ylab = "mean")
plot(sig.cb, type = 'l', ylim = c(0, 0.1), xlab = "position", ylab = "standard deviation")
plot(x.sim, ylim = c(-0.1, 1), xlab = "position", ylab = "observed value")
plot(mu.sp, type = "l", lty = 2, ylim = c(0.15, 0.85), xlab = "position", ylab = "mean")
lines(mu.ti, col = "blue")
lines(mu.smash, col = "red")
legend(x = 750, y = 0.7, legend = c("SMASH", "TI-thresh, RMAD"), col = c("red", "blue"), lty = c(1, 1), cex = 1.2, pt.cex = 0.5, bty = "n")
dev.off()