# TO DO: Explain here what this script does, and how to use it.
#
# TO DO: Change the name of this file.
#
# SET UP ENVIROMENT
# -----------------
# Load the sm and vioplot packages, and load the custom violin plot
# function.
library(ggplot2)
library(cowplot)
# library(sm)
# library(vioplot)
source("../code/vioplot_col.R")

# LOAD RESULTS
# ------------
# Load the results of the simulation experiments.
cat("Loading simulation results.\n")
load("../output/gaus-dscr.RData")

homo.data.smash <-
  res[res$.id    == "sp.3.v1" &
      res$method == "smash.s8",]
homo.data.smash.homo <-
  res[res$.id    == "sp.3.v1" &
      res$method == "smash.homo.s8",]
homo.data.tithresh <-
  res[res$.id == "sp.3.v1" &
      res$method == "tithresh.homo.s8",]
homo.data.ebayes <-
  res[res$.id    == "sp.3.v1" &
      res$method == "ebayesthresh",]
homo.data.smash.true <-
  res[res$.id == "sp.3.v1" &
  res$method  == "smash.true.s8",]
homo.data <-
  res[res$.id == "sp.3.v1" &
  (res$method == "smash.s8" |
   res$method == "ebayesthresh" |
   res$method == "tithresh.homo.s8"),]

# PLOT #1
# -------
# This reproduces Fig. 2 of the manuscript (comparison of accuracy of
# estimated mean curves for data simulated with homoskedastic Gaussian
# errors.)
pdat <-
  rbind(data.frame(method      = "smash",
                   method.type = "est",
                   mise        = homo.data.smash$mise),
        data.frame(method      = "smash.homo",
                   method.type = "homo",
                   mise        = homo.data.smash.homo$mise),
        data.frame(method      = "tithresh",
                   method.type = "homo",
                   mise        = homo.data.tithresh$mise),
        data.frame(method      = "ebayesthresh",
                   method.type = "homo",
                   mise        = homo.data.ebayes$mise),
        data.frame(method      = "smash.true",
                   method.type = "true",
                   mise        = homo.data.smash.true$mise))
pdat <-
  transform(pdat,
            method = factor(method,
                            names(sort(tapply(pdat$mise,pdat$method,mean),
                                       decreasing = TRUE))))
p <- ggplot(pdat,aes(x = method,y = mise,fill = method.type)) +
     geom_violin(fill = "skyblue",color = "skyblue") +
     geom_boxplot(width = 0.15,outlier.shape = NA) +
     scale_y_continuous(breaks = seq(6,16,2)) +
     scale_fill_manual(values = c("darkorange","dodgerblue","gold"),
                       guide = FALSE) +
     coord_flip() +
     labs(x = "MISE",y = "") +
     theme(axis.line = element_blank(),
           axis.ticks.y = element_blank())
print(p)
        par(cex.axis = 0.75,cex.lab = 1,cex.sub = 1)
vioplot.col(homo.data.smash$mise,
            homo.data.smash.homo$mise,
            homo.data.tithresh$mise,
            homo.data.ebayes$mise,
            homo.data.smash.true$mise,horizontal = TRUE,ylim = c(5,18),
            col = c("magenta","gold","gold","gold","purple"))
axis(2,1:5,las = 2,
     labels = c("SMASH","SMASH (homo)","TI-thresholding (homo)",
                "EbayesThresh (homo)","SMASH (true variance)"))
title(xlab = "MISE", line = 3.5)
abline(v = median(homo.data.smash$mise),lty = 3,col = 3)
abline(h = 4.5, lty = 3)

stop()

cblocks.fn = function(t, type) {
  pos = c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 0.76, 0.78, 0.81)
  hgt = 2.88/5 * c(4, (-5), 3, (-4), 5, (-4.2), 2.1, 4.3, (-3.1), 2.1, (-4.2))
  fn = rep(0, length(t))
  for (j in 1:length(pos))
    fn = fn + (1 + sign(t - pos[j])) * (hgt[j]/2)
  fn[fn < 0] = 0
  if (type == "mean")
    return(NULL)
  else if (type == "var")
    return(0.01 + 1 * (fn - min(fn))/max(fn))
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
  fn = 623.87 * t^3 * (1 - 2 * t) * (t >= 0 & t <= 0.5) + 187.161 *
      (0.125 - t^3) * t^4 * (t > 0.5 & t <= 0.8) +
          3708.470441 * (t - 1)^3 * (t > 0.8 & t <= 1)
  fn = (0.6/(max(fn) - min(fn))) * fn
  if (type == "mean") {
    return(fn - min(fn) + 0.2)
  } else if (type == "var") {
    return(NULL)
  }
}

spikes.fn = function(t, type) {
  spike.f = function(x)
    (0.75  * exp(-500 * (x - 0.23)^2) +
      1.5  * exp(-2000 * (x - 0.33)^2) +
      3    * exp(-8000 * (x - 0.47)^2) +
      2.25 * exp(-16000 * (x - 0.69)^2) +
      0.5 * exp(-32000 * (x - 0.83)^2))
  fn = spike.f(t)
  if (type == "mean") {
    return((1 + fn)/5)
  } else if (type == "var") {
    return(NULL)
  }
}

t = (1:1024)/1024

hetero.data.smash = res[res$.id == "sp.3.v5" & res$method == "smash.s8" , ]
hetero.data.smash.homo = res[res$.id == "sp.3.v5" &
    res$method == "smash.homo.s8" , ]
hetero.data.tithresh.homo = res[res$.id == "sp.3.v5" &
    res$method == "tithresh.homo.s8" , ]
hetero.data.tithresh.rmad = res[res$.id == "sp.3.v5" &
    res$method == "tithresh.rmad.s8" , ]
hetero.data.tithresh.smash = res[res$.id == "sp.3.v5" &
    res$method == "tithresh.smash.s8" , ]
hetero.data.tithresh.true = res[res$.id == "sp.3.v5" &
    res$method == "tithresh.true.s8" , ]
hetero.data.ebayes = res[res$.id == "sp.3.v5" &
    res$method == "ebayesthresh" , ]
hetero.data.smash.true = res[res$.id == "sp.3.v5" &
    res$method == "smash.true.s8" , ]

vioplot(hetero.data.smash$mise, hetero.data.tithresh.rmad$mise,
        hetero.data.tithresh.smash$mise, hetero.data.smash.homo$mise,
        hetero.data.tithresh.homo$mise, hetero.data.ebayes$mise,
        hetero.data.smash.true$mise, hetero.data.tithresh.true$mise)

mu = spikes.fn(t, "mean")
sigma.ini = sqrt(cblocks.fn(t, "var"))
sd.fn = sigma.ini/mean(sigma.ini) * sd(mu)/3

#horizontal
par(cex.axis = 1.5, cex.lab = 1.5, cex.sub = 1.5,
    mar = c(6.1, 11.1, 4.1, 0.5), mgp = c(3, 1.5, 0))
vioplot.col(hetero.data.smash$mise, hetero.data.tithresh.rmad$mise,
            hetero.data.tithresh.smash$mise, hetero.data.tithresh.true$mise,
            hetero.data.smash.homo$mise, hetero.data.ebayes$mise,
            hetero.data.smash.true$mise, ylim = c(10, 80), horizontal = TRUE,
            names=NULL, col = c("magenta", rep("red", 3), rep("gold", 2),
                            "purple"))
par(lheight = 0.8)
axis(side = 2, at = 1:7,
     labels = c("SMASH", "TI-thresh\nRMAD variance",
         "TI-thresh\nSMASH variance", "TI-thresh\ntrue variance",
         "SMASH\nhomo", "EbayesThresh\nhomo", "SMASH\ntrue variance"), las = 2)
title(xlab = "MISE", line = 4)
abline(v = median(hetero.data.smash$mise), lty = 3, col = 3)
abline(h = 4.5, lty = 3)
abline(h = 6.5, lty = 3)

par(cex.axis = 1.5, cex.lab = 1.5, cex.sub = 1.5, mar = c(6.1, 3.1, 4.1, 2.1))
plot(mu, type = 'l', ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
lines(mu + 2* sd.fn, col = 2, lty = 5)
lines(mu - 2* sd.fn, col = 2, lty = 5)
axis(1, labels = FALSE, tick = FALSE)
axis(2)

hetero.data.smash.2 = res[res$.id == "cor.3.v3" & res$method == "smash.s8" , ]
hetero.data.smash.homo.2 = res[res$.id == "cor.3.v3" &
    res$method == "smash.homo.s8" , ]
hetero.data.tithresh.homo.2 = res[res$.id == "cor.3.v3" &
    res$method == "tithresh.homo.s8" , ]
hetero.data.tithresh.rmad.2 = res[res$.id == "cor.3.v3" &
    res$method == "tithresh.rmad.s8" , ]
hetero.data.tithresh.smash.2 = res[res$.id == "cor.3.v3" &
    res$method == "tithresh.smash.s8" , ]
hetero.data.tithresh.true.2 = res[res$.id == "cor.3.v3" &
    res$method == "tithresh.true.s8" , ]
hetero.data.ebayes.2 = res[res$.id == "cor.3.v3" &
    res$method == "ebayesthresh" , ]
hetero.data.smash.true.2 = res[res$.id == "cor.3.v3" &
    res$method == "smash.true.s8" , ]

#horizontal
par(cex.axis = 1.5, cex.lab = 1.5, cex.sub = 1.5,
    mar = c(6.1, 11.1, 4.1, 0.5), mgp = c(3, 1.5, 0))
vioplot.col(hetero.data.smash.2$mise, hetero.data.tithresh.rmad.2$mise,
            hetero.data.tithresh.smash.2$mise,
            hetero.data.tithresh.true.2$mise, hetero.data.smash.homo.2$mise,
            hetero.data.ebayes.2$mise, hetero.data.smash.true.2$mise,
            horizontal = TRUE, ylim = c(0.5,5),
            col = c("magenta", rep("red", 3), rep("gold", 2), "purple"))
par(lheight = 0.8)
axis(side = 2, at = 1:7,
     labels = c("SMASH", "TI-thresh\nRMAD variance",
         "TI-thresh\nSMASH variance", "TI-thresh\ntrue variance",
         "SMASH\nhomo", "EbayesThresh\nhomo", "SMASH\ntrue variance"), las = 2)
title(xlab = "MISE", line = 4)
abline(v = median(hetero.data.smash.2$mise), lty = 3, col = 3)
abline(h = 4.5, lty = 3)
abline(h = 6.5, lty = 3)

#par(cex.axis = 1.5, cex.lab = 1.5, cex.sub = 1.5, mar = c(6.1, 3.1, 4.1, 2.1))
# plot(mu, type = 'l', ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
# lines(mu + 2* sd.fn.2, col = 2, lty = 5)
# lines(mu - 2* sd.fn.2, col = 2, lty = 5)
# axis(1, labels = FALSE, tick = FALSE)
# axis(2)

# SESSIONINFO
# -----------
# Print out information on the computing environment, including the
# version of R and the attached R packages, used to generate these
# results.
print(sessionInfo())
