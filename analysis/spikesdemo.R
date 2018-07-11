# TO DO: Explain here what this code does, and how to use it.
source("../code/plotting.R")
# source("code_paper/sim/image.scale.R")

# SET UP ENVIROMENT
# -----------------
# Load the ashr and smashr package.
library(ashr)
library(smashr)
# library(scales)

# Initialize the sequence of pseudorandom numbers.
set.seed(70915)

# DEFINE "SPIKES" MEAN FUNCTION
# -----------------------------
n <- 1024
t <- 1:n/n
spike.f <-
  function (x) (0.75 * exp(-500   * (x - 0.23)^2) +
                1.5  * exp(-2000  * (x - 0.33)^2) +
                3    * exp(-8000  * (x - 0.47)^2) +
                2.25 * exp(-16000 * (x - 0.69)^2) +
                0.5  * exp(-32000 * (x - 0.83)^2))
mu.sp <- spike.f(t)
mu.sp <- (1 + mu.sp)/5

# SIMULATE DATA
# -------------
cat("Generating data.\n")
pos <- c(0.1,0.13,0.15,0.23,0.25,0.4,0.44,0.65,0.76,0.78,0.81)
hgt <- 2.88/5 * c(4,-5,3,-4,5,-4.2,2.1,4.3,-3.1,2.1,-4.2)
sig.cb <- rep(0,length(t))
for (j in 1:length(pos)) 
  sig.cb <- sig.cb + (1 + sign(t - pos[j])) * (hgt[j]/2)
sig.cb[sig.cb < 0] <- 0
sig.cb <- 0.1 + (sig.cb - min(sig.cb))/max(sig.cb)
rsnr   <- sqrt(3)
sig.cb <- sig.cb/mean(sig.cb) * sd(mu.sp)/rsnr^2
x.sim  <- rnorm(n,mu.sp,sig.cb)

# Plot the simulated data.
par(cex.axis = 1.8,cex.sub = 1.8,cex.lab = 1.8,mar = c(5.1,5.1,2.1,2.1),
    mgp = c(3,1.3,0))
plot(mu.sp,type = 'l',ylim = c(-0.05,0.85),xlab = "",ylab = "",lwd = 1.7)
lines(mu.sp + 2*sig.cb,col = "darkorange",lty = 5,lwd = 1.8)
lines(mu.sp - 2*sig.cb,col = "darkorange",lty = 5,lwd = 1.8)
points(x.sim,cex = 0.7,pch = 16,col = "darkblue")
invisible(readline(prompt = "Press [enter] to continue analysis... "))

stop()

# RUN SMASH
# ---------
# Apply smash and translation invariant (TI) thresholding to the
# "spikes" data. Here, we run the TI thresholding twice---once when
# the standard deviation (s.d.) function is provided, and once when it
# is estimated using the MAD algorithm.
cat("Running smash and TI thresholding on simulated data.\n")
sig.est  <- sqrt(2/(3 * (n - 2)) *
                 sum((1/2 * x.sim[1:(n-2)] - x.sim[2:(n-1)] + x.sim[3:n])^2/2))
mu.smash <- smash(x.sim,family = "DaubLeAsymm",filter.number = 8)
mu.ti    <- ti.thresh(x.sim,method = "rmad",family = "DaubLeAsymm",
                      filter.number = 8)
mu.ti.homo <- ti.thresh(x.sim,sigma = sig.est,family = "DaubLeAsymm",
                        filter.number = 8)

# Get the wavelet coefficients and their variances.
wc.sim     <- titable(x.sim)$difftable
wc.var.sim <- titable(sig.cb^2)$sumtable
wc.true    <- titable(mu.sp)$difftable

# Get shrunken estimates of the wavelet coefficients.
wc.sim.shrunk <- vector("list",10)
wc.pres       <- vector("list",10)
for(j in 0:(log2(n) - 1)){
  wc.sim.shrunk[[j+1]] <-
    ash(wc.sim[j+2,],sqrt(wc.var.sim[j+2,]),prior = "nullbiased",
        pointmass = TRUE,mixsd = NULL,mixcompdist = "normal",
        gridmult = 2,df = NULL)$result
  wc.pres[[j+1]] <- 1/sqrt(wc.var.sim[j+2,])
}

# SUMMARIZE RESULTS
# -----------------
rbPal <- colorRampPalette(c('#191970','#4169E1','#87CEEB'))

# TO DO: Get rid of this plot.
## col.3    <- rev(rbPal(100))[as.numeric(cut(wc.pres[[3]],breaks = 100))]
## wc.sig.3 <- 1/wc.pres[[3]]
## col.bw.3 <- wc.sig.3*0.8/(max(wc.sig.3) - min(wc.sig.3)) -
##             (0.8/(max(wc.sig.3) - min(wc.sig.3)))*min(wc.sig.3)
## plot(wc.sim[4,],wc.sim.shrunk[[3]]$PosteriorMean,xlim = c(-1,1),pch = 20,
##      ylim = c(-1,1),col = grey(col.bw.3))

# TO DO: Get rid of this plot as well.
## col.9    <- rev(rbPal(100))[as.numeric(cut(wc.pres[[9]],breaks = 100))]
## wc.sig.9 <- 1/wc.pres[[9]]
## col.bw.9 <- wc.sig.9*0.7/(max(wc.sig.9) - min(wc.sig.9)) -
##     (0.7/(max(wc.sig.9) - min(wc.sig.9)))*min(wc.sig.9)
## plot(wc.sim[10,],wc.sim.shrunk[[9]]$PosteriorMean,xlim = c(-14,14),
##      pch = 20,ylim = c(-14,14),col = grey(col.bw.9))

## breaks     <- seq(min(col.bw.3),max(col.bw.3),length.out = 49)
## breaks.ori <- seq(min(wc.sig.3),max(wc.sig.3),length.out = 50)

pdf("paper/simple_eg_2.pdf",height = 6,width = 8)
par(cex.axis = 1.8,cex.sub = 1.8,cex.lab = 1.8,mar = c(5.1,5.1,2.1,2.1),mgp = c(3,1.3,0))
hist(wc.sim[4,],breaks = 2,xlab = "observed wavelet coefficients",xlim = c(-25,25),ylim = c(0,1000),col = grey(0.3),main = "")
hist(wc.sim[10,],breaks = 40,add = TRUE,col = rgb(0,1,0,0.5))
legend('topright',legend = c("wavelet coefficients,scale 1","wavelet coefficients,scale 7"),cex = 1.5,fill = c(rgb(0,1,0,0.5),grey(0.1)))
dev.off()

pdf("paper/simple_eg_3.pdf",height = 6,width = 8)
par(cex.axis = 1.8,cex.sub = 1.8,cex.lab = 1.8,mar = c(5.1,5.1,2.1,2.1),mgp = c(3,1.3,0))
plot(wc.sim[10,],wc.sim.shrunk[[9]]$PosteriorMean,xlab = "observed wavelet coefficients",ylab = "shrunken wavelet coefficients",xlim = c(-4,4),pch = 20,cex = 0.8,ylim = c(-4,4),col = rgb(0,1,0,0.5))
points(wc.sim[4,],wc.sim.shrunk[[3]]$PosteriorMean,pch = 20,cex = 0.8,col = grey(0.1))
abline(0,1,lty = 3,col = grey(0.3))
legend('bottomright',legend = c("wavelet coefficients,scale 1","wavelet coefficients,scale 7"),cex = 1.5,pch = c(20,20),col = c(rgb(0,1,0,0.5),grey(0.1)))
dev.off()

pdf("paper/simple_eg_4.pdf",height = 5.8,width = 8)
layout(matrix(c(1,1,1,3,2,0),nrow = 3),widths = c(8,2),heights = c(2,4,2))
par(mar = c(7.1,7.1,2.1,1.1),cex.axis = 2.5,cex.sub = 2.5,cex.lab = 2.5,mgp = c(3,1.5,0))
plot(wc.sim[4,],wc.sim.shrunk[[3]]$PosteriorMean,xlim = c(-1,1),pch = 20,cex = 1.4,xlab = "",ylab = "",ylim = c(-1,1),col = grey(col.bw.3))
title(xlab = "observed wavelet coefficients",ylab = "shrunken wavelet coefficients",line = 4)
legend('bottomright',legend = c("wavelet coefficients,scale 7"),cex = 1.5,pch = c(20,20),col = grey(0.1))
par(mar = c(1,1,0.1,5.1),cex.axis = 2,cex.sub = 3,cex.lab = 2,mgp = c(3,1.3,0))
image.scale(volcano,col = grey(breaks),breaks = breaks.ori,horiz = FALSE,yaxt = "n")
axis(4,at = seq(0,0.25,0.05),las = 2)
par(mar = c(0.1,0.1,0.1 ,0.1))
plot(1,1,xlim = c(0,1),ylim = c(0,1),type = "n",axes = FALSE)
legend(x = -0.25,y = 0.17,legend = "standard error",cex = 1.5,bty = "n")
dev.off()

pdf("paper/simple_eg_5.pdf",height = 6,width = 8)
par(cex.axis = 1.8,cex.sub = 1.8,cex.lab = 1.8,mar = c(6.1,6.1,2.1,2.1),mgp = c(3,1.3,0))
plot(mu.sp,type = "l",lty = 2,ylim = c(-0.05,0.85),xlab = "",ylab = "")
title(xlab = "position",ylab = "reconstructed mean",line = 4)
lines(mu.ti,col = "blue")
lines(mu.smash,col = "#FF8C00")
legend(x = 500,y = 0.9,legend = c("SMASH","TI-thresh (RMAD)"),col = c("#FF8C00","blue"),lty = c(1,1),cex = 1.6,pt.cex = 0.5,bty = "n")
dev.off()

