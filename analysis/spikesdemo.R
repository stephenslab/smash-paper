# This script illustrates key features of the smash method on a small,
# simulated data set.
source("../code/plotting.R")

# SET UP ENVIROMENT
# -----------------
# Load the ashr, smashr, scales, ggplot2 and cowplot packages.
library(ashr)
library(smashr)
library(scales)
library(ggplot2)
suppressMessages(library(cowplot))

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

# Plot the simulated data set.
par(cex.axis = 0.8,cex.sub = 1,cex.lab = 1)
plot(mu.sp,type = 'l',ylim = c(-0.05,1),xlab = "position",
     ylab = "",lwd = 1.7,xaxp = c(0,1024,4),yaxp = c(0,1,4))
lines(mu.sp + 2*sig.cb,col = "darkorange",lty = 5,lwd = 1.8)
lines(mu.sp - 2*sig.cb,col = "darkorange",lty = 5,lwd = 1.8)
points(x.sim,cex = 0.7,pch = 16,col = "darkblue")
invisible(readline(prompt = "Press [enter] to continue analysis... "))

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
# Plot the distribution of observed wavelet coefficients.
par(cex.axis = 0.8,cex.lab = 0.8)
hist(wc.sim[4,],breaks = 2,xlab = "observed wavelet coefficients",
     xlim = c(-25,25),ylim = c(0,600),col = "darkblue",xaxp = c(-25,25,10),
     yaxp = c(0,600,6),main = "")
hist(wc.sim[10,],breaks = 40,add = TRUE,col = "darkorange")
invisible(readline(prompt = "Press [enter] to continue analysis... "))

# Plot the observed wavelet coefficients (at scales 1 and 7 only)
# vs. the "shrunken" wavelet coefficients estimated by adaptive
# shrinkage.
par(cex.axis = 0.8,cex.lab = 0.8)
plot(c(),c(),xlab = "observed wavelet coefficients",
     ylab = "shrunken wavelet coefficients",
     xlim = c(-2.5,2.5),ylim = c(-2.5,2.5))
abline(0,1,lty = 1,col = "gray",lwd = 1)
points(wc.sim[10,],wc.sim.shrunk[[9]]$PosteriorMean,pch = 20,cex = 0.6,
       col = "darkorange")
points(wc.sim[4,],wc.sim.shrunk[[3]]$PosteriorMean,pch = 20,cex = 0.6,
       col = "darkblue")
invisible(readline(prompt = "Press [enter] to continue analysis... "))

# Plot the observed wavelet coefficients (at scale 7 only) vs. the
# "shrunken" wavelet coefficients estimated by adaptive shrinkage, and
# how how the amount of shrinkage depends on the standard error (s.e.)
# in the observations.
wc.sig.3 <- 1/wc.pres[[3]]
p <- ggplot(data.frame(observed = wc.sim[4,],
                       shrunken = wc.sim.shrunk[[3]]$PosteriorMean,
                       se       = wc.sig.3),
            aes(x = observed,y = shrunken,col = se)) +
  geom_point() +
  xlim(c(-1,1)) +
  ylim(c(-1,1)) +
  scale_color_gradientn(colors = c("deepskyblue","darkblue")) +
  theme_cowplot()
print(p)
invisible(readline(prompt = "Press [enter] to continue analysis... "))

# Plot the ground-truth signal (mean function) and the signals
# recovered by the TI thresholding and smash methods.
par(cex.axis = 0.8)
plot(mu.sp,type = "l",col = "black",lwd = 2,xlab = "position",ylab = "",
     ylim = c(-0.05,1),xaxp = c(0,1024,4),yaxp = c(0,1,4))
lines(mu.ti,col = "dodgerblue",lwd = 1.5)
lines(mu.smash,col = "orangered",lwd = 1.5)
cat("Demo is over.\n")

# SESSIONINFO
# -----------
# Print out information on the computing environment, including the
# version of R and the attached R packages, used to generate these
# results.
print(sessionInfo())
