# Illustrate SMASH, TI thresholding and EbayesThresh on a Gaussian
# data set simulated with heteroskedastic noise.

# Load the packages and function definitions used in the analysis below.
library(smashr)
library(EbayesThresh)
library(wavethresh)
source("../code/signals.R")
source("../dsc/code/methods/ebayesthresh.wrapper.R")

# Script parameters.
mean.fn <- cor.fn      # spike.fn    
var.fn  <- doppler.fn  # cblocks.fn

# Initialize the sequence of pseudorandom numbers.
set.seed(2)

# Simulate the data set using the "Spikes" mean function and the
# "Clipped Blocks" variance function.
n     <- 1024
t     <- (1:n)/n
mu    <- mean.fn(t,"mean")
sigma <- sqrt(var.fn(t,"var"))
sd    <- sigma/mean(sigma) * sd(mu)/3
x     <- rnorm(n,mu,sd)

# Estimate the mean signal using SMASH, TI Thresholding and EbayesThresh.
sd1 <- sqrt(2/(3*(n - 2)) * sum((x[1:(n - 2)]/2 - x[2:(n - 1)] + x[3:n])^2/2))
mu.smash      <- smash(x,family = "DaubLeAsymm",filter.number = 8)
mu.smash.true <- smash(x,sigma = sd,family = "DaubLeAsymm",filter.number = 8)
mu.smash.homo <- smash(x,sigma = sd1,family = "DaubLeAsymm",filter.number = 8)
mu.ti.rmad    <- ti.thresh(x,method = "rmad",family = "DaubLeAsymm",
                           filter.number = 8)
mu.ti.smash   <- ti.thresh(x,method = "smash",family = "DaubLeAsymm",
                           filter.number = 8)
mu.ti.true    <- ti.thresh(x,sigma = sd,method = "rmad",family = "DaubLeAsymm",
                           filter.number = 8)
mu.ebayes     <- ebayesthresh.wrapper(list(x = x,sig.est = sd1),
                                      list(family = "DaubLeAsymm",
                                           filter.number = 8))

# Plot the simulated data, the mean and variance functions used to
# simulate the data, and the estimates of the mean function produced
# by the various methods.
par(cex.axis = 1,cex.lab = 1.25)
plot(x,type = "p", ylim = c(-0.05,1),xlab = "position",ylab = "",
     col = "gray",pch = 1,cex = 0.75,xaxp = c(0,1024,4),yaxp = c(0,1,4))
lines(mu,lwd = 1.5,col = "black")
lines(mu + 2*sd,col = "black",lty = "dashed",lwd = 1.5)
lines(mu - 2*sd,col = "black",lty = "dashed",lwd = 1.5)
lines(mu.ti.rmad,col = "royalblue",lwd = 1.5)
lines(mu.ebayes,col = "limegreen",lwd = 1.5)
lines(mu.smash,col = "gold",lwd = 1.5)
