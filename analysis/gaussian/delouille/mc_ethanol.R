# TO DO: Summarize this script.

# SET UP ENVIRONMENT
# ------------------
# Load the MASS, lattice and smashr packages. The MASS and lattice
# packages are loaded only for the Motorcycle Accident and Ethanol
# data sets.
library(MASS)
library(lattice)
library(smashr)
library(scales)
# library(EbayesThresh)
source("functions_mc.R")

# PREPARE DATA
# ------------
# Load the motorcycle data.
cat("Preparing Motorcycle Accident data.\n")
data(mcycle)
x.ini.mc  = sort(mcycle$times)
y.ini.mc  = mcycle$accel[order(mcycle$times)]

# Load the Ethanol data.
cat("Preparing Ethanol data.\n")
data(ethanol)
x.ini.eth = sort(ethanol$E)
y.ini.eth = ethanol$NOx[order(ethanol$E)]

# RUN SMASH
# ---------
# Apply smash to the Motorcycle Accident data set.
cat("Running smash on Motorcycle Accident data set.\n")
res.mc = smash.wrapper(x.ini.mc, y.ini.mc)

# Apply smash to the Ethanol data set.
cat("Running smash on Ethanol data set.\n")
res.eth = smash.wrapper(x.ini.eth, y.ini.eth)

# SUMMARIZE RESULTS
# -----------------
# Create a plot showing the Motorcycle Accident data and the smash
# estimates (with the dashed red lines showing the confidence
# intervals).
par(cex.axis = 1.5, cex.sub = 1.8, cex.lab = 1.5,
    mar = c(5.1, 5.1, 2.1, 2.1), mgp = c(3, 1.3, 0))
plot(res.mc$x, res.mc$mu.est, type = 'l',
     ylim = c(min(res.mc$y - 2 * sqrt(res.mc$var.est)),
              max(res.mc$y + 2 * sqrt(res.mc$var.est))),
     xlab = "time (ms)", ylab = "acceleration (g)", lwd = 1.7)
lines(res.mc$x, res.mc$mu.est + 2 * sqrt(res.mc$var.est), col = 2,
      lty = 5, lwd = 1.8)
lines(res.mc$x, res.mc$mu.est - 2 * sqrt(res.mc$var.est), col = 2,
      lty = 5, lwd = 1.8)
points(res.mc$x, res.mc$y, cex = 0.7, pch = 16, col = alpha("black", 0.5))

stop()

#Obtain motorcycle data
x.ini.ti.cons.mc = sort(mcycle$times)
y.ini.ti.cons.mc = mcycle$accel[order(mcycle$times)]
#Get estimates
res.ti.cons.mc = tithresh.cons.wrapper(x.ini.mc, y.ini.mc)



#Wrapper function to return estimated mean given raw data using SMASH, but assuming constant variance
smash.cons.wrapper = function(x.ini, y.ini){
  #Take the median of observations with repeated x values
  x = unique(x.ini)
  y = 0
  for(i in 1:length(x)){
    y[i] = median(y.ini[x.ini == x[i]])
  }
  
  #Mirror the data twice to make it periodic and a power of 2
  y.exp = c(y, y[length(y):(2*length(y) - 2^ceiling(log2(length(y))) + 1)])
  y.final = c(y.exp, y.exp[length(y.exp):1])
  
  y.noise = sig.est.func(y.final, length(y.final))
  
  #Run SMASH
  y.est = smash.gaus(y.final, sigma = y.noise, v.est=FALSE, family = "DaubLeAsymm", filter.number = 8)
  
  y.mu.est = y.est[1:length(y)]
  return(list(x = x, y = y, mu.est = y.mu.est))
}

#Obtain motorcycle data
x.ini.cons.mc = sort(mcycle$times)
y.ini.cons.mc = mcycle$accel[order(mcycle$times)]
#Get estimates
res.cons.mc = smash.cons.wrapper(x.ini.mc, y.ini.mc)


#Wrapper function to return estimated mean given raw data for TI thresholding with variance estimated using SMASH
tithresh.wrapper = function(x.ini, y.ini){
  #Take the median of observations with repeated x values
  x = unique(x.ini)
  y = 0
  for(i in 1:length(x)){
    y[i] = median(y.ini[x.ini == x[i]])
  }
  
  #Mirror the data twice to make it periodic and a power of 2
  y.exp = c(y, y[length(y):(2*length(y) - 2^ceiling(log2(length(y))) + 1)])
  y.final = c(y.exp, y.exp[length(y.exp):1])
  
  #Run TI thresholding with emprical bayes thresholding
  y.est = ti.thresh(y.final, method = 'smash', min.level = 2)
  
  y.mu.est = y.est[1:length(y)]
  return(list(x = x, y = y, mu.est = y.mu.est))
}

#Obtain motorcycle data
x.ini.ti.mc = sort(mcycle$times)
y.ini.ti.mc = mcycle$accel[order(mcycle$times)]
#Get estimates
res.ti.mc = tithresh.wrapper(x.ini.mc, y.ini.mc)

pdf("paper/motorcycle_2.pdf", height = 6, width = 8)
par(cex.axis = 1.5, cex.sub = 1.8, cex.lab = 1.5, mar = c(5.1, 5.1, 2.1, 2.1), mgp = c(3, 1.3, 0))
plot(res.mc$x, res.mc$mu.est, type = 'l', ylim = c(min(res.mc$y - 2 * sqrt(res.mc$var.est)), max(res.mc$y + 2 * sqrt(res.mc$var.est))), xlab = "time (ms)", ylab = "acceleration (g)", lwd = 2)
lines(res.mc$x, res.mc$mu.est + 2 * sqrt(res.mc$var.est), col = 2, lty = 5, lwd = 1.8)
lines(res.mc$x, res.mc$mu.est - 2 * sqrt(res.mc$var.est), col = 2, lty = 5, lwd = 1.8)
lines(res.ti.cons.mc$x, res.ti.cons.mc$mu.est, lwd = 2, col = 4)
lines(res.cons.mc$x, res.cons.mc$mu.est, lwd = 2, col = 'brown')
lines(res.ti.mc$x, res.ti.mc$mu.est, lwd = 2, col = 'darkgreen')
points(res.mc$x, res.mc$y, cex = 0.7, pch = 16, col = alpha("black", 0.5))
legend('bottomright', legend = c('SMASH', 'SMASH, homo', 'TI-thresh, SMASH est', 'TI-EBthresh, homo'), lty = c(1, 1, 1, 1), lwd = c(2, 2, 2, 2), col = c(1, 'brown', 'darkgreen', 4))
dev.off()
