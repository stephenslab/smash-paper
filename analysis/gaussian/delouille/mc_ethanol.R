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
library(wavethresh)
library(EbayesThresh)
source("functions_mc.R")

# PREPARE DATA FOR SMASH
# ----------------------
# Load the motorcycle data.
cat("Preparing Motorcycle Accident data for smash analysis.\n")
data(mcycle)
x.ini.mc = sort(mcycle$times)
y.ini.mc = mcycle$accel[order(mcycle$times)]

# RUN SMASH
# ---------
# Apply smash to the Motorcycle Accident data set.
cat("Running smash on Motorcycle Accident data set.\n")
res.mc = smash.wrapper(x.ini.mc, y.ini.mc)

# SUMMARIZE RESULTS OF SMASH ANALYSIS
# -----------------------------------
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

# PREPARE DATA FOR SMASH WITH HOMOGENEOUS VARIANCES
# -------------------------------------------------
# Load the motorcycle data for ...
cat("Preparing Motorcycle data for...\n")
x.ini.ti.cons.mc = sort(mcycle$times)
y.ini.ti.cons.mc = mcycle$accel[order(mcycle$times)]

# Load the motorcycle data once more for the second smash analysis.
cat("Preparing Motorcycle data for smash analysis (with equal variances).\n")
x.ini.cons.mc = sort(mcycle$times)
y.ini.cons.mc = mcycle$accel[order(mcycle$times)]

# Apply smash (with equal variances) to the Motorcycle Accident data set.
cat("Running smash (with equal variances) on Motorcycle Accident data set.\n")
res.cons.mc = smash.cons.wrapper(x.ini.mc, y.ini.mc)

# Apply ... to the Motorcycle Accident data set.
cat("Running ... on Motorcycle Accident data set.\n")
res.ti.cons.mc = tithresh.cons.wrapper(x.ini.mc, y.ini.mc)

stop()

# Obtain motorcycle data
x.ini.ti.mc = sort(mcycle$times)
y.ini.ti.mc = mcycle$accel[order(mcycle$times)]

#Get estimates
res.ti.mc = tithresh.wrapper(x.ini.mc, y.ini.mc)

par(cex.axis = 1.5, cex.sub = 1.8, cex.lab = 1.5, mar = c(5.1, 5.1, 2.1, 2.1), mgp = c(3, 1.3, 0))
plot(res.mc$x, res.mc$mu.est, type = 'l', ylim = c(min(res.mc$y - 2 * sqrt(res.mc$var.est)), max(res.mc$y + 2 * sqrt(res.mc$var.est))), xlab = "time (ms)", ylab = "acceleration (g)", lwd = 2)
lines(res.mc$x, res.mc$mu.est + 2 * sqrt(res.mc$var.est), col = 2, lty = 5, lwd = 1.8)
lines(res.mc$x, res.mc$mu.est - 2 * sqrt(res.mc$var.est), col = 2, lty = 5, lwd = 1.8)
lines(res.ti.cons.mc$x, res.ti.cons.mc$mu.est, lwd = 2, col = 4)
lines(res.cons.mc$x, res.cons.mc$mu.est, lwd = 2, col = 'brown')
lines(res.ti.mc$x, res.ti.mc$mu.est, lwd = 2, col = 'darkgreen')
points(res.mc$x, res.mc$y, cex = 0.7, pch = 16, col = alpha("black", 0.5))
legend('bottomright', legend = c('SMASH', 'SMASH, homo', 'TI-thresh, SMASH est', 'TI-EBthresh, homo'), lty = c(1, 1, 1, 1), lwd = c(2, 2, 2, 2), col = c(1, 'brown', 'darkgreen', 4))

