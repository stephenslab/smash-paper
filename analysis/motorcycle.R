# An illustration of "smoothing via adaptive shrinkage" applied to the
# Motorcycle Accident data set. Compare the first plot against the
# figure in the "Motorcycle Acceleration Data" section of the paper.

# SET UP ENVIRONMENT
# ------------------
# Load the MASS, scales, wavethresh, EbayesThresh and smashr
# packages. The MASS package is loaded only for the Motorcycle
# Accident data. Some additional functions are defined in
# functions_motorcycle.R.
library(MASS)
library(lattice)
library(smashr)
library(scales)
library(wavethresh)
library(EbayesThresh)
source("functions_motorcycle.R")

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

invisible(readline(prompt = "Press [enter] to continue analysis... "))

# PREPARE DATA FOR SMASH WITH HOMOGENEOUS VARIANCES
# -------------------------------------------------
# Load the motorcycle data for Empirical Bayes thresholding.
cat("Preparing Motorcycle data for Empirical Bayes thresholding.\n")
x.ini.ti.cons.mc = sort(mcycle$times)
y.ini.ti.cons.mc = mcycle$accel[order(mcycle$times)]

# Load the motorcycle data once more for the second smash analysis.
cat("Preparing Motorcycle data for smash analysis (with equal variances).\n")
x.ini.cons.mc = sort(mcycle$times)
y.ini.cons.mc = mcycle$accel[order(mcycle$times)]

# Apply Empirical Bayes thresholding to the Motorcycle Accident data set.
cat("Running Empirical Bayes thresholding on Motorcycle Accident data set.\n")
res.ti.cons.mc = tithresh.cons.wrapper(x.ini.mc, y.ini.mc)

# Apply smash (with equal variances) to the Motorcycle Accident data set.
cat("Running smash (with equal variances) on Motorcycle Accident data set.\n")
res.cons.mc = smash.cons.wrapper(x.ini.mc, y.ini.mc)

# PREPARE DATA FOR EMPIRICAL BAYES WITH SMASH-ESTIMATED VARIANCES
# ---------------------------------------------------------------
# Load the motorcycle data for Empirical Bayes thresholding (with
# smash-estimated variances).
cat("Preparing data for EbayesThresh (with smash-estimated variances).\n")
x.ini.ti.mc = sort(mcycle$times)
y.ini.ti.mc = mcycle$accel[order(mcycle$times)]

# Apply Empirical Bayes thresholding to the Motorcycle Accident data,
# with variances estimated by smash.
cat(paste("Running EbayesThresh (with smash-estimated variances) on",
          "Motorcycle data.\n"))
res.ti.mc = tithresh.wrapper(x.ini.mc, y.ini.mc)

# SUMMARIZE RESULTS FROM ALL METHODS
# ----------------------------------
# Create another plot showing the Motorcycle Accident data and the
# estimates provided by the four methods.
par(cex.axis = 1.5, cex.sub = 1.8, cex.lab = 1.5,
    mar = c(5.1, 5.1, 2.1, 2.1), mgp = c(3, 1.3, 0))
plot(res.mc$x, res.mc$mu.est, type = 'l',
     ylim = c(min(res.mc$y - 2 * sqrt(res.mc$var.est)),
              max(res.mc$y + 2 * sqrt(res.mc$var.est))),
     xlab = "time (ms)", ylab = "acceleration (g)", lwd = 2)
lines(res.mc$x, res.mc$mu.est + 2 * sqrt(res.mc$var.est), col = 2,
      lty = 5, lwd = 1.8)
lines(res.mc$x, res.mc$mu.est - 2 * sqrt(res.mc$var.est), col = 2,
      lty = 5, lwd = 1.8)
lines(res.ti.cons.mc$x, res.ti.cons.mc$mu.est, lwd = 2, col = 4)
lines(res.cons.mc$x, res.cons.mc$mu.est, lwd = 2, col = 'brown')
lines(res.ti.mc$x, res.ti.mc$mu.est, lwd = 2, col = 'darkgreen')
points(res.mc$x, res.mc$y, cex = 0.7, pch = 16, col = alpha("black", 0.5))
legend('bottomright',
       legend = c('SMASH', 'SMASH, homo', 'TI-thresh, SMASH est',
                  'TI-EBthresh, homo'),
       lty = c(1, 1, 1, 1), lwd = c(2, 2, 2, 2),
       col = c(1, 'brown', 'darkgreen', 4))
