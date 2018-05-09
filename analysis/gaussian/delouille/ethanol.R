# An illustration of "smoothing via adaptive shrinkage" applied to the
# Ethanol data set.

# SET UP ENVIRONMENT
# ------------------
# Load the lattice and smashr packages. The lattice package is loaded
# only for the Ethanol dat aset.
library(smashr)
library(lattice)

# PREPARE DATA
# ------------
# Prepare the Ethanol data for analysis with smash.
cat("Preparing data.\n")
data(ethanol)
x.ini = sort(ethanol$E)
y.ini = ethanol$NOx[order(ethanol$E)]

x = unique(x.ini)
y = 0
for(i in 1:length(x))
  y[i] = median(y.ini[x.ini == x[i]])

y.exp = c(y,y[length(y):(2*length(y)-128+1)])
y.final = c(y.exp, y.exp[length(y.exp):1])

# RUN SMASH
# ---------
# Apply smash to the Ethanol data set.
cat("Running smash.\n")
y.est     = smash.gaus(y.final)
y.est     = y.est[1:length(y)]
y.var.est = smash.gaus(y.final,v.est = TRUE)
y.var.est = y.var.est[1:length(y)]

# SUMMARIZE RESULTS
# -----------------
plot(x,y,ylim = c(-1,5))
lines(x,y.est,col = 2)

