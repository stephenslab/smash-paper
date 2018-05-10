# Note sure what this script is for? It seems that it is incomplete.

library(smashr)
library(MASS)
library(lattice)

# Obtain initial data.
x.ini = sort(mcycle$times)
y.ini = mcycle$accel[order(mcycle$times)]

x.ini = sort(ethanol$E)
y.ini = ethanol$NOx[order(ethanol$E)]

plot(x, y, ylim = c(-150,100))
lines(x, y.mu.est, col = 2)
lines(x, y.mu.est + 2*sqrt(y.var.est), col = 4)
lines(x, y.mu.est - 2*sqrt(y.var.est), col = 4)




