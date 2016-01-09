library(smash)
library(MASS)
library(lattice)
library(scales)

#Wrapper function to return estimated mean and variance function given raw data
smash.wrapper = function(x.ini, y.ini){
  #Take the median of observations with repeated x values
  x = unique(x.ini)
  y = 0
  for(i in 1:length(x)){
    y[i] = median(y.ini[x.ini == x[i]])
  }
  
  #Mirror the data twice to make it periodic and a power of 2
  y.exp = c(y, y[length(y):(2*length(y) - 2^ceiling(log2(length(y))) + 1)])
  y.final = c(y.exp, y.exp[length(y.exp):1])
  
  #Run SMASH
  y.est = ashsmooth.gaus(y.final, v.est=TRUE, joint = TRUE, family = "DaubLeAsymm", filter.number = 8)
  
  y.mu.est = y.est$mu.res[1:length(y)]
  y.var.est = y.est$var.res[1:length(y)]
  return(list(x = x, y = y, mu.est = y.mu.est, var.est = y.var.est))
}



#Obtain motorcycle data
x.ini.mc = sort(mcycle$times)
y.ini.mc = mcycle$accel[order(mcycle$times)]
#Get estimates
res.mc = smash.wrapper(x.ini.mc, y.ini.mc)

#Obtain ethanol data
x.ini.eth = sort(ethanol$E)
y.ini.eth = ethanol$NOx[order(ethanol$E)]
#Get estimates
res.eth = smash.wrapper(x.ini.eth, y.ini.eth)



pdf("paper/motorcycle.pdf", height = 6, width = 8)
par(cex.axis = 1.3, cex.sub = 1.8, cex.lab = 1.5, mar = c(5.1, 5.1, 2.1, 2.1), mgp = c(3, 1.3, 0))
plot(res.mc$x, res.mc$mu.est, type = 'l', ylim = c(min(res.mc$y - 2 * sqrt(res.mc$var.est)), max(res.mc$y + 2 * sqrt(res.mc$var.est))), xlab = "time (ms)", ylab = "acceleration (g)", lwd = 1.7)
lines(res.mc$x, res.mc$mu.est + 2 * sqrt(res.mc$var.est), col = 2, lty = 5, lwd = 1.8)
lines(res.mc$x, res.mc$mu.est - 2 * sqrt(res.mc$var.est), col = 2, lty = 5, lwd = 1.8)
points(res.mc$x, res.mc$y, cex = 0.7, pch = 16, col = alpha("black", 0.5))
dev.off()


pdf("paper/motorcycle.pdf", height = 6, width = 8)
par(cex.axis = 1.5, cex.sub = 1.8, cex.lab = 1.5, mar = c(5.1, 5.1, 2.1, 2.1), mgp = c(3, 1.3, 0))
plot(res.mc$x, res.mc$mu.est, type = 'l', ylim = c(min(res.mc$y - 2 * sqrt(res.mc$var.est)), max(res.mc$y + 2 * sqrt(res.mc$var.est))), xlab = "time (ms)", ylab = "acceleration (g)", lwd = 1.7)
lines(res.mc$x, res.mc$mu.est + 2 * sqrt(res.mc$var.est), col = 2, lty = 5, lwd = 1.8)
lines(res.mc$x, res.mc$mu.est - 2 * sqrt(res.mc$var.est), col = 2, lty = 5, lwd = 1.8)
points(res.mc$x, res.mc$y, cex = 0.7, pch = 16, col = alpha("black", 0.5))
dev.off()



# pdf("paper/ethanol.pdf", height = 6, width = 8)
# par(cex.axis = 1.5, cex.sub = 1.8, cex.lab = 1.5, mar = c(5.1, 5.1, 2.1, 2.1), mgp = c(3, 1.3, 0))
# plot(res.eth$x, res.eth$mu.est, type = 'l', ylim = c(min(res.eth$y - 2 * sqrt(res.eth$var.est)), max(res.eth$y + 2 * sqrt(res.eth$var.est))), xlab = "Equivalence Ratio", ylab = "NOx", lwd = 1.7)
# lines(res.eth$x, res.eth$mu.est + 2 * sqrt(res.eth$var.est), col = 2, lty = 5, lwd = 1.8)
# lines(res.eth$x, res.eth$mu.est - 2 * sqrt(res.eth$var.est), col = 2, lty = 5, lwd = 1.8)
# points(res.eth$x, res.eth$y, cex = 0.7, pch = 16, col = alpha("black", 0.5))
# dev.off()
