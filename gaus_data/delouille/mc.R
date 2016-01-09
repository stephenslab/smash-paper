library(smash)
library(MASS)

#wrapper function to return estimated mean and variance function given raw data
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



#Obtain initial data
x.ini = sort(mcycle$times)
y.ini = mcycle$accel[order(mcycle$times)]

x.ini = sort(ethanol$E)
y.ini = ethanol$NOx[order(ethanol$E)]


plot(x, y, ylim = c(-150,100))
lines(x, y.mu.est, col = 2)
lines(x, y.mu.est + 2*sqrt(y.var.est), col = 4)
lines(x, y.mu.est - 2*sqrt(y.var.est), col = 4)




