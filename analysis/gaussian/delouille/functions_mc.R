# This file contains functions used by the mc.R and mc_ethanol.R scripts.

# Wrapper function to return estimated mean and variance function
# given raw data.
smash.wrapper = function (x.ini, y.ini) {
    
  # Take the median of the observations with repeated x values.
  x = unique(x.ini)
  y = 0
  for(i in 1:length(x))
    y[i] = median(y.ini[x.ini == x[i]])
  
  # Mirror the data twice to make it periodic and a power of 2.
  y.exp = c(y, y[length(y):(2*length(y) - 2^ceiling(log2(length(y))) + 1)])
  y.final = c(y.exp, y.exp[length(y.exp):1])
  
  # Run SMASH.
  y.est = smash.gaus(y.final, v.est = TRUE, joint = TRUE,
                     family = "DaubLeAsymm", filter.number = 8)
  
  y.mu.est = y.est$mu.res[1:length(y)]
  y.var.est = y.est$var.res[1:length(y)]
  return(list(x = x, y = y, mu.est = y.mu.est, var.est = y.var.est))
}

waveti.ebayes = function (x, filter.number = 10, family = "DaubLeAsymm",
                          min.level = 3, noise.level) {
  n = length(x)
  J = log2(n)
  x.w <- wd(x, filter.number, family, type = "station")
  for (j in min.level:(J - 1)) {
    x.pm = ebayesthresh(accessD(x.w, j), sdev=noise.level)
    x.w = putD(x.w, j, x.pm)
  }
  mu.est = AvBasis(convert(x.w))
  return(mu.est)
}

sig.est.func = function (x, n)
  sqrt(2/(3*(n - 2)) *
  sum((1/2 * x[1:(n - 2)]-x[2:(n - 1)] + 1/2 * x[3:n])^2))

# Wrapper function to return estimated mean given raw data for TI
# thresholding with Empirical Bayes shrinkage.
tithresh.cons.wrapper = function (x.ini, y.ini) {
    
  # Take the median of observations with repeated x values.
  x = unique(x.ini)
  y = 0
  for(i in 1:length(x)){
    y[i] = median(y.ini[x.ini == x[i]])
  }
  
  # Mirror the data twice to make it periodic and a power of 2.
  y.exp = c(y, y[length(y):(2*length(y) - 2^ceiling(log2(length(y))) + 1)])
  y.final = c(y.exp, y.exp[length(y.exp):1])
  
  y.noise = sig.est.func(y.final, length(y.final))
  
  # Run TI thresholding with emprical bayes thresholding.
  y.est = waveti.ebayes(y.final, min.level = 2, noise.level = y.noise)
  
  y.mu.est = y.est[1:length(y)]
  return(list(x = x, y = y, mu.est = y.mu.est))
}

