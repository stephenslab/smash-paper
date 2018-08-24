# This file defines some of the mean and variance functions used to
# simulate data sets.

cblocks.fn <- function(t, type) {
  pos = c(0.1,0.13,0.15,0.23,0.25,0.4,0.44,0.65,0.76,0.78,0.81)
  hgt = 2.88/5 * c(4,-5,3,-4,5,-4.2,2.1,4.3,-3.1,2.1,-4.2)
  fn  = rep(0,length(t))
  for (j in 1:length(pos))
    fn = fn + (1 + sign(t - pos[j])) * (hgt[j]/2)
  fn[fn < 0] = 0
  if (type == "mean")
    return(NULL)
  else if (type == "var")
    return(0.01 + 1 * (fn - min(fn))/max(fn))
}

doppler.fn <- function(t, type) {
  dop.f = function(x) sqrt(x * (1 - x)) * sin((2 * pi * 1.05)/(x + 0.05))
  fn    = dop.f(t)
  if (type == "mean") {
    fn = 3/(max(fn) - min(fn)) * (fn - min(fn))
    return((1 + fn)/5)
  } else if (type == "var") {
    fn = 10 * fn
    fn = fn - min(fn)
    return(1e-05 + 2 * fn)
  }
}

cor.fn <- function(t, type) {
  fn = 623.87 * t^3 * (1 - 2 * t) * (t >= 0 & t <= 0.5) +
       187.161 * (0.125 - t^3) * t^4 * (t > 0.5 & t <= 0.8) +
       3708.470441 * (t - 1)^3 * (t > 0.8 & t <= 1)
  fn = (0.6/(max(fn) - min(fn))) * fn
  if (type == "mean")
    return(fn - min(fn) + 0.2)
  else if (type == "var")
    return(NULL)
}

spikes.fn <- function(t, type) {
  spike.f = function(x)
    (0.75 * exp(-500   * (x - 0.23)^2) +
     1.5  * exp(-2000  * (x - 0.33)^2) +
     3    * exp(-8000  * (x - 0.47)^2) +
     2.25 * exp(-16000 * (x - 0.69)^2) +
     0.5  * exp(-32000 * (x - 0.83)^2))
  fn = spike.f(t)
  if (type == "mean")
    return((1 + fn)/5)
  else if (type == "var")
    return(NULL)
}

