# This file defines some of the mean, variance and intensity functions
# used to simulate Gaussian and Poisson data sets.

# This defines the "Spikes" mean function.
spike.fn <- function (t, type = c("mean","var")) {
  type <- match.arg(type)
  f    <- spike.f(t)
  if (type == "mean")
    return((1 + f)/5)
  else if (type == "var")
    return(NULL)
}

# This is used by spikes.fn to define the "Spikes" mean and variance
# functions.
spike.f <- function (t)
  0.75 * exp(-500   * (t - 0.23)^2) +
  1.5  * exp(-2000  * (t - 0.33)^2) +
  3    * exp(-8000  * (t - 0.47)^2) + 
  2.25 * exp(-16000 * (t - 0.69)^2) +
  0.5  * exp(-32000 * (t - 0.83)^2)

# This defines the "Bumps" mean and variance functions.
bumps.fn <- function (t, type = c("mean","var")) {
  type <- match.arg(type)
  f    <- bumps.f(t)
  if (type == "mean")
    return((1 + f)/5)
  else if (type == "var")
    return(1e-05 + f)
}

# This is used by bumps.fn to define the "Bumps" mean and variance
# functions.
bumps.f <- function (t) {
  pos  <- c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 0.76, 0.78, 0.81)
  hgt  <- 2.97/5*c(4, 5, 3, 4, 5, 4.2, 2.1, 4.3, 3.1, 5.1, 4.2)
  wth  <- c(0.005, 0.005, 0.006, 0.01, 0.01, 0.03, 0.01, 0.01, 0.005,
            0.008, 0.005)
  f <- rep(0,length(t))
  for (i in 1:length(pos))
    f <- f + hgt[i]/((1 + (abs(t - pos[i])/wth[i]))^4)
  return(f)
}

# This defines the "Blocks" mean function.
blocks.fn <- function (t, type = c("mean","var")) {
  type <- match.arg(type)
  f    <- blocks.f(t)  
  if (type == "mean")
    return(0.2 + 0.6 * (f - min(f))/max(f - min(f)))
  else if (type == "var")
    return(NULL)
}

# This is used by blocks.fn tto define the "Blocks" mean function.
blocks.f <- function (t) {
  pos  <- c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 0.76, 0.78, 0.81)
  hgt  <- 2.88/5 * c(4, -5, 3, -4, 5, -4.2, 2.1, 4.3, -3.1, 2.1, -4.2)
  f    <- rep(0,length(t))
  for (i in 1:length(pos))
    f <- f + (1 + sign(t - pos[i])) * hgt[i]/2
  return(f)
}

# This defines the "Angles" mean function.
angles.fn <- function (t, type = c("mean","var")) {
  f <- angles.f(t)
  if (type == "mean")
    return((1 + f)/5)
  else if (type == "var")
    return(NULL)
}

# This is used by angles.fn to define the "Angles" mean function.
angles.f <- function (t) {
  s <- ((2 * t + 0.5) * (t <= 0.15)) +
       ((-12 * (t - 0.15) + 0.8) * (t > 0.15 & t <= 0.2)) +
       0.2 * (t > 0.2 & t <= 0.5) +
       ((6 * (t - 0.5) + 0.2) * (t > 0.5 & t <= 0.6)) +
       ((-10 * (t - 0.6) + 0.8) * (t > 0.6 & t <= 0.65)) +
       ((-0.5 * (t - 0.65) + 0.3) * (t > 0.65 & t <= 0.85)) +
       ((2 * (t - 0.85) + 0.2) * (t > 0.85))
  f <- 3/5 * ((5/(max(s) - min(s))) * s - 1.6) - 0.0419569
  return(f)
}

# This defines the "Doppler" mean and variance functions.
doppler.fn <- function (t, type = c("mean","var")) {
  type <- match.arg(type)
  f    <- dop.f(t)
  if (type == "mean") {
    f <- 3/(max(f) - min(f)) * (f - min(f))
    return((1 + f)/5)
  } else if (type == "var") {
    f <- 10 * f
    f <- f - min(f)
    return(1e-05 + 2*f)
  }
}

# This is used by doppler.fn to define the "Doppler" mean and variance
# functions.
dop.f <- function (t)
  sqrt(t * (1 - t)) * sin((2 * pi * 1.05)/(t + 0.05))

# This defines the "Blip" mean function.
blip.fn <- function (t, type = c("mean","var")) {
  type <- match.arg(type)
  f    <- blip.fn(t)
  if (type == "mean")
    return(f)
  else if (type == "var")
    return(NULL)
}

# This is used by blip.fn to define the "Blip" mean function.
blip.f <- function (t)
  (0.32 + 0.6 * t + 0.3 * exp(-100 * (t - 0.3)^2)) * (t >= 0 & t <= 0.8) +
  (-0.28 + 0.6 * t + 0.3 * exp(-100 * (t - 1.3)^2)) * (t > 0.8 & t <= 1)

# This defines the "Corner" mean function.
cor.fn <- function (t, type = c("mean","var")) {
  type <- match.arg(type)
  f    <- cor.f(t)
  if (type == "mean")
    return(f - min(f) + 0.2)
  else if (type == "var")
    return(NULL)
}

# This is used by cor.fn to define the "Corner" mean function.
cor.f <- function (t) {
  f    <- 623.87 * t^3 * (1 - 2 * t) * (t >= 0 & t <= 0.5) +
          187.161 * (0.125 - t^3) * t^4 * (t > 0.5 & t <= 0.8) +
          3708.470441 * (t - 1)^3 * (t > 0.8 & t <= 1)
  f    <- (0.6/(max(f) - min(f))) * f
  return(f)
}

# This defines the constant variance function.
cons.fn = function (t, type = c("var","mean")) {
  type <- match.arg(type)
  f    <- rep(1,length(t))
  if (type == "mean") {
    return(NULL)
  } else if (type == "var") {
    return(f)
  }
}
  
# This defines the "Clipped Blocks" variance function.
cblocks.fn <- function (t, type = c("var","mean")) {
  type <- match.arg(type)
  f    <- cblocks.f(t)
  if (type == "mean")
    return(NULL)
  else if (type == "var")
    return(0.01 + 1 * (f - min(f))/max(f))
}

# This is used by cblocks.fn to define the "Clipped Blocks" variance
# function.
cblocks.f <- function (t) {
  pos <- c(0.1,0.13,0.15,0.23,0.25,0.4,0.44,0.65,0.76,0.78,0.81)
  hgt <- 2.88/5 * c(4,-5,3,-4,5,-4.2,2.1,4.3,-3.1,2.1,-4.2)
  f   <- rep(0,length(t))
  for (i in 1:length(pos))
    f <- f + (1 + sign(t - pos[i])) * (hgt[i]/2)
  f[f < 0] <- 0
  return(f)
}

# This defines the "triple exponential" variance function.
texp.fn <- function (t, type = c("var","mean")) {
  type <- match.arg(type)
  f    <- fexp.f(t)
  if (type == "mean")
    return(NULL)
  else if (type == "var")
    return(f)
}

# This is used by texp.fn to define the "triple exponential" variance
# function.
texp.f <- function (t)
  1e-04 + 4*(exp(-550 * (t - 0.2)^2) +
             exp(-200 * (t - 0.5)^2) +
             exp(-950 * (t - 0.8)^2))
