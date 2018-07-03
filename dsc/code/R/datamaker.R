spikes.fn <- function (t, type) {
  spike.f <- function (x) (0.75*exp(-500*(x - 0.23)^2) +
                          1.5*exp(-2000*(x - 0.33)^2) +
                          3*exp(-8000*(x - 0.47)^2) +
                          2.25*exp(-16000*(x - 0.69)^2) +
                          0.5*exp(-32000*(x - 0.83)^2))
  fn <- spike.f(t)
  if (type == "mean")
    return((1 + fn)/5)
  else if (type == "var")
      return(NULL)
  return(NULL)
}

bumps.fn <- function (t, type) {
  pos <- c(0.1,0.13,0.15,0.23,0.25,0.4,0.44,0.65,0.76,0.78,0.81)
  hgt <- 2.97/5 * c(4,5,3,4,5,4.2,2.1,4.3,3.1,5.1,4.2)
  wth <- c(0.005,0.005,0.006,0.01,0.01,0.03,0.01,0.01,0.005,0.008,0.005)
  fn  <- rep(0,length(t))
  for (j in 1:length(pos))
    fn <- fn + hgt[j]/((1 + (abs(t - pos[j])/wth[j]))^4)
  if (type == "mean")
    return((1 + fn)/5)
  else if (type == "var")
    return(1e-05 + fn)
  return(NULL)
}

blocks.fn <- function (t, type) {
  pos <- c(0.1,0.13,0.15,0.23,0.25,0.4,0.44,0.65,0.76,0.78,0.81)
  hgt <- 2.88/5 * c(4,-5,3,-4,5,-4.2,2.1,4.3,-3.1,2.1,-4.2)
  fn  <- rep(0,length(t))
  for (j in 1:length(pos))
    fn <- fn + (1 + sign(t - pos[j])) * (hgt[j]/2)
  if (type == "mean")
    return(0.2 + 0.6*(fn - min(fn))/max(fn - min(fn)))
  else if (type == "var")
    return(NULL)
  return(NULL)
}

angles.fn <- function (t, type) {
  sig <- ((2 * t + 0.5) * (t <= 0.15)) +
         ((-12 * (t - 0.15) + 0.8) * (t > 0.15 & t <= 0.2)) +
         0.2 * (t > 0.2 & t <= 0.5) +
         ((6 * (t - 0.5) + 0.2) * (t > 0.5 & t <= 0.6)) +
         ((-10 * (t - 0.6) + 0.8) * (t > 0.6 & t <= 0.65)) +
         ((-0.5 * (t - 0.65) + 0.3) * (t > 0.65 & t <= 0.85)) +
         ((2 * (t - 0.85) + 0.2) * (t > 0.85))
  fn <- 3/5 * ((5/(max(sig) - min(sig))) * sig - 1.6) - 0.0419569
  if (type == "mean")
    return((1 + fn)/5)
  else if (type == "var")
    return(NULL)
  return(NULL)
}

doppler.fn <- function (t, type) {
  dop.f <- function (x) sqrt(x * (1 - x)) * sin((2 * pi * 1.05)/(x + 0.05))
  fn    <- dop.f(t)
  if (type == "mean") {
    fn <- 3/(max(fn) - min(fn)) * (fn - min(fn))
    return((1 + fn)/5)
  } else if (type == "var") {
    fn <- 10 * fn
    fn <- fn - min(fn)
    return(1e-05 + 2 * fn)
  }
  return(NULL)
}

blip.fn <- function (t, type) {
  fn <- (0.32 + 0.6 * t + 0.3 * exp(-100 * (t - 0.3)^2)) * (t >= 0 & t <= 0.8) + (-0.28 + 0.6 * t + 0.3 * exp(-100 * (t - 1.3)^2)) * (t > 0.8 & t <= 1)
  if (type == "mean")
    return(fn)
  else if (type == "var")
    return(NULL)
  return(NULL)
}

cor.fn <- function(t, type) {
  fn <- 623.87 * t^3 * (1 - 2 * t) * (t >= 0 & t <= 0.5) +
        187.161 * (0.125 - t^3) * t^4 * (t > 0.5 & t <= 0.8) +
        3708.470441 * (t - 1)^3 * (t > 0.8 & t <= 1)
  fn <- (0.6/(max(fn) - min(fn))) * fn
  if (type == "mean") {
    return(fn - min(fn) + 0.2)
  } else if (type == "var") {
    return(NULL)
  }
  return(NULL)
}

cblocks.fn <- function (t, type) {
  pos <- c(0.1,0.13,0.15,0.23,0.25,0.4,0.44,0.65,0.76,0.78,0.81)
  hgt <- 2.88/5 * c(4,(-5),3,(-4),5,(-4.2),2.1,4.3,(-3.1),2.1,(-4.2))
  fn  <- rep(0,length(t))
  for (j in 1:length(pos))
    fn <- fn + (1 + sign(t - pos[j])) * (hgt[j]/2)
  fn[fn < 0] <- 0
  if (type == "mean") {
    return(NULL)
  } else if (type == "var") {
    return(1e-05 + 1 * (fn - min(fn))/max(fn))
  }
  return(NULL)
}

texp.fn <- function (t, type) {
  fn <- 1e-04 + 4 * (exp(-550 * (t - 0.2)^2) +
        exp(-200 * (t - 0.5)^2) + exp(-950 * (t - 0.8)^2))
  if (type == "mean")
    return(NULL)
  else if (type == "var")
    return(fn)
  return(NULL)
}

cons.fn <- function (t, type) {
  fn <- rep(1,length(t))
  if (type == "mean")
    return(NULL)
  else if (type == "var")
    return(fn)
  return(NULL)
}

gaussian.1d <- function (args) {
  n      <- args$n
  rsnr   <- args$rsnr
  meanfn <- args$meanfn
  varfn  <- args$varfn
  
  t <- 1:n/n
  
  mu    <- meanfn(t,"mean")
  sigma <- sqrt(varfn(t,"var"))
  
  sig.true <- sigma/mean(sigma) * sd(mu)/rsnr^2
  x.data   <- rnorm(n,mu,sig.true)
  sig.est  <- sqrt(2/(3 * (n - 2)) *
                   sum((1/2 * x.data[1:(n - 2)] -
                        x.data[2:(n - 1)] +
                        1/2 * x.data[3:n])^2))
  
  meta  <- list(mu = mu)
  input <- list(x = x.data,sig.true = sig.true,sig.est = sig.est)
  return(list(meta = meta,input = input))
} 
