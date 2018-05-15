# Specify true nonparametric mean and variance functions.
fTrue    <- function(x) return(sin(3*pi*x^2))
gTrue    <- function(x) return(exp(0.1 + cos(4*pi*x)))
loggTrue <- function(x) return(0.1 + cos(4*pi*x))

