# Runs wavelet shrinkage method SMASH and jointly estimates the
# variance function. The inputs are:
#
#   input, a list containing the data set (x), the true sigma values
#   (sig.true) and the estimated sigma values (sig.est);
#
#   args, a list containing family and filter.number, which
#   determine the wavelet basis used.
#
# The return value is the estimated (posterior) mean function.
smash.wrapper <- function (input, args)
  smash.gaus(input$x,filter.number = args$filter.number,family = args$family)

