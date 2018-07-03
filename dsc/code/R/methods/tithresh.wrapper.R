# Runs wavelet denoising with TI-thresholding, and jointly estimates
# the variance function. The inputs are:
# 
#   input, a list containing the data (x),the ground-truth sigma
#   values (sig.true), and the estimated sigma values (sig.est);
#
#   args, a list containing family and filter.number, which determine
#   the wavelet basis used.
#
# The return value is the estimated mean function.
tithresh.wrapper <- function (input, args)
  ti.thresh(input$x,method = args$method,filter.number = args$filter.number,
            family = args$family)


