# Runs wavelet denoising with TI-thresholding with the true
# variance known. The inputs are:
#
#   input, a list containing the data (x), the true sigma values
#   (sig.true), and the estimated sigma values (sig.est);
#
#   args, a list containing family and filter.number, which determine
#   the wavelet basis used.
#
# The return value is the estimated mean function.
tithresh.true.wrapper <- function(input, args)
  ti.thresh(input$x,sigma = input$sig.true,filter.number = args$filter.number,
            family = args$family)

