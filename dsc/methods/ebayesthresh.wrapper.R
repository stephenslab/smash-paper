# runs wavelet denoising using the Bayesian shrinkage/thresholding procedure EbayesThresh, assuming constant variance
# inputs: 
# input: a list containing x: the data, sig.true: the true sigma values, and sig.est: the estimated sigma values
# args: a list containing family and filter.number, which determine the wavelet basis used 
#
# returns the estimated (posterior mean) mean function
ebayesthresh.wrapper = function(input, args) {
  if (is.null(args$filter.number)) 
    args$filter.number = 8
  if (is.null(args$family)) 
    args$family = "DaubLeAsymm"
  
  n = length(input$x)
  J = log2(n)
  x.w <- wd(input$x, args$filter.number, args$family, type = "station")
  x.w <- ebayesthresh.wavelet(x.w, vscale = input$sig.est, prior = "laplace", threshrule = "mean")
  mu.est = AvBasis(convert(x.w))
  return(mu.est)
} 
