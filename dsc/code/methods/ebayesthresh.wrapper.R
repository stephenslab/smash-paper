# Runs wavelet denoising using the Bayesian shrinkage/thresholding
# procedure EbayesThresh, assuming constant variance. The inputs are:
#
#   input, a list containing the data (x), the ground-truth sigma
#   values (sig.true), and the estimated sigma values (sig.est).
#
#   args, a list containing family and filter.number, which determine
#   the wavelet basis used.
#
# The return value is the estimated (posterior) mean function.
ebayesthresh.wrapper <- function (input, args) {
  if (is.null(args$filter.number)) 
    args$filter.number <- 8
  if (is.null(args$family)) 
    args$family <- "DaubLeAsymm"
  
  n   <- length(input$x)
  J   <- log2(n)
  cat("     - Computing discrete wavelet transform.\n")
  x.w <- wd(input$x,args$filter.number,args$family,type = "station")
  cat("     - Running Empirical Bayes thresholding.\n")
  x.w <- ebayesthresh.wavelet(x.w,vscale = input$sig.est,prior = "laplace",
                              threshrule = "mean")
  mu.est <- AvBasis(convert(x.w))
  cat("     - Outputting posterior mean estimate.\n")
  return(mu.est)
} 
