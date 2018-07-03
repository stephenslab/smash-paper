#runs wavelet denoising with TI-thresholding in matlab with the true variance known
#inputs:
#input: a list containing x: the data, sig.true the true sigma values, and sig.est: the estimated sigma values
#args: a list containing family and filter.number, which determine the wavelet basis used
#
#returns the estimated mean function
tithresh.true.wrapper = function(input, args) {
  mu.est = ti.thresh(input$x, sigma = input$sig.true, filter.number = args$filter.number, family = args$family)
  return(mu.est)
} 
