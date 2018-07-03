#runs wavelet shrinkage method SMASH assuming a constant variance
#inputs:
#input: a list containing x: the data, sig.true the true sigma values, and sig.est: the estimated sigma values
#args: a list containing family and filter.number, which determine the wavelet basis used
#
#returns the estimated (posterior mean) mean function
smash.homo.wrapper = function(input, args) {
  mu.est = smash.gaus(input$x, input$sig.est, filter.number = args$filter.number, family = args$family)
  return(mu.est)
} 
