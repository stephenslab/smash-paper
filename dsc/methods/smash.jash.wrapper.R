#runs wavelet shrinkage method SMASH and jointly estimates the variance function, using JASH for variance shrinkage
#inputs:
#input: a list containing x: the data, sig.true the true sigma values, and sig.est: the estimated sigma values
#args: a list containing family and filter.number, which determine the wavelet basis used
#
#returns the estimated (posterior mean) mean function
smash.jash.wrapper = function(input, args) {
  mu.est = smash.gaus(input$x, v.basis = TRUE, filter.number = args$filter.number, family = args$family, jash = TRUE)
  return(mu.est)
} 
