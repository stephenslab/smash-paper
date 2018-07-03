# Runs wavelet shrinkage method SMASH, and jointly estimates the
# variance function, using JASH for variance shrinkage. The inputs
# are:
#
#   input, a list containing the data (x), the ground-truth sigma
#   values (sig.true), and the estimated sigma values (sig.est).
#
#   args, a list containing family and filter.number, which determine
#   the wavelet basis used.
#
# The return value is the estimated (posterior) mean function.
smash.jash.wrapper <- function(input, args)
  smash.gaus(input$x,v.basis = TRUE,filter.number = args$filter.number,
             family = args$family,jash = TRUE)
