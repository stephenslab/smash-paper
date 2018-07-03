score <- function (data, output)
  list(mise = 10000 * sum((output - data$meta$mu)^2)/sum(data$meta$mu^2))
add_score(dsc_smash,score,"mise") 
