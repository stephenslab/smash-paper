score = function(data, output) {
  mise = 10000 * sum((output - data$meta$mu)^2)/sum(data$meta$mu^2)
  return(list(mise = mise))
}

add_score(dsc_smash, score, "mise") 
