# This file defines a few functions used in the "gaussmeanest" analysis.

# This function extracts the DSC results from a simulation scenario
# with homoskedastic (constant) variance, and transforms the reseults
# into a data frame suitable for plotting using ggplot2.
get.results.homosked <- function (res, scenario) {
  dat.smash     <- res[res$.id == scenario & res$method == "smash.s8",]
  dat.smash.homo<- res[res$.id == scenario & res$method == "smash.homo.s8",]
  dat.tithresh  <- res[res$.id == scenario & res$method == "tithresh.homo.s8",]
  dat.ebayes    <- res[res$.id == scenario & res$method == "ebayesthresh",]
  dat.smash.true<- res[res$.id == scenario & res$method == "smash.true.s8",]
  return(rbind(data.frame(method      = "smash",
                          method.type = "est",
                          mise        = dat.smash$mise),
               data.frame(method      = "smash.homo",
                          method.type = "homo",
                          mise        = dat.smash.homo$mise),
               data.frame(method      = "tithresh",
                          method.type = "homo",
                          mise        = dat.tithresh$mise),
               data.frame(method      = "ebayesthresh",
                          method.type = "homo",
                          mise        = dat.ebayes$mise),
               data.frame(method      = "smash.true",
                          method.type = "true",
                          mise        = dat.smash.true$mise)))
}

    
