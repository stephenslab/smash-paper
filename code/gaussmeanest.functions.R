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

# This function is used to draw the bar charts summarizing the results
# of the Gaussian data simulations.
create.bar.plots <- function (pdat) {
    
  # Compute the mean, upper and lower quantiles, which will be used to
  # plot the barchart bars and the error bars.
  pdat <- with(pdat,tapply(mise,list(sim,snr,method),
                         function (x) c(mean(x),quantile(x,c(0.1,0.9)))))
  pdat <- adply(pdat,1:3)
  names(pdat) <- c("simulation","snr","method","mean","low","high")

  # Re-order the factor levels in the "method" column.
  pdat <- transform(pdat,
                    method = factor(method,c("smash.true","smash.homo","smash",
                                             "tithresh","ebayesthresh")))

  # Create the bar charts.
  return(ggplot(pdat,aes_string(x = "method",y = "mean",fill = "method")) +
         geom_col(width = 0.65) +
         geom_errorbar(aes_string(x = "method",ymin = "low",ymax = "high"),
                       color = "black",inherit.aes = FALSE,width = 0.3,
                       size = 0.4) +
         scale_x_discrete(breaks = NULL) +
         scale_fill_manual(values = c("gold","darkorange","indianred",
                                      "darkblue","royalblue")) +
         facet_grid(snr ~ simulation,scales = "free") +
         labs(x = "",y = "MISE") +
         theme_cowplot(font_size = 12) +
         theme(strip.background = element_blank()))
}
         
