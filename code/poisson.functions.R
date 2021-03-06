# This file defines some functions used by the "poisson" analysis.

# This is the Gaussian denoising routine used in the Haar-Fisz method;
# this specifies the "meth.1" argument in the call to "denoise.poisson"
# from the haarfisz package.
hf.la10.ti4 <- function (x) {
  TT      <- length(x)
  thresh  <- sqrt(2*log(TT))
  x.w     <- wd(x,10,"DaubLeAsymm",type = "station")
  x.w.t   <- threshold(x.w,levels = 4:(x.w$nlevels-1),policy = "manual",
                       value = thresh,type = "hard")
  x.w.t.r <- AvBasis(convert(x.w.t))
  return(x.w.t.r)
}

# Create a plot showing the test function (mu) and an example data set
# simulated using that test function (x).
create.signal.plot <- function (t, x, mu)
  ggplot(data.frame(x = t,y = x),aes_string(x = "x",y = "y")) +
    geom_point(color = "gray",shape = 1) +
    geom_line(dat = data.frame(x = t,y = mu),color = "black",size = 1) +
    labs(x = "",y = "")

# This function is used to create the tables in the "poisson"
# analysis.
create.results.table <- function (dat) {
  table.row.names <- c("SMASH","BMSM","Haar-Fisz")
  table.col.names <- c("(1/100, 3)","(1/8, 8)","(1/128, 128)")
  rownames(dat) <- table.row.names
  colnames(dat) <- table.col.names
  return(kable(dat,format = "html",
    table.attr = paste("class=\"table\";","style=\"font-family: sans-serif;",
                   "width: auto; margin-left: auto; margin-right: auto\"")))
}

# This function is used to prepare the simulation results for creating
# the violin plots in the "poisson" analysis.
get.violin.plot.data <- function (mise) {
  out <- NULL
  n   <- length(methods)
  for (i in 1:n)
    out <- rbind(out,data.frame(method = names(mise)[i],mise = mise[[i]]))
  out <- transform(out,method = factor(method))
  return(out)
}

# This function is used to create the violin plots in the "poisson"
# analysis.
create.violin.plots <- function (pdat1, pdat8, pdat128) {
  p1 <- ggplot(pdat1,aes(x = method,y = mise)) +
        geom_violin(fill = "skyblue",color = "skyblue") +
        geom_boxplot(width = 0.15,outlier.shape = NA) +
        coord_flip() +
        labs(x = "",y = "MISE",title = "intensities (1/100,3)") +
		theme_cowplot() +
		theme(plot.title = element_text(size = 14,face = "plain"),
		      axis.line = element_blank())
  p8 <- ggplot(pdat8,aes(x = method,y = mise)) +
        geom_violin(fill = "skyblue",color = "skyblue") +
        geom_boxplot(width = 0.15,outlier.shape = NA) +
        coord_flip() +
        labs(x = "",y = "MISE",title = "intensities (1/8,8)") +
		theme_cowplot() +
		theme(plot.title = element_text(size = 14,face = "plain"),
		      axis.line = element_blank())
  p128 <- ggplot(pdat128,aes(x = method,y = mise)) +
        geom_violin(fill = "skyblue",color = "skyblue") +
        geom_boxplot(width = 0.15,outlier.shape = NA) +
        coord_flip() +
        labs(x = "",y = "MISE",title = "intensities (1/128,128)") +
		theme_cowplot() +
		theme(plot.title = element_text(size = 14,face = "plain"),
		      axis.line = element_blank())
  return(plot_grid(p1,p8,p128,nrow = 1,ncol = 3))
}

# This function is used to draw the bar charts summarizing the results
# of the Poisson simulations.
create.bar.plots <- function (pdat, show.legend = TRUE) {
  intensities <- c("(1/100, 3)","(1/8, 8)","(1/128, 128)")
  
  # Prepare the data to be plotted.
  pdat <- transform(pdat,
                    method     = factor(method),
                    simulation = factor(simulation),
                    intensity  = factor(intensity,intensities))

  # Compute the mean, upper and lower quantiles, which will be used to
  # plot the barchart bars and the error bars.
  pdat <- with(pdat,tapply(mise,list(simulation,intensity,method),
                           function (x) c(mean(x),quantile(x,c(0.1,0.9)))))
  pdat <- adply(pdat,1:3)
  names(pdat) <- c("simulation","intensity","method","mean","low","high")
  
  # Create the bar charts.
  return(ggplot(pdat,aes_string(x = "method",y = "mean",fill = "method")) +
    geom_col(width = 0.65,show.legend = show.legend) +
    geom_errorbar(aes_string(x = "method",ymin = "low",ymax = "high"),
                  color = "black",inherit.aes = FALSE,width = 0.2,size = 0.3) +
    scale_x_discrete(breaks = NULL) +
    scale_fill_manual(values = c("dodgerblue","darkblue","darkorange")) +
    facet_grid(intensity ~ simulation,scales = "free") +
    labs(x = "",y = "MISE") +
    theme_cowplot(font_size = 14) +
    theme(strip.background = element_blank()))
}
