# This file defines some functions used by the "poisson" analysis.

create.results.table <- function (dat) {
  table.row.names <- c("SMASH","BMSM","Haar-Fisz")
  table.col.names <- c("(1/100, 3)","(1/8, 8)","(1/128, 128)")
  rownames(dat) <- table.row.names
  colnames(dat) <- table.col.names
  return(kable(dat,format = "html",
    table.attr = paste("class=\"table\";","style=\"font-family: sans-serif;",
                   "width: auto; margin-left: auto; margin-right: auto\"")))
}

get.violin.plot.data <- function (mise) {
  out <- NULL
  n   <- length(methods)
  for (i in 1:n)
    out <- rbind(out,data.frame(method = names(mise)[i],mise = mise[[i]]))
  out <- transform(out,method = factor(method))
  return(out)
}

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
