# TO DO: Explain here what this script does, and how to use it.
#
# TO DO: (Perhaps) change the name of this file.
#
# SET UP ENVIROMENT
# -----------------
# Load the ggplot2 and cowplot packages, and the functions definining
# the mean and variances used to simulate the data.
library(ggplot2)
library(cowplot)
source("../code/signals.R")

# LOAD RESULTS
# ------------
# Load the results of the simulation experiments.
cat("Loading simulation results.\n")
load("../output/gaus-dscr.RData")

# PLOT #1
# -------
# This reproduces the plot in Fig. 2 of the manuscript comparing the
# accuracy of estimated mean curves in the data sets simulated from
# the "spikes" mean function with constant variance.  function.
cat("Generating the first plot.\n")

# Extract the data used to generate the first plot.
homo.data.smash <-
  res[res$.id    == "sp.3.v1" &
      res$method == "smash.s8",]
homo.data.smash.homo <-
  res[res$.id    == "sp.3.v1" &
      res$method == "smash.homo.s8",]
homo.data.tithresh <-
  res[res$.id == "sp.3.v1" &
      res$method == "tithresh.homo.s8",]
homo.data.ebayes <-
  res[res$.id    == "sp.3.v1" &
      res$method == "ebayesthresh",]
homo.data.smash.true <-
  res[res$.id == "sp.3.v1" &
  res$method  == "smash.true.s8",]
homo.data <-
  res[res$.id == "sp.3.v1" &
  (res$method == "smash.s8" |
   res$method == "ebayesthresh" |
   res$method == "tithresh.homo.s8"),]

# Transform these data into a data frame suitable for ggplot2.
pdat <-
  rbind(data.frame(method      = "smash",
                   method.type = "est",
                   mise        = homo.data.smash$mise),
        data.frame(method      = "smash.homo",
                   method.type = "homo",
                   mise        = homo.data.smash.homo$mise),
        data.frame(method      = "tithresh",
                   method.type = "homo",
                   mise        = homo.data.tithresh$mise),
        data.frame(method      = "ebayesthresh",
                   method.type = "homo",
                   mise        = homo.data.ebayes$mise),
        data.frame(method      = "smash.true",
                   method.type = "true",
                   mise        = homo.data.smash.true$mise))
pdat <-
  transform(pdat,
            method = factor(method,
                            names(sort(tapply(pdat$mise,pdat$method,mean),
                                       decreasing = TRUE))))

# Create the combined boxplot and violin plot using ggplot2.
p <- ggplot(pdat,aes(x = method,y = mise,fill = method.type)) +
     geom_violin(fill = "skyblue",color = "skyblue") +
     geom_boxplot(width = 0.15,outlier.shape = NA) +
     scale_y_continuous(breaks = seq(6,16,2)) +
     scale_fill_manual(values = c("darkorange","dodgerblue","gold"),
                       guide = FALSE) +
     coord_flip() +
     labs(x = "",y = "MISE") +
     theme(axis.line = element_blank(),
           axis.ticks.y = element_blank())
print(p)
invisible(readline(prompt = "Press [enter] to continue plotting... "))

# PLOT #2
# -------
# TO DO: Explain here what is shown in the second plot.
cat("Generating the second plot.\n")

# Extract the data used to generate the second plot.
hetero.data.smash <-
  res[res$.id == "sp.3.v5" & res$method == "smash.s8",]
hetero.data.smash.homo <-
  res[res$.id == "sp.3.v5" & res$method == "smash.homo.s8",]
hetero.data.tithresh.homo <-
  res[res$.id == "sp.3.v5" & res$method == "tithresh.homo.s8",]
hetero.data.tithresh.rmad <-
  res[res$.id == "sp.3.v5" & res$method == "tithresh.rmad.s8",]
hetero.data.tithresh.smash <-
  res[res$.id == "sp.3.v5" & res$method == "tithresh.smash.s8",]
hetero.data.tithresh.true <-
  res[res$.id == "sp.3.v5" & res$method == "tithresh.true.s8",]
hetero.data.ebayes <-
  res[res$.id == "sp.3.v5" & res$method == "ebayesthresh",]
hetero.data.smash.true <-
  res[res$.id == "sp.3.v5" & res$method == "smash.true.s8",]

# Transform these data into a data frame suitable for ggplot2.
pdat <-
  rbind(data.frame(method      = "smash",
                   method.type = "est",
                   mise        = hetero.data.smash$mise),
        data.frame(method      = "smash.homo",
                   method.type = "homo",
                   mise        = hetero.data.smash.homo$mise),
        data.frame(method      = "tithresh.rmad",
                   method.type = "tithresh",
                   mise        = hetero.data.tithresh.rmad$mise),
        data.frame(method      = "tithresh.smash",
                   method.type = "tithresh",
                   mise        = hetero.data.tithresh.smash$mise),
        data.frame(method      = "tithresh.true",
                   method.type = "tithresh",
                   mise        = hetero.data.tithresh.true$mise),
        data.frame(method      = "ebayesthresh",
                   method.type = "homo",
                   mise        = hetero.data.ebayes$mise),
        data.frame(method      = "smash.true",
                   method.type = "true",
                   mise        = hetero.data.smash.true$mise))
pdat <-
  transform(pdat,
            method = factor(method,
                            names(sort(tapply(pdat$mise,pdat$method,mean),
                                       decreasing = TRUE))))

# Create the combined boxplot and violin plot using ggplot2.
p <- ggplot(pdat,aes(x = method,y = mise,fill = method.type)) +
     geom_violin(fill = "skyblue",color = "skyblue") +
     geom_boxplot(width = 0.15,outlier.shape = NA) +
     scale_fill_manual(values=c("darkorange","dodgerblue","limegreen","gold"),
                       guide = FALSE) +
     coord_flip() +
     scale_y_continuous(breaks = seq(10,70,10)) +
     labs(x = "",y = "MISE") +
     theme(axis.line = element_blank(),
           axis.ticks.y = element_blank())
print(p)
invisible(readline(prompt = "Press [enter] to continue plotting... "))

# PLOT #3
# -------
# TO DO: Explain here what is shown in the third plot.
hetero.data.smash.2 <-
  res[res$.id == "cor.3.v3" & res$method == "smash.s8",]
hetero.data.smash.homo.2 <-
  res[res$.id == "cor.3.v3" & res$method == "smash.homo.s8",]
hetero.data.tithresh.homo.2 <-
  res[res$.id == "cor.3.v3" & res$method == "tithresh.homo.s8",]
hetero.data.tithresh.rmad.2 <-
  res[res$.id == "cor.3.v3" & res$method == "tithresh.rmad.s8",]
hetero.data.tithresh.smash.2 <-
  res[res$.id == "cor.3.v3" & res$method == "tithresh.smash.s8",]
hetero.data.tithresh.true.2 <-
  res[res$.id == "cor.3.v3" & res$method == "tithresh.true.s8",]
hetero.data.ebayes.2 <-
  res[res$.id == "cor.3.v3" & res$method == "ebayesthresh",]
hetero.data.smash.true.2 <-
  res[res$.id == "cor.3.v3" & res$method == "smash.true.s8",]

# Transform these data into a data frame suitable for ggplot2.
pdat <-
  rbind(data.frame(method      = "smash",
                   method.type = "est",
                   mise        = hetero.data.smash.2$mise),
        data.frame(method      = "smash.homo",
                   method.type = "homo",
                   mise        = hetero.data.smash.homo.2$mise),
        data.frame(method      = "tithresh.rmad",
                   method.type = "tithresh",
                   mise        = hetero.data.tithresh.rmad.2$mise),
        data.frame(method      = "tithresh.smash",
                   method.type = "tithresh",
                   mise        = hetero.data.tithresh.smash.2$mise),
        data.frame(method      = "tithresh.true",
                   method.type = "tithresh",
                   mise        = hetero.data.tithresh.true.2$mise),
        data.frame(method      = "ebayesthresh",
                   method.type = "homo",
                   mise        = hetero.data.ebayes.2$mise),
        data.frame(method      = "smash.true",
                   method.type = "true",
                   mise        = hetero.data.smash.true.2$mise))
pdat <-
  transform(pdat,
            method = factor(method,
                            names(sort(tapply(pdat$mise,pdat$method,mean),
                                       decreasing = TRUE))))

# Create the combined boxplot and violin plot using ggplot2.
p <- ggplot(pdat,aes(x = method,y = mise,fill = method.type)) +
     geom_violin(fill = "skyblue",color = "skyblue") +
     geom_boxplot(width = 0.15,outlier.shape = NA) +
     scale_fill_manual(values=c("darkorange","dodgerblue","limegreen","gold"),
                       guide = FALSE) +
     coord_flip() +
     scale_y_continuous(breaks = seq(10,70,10)) +
     labs(x = "",y = "MISE") +
     theme(axis.line = element_blank(),
           axis.ticks.y = element_blank())
print(p)

# SESSIONINFO
# -----------
# Print out information on the computing environment, including the
# version of R and the attached R packages, used to generate these
# results.
cat("Retrieving session information.\n")
print(sessionInfo())
