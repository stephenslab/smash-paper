# This R script implements the "Gaussian mean estimation" simulations
# for the SMASH paper. For instructions on how to run this code, see
# the accompanying README.
library(AlgDesign)
library(EbayesThresh)
library(wavethresh)
library(smashr)

# SCRIPT PARAMETERS
# -----------------
# This specifies the command to start MATLAB. Note that the -r option
# is added to this command to run the methods implemented in MATLAB.
matlab.exec <- "matlab -nodisplay -nodesktop"

# DEFINE THE DSC
# --------------
cat("Creating DSC.\n")
if (!dir.exists("results"))
  dir.create("results")
if (!dir.exists("results/temp"))
  dir.create("results/temp")
dscr.path <- getwd()
dsc_smash <- new_dsc("mean_function_estimation","results")
source_dir("code/methods")
source("code/datamaker.R")
source("code/scenarios.R")
source("code/methods.R")
source("code/score.R")

# RUN THE DSC
# -----------
cat("Running DSC.\n")
timing <- system.time(res <- run_dsc(dsc_smash))
cat(sprintf("Computation took %d seconds.\n",round(timing["elapsed"])))

# SAVE RESULTS TO FILE
# --------------------
cat("Saving DSC results to file.\n")
save(list = c("dsc_smash","res"),file = "dscr.RData")

# SESSION INFORMATION
# -------------------
# This is the version of R and the packages that were used to generate
# these results.
sessionInfo()
