# TO DO: Explain what this script does, and how to use it.
library(dscr)
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
dscr.path <- getwd()
dsc_smash <- new_dsc("mean_function_estimation","output")
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

stop()

# SAVE RESULTS TO FILE
# --------------------
cat("Saving DSC results to file.\n")
save(list = c("dsc_smash","res"),file = "res.RData")

# SESSION INFORMATION
# -------------------
# This is the version of R and the packages that were used to generate
# these results.
sessionInfo()
