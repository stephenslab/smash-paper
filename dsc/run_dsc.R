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
cat("Constructing DSC.\n")
dscr.path <- getwd()
dsc_smash <- new_dsc("mean_function_estimation","dscr-smash-output")
source("scenarios.R")
source("methods.R")
source("score.R")

# RUN THE DSC
# -----------
cat("Running DSC.\n")
timing <- system.time(res <- run_dsc(dsc_smash))
cat(sprintf("Computation took %d seconds.\n",round(timing["elapsed"])))

# SAVE RESULTS TO FILE
# --------------------
cat("Saving DSC results to file.\n")
save(list = c("dsc_smash","res"),file = "res.RData")
