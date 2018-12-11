This directory contains scripts used to implement the "Gaussian mean
estimation" simulations in the SMASH paper. In these simulations, we
assess the ability of different signal denoising methods to recover
the true signal after being provided with Gaussian-distributed
observations of the signal. We consider scenarios in which the data
have homoskedastic errors (constant variance) and heteroskedastic
errors (non-constant variance).

To run these simulations, please follow the instructions below. These
instructions assume that you have R and MATLAB installed on your
computer.

File dscr.RData in the "output" directory of this repository contains
previously generated results.

## Instructions

1. Download or clone the [git repository][smash-github] on your computer.

2. Install the following R packages: [AlgDesign][algdesign],
   [wavethresh][wavethresh], [EbayesThresh][EbayesThresh], 
   [dscr][dscr] and [smashr][smash].

3. Run the `InstallMEX.m` script in the Wavelab850 subdirectory to
   build the MEX files from their C source.

4. Run the R script [run_dsc.R](run_dsc.R) from the "dsc" directory of
   the git repository. This can be done in batch mode (e.g., using
   Rscript), or interactively in R or RStudio. When running
   interactively, make sure your working directory is the "dsc"
   directory of this git repository. Modify the script parameters as
   needed.

5. Upon completion, results will be saved as both as an R object
   (res.Robj) and an R image file (res.RData).

[smash-github]: https://github.com/stephenslab/smash-paper
[smashr]: https://github.com/stephenslab/smashr
[dscr]: https://github.com/stephens999/dscr
[ebayesthresh]: https://github.com/stephenslab/EbayesThresh
[wavethresh]: https://cran.r-project.org/package=wavethresh
[algdesign]: https://cran.r-project.org/package=AlgDesign
