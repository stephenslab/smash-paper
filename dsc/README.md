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

Download or clone the [git repository][smash-github] on your computer.




* Clone or download this repository.

* Install these R packages: smashr (from ashwave repo; see below),
  dscr, ashr (see below), dependencies in ashr & smash (e.g. wavethresh,
  EbayesThresh).

* Matlab libraries:
  WaveLab (http://statweb.stanford.edu/~wavelab),
  WavDen (http://www-ljk.imag.fr/SMS/software/GaussianWaveDen/down.html).
  Follow the instructions provided with these software distributions
  to install these.
  
* Add path to the Matlab bin directory to R's path variable using add_path().

* Add path to dscr-smash/methods folder in Matlab.

* Run the R script run_dsc.R, either in batch mode or interactively.

* Upon completion results will be saved as both as an R object
  (res.Robj) and an R image file (res.RData).

[smash-github]: https://github.com/stephenslab/smash-paper
