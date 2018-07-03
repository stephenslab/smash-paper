# dscr-smash

A simulation study comparing the wavelet shrinkage procedure SMASH
(Xing, Carbonetto and Stephens (2018)) and other popular wavelet
denoising procedures for observations with Gaussian noise.

# Background 

For a general introduction to DSCs, see
[here](https://github.com/stephens999/dscr/blob/master/intro.md).

In this simulation study, we aim to compare several methods for
performing nonparametric regression when the observations follow a
Gaussian distribution. The problem is of the form
Y_i=\mu_i+\epsilon_i, i=1,...,n, where \mu is the underlying mean
function and assumed to be "smooth", and \epsilon_i's are independent
Gaussian noise with mean 0 and varying standard deviation
\sigma_i. The goal is to recover \mu as accurately as possible given
Y.

The simulation schemes cover a wide range of different mean and
variance functions (\mu and \sigma), as well as different signal to
noise ratios (SNRs). The various methods are run and the resulting
estimate of \mu is scored using the mean integrated squared error
(MISE), which is simply the standard mean squared error rescaled
appropriately.

# Installation

```
install.packages("devtools")
library(devtools)
install_github("zrxing/dscr-smash")
```

# To reproduce results from Xing & Stephens (2016)

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

* To visualize the simulation results using the interactive RShiny
  plot, open graphs.Rmd in RStudio and click on "Run
  Document". Additional description is contained within graphs.Rmd.

# Input, meta and output formats

This DSC uses the following formats:

`input: list(x [vector], sig.true [vector], sig.est [vector])` #x is the vector of observations. sig.true contains the true values of \sigma_i, and sig.est is an estimate of \sigma under the assumption that all \sigma_i's are equal.

`meta: list(mu [vector])` #mu contains the true values of \mu_i as defined above


`output: an estimate of \mu [vector]` 

# Scores

The performance of a method is scored by the quantity 10000*sum((\mu_i-\mu-hat_i)^2)/sum(\mu_i^2), where \mu is the true mean function and \mu-hat is the estimated mean function

See [score.R](score.R).

(NB: results run with the following versions of the softwares:
ashr commit: bf9cca351d7cf804b67a56552f8057b3af7c11cb;
ashwave commit: e59bdebb61ae06a6459c9248a5167178fe0150df)

