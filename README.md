# Source to reproduce results from Xing & Stephens (2016)

This repository contains source code and data to reproduce some of the
empirical results presented in
[Xing & Stephens (2016)][smashr-preprint]. For additional source code,
see the [dscr-smash][dscr-smash]
repository. The new methods can be found in the
[smashr package][smashr].

If you find a bug, please post an [issue][issues].

## Citing this work

If you find any of the source code in this repository useful for your
work, please cite our paper:

> Zhengrong Xing and Matthew Stephens (2016). *Smoothing via Adaptive
Shrinkage (smash): denoising Poisson and heteroskedastic Gaussian
signals.* [arXiv:1709.10066](http://arxiv.org/abs/1605.07787).

## License

Copyright (c) 2016-2018, Zhengrong Xing.

All source code and software in this repository are made available
under the terms of the
[MIT license](https://opensource.org/licenses/mit-license.html). See
the [LICENSE](LICENSE) file for the full text of the license.

## Instructions

To reproduce the results of Xing & Stephens (2016), you need
to...

1. Install the devtools, scales, wavethresh and EbayesThresh packages
   from CRAN:

   ```R
   install.packages(c("devtools","scales","wavethresh","EbayesThresh"))
   ```

2. Install the [smashr](https://github.com/stephenslab/smashr) package
   from GitHub:

   ```R
   library(devtools)
   install_github("stephenslab/smashr")
   ```

   The `install_github` command should automatically attempt to
   download and install from CRAN any additional packages that are
   needed for smashr.

## What's included

This repo contains the following directories:

`res_paper`: .RData objects which contain simulation results for all the simulation studies conducted. Specifically:

* `res_gaus_dscr.RData` : Gaussian simulations with n = 1024 (run in dsc framework; code at https://github.com/zrxing/dscr-smash) 
* `res_gaus_256.RData`, `res_gaus_512.RData` : Gaussian simulations with n = 256 and n = 512 respectively. 
* `res_pois.RData` contains the simulation results for Poisson observations, 
* `res_pois_hf.RData` contains results from experimenting with different Gaussian denoising techniques in the second stage of the Haar-Fisz Poisson de-noising algorithm.

`code`: contains three sub-directories. 

* `gaus_data` contains the code for several analyses for Gaussian data
  using SMASH, including examples not covered in the main paper.

* `pois_data` contains the code for running SMASH on the ChIP-seq
  experiment from ENCODE.

* `sim` contains several R scripts for generating figures and tables
  used in the paper and the supplementary materials.

  + `pois` contains R and Matlab code for generating the simulation results for Poisson observations (`sim_pois.R` and `sim_pois.m`), as well as an R script to test out different Haar-Fisz options (`sim_hf.R`).    
  + `plot_gaus.R` and `plot_pois.R` contain all the code for generating all the figures in the paper; 
  + `table_gaus.R` generates Tex tables containing the MISE for each simulation scenario, as well as an additional table detailing additional information about the various methods used; 
  + `table_pois` generates the full results for the Poisson simulation study which are shown in supplementary materials.

(NB: stable ashr commit: bf9cca351d7cf804b67a56552f8057b3af7c11cb;
stable smash-paper commit: e59bdebb61ae06a6459c9248a5167178fe0150df)

## Credits

This project was developed by [Zhengrong Xing][zhengrong] at the
University of Chicago, with support and contributions from
[Peter Carbonetto][peter] and [Matthew Stephens][matthew].

[smashr]: https://github.com/stephenslab/smashr
[smashr-preprint]: http://arxiv.org/abs/1605.07787
[dscr-smash]: https://github.com/stephenslab/dscr-smash
[issues]: https://github.com/stephenslab/smash-paper/issues
[zhengrong]: https://github.com/zrxing
[peter]: http://pcarbo.github.io
[matthew]: http://stephenslab.uchicago.edu
