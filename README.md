This repository contains the companion paper and supplementary materials for the R package ```smashr```. For more information regarding the usage of the package see https://github.com/stephenslab/smashr. For details on the Gaussian simulations performed in dscr refer to https://github.com/zrxing/dscr-smash.


#Description

This repo contains the following directories:

`paper`: latex and pdf files for various drafts of the paper and supplementary materials. Figures from the paper are all within this directory.

`intro`: contains an .rmd file intro illustrating some simple simulations and data analyses using SMASH. More detailed than the readme in the ```smashr``` repo.

`res_paper`: .RData objects which contain simulation results for all the simulation studies conducted. Specifically:     

* `res_gaus_dscr.RData` : Gaussian simulations with n = 1024 (run in dsc framework; code at https://github.com/zrxing/dscr-smash) 
* `res_gaus_256.RData`, `res_gaus_512.RData` : Gaussian simulations with n = 256 and n = 512 respectively. 
* `res_pois.RData` contains the simulation results for Poisson observations, 
* `res_pois_hf.RData` contains results from experimenting with different Gaussian denoising techniques in the second stage of the Haar-Fisz algorithm.

`code_paper`: contains three sub-directories. 
* `gaus_data` contains the code for several analyses for Gaussian data using SMASH. 
* `pois_data` contains the code for running SMASH on the ChIP-seq experiment from ENCODE. 
* `sim` contains several R scripts for generating figures and tables used in the paper and the supplementary materials. 
  + `sim` contains R and Matlab code for generating the simulation results for Poisson observations (`sim_pois.R` and `sim_pois.m`), as well as an R script to test out different Haar-Fisz options (`sim_hf.R`).    
  + `plot_gaus.R` and `plot_pois.R` contain all the code for generating all the figures in the paper; 
  + `table_gaus.R` generates Tex tables containing the MISE for each simulation scenario, as well as an additional table detailing additional information about the various methods used; 
  + `table_pois` generates the full results for the Poisson simulation study which are shown in supplementary materials.



(NB: stable ashr commit: bf9cca351d7cf804b67a56552f8057b3af7c11cb;
stable smash-paper commit: e59bdebb61ae06a6459c9248a5167178fe0150df)
