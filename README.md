This repository contains the companion paper and supplementary materials for the R package SMASH. For more information regarding the usage of the package see https://github.com/stephenslab/ashwave. For details on the simulation performed in dscr refer to https://github.com/zrxing/dscr-smash


#Description
This repo contains the following directories:

paper: contains the latex and pdf files for various drafts of the paper and supplementary materials. Figures used in the paper are all within this directory.

package: contains the SMASH package; for the newest version use https://github.com/zrxing/dscr-smash

intro: contains an .rmd file intro illustrating some simple simulations and data analyses using SMASH. More detailed than the readme in the SMASH repo.

res_paper: contains .RData objects which contain simulation results for (almost) all the simulation studies conducted. Specifically, res_gaus_dscr.RData contains the results from Gaussian simulations with n = 1024 based on the dsc framework, part of which is shown in the main paper. Full results are available as an RShiny app from the dscr-smash repo (https://github.com/zrxing/dscr-smash). res_gaus_256.RData and res_gaus_512.RData contains results from Gaussian simulations with n = 256 and n = 512 respectively. These results are not shown in the main paper. res_pois.RData contains the simulation results for Poisson noise, and results are partially presented in the main paper, while full results are shown in supplementary materials. res_pois_hf.RData contains results from experimenting with different Gaussian denoising techniques in the second stage of the Haar-Fisz algorithm.

code_paper: contains three sub-directories. gaus_data contains the code for several analyses for Gaussian data using SMASH. pois_data contains the code for running SMASH on the ChIP-seq experiment from ENCODE. sim contains several R scripts for generating figures and tables used in the paper and the supplementary materials. Specifically, plot_gaus.R and plot_pois.R contain all the code for generating all the figures in the paper; table_gaus.R generates Tex tables containing the MISE for each simulation scenario, as well as an additional table detailing additional information about the various methods used; and table_pois generates the full results for the Poisson simulation study which are shown in supplementary materials.



(NB: stable ashr commit: bf9cca351d7cf804b67a56552f8057b3af7c11cb;
stable ashwave commit: e59bdebb61ae06a6459c9248a5167178fe0150df)