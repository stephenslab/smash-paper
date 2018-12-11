meanfn       <- c(spikes.fn,bumps.fn,blocks.fn,angles.fn,
                  doppler.fn,blip.fn,cor.fn)
meanfn.short <- c("sp","bump","blk","ang","dop","blip","cor")
rsnr         <- sqrt(c(1,3))
varfn        <- c(cons.fn,texp.fn,doppler.fn,bumps.fn,cblocks.fn)
varfn.short  <- c("cons","texp","dop","bump","cblk")
seeds        <- 1:100

# TESTING
# -------
meanfn       <- meanfn[1:3]
meanfn.short <- meanfn.short[1:3]
varfn        <- varfn[c(1,3)]
varfn.short  <- varfn.short[c(1,3)]
seeds        <- 1:2

design <- gen.factorial(c(length(varfn),length(rsnr),length(meanfn)),
                        center = FALSE)

for (i in 1:nrow(design)) {
  scenario.name <- paste(meanfn.short[design[i,3]],
                         rsnr[design[i,2]]^2,
                         varfn.short[design[i,1]],sep = ".")
  scenario.args <- list(n = 1024,rsnr = rsnr[design[i,2]],
                        meanfn = meanfn[[design[i,3]]],
                        varfn = varfn[[design[i,1]]])
  add_scenario(dsc_smash,name = scenario.name,fn = gaussian.1d,
               args = scenario.args,seed = seeds)
}
