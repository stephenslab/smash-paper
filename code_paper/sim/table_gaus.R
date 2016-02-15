#This script produces supplementary tables for Gaussian simulations


library(xtable)

load("res_paper/res_gaus_dscr.RData")

scale.baseline = function(x) x/x[1]

mise.res = aggregate(mise ~ method + scenario, res, mean)
table.order = c(6, 10, 16, 18, 8, 3, 5, 13, 14, 12, 20)
colname.order = c(7, 4, 3, 1, 6, 2, 5)

tex.row.names = c("SMASH, joint, Haar",
                  "SMASH, joint, Symm8",
                  "TI-thresh, RMAD, Symm8",
                  "TI-thresh, SMASH, Symm8",
                  "SMASH, homo, Symm8",
                  "Ebayes, homo, Symm8",
                  "PostMean, homo, Symm8",
                  "SURE, homo, Symm8",
                  "TI-thresh, homo, Symm8",
                  "SMASH, truth, Symm8",
                  "TI-thresh, truth, Symm8"
)

tex.col.names = c("Spikes", "Bumps", "Block", "Angles", "Doppler", "Blip", "Corner")

library(AlgDesign)

rsnr=sqrt(c(1,3))
varfn=1:5

design=gen.factorial(c(length(rsnr),length(varfn)),center=FALSE)

for(i in 1:dim(design)[1]){
  table.suffix = paste0(".",rsnr[design[i,1]]^2,".","v",varfn[design[i,2]])  
  
  name.ind = grep(paste0("*",table.suffix), mise.res$scenario)
  table.temp = matrix(mise.res[name.ind, 3], 20, 7)
  table.temp = table.temp[table.order, colname.order]
  table.temp = apply(table.temp, 2, scale.baseline)
  
  rownames(table.temp) = tex.row.names
  colnames(table.temp) = tex.col.names
  
  assign(paste0("table",table.suffix),table.temp)
}

ndigits=2

edit.smallest = function(x){
  x[which(x == min(x))] = paste0(min(x), "(*)")
  return(x)
}


mark.smallest = function(table, ndigits){
  table = apply(table, 2, sprintf, fmt = paste0("%.", ndigits, "f"))
  table = apply(table, 2, edit.smallest)
  return(table)
}


print(xtable(mark.smallest(table.1.v1,2),caption="Relative MSE of various methods for iid Gaussian errors, SNR=1",digits=ndigits))
print(xtable(mark.smallest(table.3.v1,2),caption="Relative MSE of various methods for iid Gaussian errors, SNR=3",digits=ndigits))


####

print(xtable(mark.smallest(table.1.v2,2),caption="Relative MSE of various methods for variance function V2 as in Cai & Wang (2008), SNR=1, Gaussian errors",digits=ndigits))
print(xtable(mark.smallest(table.3.v2,2),caption="Relative MSE of various methods for variance function V2 as in Cai & Wang (2008), SNR=3, Gaussian errors",digits=ndigits))
print(xtable(mark.smallest(table.1.v3,2),caption="Relative MSE of various methods for the Doppler variance function, SNR=1, Gaussian errors",digits=ndigits))
print(xtable(mark.smallest(table.3.v3,2),caption="Relative MSE of various methods for the Doppler variance function, SNR=3, Gaussian errors",digits=ndigits))
print(xtable(mark.smallest(table.1.v4,2),caption="Relative MSE of various methods for the Bumps variance function, SNR=1, Gaussian errors",digits=ndigits))
print(xtable(mark.smallest(table.3.v4,2),caption="Relative MSE of various methods for the Bumps variance function, SNR=3, Gaussian errors",digits=ndigits))
print(xtable(mark.smallest(table.1.v5,2),caption="Relative MSE of various methods for the Clipped Blocks variance function, SNR=1, Gaussian errors",digits=ndigits))
print(xtable(mark.smallest(table.3.v5,2),caption="Relative MSE of various methods for the Clipped Blocks variance function, SNR=3, Gaussian errors",digits=ndigits))
#####


table.1.v5.short=table.1.v5[c(-6,-7,-8),]
table.3.v5.short=table.3.v5[c(-6,-7,-8),]
table.1.v4.short=table.1.v4[c(-6,-7,-8),]
table.3.v4.short=table.3.v4[c(-6,-7,-8),]



print(xtable(table.1.v5.short,caption="Relative MSE of various methods for the Clipped Blocks variance function, SNR=1, Gaussian errors",digits=ndigits))
print(xtable(table.3.v5.short,caption="Relative MSE of various methods for the Clipped Blocks variance function, SNR=3, Gaussian errors",digits=ndigits))
print(xtable(table.1.v4.short,caption="Relative MSE of various methods for the Bumps variance function, SNR=1, Gaussian errors",digits=ndigits))
print(xtable(table.3.v4.short,caption="Relative MSE of various methods for the Bumps variance function, SNR=3, Gaussian errors",digits=ndigits))

#####
row.names.all = c("BAMS, homo, Symm8",
                  "BlockJS, homo, Symm8",
                  "Ebayes, homo, Symm8",
                  "Neighblock, homo, Symm8",
                  "PostMean, homo, Symm8",
                  "SMASH, joint, Haar",
                  "SMASH, homo, Haar",
                  "SMASH, homo, Symm8",
                  "SMASH-JASH, joint, Haar",
                  "SMASH, joint, Symm8",
                  "SMASH, truth, Haar",
                  "SMASH, truth, Symm8",
                  "SURE, homo, Symm8",
                  "TI-thresh, homo, Symm8",
                  "TI-thresh, RMAD, Haar",
                  "TI-thresh, RMAD, Symm8",
                  "TI-thresh, SMASH, Haar",
                  "TI-thresh, SMASH, Symm8",
                  "TI-thresh, truth, Haar",
                  "TI-thresh, truth, Symm8"
)

software = c("Matlab - WaveDen", 
             "Matlab - WaveDen",
             "R - EbayesThresh",
             "Matlab - WaveDen", 
             "Matlab - WaveDen", 
             "R - SMASH", 
             "R - SMASH", 
             "R - SMASH", 
             "R - SMASH", 
             "R - SMASH", 
             "R - SMASH", 
             "R - SMASH", 
             "Matlab - WaveDen",
             "Matlab - WaveDen", 
             "R - SMASH",
             "R - SMASH",
             "R - SMASH",
             "R - SMASH",
             "R - SMASH", 
             "R - SMASH")

var.assumption = c("Homo", 
                   "Homo",
                   "Homo", 
                   "Homo", 
                   "Homo", 
                   "Hetero", 
                   "Homo", 
                   "Homo",
                   "Hetero",
                   "Hetero",
                   "Truth", 
                   "Truth",
                   "Homo",
                   "Homo",
                   "Hetero",
                   "Hetero",
                   "Hetero",
                   "Hetero",
                   "Truth",
                   "Truth"
                   )

basis = c("Symm8", 
          "Symm8",
          "Symm8",
          "Symm8",
          "Symm8",
          "Haar",
          "Haar",
          "Symm8",
          "Haar",
          "Symm8",
          "Haar",
          "Symm8",
          "Symm8",
          "Symm8",
          "Haar",
          "Symm8",
          "Haar",
          "Symm8",
          "Haar",
          "Symm8")

ti.status = c(FALSE,
              FALSE,
              TRUE,
              FALSE,
              FALSE,
              TRUE,
              TRUE,
              TRUE,
              TRUE,
              TRUE,
              TRUE,
              TRUE,
              FALSE,
              TRUE,
              TRUE,
              TRUE,
              TRUE,
              TRUE,
              TRUE,
              TRUE)

shrink.type = c("Bayesian",
                "Thresholding",
                "Empirical Bayes",
                "Thresholding",
                "Empirical Bayes",
                "Empirical Bayes",
                "Empirical Bayes",
                "Empirical Bayes",
                "Empirical Bayes",
                "Empirical Bayes",
                "Empirical Bayes",
                "Empirical Bayes",
                "Thresholding",
                "Thresholding",
                "Thresholding",
                "Thresholding",
                "Thresholding",
                "Thresholding",
                "Thresholding",
                "Thresholding")

add.notes = c("Exponential prior, parameters estimated via moment matching",
              "Non-overlapping block thresholding, James Stein rule, constant threshold",
              "Laplace prior, fixed hyperparameters",
              "Overlapping block thresholding, James Stein rule, constant threshold",
              "Spike and slab prior, parameters estimated via maximum likelihood",
              "Mixture Gaussian prior, mixing proportions estimated via maximum likelihood",
              "Mixture Gaussian prior, mixing proportions estimated via maximum likelihood",
              "Mixture Gaussian prior, mixing proportions estimated via maximum likelihood",
              "Mixture Gaussian prior, mixing proportions estimated via maximum likelihood",
              "Mixture Gaussian prior, mixing proportions estimated via maximum likelihood",
              "Mixture Gaussian prior, mixing proportions estimated via maximum likelihood",
              "Mixture Gaussian prior, mixing proportions estimated via maximum likelihood",
              "Level-dependent thresholding, soft thresholding rule, SureShrink threshold",
              "Level-dependent thresholding, hard thresholding rule, TI threshold",
              "Term by term thresholding, hard thresholding rule, modified TI threshold, sigma estimated via RMAD",
              "Term by term thresholding, hard thresholding rule, modified TI threshold, sigma estimated via RMAD",
              "Term by term thresholding, hard thresholding rule, TI threshold, sigma estimated via SMASH",
              "Term by term thresholding, hard thresholding rule, TI threshold, sigma estimated via SMASH",
              "Term by term thresholding, hard thresholding rule, TI threshold",
              "Term by term thresholding, hard thresholding rule, TI threshold" )


method.info = data.frame(method = mise.res$method[1:20],
                         software = software,
                         var.assumption = var.assumption,
                         basis = basis,
                         ti.status = ti.status,
                         shrink.type = shrink.type,
                         add.notes = add.notes)


res.final = merge(mise.res, method.info)
