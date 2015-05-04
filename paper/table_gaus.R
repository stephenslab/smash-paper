library(xtable)

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

print(xtable(table.1.v1,caption="Relative MSE of various methods for iid Gaussian errors, SNR=1",digits=ndigits))
print(xtable(table.3.v1,caption="Relative MSE of various methods for iid Gaussian errors, SNR=3",digits=ndigits))


####

print(xtable(table.1.v2,caption="Relative MSE of various methods for variance function V2 as in Cai & Wang (2008), SNR=1, Gaussian errors",digits=ndigits))
print(xtable(table.3.v2,caption="Relative MSE of various methods for variance function V2 as in Cai & Wang (2008), SNR=3, Gaussian errors",digits=ndigits))
print(xtable(table.1.v3,caption="Relative MSE of various methods for the Doppler variance function, SNR=1, Gaussian errors",digits=ndigits))
print(xtable(table.3.v3,caption="Relative MSE of various methods for the Doppler variance function, SNR=3, Gaussian errors",digits=ndigits))
print(xtable(table.1.v4,caption="Relative MSE of various methods for the Bumps variance function, SNR=1, Gaussian errors",digits=ndigits))
print(xtable(table.3.v4,caption="Relative MSE of various methods for the Bumps variance function, SNR=3, Gaussian errors",digits=ndigits))
print(xtable(table.1.v5,caption="Relative MSE of various methods for the Clipped Blocks variance function, SNR=1, Gaussian errors",digits=ndigits))
print(xtable(table.3.v5,caption="Relative MSE of various methods for the Clipped Blocks variance function, SNR=3, Gaussian errors",digits=ndigits))
#####


table.1.v5.short=table.1.v5[c(-6,-7,-8),]
table.3.v5.short=table.3.v5[c(-6,-7,-8),]
table.1.v4.short=table.1.v4[c(-6,-7,-8),]
table.3.v4.short=table.3.v4[c(-6,-7,-8),]



print(xtable(table.1.v5.short,caption="Relative MSE of various methods for the Clipped Blocks variance function, SNR=1, Gaussian errors",digits=ndigits))
print(xtable(table.3.v5.short,caption="Relative MSE of various methods for the Clipped Blocks variance function, SNR=3, Gaussian errors",digits=ndigits))
print(xtable(table.1.v4.short,caption="Relative MSE of various methods for the Bumps variance function, SNR=1, Gaussian errors",digits=ndigits))
print(xtable(table.3.v4.short,caption="Relative MSE of various methods for the Bumps variance function, SNR=3, Gaussian errors",digits=ndigits))
