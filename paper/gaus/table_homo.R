library(xtable)

tex.row.names=c("SMASH, joint, Haar",
"SMASH, joint, Symm8",
"TI, RMAD, Symm8",
"TI, SMASH, Symm8",
"SMASH, homo, Symm8",
"TI, homo, Symm8",
"Ebayes, homo, Symm8",
"SURE, homo, Symm8",
"PostMean, homo, Symm8",
"SMASH, truth, Symm8",
"TI, truth, Symm8"
)

tex.col.names=c("Spikes","Bumps","Block","Angles","Doppler","Blip","Corner")


mise.sp.1.v1.table=mise.sp.1.v1[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.sp.3.v1.table=mise.sp.3.v1[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.bump.1.v1.table=mise.bump.1.v1[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.bump.3.v1.table=mise.bump.3.v1[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.blk.1.v1.table=mise.blk.1.v1[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.blk.3.v1.table=mise.blk.3.v1[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.ang.1.v1.table=mise.ang.1.v1[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.ang.3.v1.table=mise.ang.3.v1[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.dop.1.v1.table=mise.dop.1.v1[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.dop.3.v1.table=mise.dop.3.v1[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.blip.1.v1.table=mise.blip.1.v1[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.blip.3.v1.table=mise.blip.3.v1[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.cor.1.v1.table=mise.cor.1.v1[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.cor.3.v1.table=mise.cor.3.v1[c(10,11,14,16,9,1,7,5,6,18,20)]

mise.1.v1.table=cbind(mise.sp.1.v1.table,
mise.bump.1.v1.table,
mise.blk.1.v1.table,
mise.ang.1.v1.table,
mise.dop.1.v1.table,
mise.blip.1.v1.table,
mise.cor.1.v1.table
)

mise.3.v1.table=cbind(mise.sp.3.v1.table,
mise.bump.3.v1.table,
mise.blk.3.v1.table,
mise.ang.3.v1.table,
mise.dop.3.v1.table,
mise.blip.3.v1.table,
mise.cor.3.v1.table
)

rownames(mise.1.v1.table)=tex.row.names
colnames(mise.1.v1.table)=tex.col.names


rownames(mise.3.v1.table)=tex.row.names
colnames(mise.3.v1.table)=tex.col.names


print(xtable(mise.1.v1.table,caption="MISE of various methods for iid Gaussian errors, RSNR=1",digits=1))
print(xtable(mise.3.v1.table,caption="MISE of various methods for iid Gaussian errors, RSNR=3",digits=1))
