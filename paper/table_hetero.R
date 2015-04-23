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


mise.sp.1.v2.table=mise.sp.1.v2[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.sp.3.v2.table=mise.sp.3.v2[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.bump.1.v2.table=mise.bump.1.v2[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.bump.3.v2.table=mise.bump.3.v2[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.blk.1.v2.table=mise.blk.1.v2[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.blk.3.v2.table=mise.blk.3.v2[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.ang.1.v2.table=mise.ang.1.v2[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.ang.3.v2.table=mise.ang.3.v2[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.dop.1.v2.table=mise.dop.1.v2[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.dop.3.v2.table=mise.dop.3.v2[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.blip.1.v2.table=mise.blip.1.v2[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.blip.3.v2.table=mise.blip.3.v2[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.cor.1.v2.table=mise.cor.1.v2[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.cor.3.v2.table=mise.cor.3.v2[c(10,11,14,16,9,1,7,5,6,18,20)]

mise.1.v2.table=cbind(mise.sp.1.v2.table,
mise.bump.1.v2.table,
mise.blk.1.v2.table,
mise.ang.1.v2.table,
mise.dop.1.v2.table,
mise.blip.1.v2.table,
mise.cor.1.v2.table
)

mise.3.v2.table=cbind(mise.sp.3.v2.table,
mise.bump.3.v2.table,
mise.blk.3.v2.table,
mise.ang.3.v2.table,
mise.dop.3.v2.table,
mise.blip.3.v2.table,
mise.cor.3.v2.table
)

rownames(mise.1.v2.table)=tex.row.names
colnames(mise.1.v2.table)=tex.col.names



rownames(mise.3.v2.table)=tex.row.names
colnames(mise.3.v2.table)=tex.col.names




mise.sp.1.v3.table=mise.sp.1.v3[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.sp.3.v3.table=mise.sp.3.v3[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.bump.1.v3.table=mise.bump.1.v3[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.bump.3.v3.table=mise.bump.3.v3[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.blk.1.v3.table=mise.blk.1.v3[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.blk.3.v3.table=mise.blk.3.v3[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.ang.1.v3.table=mise.ang.1.v3[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.ang.3.v3.table=mise.ang.3.v3[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.dop.1.v3.table=mise.dop.1.v3[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.dop.3.v3.table=mise.dop.3.v3[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.blip.1.v3.table=mise.blip.1.v3[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.blip.3.v3.table=mise.blip.3.v3[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.cor.1.v3.table=mise.cor.1.v3[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.cor.3.v3.table=mise.cor.3.v3[c(10,11,14,16,9,1,7,5,6,18,20)]

mise.1.v3.table=cbind(mise.sp.1.v3.table,
mise.bump.1.v3.table,
mise.blk.1.v3.table,
mise.ang.1.v3.table,
mise.dop.1.v3.table,
mise.blip.1.v3.table,
mise.cor.1.v3.table
)

mise.3.v3.table=cbind(mise.sp.3.v3.table,
mise.bump.3.v3.table,
mise.blk.3.v3.table,
mise.ang.3.v3.table,
mise.dop.3.v3.table,
mise.blip.3.v3.table,
mise.cor.3.v3.table
)

rownames(mise.1.v3.table)=tex.row.names
colnames(mise.1.v3.table)=tex.col.names


rownames(mise.3.v3.table)=tex.row.names
colnames(mise.3.v3.table)=tex.col.names






mise.sp.1.v4.table=mise.sp.1.v4[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.sp.3.v4.table=mise.sp.3.v4[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.bump.1.v4.table=mise.bump.1.v4[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.bump.3.v4.table=mise.bump.3.v4[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.blk.1.v4.table=mise.blk.1.v4[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.blk.3.v4.table=mise.blk.3.v4[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.ang.1.v4.table=mise.ang.1.v4[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.ang.3.v4.table=mise.ang.3.v4[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.dop.1.v4.table=mise.dop.1.v4[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.dop.3.v4.table=mise.dop.3.v4[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.blip.1.v4.table=mise.blip.1.v4[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.blip.3.v4.table=mise.blip.3.v4[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.cor.1.v4.table=mise.cor.1.v4[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.cor.3.v4.table=mise.cor.3.v4[c(10,11,14,16,9,1,7,5,6,18,20)]

mise.1.v4.table=cbind(mise.sp.1.v4.table,
mise.bump.1.v4.table,
mise.blk.1.v4.table,
mise.ang.1.v4.table,
mise.dop.1.v4.table,
mise.blip.1.v4.table,
mise.cor.1.v4.table
)

mise.3.v4.table=cbind(mise.sp.3.v4.table,
mise.bump.3.v4.table,
mise.blk.3.v4.table,
mise.ang.3.v4.table,
mise.dop.3.v4.table,
mise.blip.3.v4.table,
mise.cor.3.v4.table
)

rownames(mise.1.v4.table)=tex.row.names
colnames(mise.1.v4.table)=tex.col.names


rownames(mise.3.v4.table)=tex.row.names
colnames(mise.3.v4.table)=tex.col.names





mise.sp.1.v5.table=mise.sp.1.v5[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.sp.3.v5.table=mise.sp.3.v5[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.bump.1.v5.table=mise.bump.1.v5[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.bump.3.v5.table=mise.bump.3.v5[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.blk.1.v5.table=mise.blk.1.v5[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.blk.3.v5.table=mise.blk.3.v5[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.ang.1.v5.table=mise.ang.1.v5[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.ang.3.v5.table=mise.ang.3.v5[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.dop.1.v5.table=mise.dop.1.v5[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.dop.3.v5.table=mise.dop.3.v5[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.blip.1.v5.table=mise.blip.1.v5[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.blip.3.v5.table=mise.blip.3.v5[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.cor.1.v5.table=mise.cor.1.v5[c(10,11,14,16,9,1,7,5,6,18,20)]
mise.cor.3.v5.table=mise.cor.3.v5[c(10,11,14,16,9,1,7,5,6,18,20)]

mise.1.v5.table=cbind(mise.sp.1.v5.table,
mise.bump.1.v5.table,
mise.blk.1.v5.table,
mise.ang.1.v5.table,
mise.dop.1.v5.table,
mise.blip.1.v5.table,
mise.cor.1.v5.table
)

mise.3.v5.table=cbind(mise.sp.3.v5.table,
mise.bump.3.v5.table,
mise.blk.3.v5.table,
mise.ang.3.v5.table,
mise.dop.3.v5.table,
mise.blip.3.v5.table,
mise.cor.3.v5.table
)

rownames(mise.1.v5.table)=tex.row.names
colnames(mise.1.v5.table)=tex.col.names


rownames(mise.3.v5.table)=tex.row.names
colnames(mise.3.v5.table)=tex.col.names


print(xtable(mise.1.v2.table,caption="MISE of various methods for variance function V2 as in Cai & Wang (2008), RSNR=1",digits=1))
print(xtable(mise.3.v2.table,caption="MISE of various methods for variance function V2 as in Cai & Wang (2008), RSNR=3",digits=1))
print(xtable(mise.1.v3.table,caption="MISE of various methods for the Doppler variance function, RSNR=1",digits=1))
print(xtable(mise.3.v3.table,caption="MISE of various methods for the Doppler variance function, RSNR=3",digits=1))
print(xtable(mise.1.v4.table,caption="MISE of various methods for the Bumps variance function, RSNR=1",digits=1))
print(xtable(mise.3.v4.table,caption="MISE of various methods for the Bumps variance function, RSNR=3",digits=1))
print(xtable(mise.1.v5.table,caption="MISE of various methods for the Clipped Blocks variance function, RSNR=1",digits=1))
print(xtable(mise.3.v5.table,caption="MISE of various methods for the Clipped Blocks variance function, RSNR=3",digits=1))






mise.1.v2.table.short=mise.1.v2.table[c(-7,-8,-9),]
mise.3.v2.table.short=mise.3.v2.table[c(-7,-8,-9),]
mise.1.v4.table.short=mise.1.v4.table[c(-7,-8,-9),]
mise.3.v4.table.short=mise.3.v4.table[c(-7,-8,-9),]



print(xtable(mise.1.v2.table.short,caption="MISE of various methods for variance function V2 as in Cai & Wang (2008), RSNR=1",digits=1))
print(xtable(mise.3.v2.table.short,caption="MISE of various methods for variance function V2 as in Cai & Wang (2008), RSNR=3",digits=1))
print(xtable(mise.1.v4.table.short,caption="MISE of various methods for the Bumps variance function, RSNR=1",digits=1))
print(xtable(mise.3.v4.table.short,caption="MISE of various methods for the Bumps variance function, RSNR=3",digits=1))




