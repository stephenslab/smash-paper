#This script produces supplementary tables for Poisson simulations

library(xtable)

load("res_paper/res_pois.RData")

scale.baseline = function(x) x/x[1]

tex.row.names=c("SMASH",
                "BMSM",
                "BMMIM",
                "Haar-Fisz",
                "Anscombe",
                "Haar thresholds",
                "L1_penalty"
)

tex.col.names=c("intensity (0.01,3)","intensity (1/8,8)","intensity (1/128,128)")


# mise.s.1.table=mise.s.1[c(1,4,8,7,9,5,3,2)]
# mise.s.8.table=mise.s.8[c(1,4,8,7,9,5,3,2)]
# mise.s.128.table=mise.s.128[c(1,4,8,7,9,5,3,2)]
# 
# mise.ang.1.table=mise.ang.1[c(1,4,8,7,9,5,3,2)]
# mise.ang.8.table=mise.ang.8[c(1,4,8,7,9,5,3,2)]
# mise.ang.128.table=mise.ang.128[c(1,4,8,7,9,5,3,2)]
# 
# mise.bur.1.table=mise.bur.1[c(1,4,8,7,9,5,3,2)]
# mise.bur.8.table=mise.bur.8[c(1,4,8,7,9,5,3,2)]
# mise.bur.128.table=mise.bur.128[c(1,4,8,7,9,5,3,2)]
# 
# mise.cb.1.table=mise.cb.1[c(1,4,8,7,9,5,3,2)]
# mise.cb.8.table=mise.cb.8[c(1,4,8,7,9,5,3,2)]
# mise.cb.128.table=mise.cb.128[c(1,4,8,7,9,5,3,2)]
# 
# mise.b.1.table=mise.b.1[c(1,4,8,7,9,5,3,2)]
# mise.b.8.table=mise.b.8[c(1,4,8,7,9,5,3,2)]
# mise.b.128.table=mise.b.128[c(1,4,8,7,9,5,3,2)]
# 
# mise.hs.1.table=mise.hs.1[c(1,4,8,7,9,5,3,2)]
# mise.hs.8.table=mise.hs.8[c(1,4,8,7,9,5,3,2)]
# mise.hs.128.table=mise.hs.128[c(1,4,8,7,9,5,3,2)]

mise.s.1.table=mise.s.1[c(1,2,3,4,6,7,8)]
mise.s.8.table=mise.s.8[c(1,2,3,4,6,7,8)]
mise.s.128.table=mise.s.128[c(1,2,3,4,6,7,8)]

mise.ang.1.table=mise.ang.1[c(1,2,3,4,6,7,8)]
mise.ang.8.table=mise.ang.8[c(1,2,3,4,6,7,8)]
mise.ang.128.table=mise.ang.128[c(1,2,3,4,6,7,8)]

mise.bur.1.table=mise.bur.1[c(1,2,3,4,6,7,8)]
mise.bur.8.table=mise.bur.8[c(1,2,3,4,6,7,8)]
mise.bur.128.table=mise.bur.128[c(1,2,3,4,6,7,8)]

mise.cb.1.table=mise.cb.1[c(1,2,3,4,6,7,8)]
mise.cb.8.table=mise.cb.8[c(1,2,3,4,6,7,8)]
mise.cb.128.table=mise.cb.128[c(1,2,3,4,6,7,8)]

mise.b.1.table=mise.b.1[c(1,2,3,4,6,7,8)]
mise.b.8.table=mise.b.8[c(1,2,3,4,6,7,8)]
mise.b.128.table=mise.b.128[c(1,2,3,4,6,7,8)]

mise.hs.1.table=mise.hs.1[c(1,2,3,4,6,7,8)]
mise.hs.8.table=mise.hs.8[c(1,2,3,4,6,7,8)]
mise.hs.128.table=mise.hs.128[c(1,2,3,4,6,7,8)]


mise.s.table=cbind(mise.s.1.table,
mise.s.8.table,
mise.s.128.table
)

mise.ang.table=cbind(mise.ang.1.table,
mise.ang.8.table,
mise.ang.128.table
)

mise.bur.table=cbind(mise.bur.1.table,
mise.bur.8.table,
mise.bur.128.table
)

mise.cb.table=cbind(mise.cb.1.table,
mise.cb.8.table,
mise.cb.128.table
)

mise.b.table=cbind(mise.b.1.table,
mise.b.8.table,
mise.b.128.table
)

mise.hs.table=cbind(mise.hs.1.table,
mise.hs.8.table,
mise.hs.128.table
)



rownames(mise.ang.table)=tex.row.names
colnames(mise.ang.table)=tex.col.names
rownames(mise.b.table)=tex.row.names
colnames(mise.b.table)=tex.col.names
rownames(mise.cb.table)=tex.row.names
colnames(mise.cb.table)=tex.col.names
rownames(mise.s.table)=tex.row.names
colnames(mise.s.table)=tex.col.names
rownames(mise.hs.table)=tex.row.names
colnames(mise.hs.table)=tex.col.names
rownames(mise.bur.table)=tex.row.names
colnames(mise.bur.table)=tex.col.names


# mise.ang.table=apply(mise.ang.table,2,scale.baseline)
# mise.b.table=apply(mise.b.table,2,scale.baseline)
# mise.cb.table=apply(mise.cb.table,2,scale.baseline)
# mise.s.table=apply(mise.s.table,2,scale.baseline)
# mise.hs.table=apply(mise.hs.table,2,scale.baseline)
# mise.bur.table=apply(mise.bur.table,2,scale.baseline)

print(xtable(mise.s.table,caption="Comparison of methods for denoising Poisson data for the ``Spikes'' test function for 3 different (min,max) intensities ((0.01,3), (1/8,8), (1/128,128)). Performance is measured using MISE over 100 independent datasets, with smaller values indicating better performance. Values colored in red indicates the smallest MISE amongst all methods for a given (min, max) intensity.",label="table:pois_sp",digits=2))
print(xtable(mise.ang.table,caption="Comparison of methods for denoising Poisson data for the ``Angles'' test function for 3 different (min,max) intensities ((0.01,3), (1/8,8), (1/128,128)). Performance is measured using MISE over 100 independent datasets, with smaller values indicating better performance. Values colored in red indicates the smallest MISE amongst all methods for a given (min, max) intensity.",label="table:pois_ang",digits=2))
print(xtable(mise.hs.table,caption="Comparison of methods for denoising Poisson data for the ``Heavisine'' test function for 3 different (min,max) intensities ((0.01,3), (1/8,8), (1/128,128)). Performance is measured using MISE over 100 independent datasets, with smaller values indicating better performance. Values colored in red indicates the smallest MISE amongst all methods for a given (min, max) intensity.",label="table:pois_hs",digits=2))
print(xtable(mise.bur.table,caption="Comparison of methods for denoising Poisson data for the ``Bursts'' test function for 3 different (min,max) intensities ((0.01,3), (1/8,8), (1/128,128)). Performance is measured using MISE over 100 independent datasets, with smaller values indicating better performance. Values colored in red indicates the smallest MISE amongst all methods for a given (min, max) intensity.",label="table:pois_bur",digits=2))
print(xtable(mise.cb.table,caption="Comparison of methods for denoising Poisson data for the ``Clipped Blocks'' test function for 3 different (min,max) intensities ((0.01,3), (1/8,8), (1/128,128)). Performance is measured using MISE over 100 independent datasets, with smaller values indicating better performance. Values colored in red indicates the smallest MISE amongst all methods for a given (min, max) intensity.",label="table:pois_cb",digits=2))
print(xtable(mise.b.table,caption="Comparison of methods for denoising Poisson data for the ``Bumps'' test function for 3 different (min,max) intensities ((0.01,3), (1/8,8), (1/128,128)). Performance is measured using MISE over 100 independent datasets, with smaller values indicating better performance. Values colored in red indicates the smallest MISE amongst all methods for a given (min, max) intensity.",label="table:pois_b",digits=2))

