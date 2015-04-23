library(xtable)

tex.row.names=c("SMASH",
"Haar-Fisz",
"Anscombe, universal",
"Anscombe, CV",
"Haar thresholds",
"CVS",
"BMMIM",
"BMSM"
)

tex.col.names=c("(0.01,3)","(1/8,8)","(1/128,128)")


mise.s.1.table=mise.s.1[c(1,4,8,7,9,5,3,2)]
mise.s.8.table=mise.s.8[c(1,4,8,7,9,5,3,2)]
mise.s.128.table=mise.s.128[c(1,4,8,7,9,5,3,2)]

mise.ang.1.table=mise.ang.1[c(1,4,8,7,9,5,3,2)]
mise.ang.8.table=mise.ang.8[c(1,4,8,7,9,5,3,2)]
mise.ang.128.table=mise.ang.128[c(1,4,8,7,9,5,3,2)]

mise.bur.1.table=mise.bur.1[c(1,4,8,7,9,5,3,2)]
mise.bur.8.table=mise.bur.8[c(1,4,8,7,9,5,3,2)]
mise.bur.128.table=mise.bur.128[c(1,4,8,7,9,5,3,2)]

mise.cb.1.table=mise.cb.1[c(1,4,8,7,9,5,3,2)]
mise.cb.8.table=mise.cb.8[c(1,4,8,7,9,5,3,2)]
mise.cb.128.table=mise.cb.128[c(1,4,8,7,9,5,3,2)]

mise.b.1.table=mise.b.1[c(1,4,8,7,9,5,3,2)]
mise.b.8.table=mise.b.8[c(1,4,8,7,9,5,3,2)]
mise.b.128.table=mise.b.128[c(1,4,8,7,9,5,3,2)]

mise.hs.1.table=mise.hs.1[c(1,4,8,7,9,5,3,2)]
mise.hs.8.table=mise.hs.8[c(1,4,8,7,9,5,3,2)]
mise.hs.128.table=mise.hs.128[c(1,4,8,7,9,5,3,2)]


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


print(xtable(mise.ang.table,caption="MISE of various methods for the Angles function",digits=1))
print(xtable(mise.b.table,caption="MISE of various methods for the Bumps function",digits=1))
