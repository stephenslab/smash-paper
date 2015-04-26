mise.res = aggregate(mise ~ method + scenario, res, mean)
table.order = c(6, 10, 16, 17, 8, 3, 5, 13, 14, 12, 20)


#table 1
#scenario.name = sapply(mise.res$scenario, strsplit, "[.]")
name.ind = grep("*.1.v1", mise.res$scenario)
table.1.v1 = matrix(mise.res[name.ind, 3], 20, 7)
table.1.v1 = table.1.v1[table.order, ]

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

rownames(table.1.v1) = tex.row.names
colnames(table.1.v1) = tex.col.names


