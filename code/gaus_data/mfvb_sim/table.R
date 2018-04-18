library(xtable)

mse.table=rbind(c(mean(mse.mu.uneven.mfvb),mean(mse.sd.uneven.mfvb),mean(mse.mu.even.mfvb),mean(mse.sd.even.mfvb)),
c(mean(mse.mu.uneven.s8),mean(mse.sd.uneven.s8),mean(mse.mu.even.s8),mean(mse.sd.even.s8)))

rownames(mse.table)=c("MFVB","SMASH")
colnames(mse.table)=c("mean","sd","mean","sd")

print(xtable(mse.table,caption="MSEs of MFVB and SMASH for two simulation scenarios",digits=4))




