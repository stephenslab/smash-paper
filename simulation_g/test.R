wfmm.res=read.table("test_eg_beta.txt")
wfmm.res.quan=read.table("test_eg_beta_quan.txt")
wfmm.res=data.matrix(wfmm.res)
wfmm.res.quan=data.matrix(wfmm.res.quan)
load("sim_res.Robj")



pdf("test_effect.pdf")
par(mfrow=c(1,1))
plot(mu.spm-mu.sp,type='l',ylim=c(-1.2,1.6),main="1024 points, 10 vs 10 samples, snr=3",ylab="effect")

lines(res1$effect.mean,col=2)
lines(wfmm.res[2,],col=4)


lines(res1$effect.mean-2*sqrt(res1$effect.var),col=6)
lines(res1$effect.mean+2*sqrt(res1$effect.var),col=6)

lines(wfmm.res.quan[1,1025:2048],col=3)
lines(wfmm.res.quan[8,1025:2048],col=3)
dev.off()



pdf("test_baseline.pdf")
plot(mu.sp,type='l',ylim=c(-0.5,3.5),main="1024 points, 10 vs 10 samples, snr=3",ylab="baseline")

lines(res1$baseline.mean,col=2)
lines(wfmm.res[1,],col=4)

lines(res1$baseline.mean-2*sqrt(res1$baseline.var),col=6)
lines(res1$baseline.mean+2*sqrt(res1$baseline.var),col=6)

lines(wfmm.res.quan[1,1:1024],col=3)
lines(wfmm.res.quan[8,1:1024],col=3)
dev.off()
