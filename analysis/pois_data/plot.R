# TO DO: Explain here what this script does.
library(smashr)
library(scales)

load("../../data/reg_880000_1011072.RData")

bppos=880001:1011072

res1=smash.poiss(M[1,])
res2=smash.poiss(M[2,])
res3=smash.poiss(M[3,])
res4=smash.poiss(M[4,])

par(mfrow=c(4,1))
plot(bppos,res1,type='l',xlab="position",ylab="intensity",
     ylim=c(0,5.5),main="Cell GM12878, Replicate 1")
plot(bppos,res2,type='l',xlab="position",ylab="intensity",
     ylim=c(0,5.5),main="Cell GM12878, Replicate 2")
plot(bppos,res3,type='l',xlab="position",ylab="intensity",
     ylim=c(0,5.5),main="Cell H1hesc, Replicate 1")
plot(bppos,res4,type='l',xlab="position",ylab="intensity",
     ylim=c(0,5.5),main="Cell H1hesc, Replicate 2")

peaks=read.table("../../data/Gm1287peaks_chr1_sorted.txt")
peaks=peaks[order(peaks[,2]),]
all.base=1:peaks[dim(peaks)[1],3]
peaks.start=peaks[,2][peaks[,2]>=880000&peaks[,2]<=1011072]
peaks.end=peaks[,3][peaks[,3]>=880000&peaks[,3]<=1011072]

res=smash.poiss(M[1,]+M[2,],post.var=TRUE)

plot(bppos,M[1,]+M[2,],xlab="position",ylab="counts",pch=16,cex=0.5,
     col=alpha("black",0.04))

plot(bppos,res$est,type='l',xlab="position",ylab="intensity")
for(i in 1:length(peaks.start))
  segments(peaks.start[i],-0.02,peaks.end[i],-0.02,col=2,lwd=4)


