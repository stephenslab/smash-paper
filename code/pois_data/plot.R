library(smash)
library(scales)

load("reg_880000_1011072.Robj")

bppos=880001:1011072

res1=ashsmooth.pois(M[1,])
res2=ashsmooth.pois(M[2,])
res3=ashsmooth.pois(M[3,])
res4=ashsmooth.pois(M[4,])


pdf("smoothing_all.pdf")
par(mfrow=c(4,1))
plot(bppos,res1,type='l',xlab="position",ylab="intensity",ylim=c(0,5.5),main="Cell GM12878, Replicate 1")
plot(bppos,res2,type='l',xlab="position",ylab="intensity",ylim=c(0,5.5),main="Cell GM12878, Replicate 2")
plot(bppos,res3,type='l',xlab="position",ylab="intensity",ylim=c(0,5.5),main="Cell H1hesc, Replicate 1")
plot(bppos,res4,type='l',xlab="position",ylab="intensity",ylim=c(0,5.5),main="Cell H1hesc, Replicate 2")
dev.off()


peaks=read.table("pois_data/Gm1287peaks_chr1_sorted")
peaks=peaks[order(peaks[,2]),]
all.base=1:peaks[dim(peaks)[1],3]
peaks.start=peaks[,2][peaks[,2]>=880000&peaks[,2]<=1011072]
peaks.end=peaks[,3][peaks[,3]>=880000&peaks[,3]<=1011072]

res=ashsmooth.pois(M[1,]+M[2,],post.var=TRUE)

pdf("peaks_comp_a.pdf",height=5,width=10)
plot(bppos,M[1,]+M[2,],xlab="position",ylab="counts",pch=16,cex=0.5,col=alpha("black",0.04))
dev.off()

pdf("peaks_comp_b.pdf",height=5,width=10)
plot(bppos,res$est,type='l',xlab="position",ylab="intensity")
for(i in 1:length(peaks.start)){
  segments(peaks.start[i],-0.02,peaks.end[i],-0.02,col=2,lwd=4)
}
dev.off()



#lines((res$est+3*sqrt(res$var)),col=2)
#lines((res$est-3*sqrt(res$var)),col=2)

#res.thresh=res$est
#res.thresh[(res$est+3*sqrt(res$var)>=0)&(res$est-3*sqrt(res$var)<=0)]=0
#plot(res.thresh,type='l')

#res.thresh[res$est<=0.2]=0
#plot(res.thresh,type='l')
