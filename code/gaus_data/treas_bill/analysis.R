data=read.csv("D:/Grad School/Spring 2013/multiscale_ash/gaus_data/treas_bill/treas_bill_ave.csv",header=FALSE)
data=as.numeric(data[,1])


library(smash)


y=ar(data,FALSE,5)$res
y=y[!is.na(y)]
x=sort(data[5:(length(data)-1)])
y=y[order(data[5:(length(data)-1)])]


#coef=c(1.0733,-0.0423,0.0165,0.0228,-0.0773)

resid.comp=function(data,coef){
  smean=mean(data)
  ord=length(coef)
  res=NULL
  for(i in (length(coef)+1):length(data)){
    res[i-ord]=(data[i]-smean)-(sum(coef*(data[(i-1):(i-ord)]-smean)))
  }
  return(res)
}

#y=resid.comp(data,coef)
#x=sort(data[5:(length(data)-1)])
#y=y[order(data[5:(length(data)-1)])]


plot(x,y)


####
xmin=min(x)
xmax=max(x)

bin.size=(xmax-xmin)/32
x.seq=seq(xmin,xmax,bin.size)
x.seq.l=x.seq[1:32]
x.seq.u=x.seq[2:33]
y.agg=0
for(i in 1:32){
  y.agg[i]=y[x>=x.seq.l[i]&x<x.seq.u[i]]
}
y.agg[1024]=y[x>=x.seq.l[1024]&x<=x.seq.u[1024]]

y.approx=approx(x,y,n=1024,method='constant',rule=2)$y
x.approx=approx(x,y,n=1024,method='constant',rule=2)$x

lines(x.approx,y.approx,col=2)
y.est.approx=bayesmooth(y.approx)
y.est.var.approx=bayesmooth(y.approx,v.est=TRUE)

plot(x,y)
lines(x.approx,y.est.approx,col=2)
####

x.mod=unique(x)
y.mod=0
for(i in 1:length(x.mod)){
  y.mod[i]=median(y[x==x.mod[i]])
}


y.exp=c(y.mod,y.mod[length(y.mod):(2*length(y.mod)-2^10+1)])
y.final=c(y.exp,y.exp[length(y.exp):1])

y.est=ashsmooth.gaus(y.final,v.est=TRUE,joint=TRUE,post.var=TRUE)



mu.est.final=y.est$mu.res$mu.est[1:(length(y.mod))]
mu.est.ci=list(y.est$mu.res$mu.est[1:(length(y.mod))]+2*sqrt(y.est$mu.res$mu.est.var[1:(length(y.mod))]),y.est$mu.res$mu.est[1:(length(y.mod))]-2*sqrt(y.est$mu.res$mu.est.var[1:(length(y.mod))]))
var.est.final=y.est$var.res$var.est[1:(length(y.mod))]
var.est.ci=list(pmax(0,y.est$var.res$var.est[1:(length(y.mod))]-2*sqrt(y.est$var.res$var.est.var[1:(length(y.mod))])),y.est$var.res$var.est[1:(length(y.mod))]+2*sqrt(y.est$var.res$var.est.var[1:(length(y.mod))]))


pdf("figures_treasury_a.pdf")
plot(data,type='l',axes=FALSE,xlab="Year",ylab="Interest Rate (%)")
axis(1,at=c(417,939,1461),labels=c("1970","1980","1990"))
axis(2)
dev.off()

pdf("figures_treasury_b.pdf")
plot(x,y,xlab="X(t)",ylab="Y(t)",cex=0.5)
lines(x.mod,mu.est.final,col=2)
dev.off()

pdf("figures_treasury_c.pdf")
plot(x.mod,mu.est.final,type='l',ylim=c(-0.4,0.4),xlab="X",ylab="Esimated mean")
lines(x.mod,mu.est.ci[[1]],lty=2)
lines(x.mod,mu.est.ci[[2]],lty=2)
dev.off()

pdf("figures_treasury_d.pdf")
plot(x.mod,var.est.final,type='l',ylim=c(0,0.6),xlab="X",ylab="Estimated variance")
lines(x.mod,var.est.ci[[1]],lty=2)
lines(x.mod,var.est.ci[[2]],lty=2)
dev.off()


lv=log(sqrt(var.est.final))
lx=log(x.mod)
cor(lv,lx)
lm(lv~lx)


