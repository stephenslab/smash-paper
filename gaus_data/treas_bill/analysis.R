data=read.csv("gaus_data/treas_bill/treas_bill.csv",header=FALSE)
data=as.numeric(data[,1])

library(tseries)
source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/bayesmooth.R"))


y=ar(data,FALSE,5)$res
y=y[!is.na(y)]
x=sort(data[5:(length(data)-1)])
y=y[order(data[5:(length(data)-1)])]

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

y.est=bayesmooth(y.final)
y.est.var=bayesmooth(y.final,v.est=TRUE)

plot(x,y)
lines(x.mod,y.est[1:(length(y.mod))],col=2)

plot(x.mod,sqrt(y.est.var[1:(length(y.mod))]),type='l')



