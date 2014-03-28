source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/bayesmooth.R"))
source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/bayesmooth_alt2.R"))

library(MASS)

x.ini=sort(mcycle$times)
y.ini=mcycle$accel[order(mcycle$times)]

x=unique(x.ini)
y=0
for(i in 1:length(x)){
  y[i]=median(y.ini[x.ini==x[i]])
}


y.exp=c(y,y[length(y):(2*length(y)-128+1)])
y.final=c(y.exp,y.exp[length(y.exp):1])

y.est=bayesmooth(y.final)
y.est.alt2=bayesmooth.alt2(y.final)

y.est=y.est[1:length(y)]
y.est.alt2=y.est.alt2[1:length(y)]

y.var.est=bayesmooth(y.final,v.est=TRUE)
y.var.est.alt2=bayesmooth.alt2(y.final,v.est=TRUE)

y.var.est=y.var.est[1:length(y)]
y.var.est.alt2=y.var.est.alt2[1:length(y)]

plot(x,y,ylim=c(-150,100))
lines(x,y.est,col=2)
plot(x,y.var.est,type='l')




