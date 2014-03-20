source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/bayesmooth.R"))
source("simulation_1d_g/threshold_var.R")
source("simulation_1d_g/threshold_mod.R")
source("simulation_1d_g/wd_var.R")
library(wavethresh)
library(caTools)


la10.sure3=function (x) 
{
    x.w <- wd(x)
    x.w.t <- threshold.wd.mod(x.w)
    x.w.t.r <- wr(x.w.t)
    return(x.w.t.r)
}

la10.ti3u=function (x, filter.number = 10, family = "DaubLeAsymm", min.level = 3, noise.level) 
{
    TT <- length(x)
    thresh <- noise.level * sqrt(2 * log(TT))    
    x.w <- wd(x, filter.number, family, type = "station")
    x.w.t <- threshold(x.w, levels = (min.level):(x.w$nlevels - 1), policy = "manual", value = thresh, type = "hard")
    x.w.t.r <- AvBasis(convert(x.w.t))
    return(x.w.t.r)
}


waveti.var=function (x, sigma=NULL, filter.number = 8, family = "DaubLeAsymm", min.level = 3) 
{
    n <- length(x)
    J <- log2(n)
    if(length(sigma)==1) sigma <- rep(sigma,n)
    x.w <- wd(x, filter.number, family, type = "station")
    if(is.null(sigma)){
      win.size <- round(n/10)
      odd.boo <- (win.size%%2==1)
      win.size <- win.size+(1-odd.boo)
      sigma <- runmad(accessD(x.w,J-1),win.size,endrule="func")
    }
    x.w.v <- wd.var(sigma^2, filter.number, family,type='station')
    x.w.t <- threshold.wd.var(x.w, x.w.v, levels = (min.level):(J-1))
    x.w.t.r <- AvBasis(convert(x.w.t))
    return(x.w.t.r)
}


mse=function(x,y) mean((x-y)^2)
l2norm=function(x) sum(x^2)
mise=function(x,y) 10000*mean(apply(x-rep(1,100)%o%y,1,l2norm)/l2norm(y))

spike.f=function(x) (0.75*exp(-500*(x-0.23)^2)+1.5*exp(-2000*(x-0.33)^2)+3*exp(-8000*(x-0.47)^2)+2.25*exp(-16000*(x-0.69)^2)+0.5*exp(-32000*(x-0.83)^2))
cor.f=function(t) 623.87*t^3*(1-2*t)*(t>=0&t<=0.5)+187.161*(0.125-t^3)*t^4*(t>0.5&t<=0.8)+3708.470441*(t-1)^3*(t>0.8&t<=1)
n=1024
t=sort(rnorm(1024,0,1))
t=(t-min(t))/(max(t)-min(t))
t=1:n/n
mu.s=spike.f(t)
mu.hsm=4 * sin(4 * pi * t) - 2*sign(t - 0.3) - 2*sign(0.72 - t)

mu.t=(1+mu.s)/5
var1=rep(1,n)
var2=(0.0001+4*(exp(-550*(t-0.2)^2)+exp(-200*(t-0.5)^2)+exp(-950*(t-0.8)^2)))/1.35
sigma.ini=sqrt(var1)
rsnr=sqrt(1)

sigma.t=sigma.ini/mean(sigma.ini)*sd(mu.t)/rsnr^2

plot(t,mu.t)

set.seed(1025)
X.s=rnorm(n,mu.t,sigma.t)
set.seed(1107)
X.s=matrix(rnorm(100*n,mu.t,sigma.t),nrow=100,byrow=TRUE)


###single
system.time(mu.est<-bayesmooth(X.s))
mu.est.sure3=la10.sure3(X.s)
mu.est.ti3=la10.ti3u(X.s,noise.level=sigma.t[1])
mu.est.waveti=waveti.var(X.s)
X.s.new=inter.samp(t,X.s,1024)$y
t.new=inter.samp(t,X.s,1024)$x
mu.est.samp=bayesmooth(X.s.new)
mu.est.ti3.samp=la10.ti3u(X.s.new,noise.level=sigma.t[1])


xgrid=1:n/n
mu.est.inter=approx(t,mu.est,xgrid)$y
mu.est.ti3.inter=approx(t,mu.est.ti3,xgrid)$y
mu.est.waveti.inter=approx(t,mu.est.waveti,xgrid)$y
mu.t.inter=(spike.f(xgrid)+1)/5

mse(mu.est.inter,mu.t.inter)
mse(mu.est.samp,mu.t.inter)
mse(mu.est.ti3.inter,mu.t.inter)
mse(mu.est.ti3.samp,mu.t.inter)
mse(mu.est.waveti.inter,mu.t.inter)

plot(xgrid,mu.t.inter,type='l')
lines(xgrid,mu.est.inter,col=2)
lines(xgrid,mu.est.samp,col=4)

lines(xgrid,mu.est.ti3.inter,col=3)
lines(xgrid,mu.est.ti3.samp,col=4)


###multiple
mu.est=apply(X.s,1,bayesmooth)
mu.est.ti3=apply(X.s,1,la10.ti3u,noise.level=sigma.t[1])
mu.est.waveti=apply(X.s,1,waveti.var)


xgrid=1:n/n
mu.est.inter=matrix(0,100,n)
mu.est.ti3.inter=matrix(0,100,n)
mu.est.waveti.inter=matrix(0,100,n)
for(i in 1:100){
  mu.est.inter[i,]=approx(t,mu.est[,i],xgrid)$y
  mu.est.ti3.inter[i,]=approx(t,mu.est.ti3[,i],xgrid)$y
  mu.est.waveti.inter[i,]=approx(t,mu.est.waveti[,i],xgrid)$y
}
mu.t.inter=(spike.f(xgrid)+1)/5

mise(mu.est.inter,mu.t.inter)
mise(mu.est.ti3.inter,mu.t.inter)
mise(mu.est.waveti.inter,mu.t.inter)
