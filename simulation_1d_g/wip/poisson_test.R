spike.f=function(x) (0.75*exp(-500*(x-0.23)^2)+1.5*exp(-2000*(x-0.33)^2)+3*exp(-8000*(x-0.47)^2)+2.25*exp(-16000*(x-0.69)^2)+0.5*exp(-32000*(x-0.83)^2))
#spike.f=function(x) (0.75*exp(-500*(x-0.23)^2)+1.5*exp(-1000*(x-0.33)^2)+3*exp(-500*(x-0.47)^2)+2.25*exp(-500*(x-0.69)^2)+0.5*exp(-500*(x-0.83)^2))

n=4096
t=1:n/n
mu.s=spike.f(1:n/n)
var2=(1+4*(exp(-550*(t-0.2)^2)+exp(-200*(t-0.5)^2)+exp(-950*(t-0.8)^2)))/1.35

mu.s0=100
set.seed(415)
X.s=rpois(n,mu.s0*(1+var2))
X.ss=rnorm(n,mu.s0*(var2),sqrt(mu.s0*(var2)))

library(multiseq)
source("D:/Grad School/Spring 2013/multiscale_ash/ash/bayesmooth.R")

mu.est.p=binnorm.smooth(X.s)
mu.est.n=bayesmooth(X.s)
var.est.n=bayesmooth(X.s,v.est=TRUE)


mu.est1.p=binnorm.smooth(X.ss)
mu.est1.n=bayesmooth(X.ss)
var.est1.n=bayesmooth(X.ss,v.est=TRUE)


plot(mu.s0*(1+var2),type='l')
lines(mu.est.p,col=2)
lines(mu.est.n,col=4)
lines(var.est.n,col=3)


plot(mu.s0*var2,type='l')
lines(mu.est1.p,col=2)
lines(mu.est1.n,col=4)
lines(var.est1.n,col=2)