spike.f=function(x) (0.75*exp(-500*(x-0.23)^2)+1.5*exp(-2000*(x-0.33)^2)+3*exp(-8000*(x-0.47)^2)+2.25*exp(-16000*(x-0.69)^2)+0.5*exp(-32000*(x-0.83)^2))
n=4096
t=1:n/n
mu.s=spike.f(t)
mu.sine=sin(20*t)
var1=1.42*((3-20*t)*(t>=0&t<0.1)+(20*t-1)*(t>=0.1&t<0.25)+(4+(1-4*t)*18/19)*(t>=0.25&t<0.725)+(2.2+10*(t-0.725))*(t>=0.725&t<0.89)+(3.85-85*(t-0.89)/11)*(t>=0.89&t<=1))
var2=(1+4*(exp(-550*(t-0.2)^2)+exp(-200*(t-0.5)^2)+exp(-950*(t-0.8)^2)))/1.35
var3=3.4*(2+sqrt(t*(1-t))*sin((2*pi*(1-0.05))/(t+0.05)))
var5=(t-0.5)^2+0.5
var6=exp(-100*(t-0.5)^2)
sigma.t=sqrt(var3)
#sigma.t=seq(2,10,length.out=n)
#sigma.t=rep(5,n)
l2norm=function(x) sum(x^2)
mise=function(x,y) l2norm(x-y)/l2norm(y)
mse=function(x,y) mean((x-y)^2)
mu.s0=100
set.seed(444)
X.s=rnorm(n,var3,sigma.t)
#X.s=matrix(rnorm(500*n,0,sigma.t),nrow=500,byrow=TRUE)


source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/cs.smooth.R"))
source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/var.smooth.R"))
source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/var_test.R"))
mu.est=cs.smooth(X.s,1)

system.time(var.est.1<-var.test(X.s,1))
system.time(var.est.2<-var.test(X.s,2))
system.time(var.est.3<-var.smooth(X.s,3))
system.time(var.est.4<-var.test(X.s,4,mu.est))
#system.time(var.est.5<-var.smooth(X.s,5))
#system.time(var.est.6<-var.smooth(X.s,6))


plot(sigma.t^2,type='l',ylim=c(4,9))
lines(var.est.1,col=2)
lines(var.est.2,col=2)
lines(var.est.3,col=4)
lines(var.est.4,col=3)
#lines(var.est.5,col=4)
#lines(var.est.6,col=3)

mse(var.est.1,sigma.t^2)
mse(var.est.2,sigma.t^2)
mse(var.est.3,sigma.t^2)
mse(var.est.4,sigma.t^2)
#mse(var.est.5,sigma.t^2)
#mse(var.est.6,sigma.t^2)


var.est=apply(X.s,1,var.smooth,2)
mean(apply((t(var.est)-rep(1,500)%o%sigma.t^2)^2,1,mean))
plot(sigma.t^2,type='l',ylim=c(0,4))
lines(apply(t(var.est),2,mean),col=2)



5:0.05935739



plot(((X.s-lshift(X.s))^2/2))
plot(mu.s0*(mu.s-1))
plot(X.s)