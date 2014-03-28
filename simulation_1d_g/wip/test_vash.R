spike.f=function(x) (0.75*exp(-500*(x-0.23)^2)+1.5*exp(-2000*(x-0.33)^2)+3*exp(-8000*(x-0.47)^2)+2.25*exp(-16000*(x-0.69)^2)+0.5*exp(-32000*(x-0.83)^2))
n=2^12
t=1:n/n
mu.s=spike.f(t)
pos = c(.1, .13, .15, .23, .25, .40, .44, .65, .76, .78, .81)
hgt = c(2, 3.6, 2, 5, 7.5, 6.9, 2, 4.8, 2, 4.1, 2.3)
wth = c(.005, .005, .006, .01, .01, .03, .01, .01, .005, .008, .005)
mu.b = rep(0,n)
for(j in 1:length(pos)){
  mu.b = mu.b + hgt[j]/(( 1 + (abs(t - pos[j])/wth[j]))^4)
}
dop.f=function(x) sqrt(x*(1-x))*sin((2*pi*1.05)/(x+0.05))
mu.dop=dop.f(t)
#mu.sine=3/4*sin(40*pi*t)
mu.sine=sin(20*t)
var1=1.42*((3-20*t)*(t>=0&t<0.1)+(20*t-1)*(t>=0.1&t<0.25)+(4+(1-4*t)*18/19)*(t>=0.25&t<0.725)+(2.2+10*(t-0.725))*(t>=0.725&t<0.89)+(3.85-85*(t-0.89)/11)*(t>=0.89&t<=1))
var2=(1+4*(exp(-550*(t-0.2)^2)+exp(-200*(t-0.5)^2)+exp(-950*(t-0.8)^2)))/1.35
var4=3.4*(2+mu.dop)
var3=mu.b
var5=(t-0.5)^2+0.5
var6=exp(-100*(t-0.5)^2)
var7=0.5+1.5*mu.s
sigma.t=sqrt(var3)
#sigma.t=seq(2,10,length.out=n)
#sigma.t=rep(2,n)
l2norm=function(x) sum(x^2)
mise=function(x,y) l2norm(x-y)/l2norm(y)
mse=function(x,y) mean((x-y)^2)
mu.s0=100
set.seed(1121)
mu.t=rep(0,n)
X.s=rnorm(n,0,sigma.t)

source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/bayesmooth.R"))
source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/bayesmooth_var.R"))
source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/bayesmooth_vash.R"))
source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/bayesmooth_var_vash.R"))

system.time(var.est<-bayesmooth(X.s,v.est=TRUE))
system.time(var.est.v<-bayesmooth.vash(X.s,v.est=TRUE))
system.time(var.est.t<-bayesmooth.var(X.s,mu.t))
system.time(var.est.v.t<-bayesmooth.var.vash(X.s,mu.t,gridmult=2))


mse(var.est,sigma.t^2)
mse(var.est.v,sigma.t^2)
mse(var.est.t,sigma.t^2)
mse(var.est.v.t,sigma.t^2)

plot(sigma.t^2,type='l')
lines(var.est,col=2)
lines(var.est.v,col=4)