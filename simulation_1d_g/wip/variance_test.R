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
X.s=rnorm(n,0,sigma.t)
set.seed(1121)
X.s=matrix(rnorm(100*n,5*mu.sine,sigma.t),nrow=100,byrow=TRUE)


source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/bayesmooth.R"))
source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/bayesmooth_alt.R"))
source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/bayesmooth_alt1.R"))
source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/bayesmooth_alt2.R"))

###single
system.time(var.est<-bayesmooth(X.s,v.est=TRUE))
system.time(var.est.alt<-bayesmooth.alt(X.s,v.est=TRUE))
system.time(var.est.alt1<-bayesmooth.alt1(X.s,v.est=TRUE))
system.time(var.est.alt2<-bayesmooth.alt2(X.s,v.est=TRUE))

system.time(mu.est<-bayesmooth(X.s))
system.time(mu.est.alt<-bayesmooth.alt(X.s))
system.time(mu.est.alt1<-bayesmooth.alt1(X.s))
system.time(mu.est.alt2<-bayesmooth.alt2(X.s))


mse(var.est,sigma.t^2)
mse(var.est.alt,sigma.t^2)
mse(var.est.alt1,sigma.t^2)
mse(var.est.alt2,sigma.t^2)

plot(sigma.t^2,type='l')
lines(var.est,col=2)
lines(var.est.alt,col=3)
lines(var.est.alt2,col=4)


mse(mu.est,0)
mse(mu.est.alt,0)
mse(mu.est.alt1,0)
mse(mu.est.alt2,0)

plot(rep(0,n),type='l')
lines(mu.est,col=2)
lines(mu.est.alt,col=3)
lines(mu.est.alt2,col=4)


system.time(var.est.2<-var.smooth(X.s,2))
system.time(var.est.3<-var.smooth(X.s,3))

system.time(var.est.4<-var.smooth(X.s,4,mu.est))
system.time(var.est.4.old<-var.smooth_old(X.s,4,mu.est))

#system.time(var.est.5<-var.smooth(X.s,5))
#system.time(var.est.6<-var.smooth(X.s,6))


plot(sigma.t^2,type='l',ylim=c(0,7))
lines(var.est,col=2)
lines(var.est.2,col=2)
lines(var.est.3,col=4)
lines(var.est.4,col=3)
#lines(var.est.5,col=4)
#lines(var.est.6,col=3)

mse(var.est.2,sigma.t^2)
mse(var.est.3,sigma.t^2)
mse(var.est.4,sigma.t^2)
#mse(var.est.5,sigma.t^2)
#mse(var.est.6,sigma.t^2)

###multiple
mvar.est=apply(X.s,1,bayesmooth,v.est=TRUE)
mvar.est.alt=apply(X.s,1,bayesmooth.alt,v.est=TRUE)
mvar.est.alt1=apply(X.s,1,bayesmooth.alt1,v.est=TRUE)
mvar.est.alt2=apply(X.s,1,bayesmooth.alt2,v.est=TRUE)

mean(apply((t(mvar.est)-rep(1,100)%o%sigma.t^2)^2,1,mean))
mean(apply((t(mvar.est.alt)-rep(1,100)%o%sigma.t^2)^2,1,mean))
mean(apply((t(mvar.est.alt1)-rep(1,100)%o%sigma.t^2)^2,1,mean))
mean(apply((t(mvar.est.alt2)-rep(1,100)%o%sigma.t^2)^2,1,mean))


plot(sigma.t^2,type='l',ylim=c(0,20))
lines(apply(t(var.est),2,mean),col=2)



#2:0.0703
#5:0.05935739

#ash2, ash, 0.0553797
#ash, ash2  0.0516653
#ash, ash   0.05578719
#ash2, ash2 0.05120646

var.est=matrix(0,500,n)
for(i in 1:500){
  #source("D:/Grad School/Spring 2013/ash_test/Rcode/ash.R")
  #source("D:/Grad School/Spring 2013/ash_test/Rcode/mix.R")
  #require(ashr)
  mu.est=cs.smooth(X.s[i,],1)
  #source("D:/Grad School/Spring 2013/ash_test/Rcode/ash.R")
  #source("D:/Grad School/Spring 2013/ash_test/Rcode/mix.R")
  #require(ashr)  
  var.est[i,]=var.smooth(X.s[i,],4,mu.est)
  print(i)
}
mean(apply((var.est-rep(1,500)%o%sigma.t^2)^2,1,mean))
plot(sigma.t^2,type='l',ylim=c(0,7))
lines(apply(var.est,2,mean),col=2)


plot(((X.s-lshift(X.s))^2/2))
plot(mu.s0*(mu.s-1))
plot(X.s)