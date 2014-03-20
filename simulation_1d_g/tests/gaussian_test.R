spike.f=function(x) (0.75*exp(-500*(x-0.23)^2)+1.5*exp(-2000*(x-0.33)^2)+3*exp(-8000*(x-0.47)^2)+2.25*exp(-16000*(x-0.69)^2)+0.5*exp(-32000*(x-0.83)^2))
n=1024
t=1:n/n
mu.s=spike.f(1:n/n)
#sigma.t=seq(2,10,length.out=n)
sigma.t=rep(0.17,n)

mu.s0=2
set.seed(425)
X.s=rnorm(n,mu.s0*(mu.s),sigma.t)

pos = c(.1, .13, .15, .23, .25, .40, .44, .65, .76, .78, .81)
hgt = 2.97/5*c(4, 5, 3, 4, 5, 4.2, 2.1, 4.3, 3.1, 5.1, 4.2)
wth = c(.005, .005, .006, .01, .01, .03, .01, .01, .005, .008, .005)
mu.b = rep(0,n)
for(j in 1:length(pos)){
  mu.b = mu.b + hgt[j]/(( 1 + (abs(t - pos[j])/wth[j]))^4)
}

mu.b0=2
set.seed(415)
#X.b=rnorm(n,mu.b0*(mu.b),sigma.t)
X.b=rnorm(n,(mu.b/max(mu.b)),sigma.t)


dop.f=function(x) sqrt(x*(1-x))*sin((2*pi*1.05)/(x+0.05))
mu.dop=dop.f(t)
mu.dop=3/(max(mu.dop)-min(mu.dop))*(mu.dop-min(mu.dop))

sig=((2*t + 0.5)*(t <= 0.15)) + 
       ((-12*(t-0.15) + 0.8)*(t > 0.15 & t <= 0.2)) + 
       0.2*(t > 0.2 & t <= 0.5) + 
       ((6*(t - 0.5) + 0.2)*(t > 0.5 & t <= 0.6)) + 
       ((-10*(t - 0.6) + 0.8)*(t > 0.6 & t <= 0.65)) + 
       ((-0.5*(t - 0.65) + 0.3)*(t > 0.65 & t <= 0.85)) + 
       ((2*(t - 0.85) + 0.2)*(t > 0.85))
mu.ang = 3/5*((5/(max(sig)-min(sig)))*sig - 1.6)-0.0419569



heavi = 4 * sin(4 * pi * t) - sign(t - 0.3) - sign(0.72 - t)
mu.hs = heavi/sqrt(var(heavi)) * 1 * 2.99 / 3.366185
mu.hs = mu.hs - min(mu.hs)



pos=c(.1, .13, .15, .23, .25, .40, .44, .65, .76, .78, .81)
hgt = 2.88/5*c(4, (-5), 3, (-4), 5, (-4.2), 2.1, 4.3, (-3.1), 2.1, (-4.2))
mu.blk = rep(0,n)
for(j in 1:length(pos)){
  mu.blk = mu.blk + (1 + sign(t-pos[j]))*(hgt[j]/2)
}


mse=function(x,y) mean((x-y)^2)
l2norm=function(x) sum(x^2)
mise=function(x,y) 10000*mean(apply(x-rep(1,100)%o%y,1,l2norm)/l2norm(y))




source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/bayesmooth_alt.R"))
source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/bayesmooth.R"))
source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/bayesmooth_test.R"))





n=1024
t=1:n/n
mu.t=(1+mu.s)/5
#mu.t=(1+mu.b)/5
#mu.t=(1+mu.blk)/5
#mu.t=(1+mu.ang)/5
#mu.t=(1+mu.dop)/5

#mu.t=mu.blip
#mu.t=mu.cor

rsnr=sqrt(1)

var1=rep(1,n)
var2=(0.0001+4*(exp(-550*(t-0.2)^2)+exp(-200*(t-0.5)^2)+exp(-950*(t-0.8)^2)))/1.35
var4=3.4*(2+mu.dop)
sigma.ini=sqrt(var2)
sigma.t=sigma.ini/mean(sigma.ini)*sd(mu.t)/rsnr^2

#sigma.t=sd(mu.t)/rsnr^2

set.seed(1025)
X.s=rnorm(n,mu.t,sigma.t)
set.seed(1107)
X.s=matrix(rnorm(100*n,mu.t,sigma.t),nrow=100,byrow=TRUE)

###single
system.time(mu.est<-bayesmooth(X.s))
system.time(mu.est.sd<-bayesmooth(X.s/sigma.t,sigma=rep(1,n)))
system.time(mu.est.alt<-bayesmooth.alt(X.s))
system.time(mu.est.test<-bayesmooth.test(X.s,sigma=sigma.t))

mse(mu.est,mu.t)
mse(mu.est.sd*sigma.t,mu.t)
mse(mu.est.alt,mu.t)
mse(mu.est.test,mu.t)

###multiple
mu.est=apply(X.s,1,bayesmooth)
mu.est.alt=apply(X.s,1,bayesmooth.alt)
mu.est.test=apply(X.s,1,bayesmooth.test,sigma=sigma.t)

mise(t(mu.est),mu.t)
mise(t(mu.est.alt),mu.t)
mise(t(mu.est.test),mu.t)

