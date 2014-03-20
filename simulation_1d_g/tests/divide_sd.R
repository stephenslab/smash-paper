source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/bayesmooth.R"))

mse=function(x,y) mean((x-y)^2)
l2norm=function(x) sum(x^2)
mise=function(x,y) 10000*mean(apply(x-rep(1,100)%o%y,1,l2norm)/l2norm(y))

spike.f=function(x) (0.75*exp(-500*(x-0.23)^2)+1.5*exp(-2000*(x-0.33)^2)+3*exp(-8000*(x-0.47)^2)+2.25*exp(-16000*(x-0.69)^2)+0.5*exp(-32000*(x-0.83)^2))

n=1024
t=1:n/n

mu.s=spike.f(t)

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
sigma.ini=sqrt(var2)
sigma.t=sigma.ini/mean(sigma.ini)*sd(mu.t)/rsnr^2


###single
set.seed(1025)
X.s=rnorm(n,mu.t,sigma.t)

mu.est<-bayesmooth(X.s,sigma=sigma.t)
mu.est.sig=sigma.t*bayesmooth(X.s/sigma.t,sigma=rep(1,n))

mse(mu.est,mu.t)
mse(mu.est.sig,mu.t)

plot(mu.t,type='l')
lines(mu.est,col=2)
lines(mu.est.sig,col=4)

###multiple
set.seed(1107)
X.s=matrix(rnorm(100*n,mu.t,sigma.t),nrow=100,byrow=TRUE)

mu.est=apply(X.s,1,bayesmooth,sigma=sigma.t)
mu.est.sig=sigma.t%o%rep(1,100)*apply(X.s/(rep(1,100)%o%sigma.t),1,bayesmooth,sigma=rep(1,n))

mise(t(mu.est),mu.t)
mise(t(mu.est.sig),mu.t)
