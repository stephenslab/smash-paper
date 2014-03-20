source("D:/Grad School/Spring 2013/multiscale_ash/simulation_1d_g/wip/s8_comp.R")
source("D:/Grad School/Spring 2013/multiscale_ash/ash/bayesmooth.R")


spike.f=function(x) (0.75*exp(-500*(x-0.23)^2)+1.5*exp(-2000*(x-0.33)^2)+3*exp(-8000*(x-0.47)^2)+2.25*exp(-16000*(x-0.69)^2)+0.5*exp(-32000*(x-0.83)^2))
n=1024
t=1:n/n

mu.s=spike.f(t)

pos = c(.1, .13, .15, .23, .25, .40, .44, .65, .76, .78, .81)
hgt = 2.97/5*c(4, 5, 3, 4, 5, 4.2, 2.1, 4.3, 3.1, 5.1, 4.2)
wth = c(.005, .005, .006, .01, .01, .03, .01, .01, .005, .008, .005)
mu.b = rep(0,n)
for(j in 1:length(pos)){
  mu.b = mu.b + hgt[j]/(( 1 + (abs(t - pos[j])/wth[j]))^4)
}


rsnr=sqrt(3)
mu.t=(1+mu.b)/5
var2=(0.01+4*(exp(-550*(t-0.2)^2)+exp(-200*(t-0.5)^2)+exp(-950*(t-0.8)^2)))/1.35
sigma.ini=sqrt(var2)
sigma.t=sigma.ini/mean(sigma.ini)*sd(mu.t)/rsnr^2

x=rnorm(n,mu.t,sigma.t)

mse=function(x,y) mean((x-y)^2)
l2norm=function(x) sum(x^2)
mise=function(x,y) 10000*mean(apply(x-rep(1,100)%o%y,1,l2norm)/l2norm(y))


system.time(mu.true.s8<-bayesmooth(x,sigma.t,basis="symm8",gridmult=0))
system.time(mu.true.haar<-bayesmooth(x,sigma.t,gridmult=0))

system.time(mu.est.s8<-bayesmooth(x,basis="symm8",gridmult=0))
system.time(mu.est.haar<-bayesmooth(x,gridmult=0))


mse(mu.true.s8,mu.t)
mse(mu.true.haar,mu.t)
mse(mu.est.s8,mu.t)
mse(mu.est.haar,mu.t)

par(mfrow=c(1,1))
plot(mu.t,type='l')
lines(mu.est.s8,col=2)
lines(mu.est.haar,col=4)

plot(mu.t,type='l')
lines(mu.true.haar,col=4)
lines(mu.true.s8,col=2)




