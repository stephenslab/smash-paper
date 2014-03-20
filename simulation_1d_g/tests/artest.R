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


mu.blip=(0.32+0.6*t+0.3*exp(-100*(t-0.3)^2))*(t>=0&t<=0.8)+(-0.28+0.6*t+0.3*exp(-100*(t-1.3)^2))*(t>0.8&t<=1)

para=function(x,c) (x-c)^2*(x>c&x<=1)
mu.para=0.8-30*para(t,0.1)+60*para(t,0.2)-30*para(t,0.3)+500*para(t,0.35)-1000*para(t,0.37)+1000*para(t,0.41)-500*para(t,0.43)+7.5*para(t,0.5)-15*para(t,0.7)+7.5*para(t,0.9)

tssine=function(x) (1-cos(pi*x))/2
mu.tss=0.3*sin(3*pi*(tssine(tssine(tssine(tssine(t))))+t))+0.5


mse=function(x,y) mean((x-y)^2)
l2norm=function(x) sum(x^2)
mise=function(x,y) 10000*mean(apply(x-rep(1,100)%o%y,1,l2norm)/l2norm(y))

mu.cor=623.87*t^3*(1-2*t)*(t>=0&t<=0.5)+187.161*(0.125-t^3)*t^4*(t>0.5&t<=0.8)+3708.470441*(t-1)^3*(t>0.8&t<=1)
mu.cor=(0.6/(max(mu.cor)-min(mu.cor)))*mu.cor
mu.cor=mu.cor-min(mu.cor)+0.2


mu.wave=0.5+0.2*cos(4*pi*t)+0.1*cos(24*pi*t)



source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/bayesmooth.R"))
source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/bayesmooth_mc.R"))
source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/bayesmooth_test.R"))
source("D:/Grad School/Spring 2013/multiscale_ash/simulation_1d_g/threshold_var.R")
source("D:/Grad School/Spring 2013/multiscale_ash/simulation_1d_g/wd_var.R")
library(wavethresh)
library(EbayesThresh)
library(caTools)



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


waveti=function (x, filter.number = 8, family = "DaubLeAsymm", min.level = 3) 
{
    x.w <- wd(x, filter.number, family, type = "station")
    x.w.t <- threshold(x.w, by.level=TRUE, levels = (min.level):(x.w$nlevels - 1), policy = "universal", type = "hard")
    x.w.t.r <- AvBasis(convert(x.w.t))
    return(x.w.t.r)
}

s.est=function(x) {
  n=length(x)
  if(n<3){
    sqrt(1/(2*(n-1))*sum((x[1:(n-1)]-x[2:n])^2))
  }else{
    sqrt(2/(3*(n-2))*sum((1/2*x[1:(n-2)]-x[2:(n-1)]+1/2*x[3:n])^2))
  }
}

n=1024
t=1:n/n
mu.t=(1+mu.s)/5
#mu.t=(1+mu.b)/5
#mu.t=(1+mu.blk)/5
#mu.t=(1+mu.ang)/5
#mu.t=(1+mu.dop)/5

#mu.t=mu.blip
#mu.t=mu.cor

rsnr=sqrt(3)

#var1=rep(1,n)
var2=(0.1+4*(exp(-550*(t-0.2)^2)+exp(-200*(t-0.5)^2)+exp(-950*(t-0.8)^2)))/1.35
var4=3.4*(2+mu.dop)
sigma.ini=sqrt(var2)
sigma.t=sigma.ini/mean(sigma.ini)*sd(mu.t)/rsnr^2

sigma.t=sd(mu.t)/rsnr^2

set.seed(1025)
X1=rnorm(n,mu.t,sigma.t)
set.seed(1025)
X2=arima.sim(list(order = c(1,0,0),ar=0.9),n=1024,sd=sigma.t)
X2=mu.t+X2

plot(X1,type='l')
lines(X2,col=2)


plot(arima.sim(list(order = c(1,0,0),ar=0.9),n=1024,sd=sigma.t))
plot(rnorm(n,0,sigma.t),type='l')

mu.est1=bayesmooth(X1,sigma.t)
mu.est1.test=bayesmooth.test(X1,sigma.t)
mu.est2=bayesmooth(as.numeric(X2),sigma.t)
mu.est2.test=bayesmooth.test(as.numeric(X2),sigma.t)

mu.est2.ti=waveti(as.numeric(X2))

mse(mu.est1,mu.t)
mse(mu.est1.test,mu.t)

mse(mu.est2,mu.t)
mse(mu.est2.test,mu.t)
mse(mu.est2.ti,mu.t)

plot(mu.t,type='l')
lines(mu.est2,col=2)
lines(mu.est2.test,col=4)





