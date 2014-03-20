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
source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/bayesmooth_test.R"))
source("simulation_1d_g/threshold_var.R")
source("simulation_1d_g/wd_var.R")
library(wavethresh)
library(EbayesThresh)
library(caTools)


la10.sure0=function (x, min.level = 0) 
{
    x.w <- wd(x)
    x.w.t <- threshold.wd.mod(x.w, levels = (min.level):(x.w$nlevels - 1))
    x.w.t.r <- wr(x.w.t)
    return(x.w.t.r)
}

la10.sure3=function (x) 
{
    x.w <- wd(x)
    x.w.t <- threshold.wd.mod(x.w)
    x.w.t.r <- wr(x.w.t)
    return(x.w.t.r)
}


la10.ti0u=function (x, filter.number = 10, family = "DaubLeAsymm", min.level = 0, noise.level) 
{
    TT <- length(x)
    thresh <- noise.level * sqrt(2 * log(TT))
    x.w <- wd(x, filter.number, family, type = "station")
    x.w.t <- threshold(x.w, levels = (min.level):(x.w$nlevels - 1), policy = "manual", value = thresh, type = "hard")
    x.w.t.r <- AvBasis(convert(x.w.t))
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


la10.ti0=function (x, filter.number = 10, family = "DaubLeAsymm", min.level = 0, noise.level) 
{
    TT <- length(x)
    thresh <- noise.level * sqrt(2 * log(TT * log2(TT)))    
    x.w <- wd(x, filter.number, family, type = "station")
    x.w.t <- threshold(x.w, levels = (min.level):(x.w$nlevels - 1), policy = "manual", value = thresh, type = "hard")
    x.w.t.r <- AvBasis(convert(x.w.t))
    return(x.w.t.r)
}


la10.ti3=function (x, filter.number = 10, family = "DaubLeAsymm", min.level = 3, noise.level) 
{
    TT <- length(x)
    thresh <- noise.level * sqrt(2 * log(TT * log2(TT)))    
    x.w <- wd(x, filter.number, family, type = "station")
    x.w.t <- threshold(x.w, levels = (min.level):(x.w$nlevels - 1), policy = "manual", value = thresh, type = "hard")
    x.w.t.r <- AvBasis(convert(x.w.t))
    return(x.w.t.r)
}

la10.ti0.mod=function (x, filter.number = 10, family = "DaubLeAsymm", min.level = 0) 
{
    x.w <- wd(x, filter.number, family, type = "station")
    x.w.t <- threshold.wd.mod(x.w, levels = (min.level):(x.w$nlevels - 1), policy = "universal", type = "hard")
    x.w.t.r <- AvBasis(convert(x.w.t))
    return(x.w.t.r)
}

la10.ti3.mod=function (x, filter.number = 10, family = "DaubLeAsymm", min.level = 3) 
{
    x.w <- wd(x, filter.number, family, type = "station")
    x.w.t <- threshold.wd.mod(x.w, levels = (min.level):(x.w$nlevels - 1), policy = "universal", type = "hard")
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

var2=(1+4*(exp(-550*(t-0.2)^2)+exp(-200*(t-0.5)^2)+exp(-950*(t-0.8)^2)))/1.35
var4=3.4*(2+mu.dop)
sigma.ini=sqrt(var2)
sigma.t=sigma.ini/mean(sigma.ini)*sd(mu.t)/rsnr^2

#sigma.t=sd(mu.t)/rsnr^2

set.seed(1025)
X.s=matrix(rnorm(100*n,mu.t,sigma.t),nrow=100,byrow=TRUE)
#write(t(X.s),"simulation_1d_g/sim.txt",ncol=n)

mu.est.asht.haar=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8=matrix(0,nrow=100,ncol=n)
mu.est.ashte.haar=matrix(0,nrow=100,ncol=n)
mu.est.ashte.s8=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar=matrix(0,nrow=100,ncol=n)
mu.est.ti0=matrix(0,nrow=100,ncol=n)
mu.est.ti3=matrix(0,nrow=100,ncol=n)
mu.est.ti0u=matrix(0,nrow=100,ncol=n)
mu.est.ti3u=matrix(0,nrow=100,ncol=n)
mu.est.sure0=matrix(0,nrow=100,ncol=n)
mu.est.sure3=matrix(0,nrow=100,ncol=n)
mu.est.ti0.mod=matrix(0,nrow=100,ncol=n)
mu.est.ti3.mod=matrix(0,nrow=100,ncol=n)
mu.est.ti0.mod=matrix(0,nrow=100,ncol=n)
mu.est.ti3.mod=matrix(0,nrow=100,ncol=n)
mu.est.ebcmn=matrix(0,nrow=100,ncol=n)
mu.est.ebcmd=matrix(0,nrow=100,ncol=n)
mu.est.eblmn=matrix(0,nrow=100,ncol=n)
mu.est.eblmd=matrix(0,nrow=100,ncol=n)


for(i in 1:100){
  sigma=sqrt(2/(3*(n-2))*sum((1/2*X.s[i,1:(n-2)]-X.s[i,2:(n-1)]+1/2*X.s[i,3:n])^2))
  #mu.est.ash[i,]=bayesmooth(X.s[i,],sigma)
  mu.est.asht.haar[i,]=bayesmooth(X.s[i,],sigma=sigma,basis="haar",gridmult=0)
  mu.est.asht.s8[i,]=bayesmooth(X.s[i,],sigma=sigma,basis="symm8",gridmult=0)
  #mu.est.ashe.haar[i,]=bayesmooth(X.s[i,],basis="haar",gridmult=0)
  #mu.est.ashe.s8[i,]=bayesmooth(X.s[i,],basis="symm8",gridmult=0)
  #mu.est.ashte.haar[i,]=bayesmooth.test(X.s[i,],basis="haar",gridmult=0)
  #mu.est.ashte.s8[i,]=bayesmooth.test(X.s[i,],basis="symm8",gridmult=0)
  #mu.est.ashe[i,]=bayesmooth(X.s[i,],basis="haar")
  #mu.est.tie.s8[i,]=waveti.var(X.s[i,])
  #mu.est.tie.haar[i,]=waveti.var(X.s[i,],filter.number=1,family="DaubExPhase")
  #mu.est.ti0[i,]=la10.ti0(X.s[i,],noise.level=sigma)
  #mu.est.ti3[i,]=la10.ti3(X.s[i,],noise.level=sigma)
  #mu.est.sure0[i,]=la10.sure0(X.s[i,])
  #mu.est.sure3[i,]=la10.sure3(X.s[i,])
  #mu.est.ti0u[i,]=la10.ti0u(X.s[i,],noise.level=sigma)
  #mu.est.ti3u[i,]=la10.ti3u(X.s[i,],noise.level=sigma)
  #mu.est.ti0.mod[i,]=la10.ti0.mod(X.s[i,])
  #mu.est.ti3.mod[i,]=la10.ti3.mod(X.s[i,])
  #mu.est.ebcmn[i,]=ebayesthresh(X.s[i,],prior='cauchy',sdev=sigma.t,threshrule='mean',verbose=FALSE)
  #mu.est.ebcmd[i,]=ebayesthresh(X.s[i,],prior='cauchy',sdev=sigma.t,threshrule='median')
  #mu.est.eblmn[i,]=ebayesthresh(X.s[i,],prior='laplace',a=NA,sdev=sigma.t,threshrule='mean')
  #mu.est.eblmd[i,]=ebayesthresh(X.s[i,],prior='laplace',a=NA,sdev=sigma.t,threshrule='median')
  print(i)
}


mu.est.m=data.matrix(read.csv('simulation_1d_g/est.csv',header=F))
#(mise(mu.est.m,mu.t))
(mise(mu.est.asht.haar,mu.t))
(mise(mu.est.asht.s8,mu.t))
(mise(mu.est.ashe.haar,mu.t))
(mise(mu.est.ashe.s8,mu.t))
(mise(mu.est.ashte.haar,mu.t))
(mise(mu.est.ashte.s8,mu.t))

plot(mu.t,type='l')
lines(mu.est.asht.s8[1,],col=4)
lines(mu.est.ashe.s8[1,],col=2)

(mise(mu.est.ash.nh,mu.t))
(mise(mu.est.tie.s8,mu.t))
(mise(mu.est.tie.haar,mu.t))
(mise(mu.est.ti0u,mu.t))
(mise(mu.est.ti3u,mu.t))


(mise(mu.est.ti0,mu.t))
(mise(mu.est.ti3,mu.t))
(mise(mu.est.ti0.mod,mu.t))
(mise(mu.est.ti3.mod,mu.t))
(mise(mu.est.sure0,mu.t))
(mise(mu.est.sure3,mu.t))

(mise(mu.est.ebcmn,mu.t))
(mise(mu.est.ebcmd,mu.t))
(mise(mu.est.eblmn,mu.t))
(mise(mu.est.eblmd,mu.t))


sqrt(mean(apply((mu.est.m-rep(1,100)%o%mu.t)^2,1,mean)))
sqrt(mean(apply((mu.est.ash-rep(1,100)%o%mu.t)^2,1,mean)))
sqrt(mean(apply((mu.est.ashe-rep(1,100)%o%mu.t)^2,1,mean)))
sqrt(mean(apply((mu.est.ti0u-rep(1,100)%o%mu.t)^2,1,mean)))
sqrt(mean(apply((mu.est.ti3u-rep(1,100)%o%mu.t)^2,1,mean)))
sqrt(mean(apply((mu.est.ti0-rep(1,100)%o%mu.t)^2,1,mean)))
sqrt(mean(apply((mu.est.ti3-rep(1,100)%o%mu.t)^2,1,mean)))
sqrt(mean(apply((mu.est.ti0.mod-rep(1,100)%o%mu.t)^2,1,mean)))
sqrt(mean(apply((mu.est.ti3.mod-rep(1,100)%o%mu.t)^2,1,mean)))
sqrt(mean(apply((mu.est.sure0-rep(1,100)%o%mu.t)^2,1,mean)))
sqrt(mean(apply((mu.est.sure3-rep(1,100)%o%mu.t)^2,1,mean)))
sqrt(mean(apply((X.s-rep(1,100)%o%mu.t)^2,1,mean)))


sqrt(mse(apply(mu.est.m,2,mean),mu.t))
sqrt(mse(apply(mu.est.ash,2,mean),mu.t))
sqrt(mse(apply(mu.est.ashe,2,mean),mu.t))
sqrt(mse(apply(mu.est.ti0u,2,mean),mu.t))
sqrt(mse(apply(mu.est.ti3u,2,mean),mu.t))
sqrt(mse(apply(mu.est.ti0,2,mean),mu.t))
sqrt(mse(apply(mu.est.ti3,2,mean),mu.t))
sqrt(mse(apply(mu.est.ti0.mod,2,mean),mu.t))
sqrt(mse(apply(mu.est.ti3.mod,2,mean),mu.t))
sqrt(mse(apply(mu.est.sure0,2,mean),mu.t))
sqrt(mse(apply(mu.est.sure3,2,mean),mu.t))


plot(mu.t,type='l')
lines(mu.est.ash.nh[1,],col=4)
lines(mu.est.ashe[1,],col=2)
lines(mu.est.m[1,],col=4)
lines(mu.est.tie.s8[1,],col=2)
lines(mu.est.tie.haar[1,],col=2)


X.test=rnorm(n,100*mu.t,1)
ebest=ebayesthresh(X.test,sdev=1,bayesfac=FALSE,threshrule='soft')
ashest=bayesmooth(X.test)
plot(ebest,type='l')
lines(ashest,col=2)
mse(ebest,100*mu.t)
mse(ashest,100*mu.t)
mse(X.test,100*mu.t)


set.seed(1)
x <- rnorm(1000) + sample(c( runif(25,-7,7), rep(0,975)))
ebest=ebayesthresh(x, sdev=1)
plot(x,type='l')
lines(ebest,col=2)


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
sigma.t=sqrt(var2)
#sigma.t=seq(2,10,length.out=n)
#sigma.t=rep(2,n)
set.seed(1025)
X.s=rnorm(n,0,sigma.t)
#X.s=matrix(rnorm(500*n,mu.sine,sigma.t),nrow=500,byrow=TRUE)

source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/bayesmooth.R"))
source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/bayesmooth_test_var.R"))


var.est=bayesmooth(X.s,v.est=TRUE)
var.est.test=bayesmooth.var(X.s,v.est=TRUE)

mse(var.est,sigma.t^2)
mse(var.est.test,sigma.t^2)

plot(sigma.t^2,type='l')
lines(var.est,col=2)
lines(var.est.test,col=4)


var.est=matrix(0,500,n)
for(i in 1:500){
  var.est[i,]=bayesmooth(X.s[i,],v.est=TRUE)
  print(i)
}
mean(apply((var.est-rep(1,500)%o%sigma.t^2)^2,1,mean))
plot(sigma.t^2,type='l',ylim=c(0,7))
lines(apply(var.est,2,mean),col=2)

