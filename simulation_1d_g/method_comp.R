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


dop.f=function(x) sqrt(x*(1-x))*sin((2*pi*1.05)/(x+0.05))
mu.dop=dop.f(t)
mu.dop=3/(max(mu.dop)-min(mu.dop))*(mu.dop-min(mu.dop))
mu.dop.var=10*dop.f(t)
mu.dop.var=mu.dop.var-min(mu.dop.var)



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

mu.cblk=mu.blk
mu.cblk[mu.cblk<0]=0

mu.blip=(0.32+0.6*t+0.3*exp(-100*(t-0.3)^2))*(t>=0&t<=0.8)+(-0.28+0.6*t+0.3*exp(-100*(t-1.3)^2))*(t>0.8&t<=1)

para=function(x,c) (x-c)^2*(x>c&x<=1)
mu.para=0.8-30*para(t,0.1)+60*para(t,0.2)-30*para(t,0.3)+500*para(t,0.35)-1000*para(t,0.37)+1000*para(t,0.41)-500*para(t,0.43)+7.5*para(t,0.5)-15*para(t,0.7)+7.5*para(t,0.9)

tssine=function(x) (1-cos(pi*x))/2
mu.tss=0.3*sin(3*pi*(tssine(tssine(tssine(tssine(t))))+t))+0.5



mu.cor=623.87*t^3*(1-2*t)*(t>=0&t<=0.5)+187.161*(0.125-t^3)*t^4*(t>0.5&t<=0.8)+3708.470441*(t-1)^3*(t>0.8&t<=1)
mu.cor=(0.6/(max(mu.cor)-min(mu.cor)))*mu.cor
mu.cor=mu.cor-min(mu.cor)+0.2


mu.wave=0.5+0.2*cos(4*pi*t)+0.1*cos(24*pi*t)


sig.est.func=function(x,n) sqrt(2/(3*(n-2))*sum((1/2*x[1:(n-2)]-x[2:(n-1)]+1/2*x[3:n])^2))
mse=function(x,y) mean((x-y)^2)
l2norm=function(x) sum(x^2)
mise=function(x,y) 10000*mean(apply(x-rep(1,100)%o%y,1,l2norm)/l2norm(y))



source(file.path("~/ashwave/Rcode/bayesmooth.R"))
source(file.path("~/ashwave/Rcode/ti_thresh.R"))
library(wavethresh)
library(EbayesThresh)


waveti.ebayes=function (x, filter.number = 10, family = "DaubLeAsymm", min.level = 3, noise.level){
    n=length(x)
    J=log2(n)
    x.w <- wd(x, filter.number, family, type = "station")
    for(j in min.level:(J-1)){
      x.pm = ebayesthresh(accessD(x.w,j),sdev=noise.level)
      x.w = putD(x.w,j,x.pm)
    }
    mu.est=AvBasis(convert(x.w))
    return(mu.est)
}

##generate data
var1=rep(1,n)
var2=(0.0001+4*(exp(-550*(t-0.2)^2)+exp(-200*(t-0.5)^2)+exp(-950*(t-0.8)^2)))
var3=(0.00001+2*mu.dop.var)
var4=0.00001+mu.b
var5=0.00001+1*(mu.cblk-min(mu.cblk))/max(mu.cblk)
sigma.ini.v1=sqrt(var1)
sigma.ini.v2=sqrt(var2)
sigma.ini.v3=sqrt(var3)
sigma.ini.v4=sqrt(var4)
sigma.ini.v5=sqrt(var5)


mu.sp=(1+mu.s)/5

rsnr=sqrt(1)
sigma.sp.1.v1=sigma.ini.v1/mean(sigma.ini.v1)*sd(mu.sp)/rsnr^2
sigma.sp.1.v2=sigma.ini.v2/mean(sigma.ini.v2)*sd(mu.sp)/rsnr^2
sigma.sp.1.v3=sigma.ini.v3/mean(sigma.ini.v3)*sd(mu.sp)/rsnr^2
sigma.sp.1.v4=sigma.ini.v4/mean(sigma.ini.v4)*sd(mu.sp)/rsnr^2
sigma.sp.1.v5=sigma.ini.v5/mean(sigma.ini.v5)*sd(mu.sp)/rsnr^2

rsnr=sqrt(3)
sigma.sp.3.v1=sigma.ini.v1/mean(sigma.ini.v1)*sd(mu.sp)/rsnr^2
sigma.sp.3.v2=sigma.ini.v2/mean(sigma.ini.v2)*sd(mu.sp)/rsnr^2
sigma.sp.3.v3=sigma.ini.v3/mean(sigma.ini.v3)*sd(mu.sp)/rsnr^2
sigma.sp.3.v4=sigma.ini.v4/mean(sigma.ini.v4)*sd(mu.sp)/rsnr^2
sigma.sp.3.v5=sigma.ini.v5/mean(sigma.ini.v5)*sd(mu.sp)/rsnr^2

mu.bump=(1+mu.b)/5

rsnr=sqrt(1)
sigma.bump.1.v1=sigma.ini.v1/mean(sigma.ini.v1)*sd(mu.bump)/rsnr^2
sigma.bump.1.v2=sigma.ini.v2/mean(sigma.ini.v2)*sd(mu.bump)/rsnr^2
sigma.bump.1.v3=sigma.ini.v3/mean(sigma.ini.v3)*sd(mu.bump)/rsnr^2
sigma.bump.1.v4=sigma.ini.v4/mean(sigma.ini.v4)*sd(mu.bump)/rsnr^2
sigma.bump.1.v5=sigma.ini.v5/mean(sigma.ini.v5)*sd(mu.bump)/rsnr^2

rsnr=sqrt(3)
sigma.bump.3.v1=sigma.ini.v1/mean(sigma.ini.v1)*sd(mu.bump)/rsnr^2
sigma.bump.3.v2=sigma.ini.v2/mean(sigma.ini.v2)*sd(mu.bump)/rsnr^2
sigma.bump.3.v3=sigma.ini.v3/mean(sigma.ini.v3)*sd(mu.bump)/rsnr^2
sigma.bump.3.v4=sigma.ini.v4/mean(sigma.ini.v4)*sd(mu.bump)/rsnr^2
sigma.bump.3.v5=sigma.ini.v5/mean(sigma.ini.v5)*sd(mu.bump)/rsnr^2

mu.blk=0.2+0.6*(mu.blk-min(mu.blk))/max(mu.blk-min(mu.blk))

rsnr=sqrt(1)
sigma.blk.1.v1=sigma.ini.v1/mean(sigma.ini.v1)*sd(mu.blk)/rsnr^2
sigma.blk.1.v2=sigma.ini.v2/mean(sigma.ini.v2)*sd(mu.blk)/rsnr^2
sigma.blk.1.v3=sigma.ini.v3/mean(sigma.ini.v3)*sd(mu.blk)/rsnr^2
sigma.blk.1.v4=sigma.ini.v4/mean(sigma.ini.v4)*sd(mu.blk)/rsnr^2
sigma.blk.1.v5=sigma.ini.v5/mean(sigma.ini.v5)*sd(mu.blk)/rsnr^2

rsnr=sqrt(3)
sigma.blk.3.v1=sigma.ini.v1/mean(sigma.ini.v1)*sd(mu.blk)/rsnr^2
sigma.blk.3.v2=sigma.ini.v2/mean(sigma.ini.v2)*sd(mu.blk)/rsnr^2
sigma.blk.3.v3=sigma.ini.v3/mean(sigma.ini.v3)*sd(mu.blk)/rsnr^2
sigma.blk.3.v4=sigma.ini.v4/mean(sigma.ini.v4)*sd(mu.blk)/rsnr^2
sigma.blk.3.v5=sigma.ini.v5/mean(sigma.ini.v5)*sd(mu.blk)/rsnr^2


mu.ang=(1+mu.ang)/5

rsnr=sqrt(1)
sigma.ang.1.v1=sigma.ini.v1/mean(sigma.ini.v1)*sd(mu.ang)/rsnr^2
sigma.ang.1.v2=sigma.ini.v2/mean(sigma.ini.v2)*sd(mu.ang)/rsnr^2
sigma.ang.1.v3=sigma.ini.v3/mean(sigma.ini.v3)*sd(mu.ang)/rsnr^2
sigma.ang.1.v4=sigma.ini.v4/mean(sigma.ini.v4)*sd(mu.ang)/rsnr^2
sigma.ang.1.v5=sigma.ini.v5/mean(sigma.ini.v5)*sd(mu.ang)/rsnr^2

rsnr=sqrt(3)
sigma.ang.3.v1=sigma.ini.v1/mean(sigma.ini.v1)*sd(mu.ang)/rsnr^2
sigma.ang.3.v2=sigma.ini.v2/mean(sigma.ini.v2)*sd(mu.ang)/rsnr^2
sigma.ang.3.v3=sigma.ini.v3/mean(sigma.ini.v3)*sd(mu.ang)/rsnr^2
sigma.ang.3.v4=sigma.ini.v4/mean(sigma.ini.v4)*sd(mu.ang)/rsnr^2
sigma.ang.3.v5=sigma.ini.v5/mean(sigma.ini.v5)*sd(mu.ang)/rsnr^2


mu.dop=(1+mu.dop)/5

rsnr=sqrt(1)
sigma.dop.1.v1=sigma.ini.v1/mean(sigma.ini.v1)*sd(mu.dop)/rsnr^2
sigma.dop.1.v2=sigma.ini.v2/mean(sigma.ini.v2)*sd(mu.dop)/rsnr^2
sigma.dop.1.v3=sigma.ini.v3/mean(sigma.ini.v3)*sd(mu.dop)/rsnr^2
sigma.dop.1.v4=sigma.ini.v4/mean(sigma.ini.v4)*sd(mu.dop)/rsnr^2
sigma.dop.1.v5=sigma.ini.v5/mean(sigma.ini.v5)*sd(mu.dop)/rsnr^2

rsnr=sqrt(3)
sigma.dop.3.v1=sigma.ini.v1/mean(sigma.ini.v1)*sd(mu.dop)/rsnr^2
sigma.dop.3.v2=sigma.ini.v2/mean(sigma.ini.v2)*sd(mu.dop)/rsnr^2
sigma.dop.3.v3=sigma.ini.v3/mean(sigma.ini.v3)*sd(mu.dop)/rsnr^2
sigma.dop.3.v4=sigma.ini.v4/mean(sigma.ini.v4)*sd(mu.dop)/rsnr^2
sigma.dop.3.v5=sigma.ini.v5/mean(sigma.ini.v5)*sd(mu.dop)/rsnr^2


mu.blip=mu.blip

rsnr=sqrt(1)
sigma.blip.1.v1=sigma.ini.v1/mean(sigma.ini.v1)*sd(mu.blip)/rsnr^2
sigma.blip.1.v2=sigma.ini.v2/mean(sigma.ini.v2)*sd(mu.blip)/rsnr^2
sigma.blip.1.v3=sigma.ini.v3/mean(sigma.ini.v3)*sd(mu.blip)/rsnr^2
sigma.blip.1.v4=sigma.ini.v4/mean(sigma.ini.v4)*sd(mu.blip)/rsnr^2
sigma.blip.1.v5=sigma.ini.v5/mean(sigma.ini.v5)*sd(mu.blip)/rsnr^2

rsnr=sqrt(3)
sigma.blip.3.v1=sigma.ini.v1/mean(sigma.ini.v1)*sd(mu.blip)/rsnr^2
sigma.blip.3.v2=sigma.ini.v2/mean(sigma.ini.v2)*sd(mu.blip)/rsnr^2
sigma.blip.3.v3=sigma.ini.v3/mean(sigma.ini.v3)*sd(mu.blip)/rsnr^2
sigma.blip.3.v4=sigma.ini.v4/mean(sigma.ini.v4)*sd(mu.blip)/rsnr^2
sigma.blip.3.v5=sigma.ini.v5/mean(sigma.ini.v5)*sd(mu.blip)/rsnr^2


mu.cor=mu.cor

rsnr=sqrt(1)
sigma.cor.1.v1=sigma.ini.v1/mean(sigma.ini.v1)*sd(mu.cor)/rsnr^2
sigma.cor.1.v2=sigma.ini.v2/mean(sigma.ini.v2)*sd(mu.cor)/rsnr^2
sigma.cor.1.v3=sigma.ini.v3/mean(sigma.ini.v3)*sd(mu.cor)/rsnr^2
sigma.cor.1.v4=sigma.ini.v4/mean(sigma.ini.v4)*sd(mu.cor)/rsnr^2
sigma.cor.1.v5=sigma.ini.v5/mean(sigma.ini.v5)*sd(mu.cor)/rsnr^2

rsnr=sqrt(3)
sigma.cor.3.v1=sigma.ini.v1/mean(sigma.ini.v1)*sd(mu.cor)/rsnr^2
sigma.cor.3.v2=sigma.ini.v2/mean(sigma.ini.v2)*sd(mu.cor)/rsnr^2
sigma.cor.3.v3=sigma.ini.v3/mean(sigma.ini.v3)*sd(mu.cor)/rsnr^2
sigma.cor.3.v4=sigma.ini.v4/mean(sigma.ini.v4)*sd(mu.cor)/rsnr^2
sigma.cor.3.v5=sigma.ini.v5/mean(sigma.ini.v5)*sd(mu.cor)/rsnr^2



set.seed(1107)
sim.m.sp.1.v1=matrix(rnorm(100*n,mu.sp,sigma.sp.1.v1),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.sp.1.v2=matrix(rnorm(100*n,mu.sp,sigma.sp.1.v2),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.sp.1.v3=matrix(rnorm(100*n,mu.sp,sigma.sp.1.v3),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.sp.1.v4=matrix(rnorm(100*n,mu.sp,sigma.sp.1.v4),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.sp.1.v5=matrix(rnorm(100*n,mu.sp,sigma.sp.1.v5),nrow=100,byrow=TRUE)


set.seed(1107)
sim.m.sp.3.v1=matrix(rnorm(100*n,mu.sp,sigma.sp.3.v1),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.sp.3.v2=matrix(rnorm(100*n,mu.sp,sigma.sp.3.v2),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.sp.3.v3=matrix(rnorm(100*n,mu.sp,sigma.sp.3.v3),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.sp.3.v4=matrix(rnorm(100*n,mu.sp,sigma.sp.3.v4),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.sp.3.v5=matrix(rnorm(100*n,mu.sp,sigma.sp.3.v5),nrow=100,byrow=TRUE)


set.seed(1107)
sim.m.bump.1.v1=matrix(rnorm(100*n,mu.bump,sigma.bump.1.v1),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.bump.1.v2=matrix(rnorm(100*n,mu.bump,sigma.bump.1.v2),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.bump.1.v3=matrix(rnorm(100*n,mu.bump,sigma.bump.1.v3),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.bump.1.v4=matrix(rnorm(100*n,mu.bump,sigma.bump.1.v4),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.bump.1.v5=matrix(rnorm(100*n,mu.bump,sigma.bump.1.v5),nrow=100,byrow=TRUE)


set.seed(1107)
sim.m.bump.3.v1=matrix(rnorm(100*n,mu.bump,sigma.bump.3.v1),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.bump.3.v2=matrix(rnorm(100*n,mu.bump,sigma.bump.3.v2),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.bump.3.v3=matrix(rnorm(100*n,mu.bump,sigma.bump.3.v3),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.bump.3.v4=matrix(rnorm(100*n,mu.bump,sigma.bump.3.v4),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.bump.3.v5=matrix(rnorm(100*n,mu.bump,sigma.bump.3.v5),nrow=100,byrow=TRUE)



set.seed(1107)
sim.m.blk.1.v1=matrix(rnorm(100*n,mu.blk,sigma.blk.1.v1),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.blk.1.v2=matrix(rnorm(100*n,mu.blk,sigma.blk.1.v2),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.blk.1.v3=matrix(rnorm(100*n,mu.blk,sigma.blk.1.v3),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.blk.1.v4=matrix(rnorm(100*n,mu.blk,sigma.blk.1.v4),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.blk.1.v5=matrix(rnorm(100*n,mu.blk,sigma.blk.1.v5),nrow=100,byrow=TRUE)


set.seed(1107)
sim.m.blk.3.v1=matrix(rnorm(100*n,mu.blk,sigma.blk.3.v1),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.blk.3.v2=matrix(rnorm(100*n,mu.blk,sigma.blk.3.v2),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.blk.3.v3=matrix(rnorm(100*n,mu.blk,sigma.blk.3.v3),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.blk.3.v4=matrix(rnorm(100*n,mu.blk,sigma.blk.3.v4),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.blk.3.v5=matrix(rnorm(100*n,mu.blk,sigma.blk.3.v5),nrow=100,byrow=TRUE)



set.seed(1107)
sim.m.ang.1.v1=matrix(rnorm(100*n,mu.ang,sigma.ang.1.v1),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.ang.1.v2=matrix(rnorm(100*n,mu.ang,sigma.ang.1.v2),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.ang.1.v3=matrix(rnorm(100*n,mu.ang,sigma.ang.1.v3),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.ang.1.v4=matrix(rnorm(100*n,mu.ang,sigma.ang.1.v4),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.ang.1.v5=matrix(rnorm(100*n,mu.ang,sigma.ang.1.v5),nrow=100,byrow=TRUE)


set.seed(1107)
sim.m.ang.3.v1=matrix(rnorm(100*n,mu.ang,sigma.ang.3.v1),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.ang.3.v2=matrix(rnorm(100*n,mu.ang,sigma.ang.3.v2),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.ang.3.v3=matrix(rnorm(100*n,mu.ang,sigma.ang.3.v3),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.ang.3.v4=matrix(rnorm(100*n,mu.ang,sigma.ang.3.v4),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.ang.3.v5=matrix(rnorm(100*n,mu.ang,sigma.ang.3.v5),nrow=100,byrow=TRUE)


set.seed(1107)
sim.m.dop.1.v1=matrix(rnorm(100*n,mu.dop,sigma.dop.1.v1),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.dop.1.v2=matrix(rnorm(100*n,mu.dop,sigma.dop.1.v2),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.dop.1.v3=matrix(rnorm(100*n,mu.dop,sigma.dop.1.v3),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.dop.1.v4=matrix(rnorm(100*n,mu.dop,sigma.dop.1.v4),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.dop.1.v5=matrix(rnorm(100*n,mu.dop,sigma.dop.1.v5),nrow=100,byrow=TRUE)


set.seed(1107)
sim.m.dop.3.v1=matrix(rnorm(100*n,mu.dop,sigma.dop.3.v1),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.dop.3.v2=matrix(rnorm(100*n,mu.dop,sigma.dop.3.v2),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.dop.3.v3=matrix(rnorm(100*n,mu.dop,sigma.dop.3.v3),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.dop.3.v4=matrix(rnorm(100*n,mu.dop,sigma.dop.3.v4),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.dop.3.v5=matrix(rnorm(100*n,mu.dop,sigma.dop.3.v5),nrow=100,byrow=TRUE)

set.seed(1107)
sim.m.blip.1.v1=matrix(rnorm(100*n,mu.blip,sigma.blip.1.v1),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.blip.1.v2=matrix(rnorm(100*n,mu.blip,sigma.blip.1.v2),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.blip.1.v3=matrix(rnorm(100*n,mu.blip,sigma.blip.1.v3),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.blip.1.v4=matrix(rnorm(100*n,mu.blip,sigma.blip.1.v4),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.blip.1.v5=matrix(rnorm(100*n,mu.blip,sigma.blip.1.v5),nrow=100,byrow=TRUE)


set.seed(1107)
sim.m.blip.3.v1=matrix(rnorm(100*n,mu.blip,sigma.blip.3.v1),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.blip.3.v2=matrix(rnorm(100*n,mu.blip,sigma.blip.3.v2),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.blip.3.v3=matrix(rnorm(100*n,mu.blip,sigma.blip.3.v3),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.blip.3.v4=matrix(rnorm(100*n,mu.blip,sigma.blip.3.v4),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.blip.3.v5=matrix(rnorm(100*n,mu.blip,sigma.blip.3.v5),nrow=100,byrow=TRUE)



set.seed(1107)
sim.m.cor.1.v1=matrix(rnorm(100*n,mu.cor,sigma.cor.1.v1),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.cor.1.v2=matrix(rnorm(100*n,mu.cor,sigma.cor.1.v2),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.cor.1.v3=matrix(rnorm(100*n,mu.cor,sigma.cor.1.v3),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.cor.1.v4=matrix(rnorm(100*n,mu.cor,sigma.cor.1.v4),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.cor.1.v5=matrix(rnorm(100*n,mu.cor,sigma.cor.1.v5),nrow=100,byrow=TRUE)


set.seed(1107)
sim.m.cor.3.v1=matrix(rnorm(100*n,mu.cor,sigma.cor.3.v1),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.cor.3.v2=matrix(rnorm(100*n,mu.cor,sigma.cor.3.v2),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.cor.3.v3=matrix(rnorm(100*n,mu.cor,sigma.cor.3.v3),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.cor.3.v4=matrix(rnorm(100*n,mu.cor,sigma.cor.3.v4),nrow=100,byrow=TRUE)
set.seed(1107)
sim.m.cor.3.v5=matrix(rnorm(100*n,mu.cor,sigma.cor.3.v5),nrow=100,byrow=TRUE)



##############################
mu.est.ash.haar.sp.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.sp.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.sp.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.sp.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.sp.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.sp.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.sp.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.sp.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.sp.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.sp.3.v5=matrix(0,nrow=100,ncol=n)

mu.est.ash.s8.sp.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.sp.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.sp.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.sp.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.sp.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.sp.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.sp.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.sp.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.sp.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.sp.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.ashe.haar.sp.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.sp.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.sp.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.sp.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.sp.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.sp.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.sp.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.sp.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.sp.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.sp.3.v5=matrix(0,nrow=100,ncol=n)

mu.est.ashe.s8.sp.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.sp.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.sp.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.sp.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.sp.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.sp.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.sp.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.sp.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.sp.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.sp.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.tie.haar.sp.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.sp.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.sp.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.sp.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.sp.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.sp.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.sp.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.sp.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.sp.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.sp.3.v5=matrix(0,nrow=100,ncol=n)

mu.est.tie.s8.sp.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.sp.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.sp.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.sp.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.sp.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.sp.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.sp.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.sp.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.sp.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.sp.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.asht.haar.sp.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.sp.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.sp.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.sp.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.sp.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.sp.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.sp.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.sp.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.sp.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.sp.3.v5=matrix(0,nrow=100,ncol=n)

mu.est.asht.s8.sp.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.sp.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.sp.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.sp.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.sp.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.sp.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.sp.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.sp.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.sp.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.sp.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.tit.haar.sp.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.sp.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.sp.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.sp.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.sp.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.sp.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.sp.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.sp.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.sp.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.sp.3.v5=matrix(0,nrow=100,ncol=n)

mu.est.tit.s8.sp.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.sp.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.sp.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.sp.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.sp.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.sp.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.sp.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.sp.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.sp.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.sp.3.v5=matrix(0,nrow=100,ncol=n)



mu.est.ash.haar.bump.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.bump.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.bump.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.bump.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.bump.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.bump.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.bump.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.bump.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.bump.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.bump.3.v5=matrix(0,nrow=100,ncol=n)

mu.est.ash.s8.bump.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.bump.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.bump.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.bump.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.bump.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.bump.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.bump.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.bump.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.bump.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.bump.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.ashe.haar.bump.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.bump.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.bump.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.bump.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.bump.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.bump.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.bump.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.bump.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.bump.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.bump.3.v5=matrix(0,nrow=100,ncol=n)

mu.est.ashe.s8.bump.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.bump.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.bump.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.bump.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.bump.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.bump.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.bump.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.bump.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.bump.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.bump.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.tie.haar.bump.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.bump.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.bump.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.bump.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.bump.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.bump.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.bump.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.bump.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.bump.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.bump.3.v5=matrix(0,nrow=100,ncol=n)

mu.est.tie.s8.bump.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.bump.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.bump.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.bump.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.bump.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.bump.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.bump.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.bump.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.bump.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.bump.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.asht.haar.bump.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.bump.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.bump.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.bump.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.bump.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.bump.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.bump.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.bump.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.bump.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.bump.3.v5=matrix(0,nrow=100,ncol=n)

mu.est.asht.s8.bump.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.bump.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.bump.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.bump.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.bump.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.bump.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.bump.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.bump.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.bump.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.bump.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.tit.haar.bump.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.bump.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.bump.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.bump.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.bump.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.bump.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.bump.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.bump.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.bump.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.bump.3.v5=matrix(0,nrow=100,ncol=n)

mu.est.tit.s8.bump.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.bump.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.bump.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.bump.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.bump.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.bump.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.bump.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.bump.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.bump.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.bump.3.v5=matrix(0,nrow=100,ncol=n)



mu.est.ash.haar.blk.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.blk.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.blk.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.blk.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.blk.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.blk.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.blk.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.blk.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.blk.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.blk.3.v5=matrix(0,nrow=100,ncol=n)

mu.est.ash.s8.blk.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.blk.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.blk.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.blk.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.blk.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.blk.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.blk.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.blk.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.blk.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.blk.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.ashe.haar.blk.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.blk.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.blk.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.blk.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.blk.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.blk.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.blk.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.blk.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.blk.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.blk.3.v5=matrix(0,nrow=100,ncol=n)

mu.est.ashe.s8.blk.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.blk.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.blk.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.blk.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.blk.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.blk.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.blk.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.blk.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.blk.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.blk.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.tie.haar.blk.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.blk.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.blk.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.blk.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.blk.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.blk.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.blk.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.blk.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.blk.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.blk.3.v5=matrix(0,nrow=100,ncol=n)

mu.est.tie.s8.blk.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.blk.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.blk.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.blk.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.blk.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.blk.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.blk.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.blk.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.blk.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.blk.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.asht.haar.blk.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.blk.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.blk.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.blk.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.blk.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.blk.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.blk.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.blk.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.blk.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.blk.3.v5=matrix(0,nrow=100,ncol=n)

mu.est.asht.s8.blk.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.blk.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.blk.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.blk.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.blk.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.blk.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.blk.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.blk.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.blk.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.blk.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.tit.haar.blk.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.blk.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.blk.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.blk.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.blk.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.blk.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.blk.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.blk.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.blk.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.blk.3.v5=matrix(0,nrow=100,ncol=n)

mu.est.tit.s8.blk.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.blk.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.blk.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.blk.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.blk.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.blk.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.blk.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.blk.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.blk.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.blk.3.v5=matrix(0,nrow=100,ncol=n)



mu.est.ash.haar.ang.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.ang.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.ang.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.ang.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.ang.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.ang.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.ang.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.ang.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.ang.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.ang.3.v5=matrix(0,nrow=100,ncol=n)

mu.est.ash.s8.ang.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.ang.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.ang.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.ang.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.ang.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.ang.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.ang.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.ang.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.ang.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.ang.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.ashe.haar.ang.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.ang.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.ang.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.ang.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.ang.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.ang.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.ang.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.ang.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.ang.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.ang.3.v5=matrix(0,nrow=100,ncol=n)

mu.est.ashe.s8.ang.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.ang.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.ang.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.ang.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.ang.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.ang.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.ang.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.ang.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.ang.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.ang.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.tie.haar.ang.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.ang.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.ang.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.ang.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.ang.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.ang.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.ang.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.ang.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.ang.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.ang.3.v5=matrix(0,nrow=100,ncol=n)

mu.est.tie.s8.ang.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.ang.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.ang.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.ang.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.ang.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.ang.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.ang.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.ang.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.ang.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.ang.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.asht.haar.ang.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.ang.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.ang.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.ang.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.ang.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.ang.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.ang.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.ang.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.ang.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.ang.3.v5=matrix(0,nrow=100,ncol=n)

mu.est.asht.s8.ang.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.ang.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.ang.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.ang.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.ang.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.ang.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.ang.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.ang.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.ang.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.ang.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.tit.haar.ang.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.ang.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.ang.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.ang.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.ang.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.ang.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.ang.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.ang.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.ang.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.ang.3.v5=matrix(0,nrow=100,ncol=n)

mu.est.tit.s8.ang.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.ang.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.ang.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.ang.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.ang.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.ang.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.ang.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.ang.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.ang.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.ang.3.v5=matrix(0,nrow=100,ncol=n)



mu.est.ash.haar.dop.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.dop.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.dop.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.dop.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.dop.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.dop.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.dop.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.dop.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.dop.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.dop.3.v5=matrix(0,nrow=100,ncol=n)

mu.est.ash.s8.dop.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.dop.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.dop.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.dop.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.dop.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.dop.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.dop.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.dop.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.dop.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.dop.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.ashe.haar.dop.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.dop.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.dop.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.dop.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.dop.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.dop.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.dop.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.dop.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.dop.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.dop.3.v5=matrix(0,nrow=100,ncol=n)

mu.est.ashe.s8.dop.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.dop.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.dop.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.dop.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.dop.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.dop.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.dop.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.dop.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.dop.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.dop.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.tie.haar.dop.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.dop.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.dop.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.dop.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.dop.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.dop.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.dop.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.dop.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.dop.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.dop.3.v5=matrix(0,nrow=100,ncol=n)

mu.est.tie.s8.dop.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.dop.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.dop.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.dop.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.dop.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.dop.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.dop.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.dop.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.dop.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.dop.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.asht.haar.dop.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.dop.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.dop.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.dop.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.dop.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.dop.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.dop.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.dop.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.dop.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.dop.3.v5=matrix(0,nrow=100,ncol=n)

mu.est.asht.s8.dop.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.dop.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.dop.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.dop.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.dop.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.dop.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.dop.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.dop.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.dop.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.dop.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.tit.haar.dop.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.dop.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.dop.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.dop.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.dop.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.dop.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.dop.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.dop.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.dop.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.dop.3.v5=matrix(0,nrow=100,ncol=n)

mu.est.tit.s8.dop.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.dop.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.dop.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.dop.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.dop.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.dop.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.dop.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.dop.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.dop.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.dop.3.v5=matrix(0,nrow=100,ncol=n)



mu.est.ash.haar.blip.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.blip.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.blip.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.blip.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.blip.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.blip.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.blip.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.blip.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.blip.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.blip.3.v5=matrix(0,nrow=100,ncol=n)

mu.est.ash.s8.blip.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.blip.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.blip.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.blip.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.blip.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.blip.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.blip.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.blip.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.blip.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.blip.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.ashe.haar.blip.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.blip.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.blip.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.blip.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.blip.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.blip.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.blip.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.blip.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.blip.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.blip.3.v5=matrix(0,nrow=100,ncol=n)

mu.est.ashe.s8.blip.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.blip.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.blip.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.blip.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.blip.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.blip.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.blip.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.blip.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.blip.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.blip.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.tie.haar.blip.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.blip.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.blip.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.blip.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.blip.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.blip.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.blip.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.blip.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.blip.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.blip.3.v5=matrix(0,nrow=100,ncol=n)

mu.est.tie.s8.blip.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.blip.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.blip.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.blip.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.blip.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.blip.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.blip.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.blip.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.blip.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.blip.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.asht.haar.blip.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.blip.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.blip.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.blip.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.blip.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.blip.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.blip.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.blip.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.blip.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.blip.3.v5=matrix(0,nrow=100,ncol=n)

mu.est.asht.s8.blip.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.blip.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.blip.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.blip.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.blip.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.blip.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.blip.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.blip.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.blip.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.blip.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.tit.haar.blip.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.blip.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.blip.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.blip.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.blip.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.blip.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.blip.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.blip.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.blip.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.blip.3.v5=matrix(0,nrow=100,ncol=n)

mu.est.tit.s8.blip.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.blip.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.blip.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.blip.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.blip.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.blip.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.blip.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.blip.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.blip.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.blip.3.v5=matrix(0,nrow=100,ncol=n)




mu.est.ash.haar.cor.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.cor.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.cor.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.cor.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.cor.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.cor.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.cor.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.cor.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.cor.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ash.haar.cor.3.v5=matrix(0,nrow=100,ncol=n)

mu.est.ash.s8.cor.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.cor.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.cor.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.cor.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.cor.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.cor.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.cor.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.cor.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.cor.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ash.s8.cor.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.ashe.haar.cor.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.cor.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.cor.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.cor.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.cor.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.cor.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.cor.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.cor.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.cor.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.cor.3.v5=matrix(0,nrow=100,ncol=n)

mu.est.ashe.s8.cor.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.cor.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.cor.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.cor.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.cor.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.cor.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.cor.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.cor.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.cor.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.cor.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.tie.haar.cor.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.cor.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.cor.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.cor.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.cor.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.cor.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.cor.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.cor.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.cor.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.tie.haar.cor.3.v5=matrix(0,nrow=100,ncol=n)

mu.est.tie.s8.cor.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.cor.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.cor.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.cor.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.cor.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.cor.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.cor.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.cor.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.cor.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.tie.s8.cor.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.asht.haar.cor.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.cor.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.cor.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.cor.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.cor.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.cor.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.cor.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.cor.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.cor.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.asht.haar.cor.3.v5=matrix(0,nrow=100,ncol=n)

mu.est.asht.s8.cor.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.cor.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.cor.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.cor.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.cor.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.cor.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.cor.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.cor.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.cor.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.asht.s8.cor.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.tit.haar.cor.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.cor.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.cor.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.cor.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.cor.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.cor.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.cor.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.cor.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.cor.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.tit.haar.cor.3.v5=matrix(0,nrow=100,ncol=n)

mu.est.tit.s8.cor.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.cor.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.cor.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.cor.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.cor.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.cor.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.cor.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.cor.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.cor.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.tit.s8.cor.3.v5=matrix(0,nrow=100,ncol=n)





mu.est.tieb.haar.sp.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.sp.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.sp.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.sp.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.sp.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.sp.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.sp.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.sp.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.sp.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.sp.3.v5=matrix(0,nrow=100,ncol=n)

mu.est.tieb.s8.sp.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.sp.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.sp.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.sp.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.sp.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.sp.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.sp.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.sp.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.sp.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.sp.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.tieb.haar.bump.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.bump.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.bump.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.bump.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.bump.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.bump.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.bump.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.bump.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.bump.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.bump.3.v5=matrix(0,nrow=100,ncol=n)

mu.est.tieb.s8.bump.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.bump.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.bump.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.bump.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.bump.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.bump.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.bump.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.bump.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.bump.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.bump.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.tieb.haar.blk.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.blk.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.blk.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.blk.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.blk.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.blk.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.blk.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.blk.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.blk.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.blk.3.v5=matrix(0,nrow=100,ncol=n)

mu.est.tieb.s8.blk.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.blk.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.blk.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.blk.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.blk.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.blk.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.blk.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.blk.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.blk.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.blk.3.v5=matrix(0,nrow=100,ncol=n)



mu.est.tieb.haar.ang.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.ang.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.ang.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.ang.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.ang.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.ang.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.ang.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.ang.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.ang.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.ang.3.v5=matrix(0,nrow=100,ncol=n)

mu.est.tieb.s8.ang.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.ang.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.ang.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.ang.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.ang.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.ang.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.ang.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.ang.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.ang.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.ang.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.tieb.haar.dop.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.dop.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.dop.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.dop.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.dop.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.dop.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.dop.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.dop.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.dop.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.dop.3.v5=matrix(0,nrow=100,ncol=n)

mu.est.tieb.s8.dop.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.dop.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.dop.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.dop.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.dop.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.dop.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.dop.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.dop.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.dop.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.dop.3.v5=matrix(0,nrow=100,ncol=n)



mu.est.tieb.haar.blip.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.blip.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.blip.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.blip.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.blip.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.blip.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.blip.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.blip.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.blip.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.blip.3.v5=matrix(0,nrow=100,ncol=n)

mu.est.tieb.s8.blip.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.blip.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.blip.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.blip.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.blip.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.blip.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.blip.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.blip.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.blip.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.blip.3.v5=matrix(0,nrow=100,ncol=n)



mu.est.tieb.haar.cor.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.cor.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.cor.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.cor.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.cor.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.cor.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.cor.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.cor.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.cor.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.tieb.haar.cor.3.v5=matrix(0,nrow=100,ncol=n)

mu.est.tieb.s8.cor.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.cor.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.cor.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.cor.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.cor.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.cor.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.cor.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.cor.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.cor.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.tieb.s8.cor.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.ebayes.sp.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.sp.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.sp.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.sp.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.sp.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.sp.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.sp.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.sp.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.sp.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.sp.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.ebayes.bump.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.bump.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.bump.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.bump.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.bump.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.bump.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.bump.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.bump.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.bump.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.bump.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.ebayes.blk.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.blk.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.blk.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.blk.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.blk.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.blk.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.blk.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.blk.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.blk.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.blk.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.ebayes.ang.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.ang.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.ang.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.ang.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.ang.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.ang.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.ang.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.ang.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.ang.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.ang.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.ebayes.dop.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.dop.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.dop.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.dop.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.dop.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.dop.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.dop.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.dop.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.dop.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.dop.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.ebayes.blip.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.blip.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.blip.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.blip.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.blip.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.blip.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.blip.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.blip.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.blip.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.blip.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.ebayes.cor.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.cor.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.cor.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.cor.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.cor.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.cor.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.cor.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.cor.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.cor.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ebayes.cor.3.v5=matrix(0,nrow=100,ncol=n)



mu.est.ashe.haar.j.sp.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.sp.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.sp.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.sp.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.sp.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.sp.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.sp.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.sp.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.sp.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.sp.3.v5=matrix(0,nrow=100,ncol=n)

mu.est.ashe.s8.j.sp.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.sp.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.sp.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.sp.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.sp.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.sp.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.sp.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.sp.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.sp.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.sp.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.ashe.haar.j.bump.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.bump.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.bump.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.bump.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.bump.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.bump.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.bump.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.bump.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.bump.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.bump.3.v5=matrix(0,nrow=100,ncol=n)

mu.est.ashe.s8.j.bump.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.bump.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.bump.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.bump.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.bump.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.bump.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.bump.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.bump.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.bump.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.bump.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.ashe.haar.j.blk.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.blk.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.blk.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.blk.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.blk.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.blk.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.blk.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.blk.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.blk.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.blk.3.v5=matrix(0,nrow=100,ncol=n)

mu.est.ashe.s8.j.blk.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.blk.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.blk.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.blk.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.blk.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.blk.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.blk.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.blk.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.blk.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.blk.3.v5=matrix(0,nrow=100,ncol=n)



mu.est.ashe.haar.j.ang.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.ang.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.ang.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.ang.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.ang.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.ang.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.ang.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.ang.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.ang.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.ang.3.v5=matrix(0,nrow=100,ncol=n)

mu.est.ashe.s8.j.ang.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.ang.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.ang.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.ang.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.ang.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.ang.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.ang.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.ang.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.ang.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.ang.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.ashe.haar.j.dop.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.dop.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.dop.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.dop.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.dop.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.dop.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.dop.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.dop.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.dop.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.dop.3.v5=matrix(0,nrow=100,ncol=n)

mu.est.ashe.s8.j.dop.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.dop.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.dop.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.dop.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.dop.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.dop.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.dop.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.dop.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.dop.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.dop.3.v5=matrix(0,nrow=100,ncol=n)



mu.est.ashe.haar.j.blip.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.blip.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.blip.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.blip.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.blip.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.blip.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.blip.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.blip.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.blip.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.blip.3.v5=matrix(0,nrow=100,ncol=n)

mu.est.ashe.s8.j.blip.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.blip.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.blip.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.blip.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.blip.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.blip.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.blip.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.blip.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.blip.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.blip.3.v5=matrix(0,nrow=100,ncol=n)



mu.est.ashe.haar.j.cor.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.cor.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.cor.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.cor.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.cor.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.cor.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.cor.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.cor.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.cor.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.haar.j.cor.3.v5=matrix(0,nrow=100,ncol=n)

mu.est.ashe.s8.j.cor.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.cor.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.cor.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.cor.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.cor.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.cor.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.cor.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.cor.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.cor.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.cor.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.ashe.s8.j.h.sp.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.sp.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.sp.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.sp.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.sp.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.sp.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.sp.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.sp.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.sp.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.sp.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.ashe.s8.j.h.bump.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.bump.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.bump.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.bump.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.bump.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.bump.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.bump.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.bump.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.bump.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.bump.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.ashe.s8.j.h.blk.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.blk.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.blk.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.blk.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.blk.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.blk.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.blk.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.blk.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.blk.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.blk.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.ashe.s8.j.h.ang.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.ang.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.ang.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.ang.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.ang.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.ang.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.ang.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.ang.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.ang.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.ang.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.ashe.s8.j.h.dop.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.dop.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.dop.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.dop.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.dop.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.dop.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.dop.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.dop.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.dop.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.dop.3.v5=matrix(0,nrow=100,ncol=n)



mu.est.ashe.s8.j.h.blip.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.blip.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.blip.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.blip.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.blip.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.blip.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.blip.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.blip.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.blip.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.blip.3.v5=matrix(0,nrow=100,ncol=n)



mu.est.ashe.s8.j.h.cor.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.cor.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.cor.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.cor.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.cor.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.cor.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.cor.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.cor.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.cor.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashe.s8.j.h.cor.3.v5=matrix(0,nrow=100,ncol=n)





mu.est.ashef.haar.sp.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.sp.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.sp.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.sp.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.sp.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.sp.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.sp.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.sp.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.sp.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.sp.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.ashef.haar.bump.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.bump.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.bump.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.bump.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.bump.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.bump.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.bump.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.bump.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.bump.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.bump.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.ashef.haar.blk.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.blk.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.blk.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.blk.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.blk.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.blk.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.blk.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.blk.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.blk.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.blk.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.ashef.haar.ang.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.ang.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.ang.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.ang.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.ang.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.ang.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.ang.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.ang.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.ang.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.ang.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.ashef.haar.dop.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.dop.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.dop.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.dop.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.dop.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.dop.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.dop.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.dop.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.dop.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.dop.3.v5=matrix(0,nrow=100,ncol=n)



mu.est.ashef.haar.blip.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.blip.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.blip.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.blip.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.blip.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.blip.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.blip.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.blip.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.blip.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.blip.3.v5=matrix(0,nrow=100,ncol=n)



mu.est.ashef.haar.cor.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.cor.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.cor.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.cor.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.cor.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.cor.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.cor.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.cor.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.cor.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashef.haar.cor.3.v5=matrix(0,nrow=100,ncol=n)



mu.est.ashef.s8.sp.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.sp.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.sp.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.sp.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.sp.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.sp.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.sp.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.sp.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.sp.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.sp.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.ashef.s8.bump.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.bump.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.bump.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.bump.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.bump.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.bump.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.bump.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.bump.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.bump.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.bump.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.ashef.s8.blk.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.blk.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.blk.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.blk.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.blk.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.blk.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.blk.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.blk.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.blk.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.blk.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.ashef.s8.ang.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.ang.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.ang.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.ang.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.ang.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.ang.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.ang.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.ang.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.ang.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.ang.3.v5=matrix(0,nrow=100,ncol=n)


mu.est.ashef.s8.dop.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.dop.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.dop.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.dop.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.dop.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.dop.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.dop.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.dop.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.dop.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.dop.3.v5=matrix(0,nrow=100,ncol=n)



mu.est.ashef.s8.blip.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.blip.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.blip.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.blip.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.blip.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.blip.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.blip.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.blip.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.blip.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.blip.3.v5=matrix(0,nrow=100,ncol=n)



mu.est.ashef.s8.cor.1.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.cor.1.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.cor.1.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.cor.1.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.cor.1.v5=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.cor.3.v1=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.cor.3.v2=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.cor.3.v3=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.cor.3.v4=matrix(0,nrow=100,ncol=n)
mu.est.ashef.s8.cor.3.v5=matrix(0,nrow=100,ncol=n)


for(i in 1:100){
  sigma.sp.1.v1.est=sig.est.func(sim.m.sp.1.v1[i,],n)
  mu.est.ash.haar.sp.1.v1[i,]=bayesmooth(sim.m.sp.1.v1[i,],sigma.sp.1.v1.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.sp.1.v1[i,]=bayesmooth(sim.m.sp.1.v1[i,],sigma.sp.1.v1.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.sp.1.v1[i,]=bayesmooth(sim.m.sp.1.v1[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.sp.1.v1[i,]=bayesmooth(sim.m.sp.1.v1[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.sp.1.v1[i,]=ti.thresh(sim.m.sp.1.v1[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.sp.1.v1[i,]=ti.thresh(sim.m.sp.1.v1[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.sp.1.v1[i,]=bayesmooth(sim.m.sp.1.v1[i,],sigma=sigma.sp.1.v1,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.sp.1.v1[i,]=bayesmooth(sim.m.sp.1.v1[i,],sigma=sigma.sp.1.v1,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.sp.1.v1[i,]=ti.thresh(sim.m.sp.1.v1[i,],sigma=sigma.sp.1.v1,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.sp.1.v1[i,]=ti.thresh(sim.m.sp.1.v1[i,],sigma=sigma.sp.1.v1,filter.number=8,family="DaubLeAsymm")
  sigma.sp.3.v1.est=sig.est.func(sim.m.sp.3.v1[i,],n)
  mu.est.ash.haar.sp.3.v1[i,]=bayesmooth(sim.m.sp.3.v1[i,],sigma.sp.3.v1.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.sp.3.v1[i,]=bayesmooth(sim.m.sp.3.v1[i,],sigma.sp.3.v1.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.sp.3.v1[i,]=bayesmooth(sim.m.sp.3.v1[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.sp.3.v1[i,]=bayesmooth(sim.m.sp.3.v1[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.sp.3.v1[i,]=ti.thresh(sim.m.sp.3.v1[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.sp.3.v1[i,]=ti.thresh(sim.m.sp.3.v1[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.sp.3.v1[i,]=bayesmooth(sim.m.sp.3.v1[i,],sigma=sigma.sp.3.v1,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.sp.3.v1[i,]=bayesmooth(sim.m.sp.3.v1[i,],sigma=sigma.sp.3.v1,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.sp.3.v1[i,]=ti.thresh(sim.m.sp.3.v1[i,],sigma=sigma.sp.3.v1,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.sp.3.v1[i,]=ti.thresh(sim.m.sp.3.v1[i,],sigma=sigma.sp.3.v1,filter.number=8,family="DaubLeAsymm")

  sigma.bump.1.v1.est=sig.est.func(sim.m.bump.1.v1[i,],n)
  mu.est.ash.haar.bump.1.v1[i,]=bayesmooth(sim.m.bump.1.v1[i,],sigma.bump.1.v1.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.bump.1.v1[i,]=bayesmooth(sim.m.bump.1.v1[i,],sigma.bump.1.v1.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.bump.1.v1[i,]=bayesmooth(sim.m.bump.1.v1[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.bump.1.v1[i,]=bayesmooth(sim.m.bump.1.v1[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.bump.1.v1[i,]=ti.thresh(sim.m.bump.1.v1[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.bump.1.v1[i,]=ti.thresh(sim.m.bump.1.v1[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.bump.1.v1[i,]=bayesmooth(sim.m.bump.1.v1[i,],sigma=sigma.bump.1.v1,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.bump.1.v1[i,]=bayesmooth(sim.m.bump.1.v1[i,],sigma=sigma.bump.1.v1,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.bump.1.v1[i,]=ti.thresh(sim.m.bump.1.v1[i,],sigma=sigma.bump.1.v1,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.bump.1.v1[i,]=ti.thresh(sim.m.bump.1.v1[i,],sigma=sigma.bump.1.v1,filter.number=8,family="DaubLeAsymm")
  sigma.bump.3.v1.est=sig.est.func(sim.m.bump.3.v1[i,],n)
  mu.est.ash.haar.bump.3.v1[i,]=bayesmooth(sim.m.bump.3.v1[i,],sigma.bump.3.v1.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.bump.3.v1[i,]=bayesmooth(sim.m.bump.3.v1[i,],sigma.bump.3.v1.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.bump.3.v1[i,]=bayesmooth(sim.m.bump.3.v1[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.bump.3.v1[i,]=bayesmooth(sim.m.bump.3.v1[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.bump.3.v1[i,]=ti.thresh(sim.m.bump.3.v1[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.bump.3.v1[i,]=ti.thresh(sim.m.bump.3.v1[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.bump.3.v1[i,]=bayesmooth(sim.m.bump.3.v1[i,],sigma=sigma.bump.3.v1,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.bump.3.v1[i,]=bayesmooth(sim.m.bump.3.v1[i,],sigma=sigma.bump.3.v1,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.bump.3.v1[i,]=ti.thresh(sim.m.bump.3.v1[i,],sigma=sigma.bump.3.v1,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.bump.3.v1[i,]=ti.thresh(sim.m.bump.3.v1[i,],sigma=sigma.bump.3.v1,filter.number=8,family="DaubLeAsymm")
  
  sigma.blk.1.v1.est=sig.est.func(sim.m.blk.1.v1[i,],n)
  mu.est.ash.haar.blk.1.v1[i,]=bayesmooth(sim.m.blk.1.v1[i,],sigma.blk.1.v1.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.blk.1.v1[i,]=bayesmooth(sim.m.blk.1.v1[i,],sigma.blk.1.v1.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.blk.1.v1[i,]=bayesmooth(sim.m.blk.1.v1[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.blk.1.v1[i,]=bayesmooth(sim.m.blk.1.v1[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.blk.1.v1[i,]=ti.thresh(sim.m.blk.1.v1[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.blk.1.v1[i,]=ti.thresh(sim.m.blk.1.v1[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.blk.1.v1[i,]=bayesmooth(sim.m.blk.1.v1[i,],sigma=sigma.blk.1.v1,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.blk.1.v1[i,]=bayesmooth(sim.m.blk.1.v1[i,],sigma=sigma.blk.1.v1,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.blk.1.v1[i,]=ti.thresh(sim.m.blk.1.v1[i,],sigma=sigma.blk.1.v1,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.blk.1.v1[i,]=ti.thresh(sim.m.blk.1.v1[i,],sigma=sigma.blk.1.v1,filter.number=8,family="DaubLeAsymm")
  sigma.blk.3.v1.est=sig.est.func(sim.m.blk.3.v1[i,],n)
  mu.est.ash.haar.blk.3.v1[i,]=bayesmooth(sim.m.blk.3.v1[i,],sigma.blk.3.v1.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.blk.3.v1[i,]=bayesmooth(sim.m.blk.3.v1[i,],sigma.blk.3.v1.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.blk.3.v1[i,]=bayesmooth(sim.m.blk.3.v1[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.blk.3.v1[i,]=bayesmooth(sim.m.blk.3.v1[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.blk.3.v1[i,]=ti.thresh(sim.m.blk.3.v1[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.blk.3.v1[i,]=ti.thresh(sim.m.blk.3.v1[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.blk.3.v1[i,]=bayesmooth(sim.m.blk.3.v1[i,],sigma=sigma.blk.3.v1,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.blk.3.v1[i,]=bayesmooth(sim.m.blk.3.v1[i,],sigma=sigma.blk.3.v1,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.blk.3.v1[i,]=ti.thresh(sim.m.blk.3.v1[i,],sigma=sigma.blk.3.v1,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.blk.3.v1[i,]=ti.thresh(sim.m.blk.3.v1[i,],sigma=sigma.blk.3.v1,filter.number=8,family="DaubLeAsymm")

  sigma.ang.1.v1.est=sig.est.func(sim.m.ang.1.v1[i,],n)
  mu.est.ash.haar.ang.1.v1[i,]=bayesmooth(sim.m.ang.1.v1[i,],sigma.ang.1.v1.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.ang.1.v1[i,]=bayesmooth(sim.m.ang.1.v1[i,],sigma.ang.1.v1.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.ang.1.v1[i,]=bayesmooth(sim.m.ang.1.v1[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.ang.1.v1[i,]=bayesmooth(sim.m.ang.1.v1[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.ang.1.v1[i,]=ti.thresh(sim.m.ang.1.v1[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.ang.1.v1[i,]=ti.thresh(sim.m.ang.1.v1[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.ang.1.v1[i,]=bayesmooth(sim.m.ang.1.v1[i,],sigma=sigma.ang.1.v1,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.ang.1.v1[i,]=bayesmooth(sim.m.ang.1.v1[i,],sigma=sigma.ang.1.v1,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.ang.1.v1[i,]=ti.thresh(sim.m.ang.1.v1[i,],sigma=sigma.ang.1.v1,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.ang.1.v1[i,]=ti.thresh(sim.m.ang.1.v1[i,],sigma=sigma.ang.1.v1,filter.number=8,family="DaubLeAsymm")
  sigma.ang.3.v1.est=sig.est.func(sim.m.ang.3.v1[i,],n)
  mu.est.ash.haar.ang.3.v1[i,]=bayesmooth(sim.m.ang.3.v1[i,],sigma.ang.3.v1.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.ang.3.v1[i,]=bayesmooth(sim.m.ang.3.v1[i,],sigma.ang.3.v1.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.ang.3.v1[i,]=bayesmooth(sim.m.ang.3.v1[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.ang.3.v1[i,]=bayesmooth(sim.m.ang.3.v1[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.ang.3.v1[i,]=ti.thresh(sim.m.ang.3.v1[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.ang.3.v1[i,]=ti.thresh(sim.m.ang.3.v1[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.ang.3.v1[i,]=bayesmooth(sim.m.ang.3.v1[i,],sigma=sigma.ang.3.v1,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.ang.3.v1[i,]=bayesmooth(sim.m.ang.3.v1[i,],sigma=sigma.ang.3.v1,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.ang.3.v1[i,]=ti.thresh(sim.m.ang.3.v1[i,],sigma=sigma.ang.3.v1,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.ang.3.v1[i,]=ti.thresh(sim.m.ang.3.v1[i,],sigma=sigma.ang.3.v1,filter.number=8,family="DaubLeAsymm")

  sigma.dop.1.v1.est=sig.est.func(sim.m.dop.1.v1[i,],n)
  mu.est.ash.haar.dop.1.v1[i,]=bayesmooth(sim.m.dop.1.v1[i,],sigma.dop.1.v1.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.dop.1.v1[i,]=bayesmooth(sim.m.dop.1.v1[i,],sigma.dop.1.v1.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.dop.1.v1[i,]=bayesmooth(sim.m.dop.1.v1[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.dop.1.v1[i,]=bayesmooth(sim.m.dop.1.v1[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.dop.1.v1[i,]=ti.thresh(sim.m.dop.1.v1[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.dop.1.v1[i,]=ti.thresh(sim.m.dop.1.v1[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.dop.1.v1[i,]=bayesmooth(sim.m.dop.1.v1[i,],sigma=sigma.dop.1.v1,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.dop.1.v1[i,]=bayesmooth(sim.m.dop.1.v1[i,],sigma=sigma.dop.1.v1,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.dop.1.v1[i,]=ti.thresh(sim.m.dop.1.v1[i,],sigma=sigma.dop.1.v1,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.dop.1.v1[i,]=ti.thresh(sim.m.dop.1.v1[i,],sigma=sigma.dop.1.v1,filter.number=8,family="DaubLeAsymm")
  sigma.dop.3.v1.est=sig.est.func(sim.m.dop.3.v1[i,],n)
  mu.est.ash.haar.dop.3.v1[i,]=bayesmooth(sim.m.dop.3.v1[i,],sigma.dop.3.v1.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.dop.3.v1[i,]=bayesmooth(sim.m.dop.3.v1[i,],sigma.dop.3.v1.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.dop.3.v1[i,]=bayesmooth(sim.m.dop.3.v1[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.dop.3.v1[i,]=bayesmooth(sim.m.dop.3.v1[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.dop.3.v1[i,]=ti.thresh(sim.m.dop.3.v1[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.dop.3.v1[i,]=ti.thresh(sim.m.dop.3.v1[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.dop.3.v1[i,]=bayesmooth(sim.m.dop.3.v1[i,],sigma=sigma.dop.3.v1,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.dop.3.v1[i,]=bayesmooth(sim.m.dop.3.v1[i,],sigma=sigma.dop.3.v1,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.dop.3.v1[i,]=ti.thresh(sim.m.dop.3.v1[i,],sigma=sigma.dop.3.v1,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.dop.3.v1[i,]=ti.thresh(sim.m.dop.3.v1[i,],sigma=sigma.dop.3.v1,filter.number=8,family="DaubLeAsymm")

  sigma.blip.1.v1.est=sig.est.func(sim.m.blip.1.v1[i,],n)
  mu.est.ash.haar.blip.1.v1[i,]=bayesmooth(sim.m.blip.1.v1[i,],sigma.blip.1.v1.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.blip.1.v1[i,]=bayesmooth(sim.m.blip.1.v1[i,],sigma.blip.1.v1.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.blip.1.v1[i,]=bayesmooth(sim.m.blip.1.v1[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.blip.1.v1[i,]=bayesmooth(sim.m.blip.1.v1[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.blip.1.v1[i,]=ti.thresh(sim.m.blip.1.v1[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.blip.1.v1[i,]=ti.thresh(sim.m.blip.1.v1[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.blip.1.v1[i,]=bayesmooth(sim.m.blip.1.v1[i,],sigma=sigma.blip.1.v1,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.blip.1.v1[i,]=bayesmooth(sim.m.blip.1.v1[i,],sigma=sigma.blip.1.v1,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.blip.1.v1[i,]=ti.thresh(sim.m.blip.1.v1[i,],sigma=sigma.blip.1.v1,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.blip.1.v1[i,]=ti.thresh(sim.m.blip.1.v1[i,],sigma=sigma.blip.1.v1,filter.number=8,family="DaubLeAsymm")
  sigma.blip.3.v1.est=sig.est.func(sim.m.blip.3.v1[i,],n)
  mu.est.ash.haar.blip.3.v1[i,]=bayesmooth(sim.m.blip.3.v1[i,],sigma.blip.3.v1.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.blip.3.v1[i,]=bayesmooth(sim.m.blip.3.v1[i,],sigma.blip.3.v1.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.blip.3.v1[i,]=bayesmooth(sim.m.blip.3.v1[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.blip.3.v1[i,]=bayesmooth(sim.m.blip.3.v1[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.blip.3.v1[i,]=ti.thresh(sim.m.blip.3.v1[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.blip.3.v1[i,]=ti.thresh(sim.m.blip.3.v1[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.blip.3.v1[i,]=bayesmooth(sim.m.blip.3.v1[i,],sigma=sigma.blip.3.v1,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.blip.3.v1[i,]=bayesmooth(sim.m.blip.3.v1[i,],sigma=sigma.blip.3.v1,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.blip.3.v1[i,]=ti.thresh(sim.m.blip.3.v1[i,],sigma=sigma.blip.3.v1,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.blip.3.v1[i,]=ti.thresh(sim.m.blip.3.v1[i,],sigma=sigma.blip.3.v1,filter.number=8,family="DaubLeAsymm")

  sigma.cor.1.v1.est=sig.est.func(sim.m.cor.1.v1[i,],n)
  mu.est.ash.haar.cor.1.v1[i,]=bayesmooth(sim.m.cor.1.v1[i,],sigma.cor.1.v1.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.cor.1.v1[i,]=bayesmooth(sim.m.cor.1.v1[i,],sigma.cor.1.v1.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.cor.1.v1[i,]=bayesmooth(sim.m.cor.1.v1[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.cor.1.v1[i,]=bayesmooth(sim.m.cor.1.v1[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.cor.1.v1[i,]=ti.thresh(sim.m.cor.1.v1[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.cor.1.v1[i,]=ti.thresh(sim.m.cor.1.v1[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.cor.1.v1[i,]=bayesmooth(sim.m.cor.1.v1[i,],sigma=sigma.cor.1.v1,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.cor.1.v1[i,]=bayesmooth(sim.m.cor.1.v1[i,],sigma=sigma.cor.1.v1,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.cor.1.v1[i,]=ti.thresh(sim.m.cor.1.v1[i,],sigma=sigma.cor.1.v1,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.cor.1.v1[i,]=ti.thresh(sim.m.cor.1.v1[i,],sigma=sigma.cor.1.v1,filter.number=8,family="DaubLeAsymm")
  sigma.cor.3.v1.est=sig.est.func(sim.m.cor.3.v1[i,],n)
  mu.est.ash.haar.cor.3.v1[i,]=bayesmooth(sim.m.cor.3.v1[i,],sigma.cor.3.v1.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.cor.3.v1[i,]=bayesmooth(sim.m.cor.3.v1[i,],sigma.cor.3.v1.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.cor.3.v1[i,]=bayesmooth(sim.m.cor.3.v1[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.cor.3.v1[i,]=bayesmooth(sim.m.cor.3.v1[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.cor.3.v1[i,]=ti.thresh(sim.m.cor.3.v1[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.cor.3.v1[i,]=ti.thresh(sim.m.cor.3.v1[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.cor.3.v1[i,]=bayesmooth(sim.m.cor.3.v1[i,],sigma=sigma.cor.3.v1,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.cor.3.v1[i,]=bayesmooth(sim.m.cor.3.v1[i,],sigma=sigma.cor.3.v1,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.cor.3.v1[i,]=ti.thresh(sim.m.cor.3.v1[i,],sigma=sigma.cor.3.v1,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.cor.3.v1[i,]=ti.thresh(sim.m.cor.3.v1[i,],sigma=sigma.cor.3.v1,filter.number=8,family="DaubLeAsymm")



  sigma.sp.1.v2.est=sig.est.func(sim.m.sp.1.v2[i,],n)
  mu.est.ash.haar.sp.1.v2[i,]=bayesmooth(sim.m.sp.1.v2[i,],sigma.sp.1.v2.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.sp.1.v2[i,]=bayesmooth(sim.m.sp.1.v2[i,],sigma.sp.1.v2.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.sp.1.v2[i,]=bayesmooth(sim.m.sp.1.v2[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.sp.1.v2[i,]=bayesmooth(sim.m.sp.1.v2[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.sp.1.v2[i,]=ti.thresh(sim.m.sp.1.v2[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.sp.1.v2[i,]=ti.thresh(sim.m.sp.1.v2[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.sp.1.v2[i,]=bayesmooth(sim.m.sp.1.v2[i,],sigma=sigma.sp.1.v2,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.sp.1.v2[i,]=bayesmooth(sim.m.sp.1.v2[i,],sigma=sigma.sp.1.v2,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.sp.1.v2[i,]=ti.thresh(sim.m.sp.1.v2[i,],sigma=sigma.sp.1.v2,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.sp.1.v2[i,]=ti.thresh(sim.m.sp.1.v2[i,],sigma=sigma.sp.1.v2,filter.number=8,family="DaubLeAsymm")
  sigma.sp.3.v2.est=sig.est.func(sim.m.sp.3.v2[i,],n)
  mu.est.ash.haar.sp.3.v2[i,]=bayesmooth(sim.m.sp.3.v2[i,],sigma.sp.3.v2.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.sp.3.v2[i,]=bayesmooth(sim.m.sp.3.v2[i,],sigma.sp.3.v2.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.sp.3.v2[i,]=bayesmooth(sim.m.sp.3.v2[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.sp.3.v2[i,]=bayesmooth(sim.m.sp.3.v2[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.sp.3.v2[i,]=ti.thresh(sim.m.sp.3.v2[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.sp.3.v2[i,]=ti.thresh(sim.m.sp.3.v2[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.sp.3.v2[i,]=bayesmooth(sim.m.sp.3.v2[i,],sigma=sigma.sp.3.v2,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.sp.3.v2[i,]=bayesmooth(sim.m.sp.3.v2[i,],sigma=sigma.sp.3.v2,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.sp.3.v2[i,]=ti.thresh(sim.m.sp.3.v2[i,],sigma=sigma.sp.3.v2,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.sp.3.v2[i,]=ti.thresh(sim.m.sp.3.v2[i,],sigma=sigma.sp.3.v2,filter.number=8,family="DaubLeAsymm")

  sigma.bump.1.v2.est=sig.est.func(sim.m.bump.1.v2[i,],n)
  mu.est.ash.haar.bump.1.v2[i,]=bayesmooth(sim.m.bump.1.v2[i,],sigma.bump.1.v2.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.bump.1.v2[i,]=bayesmooth(sim.m.bump.1.v2[i,],sigma.bump.1.v2.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.bump.1.v2[i,]=bayesmooth(sim.m.bump.1.v2[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.bump.1.v2[i,]=bayesmooth(sim.m.bump.1.v2[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.bump.1.v2[i,]=ti.thresh(sim.m.bump.1.v2[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.bump.1.v2[i,]=ti.thresh(sim.m.bump.1.v2[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.bump.1.v2[i,]=bayesmooth(sim.m.bump.1.v2[i,],sigma=sigma.bump.1.v2,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.bump.1.v2[i,]=bayesmooth(sim.m.bump.1.v2[i,],sigma=sigma.bump.1.v2,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.bump.1.v2[i,]=ti.thresh(sim.m.bump.1.v2[i,],sigma=sigma.bump.1.v2,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.bump.1.v2[i,]=ti.thresh(sim.m.bump.1.v2[i,],sigma=sigma.bump.1.v2,filter.number=8,family="DaubLeAsymm")
  sigma.bump.3.v2.est=sig.est.func(sim.m.bump.3.v2[i,],n)
  mu.est.ash.haar.bump.3.v2[i,]=bayesmooth(sim.m.bump.3.v2[i,],sigma.bump.3.v2.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.bump.3.v2[i,]=bayesmooth(sim.m.bump.3.v2[i,],sigma.bump.3.v2.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.bump.3.v2[i,]=bayesmooth(sim.m.bump.3.v2[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.bump.3.v2[i,]=bayesmooth(sim.m.bump.3.v2[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.bump.3.v2[i,]=ti.thresh(sim.m.bump.3.v2[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.bump.3.v2[i,]=ti.thresh(sim.m.bump.3.v2[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.bump.3.v2[i,]=bayesmooth(sim.m.bump.3.v2[i,],sigma=sigma.bump.3.v2,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.bump.3.v2[i,]=bayesmooth(sim.m.bump.3.v2[i,],sigma=sigma.bump.3.v2,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.bump.3.v2[i,]=ti.thresh(sim.m.bump.3.v2[i,],sigma=sigma.bump.3.v2,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.bump.3.v2[i,]=ti.thresh(sim.m.bump.3.v2[i,],sigma=sigma.bump.3.v2,filter.number=8,family="DaubLeAsymm")
  
  sigma.blk.1.v2.est=sig.est.func(sim.m.blk.1.v2[i,],n)
  mu.est.ash.haar.blk.1.v2[i,]=bayesmooth(sim.m.blk.1.v2[i,],sigma.blk.1.v2.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.blk.1.v2[i,]=bayesmooth(sim.m.blk.1.v2[i,],sigma.blk.1.v2.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.blk.1.v2[i,]=bayesmooth(sim.m.blk.1.v2[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.blk.1.v2[i,]=bayesmooth(sim.m.blk.1.v2[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.blk.1.v2[i,]=ti.thresh(sim.m.blk.1.v2[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.blk.1.v2[i,]=ti.thresh(sim.m.blk.1.v2[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.blk.1.v2[i,]=bayesmooth(sim.m.blk.1.v2[i,],sigma=sigma.blk.1.v2,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.blk.1.v2[i,]=bayesmooth(sim.m.blk.1.v2[i,],sigma=sigma.blk.1.v2,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.blk.1.v2[i,]=ti.thresh(sim.m.blk.1.v2[i,],sigma=sigma.blk.1.v2,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.blk.1.v2[i,]=ti.thresh(sim.m.blk.1.v2[i,],sigma=sigma.blk.1.v2,filter.number=8,family="DaubLeAsymm")
  sigma.blk.3.v2.est=sig.est.func(sim.m.blk.3.v2[i,],n)
  mu.est.ash.haar.blk.3.v2[i,]=bayesmooth(sim.m.blk.3.v2[i,],sigma.blk.3.v2.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.blk.3.v2[i,]=bayesmooth(sim.m.blk.3.v2[i,],sigma.blk.3.v2.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.blk.3.v2[i,]=bayesmooth(sim.m.blk.3.v2[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.blk.3.v2[i,]=bayesmooth(sim.m.blk.3.v2[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.blk.3.v2[i,]=ti.thresh(sim.m.blk.3.v2[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.blk.3.v2[i,]=ti.thresh(sim.m.blk.3.v2[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.blk.3.v2[i,]=bayesmooth(sim.m.blk.3.v2[i,],sigma=sigma.blk.3.v2,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.blk.3.v2[i,]=bayesmooth(sim.m.blk.3.v2[i,],sigma=sigma.blk.3.v2,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.blk.3.v2[i,]=ti.thresh(sim.m.blk.3.v2[i,],sigma=sigma.blk.3.v2,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.blk.3.v2[i,]=ti.thresh(sim.m.blk.3.v2[i,],sigma=sigma.blk.3.v2,filter.number=8,family="DaubLeAsymm")

  sigma.ang.1.v2.est=sig.est.func(sim.m.ang.1.v2[i,],n)
  mu.est.ash.haar.ang.1.v2[i,]=bayesmooth(sim.m.ang.1.v2[i,],sigma.ang.1.v2.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.ang.1.v2[i,]=bayesmooth(sim.m.ang.1.v2[i,],sigma.ang.1.v2.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.ang.1.v2[i,]=bayesmooth(sim.m.ang.1.v2[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.ang.1.v2[i,]=bayesmooth(sim.m.ang.1.v2[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.ang.1.v2[i,]=ti.thresh(sim.m.ang.1.v2[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.ang.1.v2[i,]=ti.thresh(sim.m.ang.1.v2[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.ang.1.v2[i,]=bayesmooth(sim.m.ang.1.v2[i,],sigma=sigma.ang.1.v2,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.ang.1.v2[i,]=bayesmooth(sim.m.ang.1.v2[i,],sigma=sigma.ang.1.v2,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.ang.1.v2[i,]=ti.thresh(sim.m.ang.1.v2[i,],sigma=sigma.ang.1.v2,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.ang.1.v2[i,]=ti.thresh(sim.m.ang.1.v2[i,],sigma=sigma.ang.1.v2,filter.number=8,family="DaubLeAsymm")
  sigma.ang.3.v2.est=sig.est.func(sim.m.ang.3.v2[i,],n)
  mu.est.ash.haar.ang.3.v2[i,]=bayesmooth(sim.m.ang.3.v2[i,],sigma.ang.3.v2.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.ang.3.v2[i,]=bayesmooth(sim.m.ang.3.v2[i,],sigma.ang.3.v2.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.ang.3.v2[i,]=bayesmooth(sim.m.ang.3.v2[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.ang.3.v2[i,]=bayesmooth(sim.m.ang.3.v2[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.ang.3.v2[i,]=ti.thresh(sim.m.ang.3.v2[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.ang.3.v2[i,]=ti.thresh(sim.m.ang.3.v2[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.ang.3.v2[i,]=bayesmooth(sim.m.ang.3.v2[i,],sigma=sigma.ang.3.v2,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.ang.3.v2[i,]=bayesmooth(sim.m.ang.3.v2[i,],sigma=sigma.ang.3.v2,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.ang.3.v2[i,]=ti.thresh(sim.m.ang.3.v2[i,],sigma=sigma.ang.3.v2,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.ang.3.v2[i,]=ti.thresh(sim.m.ang.3.v2[i,],sigma=sigma.ang.3.v2,filter.number=8,family="DaubLeAsymm")

  sigma.dop.1.v2.est=sig.est.func(sim.m.dop.1.v2[i,],n)
  mu.est.ash.haar.dop.1.v2[i,]=bayesmooth(sim.m.dop.1.v2[i,],sigma.dop.1.v2.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.dop.1.v2[i,]=bayesmooth(sim.m.dop.1.v2[i,],sigma.dop.1.v2.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.dop.1.v2[i,]=bayesmooth(sim.m.dop.1.v2[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.dop.1.v2[i,]=bayesmooth(sim.m.dop.1.v2[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.dop.1.v2[i,]=ti.thresh(sim.m.dop.1.v2[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.dop.1.v2[i,]=ti.thresh(sim.m.dop.1.v2[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.dop.1.v2[i,]=bayesmooth(sim.m.dop.1.v2[i,],sigma=sigma.dop.1.v2,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.dop.1.v2[i,]=bayesmooth(sim.m.dop.1.v2[i,],sigma=sigma.dop.1.v2,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.dop.1.v2[i,]=ti.thresh(sim.m.dop.1.v2[i,],sigma=sigma.dop.1.v2,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.dop.1.v2[i,]=ti.thresh(sim.m.dop.1.v2[i,],sigma=sigma.dop.1.v2,filter.number=8,family="DaubLeAsymm")
  sigma.dop.3.v2.est=sig.est.func(sim.m.dop.3.v2[i,],n)
  mu.est.ash.haar.dop.3.v2[i,]=bayesmooth(sim.m.dop.3.v2[i,],sigma.dop.3.v2.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.dop.3.v2[i,]=bayesmooth(sim.m.dop.3.v2[i,],sigma.dop.3.v2.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.dop.3.v2[i,]=bayesmooth(sim.m.dop.3.v2[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.dop.3.v2[i,]=bayesmooth(sim.m.dop.3.v2[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.dop.3.v2[i,]=ti.thresh(sim.m.dop.3.v2[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.dop.3.v2[i,]=ti.thresh(sim.m.dop.3.v2[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.dop.3.v2[i,]=bayesmooth(sim.m.dop.3.v2[i,],sigma=sigma.dop.3.v2,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.dop.3.v2[i,]=bayesmooth(sim.m.dop.3.v2[i,],sigma=sigma.dop.3.v2,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.dop.3.v2[i,]=ti.thresh(sim.m.dop.3.v2[i,],sigma=sigma.dop.3.v2,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.dop.3.v2[i,]=ti.thresh(sim.m.dop.3.v2[i,],sigma=sigma.dop.3.v2,filter.number=8,family="DaubLeAsymm")

  sigma.blip.1.v2.est=sig.est.func(sim.m.blip.1.v2[i,],n)
  mu.est.ash.haar.blip.1.v2[i,]=bayesmooth(sim.m.blip.1.v2[i,],sigma.blip.1.v2.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.blip.1.v2[i,]=bayesmooth(sim.m.blip.1.v2[i,],sigma.blip.1.v2.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.blip.1.v2[i,]=bayesmooth(sim.m.blip.1.v2[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.blip.1.v2[i,]=bayesmooth(sim.m.blip.1.v2[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.blip.1.v2[i,]=ti.thresh(sim.m.blip.1.v2[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.blip.1.v2[i,]=ti.thresh(sim.m.blip.1.v2[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.blip.1.v2[i,]=bayesmooth(sim.m.blip.1.v2[i,],sigma=sigma.blip.1.v2,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.blip.1.v2[i,]=bayesmooth(sim.m.blip.1.v2[i,],sigma=sigma.blip.1.v2,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.blip.1.v2[i,]=ti.thresh(sim.m.blip.1.v2[i,],sigma=sigma.blip.1.v2,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.blip.1.v2[i,]=ti.thresh(sim.m.blip.1.v2[i,],sigma=sigma.blip.1.v2,filter.number=8,family="DaubLeAsymm")
  sigma.blip.3.v2.est=sig.est.func(sim.m.blip.3.v2[i,],n)
  mu.est.ash.haar.blip.3.v2[i,]=bayesmooth(sim.m.blip.3.v2[i,],sigma.blip.3.v2.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.blip.3.v2[i,]=bayesmooth(sim.m.blip.3.v2[i,],sigma.blip.3.v2.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.blip.3.v2[i,]=bayesmooth(sim.m.blip.3.v2[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.blip.3.v2[i,]=bayesmooth(sim.m.blip.3.v2[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.blip.3.v2[i,]=ti.thresh(sim.m.blip.3.v2[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.blip.3.v2[i,]=ti.thresh(sim.m.blip.3.v2[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.blip.3.v2[i,]=bayesmooth(sim.m.blip.3.v2[i,],sigma=sigma.blip.3.v2,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.blip.3.v2[i,]=bayesmooth(sim.m.blip.3.v2[i,],sigma=sigma.blip.3.v2,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.blip.3.v2[i,]=ti.thresh(sim.m.blip.3.v2[i,],sigma=sigma.blip.3.v2,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.blip.3.v2[i,]=ti.thresh(sim.m.blip.3.v2[i,],sigma=sigma.blip.3.v2,filter.number=8,family="DaubLeAsymm")

  sigma.cor.1.v2.est=sig.est.func(sim.m.cor.1.v2[i,],n)
  mu.est.ash.haar.cor.1.v2[i,]=bayesmooth(sim.m.cor.1.v2[i,],sigma.cor.1.v2.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.cor.1.v2[i,]=bayesmooth(sim.m.cor.1.v2[i,],sigma.cor.1.v2.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.cor.1.v2[i,]=bayesmooth(sim.m.cor.1.v2[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.cor.1.v2[i,]=bayesmooth(sim.m.cor.1.v2[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.cor.1.v2[i,]=ti.thresh(sim.m.cor.1.v2[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.cor.1.v2[i,]=ti.thresh(sim.m.cor.1.v2[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.cor.1.v2[i,]=bayesmooth(sim.m.cor.1.v2[i,],sigma=sigma.cor.1.v2,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.cor.1.v2[i,]=bayesmooth(sim.m.cor.1.v2[i,],sigma=sigma.cor.1.v2,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.cor.1.v2[i,]=ti.thresh(sim.m.cor.1.v2[i,],sigma=sigma.cor.1.v2,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.cor.1.v2[i,]=ti.thresh(sim.m.cor.1.v2[i,],sigma=sigma.cor.1.v2,filter.number=8,family="DaubLeAsymm")
  sigma.cor.3.v2.est=sig.est.func(sim.m.cor.3.v2[i,],n)
  mu.est.ash.haar.cor.3.v2[i,]=bayesmooth(sim.m.cor.3.v2[i,],sigma.cor.3.v2.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.cor.3.v2[i,]=bayesmooth(sim.m.cor.3.v2[i,],sigma.cor.3.v2.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.cor.3.v2[i,]=bayesmooth(sim.m.cor.3.v2[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.cor.3.v2[i,]=bayesmooth(sim.m.cor.3.v2[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.cor.3.v2[i,]=ti.thresh(sim.m.cor.3.v2[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.cor.3.v2[i,]=ti.thresh(sim.m.cor.3.v2[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.cor.3.v2[i,]=bayesmooth(sim.m.cor.3.v2[i,],sigma=sigma.cor.3.v2,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.cor.3.v2[i,]=bayesmooth(sim.m.cor.3.v2[i,],sigma=sigma.cor.3.v2,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.cor.3.v2[i,]=ti.thresh(sim.m.cor.3.v2[i,],sigma=sigma.cor.3.v2,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.cor.3.v2[i,]=ti.thresh(sim.m.cor.3.v2[i,],sigma=sigma.cor.3.v2,filter.number=8,family="DaubLeAsymm")




  sigma.sp.1.v3.est=sig.est.func(sim.m.sp.1.v3[i,],n)
  mu.est.ash.haar.sp.1.v3[i,]=bayesmooth(sim.m.sp.1.v3[i,],sigma.sp.1.v3.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.sp.1.v3[i,]=bayesmooth(sim.m.sp.1.v3[i,],sigma.sp.1.v3.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.sp.1.v3[i,]=bayesmooth(sim.m.sp.1.v3[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.sp.1.v3[i,]=bayesmooth(sim.m.sp.1.v3[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.sp.1.v3[i,]=ti.thresh(sim.m.sp.1.v3[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.sp.1.v3[i,]=ti.thresh(sim.m.sp.1.v3[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.sp.1.v3[i,]=bayesmooth(sim.m.sp.1.v3[i,],sigma=sigma.sp.1.v3,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.sp.1.v3[i,]=bayesmooth(sim.m.sp.1.v3[i,],sigma=sigma.sp.1.v3,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.sp.1.v3[i,]=ti.thresh(sim.m.sp.1.v3[i,],sigma=sigma.sp.1.v3,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.sp.1.v3[i,]=ti.thresh(sim.m.sp.1.v3[i,],sigma=sigma.sp.1.v3,filter.number=8,family="DaubLeAsymm")
  sigma.sp.3.v3.est=sig.est.func(sim.m.sp.3.v3[i,],n)
  mu.est.ash.haar.sp.3.v3[i,]=bayesmooth(sim.m.sp.3.v3[i,],sigma.sp.3.v3.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.sp.3.v3[i,]=bayesmooth(sim.m.sp.3.v3[i,],sigma.sp.3.v3.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.sp.3.v3[i,]=bayesmooth(sim.m.sp.3.v3[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.sp.3.v3[i,]=bayesmooth(sim.m.sp.3.v3[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.sp.3.v3[i,]=ti.thresh(sim.m.sp.3.v3[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.sp.3.v3[i,]=ti.thresh(sim.m.sp.3.v3[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.sp.3.v3[i,]=bayesmooth(sim.m.sp.3.v3[i,],sigma=sigma.sp.3.v3,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.sp.3.v3[i,]=bayesmooth(sim.m.sp.3.v3[i,],sigma=sigma.sp.3.v3,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.sp.3.v3[i,]=ti.thresh(sim.m.sp.3.v3[i,],sigma=sigma.sp.3.v3,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.sp.3.v3[i,]=ti.thresh(sim.m.sp.3.v3[i,],sigma=sigma.sp.3.v3,filter.number=8,family="DaubLeAsymm")

  sigma.bump.1.v3.est=sig.est.func(sim.m.bump.1.v3[i,],n)
  mu.est.ash.haar.bump.1.v3[i,]=bayesmooth(sim.m.bump.1.v3[i,],sigma.bump.1.v3.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.bump.1.v3[i,]=bayesmooth(sim.m.bump.1.v3[i,],sigma.bump.1.v3.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.bump.1.v3[i,]=bayesmooth(sim.m.bump.1.v3[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.bump.1.v3[i,]=bayesmooth(sim.m.bump.1.v3[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.bump.1.v3[i,]=ti.thresh(sim.m.bump.1.v3[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.bump.1.v3[i,]=ti.thresh(sim.m.bump.1.v3[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.bump.1.v3[i,]=bayesmooth(sim.m.bump.1.v3[i,],sigma=sigma.bump.1.v3,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.bump.1.v3[i,]=bayesmooth(sim.m.bump.1.v3[i,],sigma=sigma.bump.1.v3,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.bump.1.v3[i,]=ti.thresh(sim.m.bump.1.v3[i,],sigma=sigma.bump.1.v3,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.bump.1.v3[i,]=ti.thresh(sim.m.bump.1.v3[i,],sigma=sigma.bump.1.v3,filter.number=8,family="DaubLeAsymm")
  sigma.bump.3.v3.est=sig.est.func(sim.m.bump.3.v3[i,],n)
  mu.est.ash.haar.bump.3.v3[i,]=bayesmooth(sim.m.bump.3.v3[i,],sigma.bump.3.v3.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.bump.3.v3[i,]=bayesmooth(sim.m.bump.3.v3[i,],sigma.bump.3.v3.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.bump.3.v3[i,]=bayesmooth(sim.m.bump.3.v3[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.bump.3.v3[i,]=bayesmooth(sim.m.bump.3.v3[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.bump.3.v3[i,]=ti.thresh(sim.m.bump.3.v3[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.bump.3.v3[i,]=ti.thresh(sim.m.bump.3.v3[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.bump.3.v3[i,]=bayesmooth(sim.m.bump.3.v3[i,],sigma=sigma.bump.3.v3,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.bump.3.v3[i,]=bayesmooth(sim.m.bump.3.v3[i,],sigma=sigma.bump.3.v3,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.bump.3.v3[i,]=ti.thresh(sim.m.bump.3.v3[i,],sigma=sigma.bump.3.v3,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.bump.3.v3[i,]=ti.thresh(sim.m.bump.3.v3[i,],sigma=sigma.bump.3.v3,filter.number=8,family="DaubLeAsymm")
  
  sigma.blk.1.v3.est=sig.est.func(sim.m.blk.1.v3[i,],n)
  mu.est.ash.haar.blk.1.v3[i,]=bayesmooth(sim.m.blk.1.v3[i,],sigma.blk.1.v3.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.blk.1.v3[i,]=bayesmooth(sim.m.blk.1.v3[i,],sigma.blk.1.v3.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.blk.1.v3[i,]=bayesmooth(sim.m.blk.1.v3[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.blk.1.v3[i,]=bayesmooth(sim.m.blk.1.v3[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.blk.1.v3[i,]=ti.thresh(sim.m.blk.1.v3[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.blk.1.v3[i,]=ti.thresh(sim.m.blk.1.v3[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.blk.1.v3[i,]=bayesmooth(sim.m.blk.1.v3[i,],sigma=sigma.blk.1.v3,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.blk.1.v3[i,]=bayesmooth(sim.m.blk.1.v3[i,],sigma=sigma.blk.1.v3,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.blk.1.v3[i,]=ti.thresh(sim.m.blk.1.v3[i,],sigma=sigma.blk.1.v3,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.blk.1.v3[i,]=ti.thresh(sim.m.blk.1.v3[i,],sigma=sigma.blk.1.v3,filter.number=8,family="DaubLeAsymm")
  sigma.blk.3.v3.est=sig.est.func(sim.m.blk.3.v3[i,],n)
  mu.est.ash.haar.blk.3.v3[i,]=bayesmooth(sim.m.blk.3.v3[i,],sigma.blk.3.v3.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.blk.3.v3[i,]=bayesmooth(sim.m.blk.3.v3[i,],sigma.blk.3.v3.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.blk.3.v3[i,]=bayesmooth(sim.m.blk.3.v3[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.blk.3.v3[i,]=bayesmooth(sim.m.blk.3.v3[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.blk.3.v3[i,]=ti.thresh(sim.m.blk.3.v3[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.blk.3.v3[i,]=ti.thresh(sim.m.blk.3.v3[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.blk.3.v3[i,]=bayesmooth(sim.m.blk.3.v3[i,],sigma=sigma.blk.3.v3,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.blk.3.v3[i,]=bayesmooth(sim.m.blk.3.v3[i,],sigma=sigma.blk.3.v3,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.blk.3.v3[i,]=ti.thresh(sim.m.blk.3.v3[i,],sigma=sigma.blk.3.v3,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.blk.3.v3[i,]=ti.thresh(sim.m.blk.3.v3[i,],sigma=sigma.blk.3.v3,filter.number=8,family="DaubLeAsymm")

  sigma.ang.1.v3.est=sig.est.func(sim.m.ang.1.v3[i,],n)
  mu.est.ash.haar.ang.1.v3[i,]=bayesmooth(sim.m.ang.1.v3[i,],sigma.ang.1.v3.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.ang.1.v3[i,]=bayesmooth(sim.m.ang.1.v3[i,],sigma.ang.1.v3.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.ang.1.v3[i,]=bayesmooth(sim.m.ang.1.v3[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.ang.1.v3[i,]=bayesmooth(sim.m.ang.1.v3[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.ang.1.v3[i,]=ti.thresh(sim.m.ang.1.v3[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.ang.1.v3[i,]=ti.thresh(sim.m.ang.1.v3[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.ang.1.v3[i,]=bayesmooth(sim.m.ang.1.v3[i,],sigma=sigma.ang.1.v3,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.ang.1.v3[i,]=bayesmooth(sim.m.ang.1.v3[i,],sigma=sigma.ang.1.v3,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.ang.1.v3[i,]=ti.thresh(sim.m.ang.1.v3[i,],sigma=sigma.ang.1.v3,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.ang.1.v3[i,]=ti.thresh(sim.m.ang.1.v3[i,],sigma=sigma.ang.1.v3,filter.number=8,family="DaubLeAsymm")
  sigma.ang.3.v3.est=sig.est.func(sim.m.ang.3.v3[i,],n)
  mu.est.ash.haar.ang.3.v3[i,]=bayesmooth(sim.m.ang.3.v3[i,],sigma.ang.3.v3.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.ang.3.v3[i,]=bayesmooth(sim.m.ang.3.v3[i,],sigma.ang.3.v3.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.ang.3.v3[i,]=bayesmooth(sim.m.ang.3.v3[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.ang.3.v3[i,]=bayesmooth(sim.m.ang.3.v3[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.ang.3.v3[i,]=ti.thresh(sim.m.ang.3.v3[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.ang.3.v3[i,]=ti.thresh(sim.m.ang.3.v3[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.ang.3.v3[i,]=bayesmooth(sim.m.ang.3.v3[i,],sigma=sigma.ang.3.v3,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.ang.3.v3[i,]=bayesmooth(sim.m.ang.3.v3[i,],sigma=sigma.ang.3.v3,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.ang.3.v3[i,]=ti.thresh(sim.m.ang.3.v3[i,],sigma=sigma.ang.3.v3,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.ang.3.v3[i,]=ti.thresh(sim.m.ang.3.v3[i,],sigma=sigma.ang.3.v3,filter.number=8,family="DaubLeAsymm")

  sigma.dop.1.v3.est=sig.est.func(sim.m.dop.1.v3[i,],n)
  mu.est.ash.haar.dop.1.v3[i,]=bayesmooth(sim.m.dop.1.v3[i,],sigma.dop.1.v3.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.dop.1.v3[i,]=bayesmooth(sim.m.dop.1.v3[i,],sigma.dop.1.v3.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.dop.1.v3[i,]=bayesmooth(sim.m.dop.1.v3[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.dop.1.v3[i,]=bayesmooth(sim.m.dop.1.v3[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.dop.1.v3[i,]=ti.thresh(sim.m.dop.1.v3[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.dop.1.v3[i,]=ti.thresh(sim.m.dop.1.v3[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.dop.1.v3[i,]=bayesmooth(sim.m.dop.1.v3[i,],sigma=sigma.dop.1.v3,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.dop.1.v3[i,]=bayesmooth(sim.m.dop.1.v3[i,],sigma=sigma.dop.1.v3,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.dop.1.v3[i,]=ti.thresh(sim.m.dop.1.v3[i,],sigma=sigma.dop.1.v3,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.dop.1.v3[i,]=ti.thresh(sim.m.dop.1.v3[i,],sigma=sigma.dop.1.v3,filter.number=8,family="DaubLeAsymm")
  sigma.dop.3.v3.est=sig.est.func(sim.m.dop.3.v3[i,],n)
  mu.est.ash.haar.dop.3.v3[i,]=bayesmooth(sim.m.dop.3.v3[i,],sigma.dop.3.v3.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.dop.3.v3[i,]=bayesmooth(sim.m.dop.3.v3[i,],sigma.dop.3.v3.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.dop.3.v3[i,]=bayesmooth(sim.m.dop.3.v3[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.dop.3.v3[i,]=bayesmooth(sim.m.dop.3.v3[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.dop.3.v3[i,]=ti.thresh(sim.m.dop.3.v3[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.dop.3.v3[i,]=ti.thresh(sim.m.dop.3.v3[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.dop.3.v3[i,]=bayesmooth(sim.m.dop.3.v3[i,],sigma=sigma.dop.3.v3,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.dop.3.v3[i,]=bayesmooth(sim.m.dop.3.v3[i,],sigma=sigma.dop.3.v3,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.dop.3.v3[i,]=ti.thresh(sim.m.dop.3.v3[i,],sigma=sigma.dop.3.v3,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.dop.3.v3[i,]=ti.thresh(sim.m.dop.3.v3[i,],sigma=sigma.dop.3.v3,filter.number=8,family="DaubLeAsymm")

  sigma.blip.1.v3.est=sig.est.func(sim.m.blip.1.v3[i,],n)
  mu.est.ash.haar.blip.1.v3[i,]=bayesmooth(sim.m.blip.1.v3[i,],sigma.blip.1.v3.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.blip.1.v3[i,]=bayesmooth(sim.m.blip.1.v3[i,],sigma.blip.1.v3.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.blip.1.v3[i,]=bayesmooth(sim.m.blip.1.v3[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.blip.1.v3[i,]=bayesmooth(sim.m.blip.1.v3[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.blip.1.v3[i,]=ti.thresh(sim.m.blip.1.v3[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.blip.1.v3[i,]=ti.thresh(sim.m.blip.1.v3[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.blip.1.v3[i,]=bayesmooth(sim.m.blip.1.v3[i,],sigma=sigma.blip.1.v3,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.blip.1.v3[i,]=bayesmooth(sim.m.blip.1.v3[i,],sigma=sigma.blip.1.v3,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.blip.1.v3[i,]=ti.thresh(sim.m.blip.1.v3[i,],sigma=sigma.blip.1.v3,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.blip.1.v3[i,]=ti.thresh(sim.m.blip.1.v3[i,],sigma=sigma.blip.1.v3,filter.number=8,family="DaubLeAsymm")
  sigma.blip.3.v3.est=sig.est.func(sim.m.blip.3.v3[i,],n)
  mu.est.ash.haar.blip.3.v3[i,]=bayesmooth(sim.m.blip.3.v3[i,],sigma.blip.3.v3.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.blip.3.v3[i,]=bayesmooth(sim.m.blip.3.v3[i,],sigma.blip.3.v3.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.blip.3.v3[i,]=bayesmooth(sim.m.blip.3.v3[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.blip.3.v3[i,]=bayesmooth(sim.m.blip.3.v3[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.blip.3.v3[i,]=ti.thresh(sim.m.blip.3.v3[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.blip.3.v3[i,]=ti.thresh(sim.m.blip.3.v3[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.blip.3.v3[i,]=bayesmooth(sim.m.blip.3.v3[i,],sigma=sigma.blip.3.v3,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.blip.3.v3[i,]=bayesmooth(sim.m.blip.3.v3[i,],sigma=sigma.blip.3.v3,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.blip.3.v3[i,]=ti.thresh(sim.m.blip.3.v3[i,],sigma=sigma.blip.3.v3,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.blip.3.v3[i,]=ti.thresh(sim.m.blip.3.v3[i,],sigma=sigma.blip.3.v3,filter.number=8,family="DaubLeAsymm")

  sigma.cor.1.v3.est=sig.est.func(sim.m.cor.1.v3[i,],n)
  mu.est.ash.haar.cor.1.v3[i,]=bayesmooth(sim.m.cor.1.v3[i,],sigma.cor.1.v3.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.cor.1.v3[i,]=bayesmooth(sim.m.cor.1.v3[i,],sigma.cor.1.v3.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.cor.1.v3[i,]=bayesmooth(sim.m.cor.1.v3[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.cor.1.v3[i,]=bayesmooth(sim.m.cor.1.v3[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.cor.1.v3[i,]=ti.thresh(sim.m.cor.1.v3[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.cor.1.v3[i,]=ti.thresh(sim.m.cor.1.v3[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.cor.1.v3[i,]=bayesmooth(sim.m.cor.1.v3[i,],sigma=sigma.cor.1.v3,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.cor.1.v3[i,]=bayesmooth(sim.m.cor.1.v3[i,],sigma=sigma.cor.1.v3,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.cor.1.v3[i,]=ti.thresh(sim.m.cor.1.v3[i,],sigma=sigma.cor.1.v3,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.cor.1.v3[i,]=ti.thresh(sim.m.cor.1.v3[i,],sigma=sigma.cor.1.v3,filter.number=8,family="DaubLeAsymm")
  sigma.cor.3.v3.est=sig.est.func(sim.m.cor.3.v3[i,],n)
  mu.est.ash.haar.cor.3.v3[i,]=bayesmooth(sim.m.cor.3.v3[i,],sigma.cor.3.v3.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.cor.3.v3[i,]=bayesmooth(sim.m.cor.3.v3[i,],sigma.cor.3.v3.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.cor.3.v3[i,]=bayesmooth(sim.m.cor.3.v3[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.cor.3.v3[i,]=bayesmooth(sim.m.cor.3.v3[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.cor.3.v3[i,]=ti.thresh(sim.m.cor.3.v3[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.cor.3.v3[i,]=ti.thresh(sim.m.cor.3.v3[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.cor.3.v3[i,]=bayesmooth(sim.m.cor.3.v3[i,],sigma=sigma.cor.3.v3,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.cor.3.v3[i,]=bayesmooth(sim.m.cor.3.v3[i,],sigma=sigma.cor.3.v3,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.cor.3.v3[i,]=ti.thresh(sim.m.cor.3.v3[i,],sigma=sigma.cor.3.v3,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.cor.3.v3[i,]=ti.thresh(sim.m.cor.3.v3[i,],sigma=sigma.cor.3.v3,filter.number=8,family="DaubLeAsymm")


  sigma.sp.1.v4.est=sig.est.func(sim.m.sp.1.v4[i,],n)
  mu.est.ash.haar.sp.1.v4[i,]=bayesmooth(sim.m.sp.1.v4[i,],sigma.sp.1.v4.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.sp.1.v4[i,]=bayesmooth(sim.m.sp.1.v4[i,],sigma.sp.1.v4.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.sp.1.v4[i,]=bayesmooth(sim.m.sp.1.v4[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.sp.1.v4[i,]=bayesmooth(sim.m.sp.1.v4[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.sp.1.v4[i,]=ti.thresh(sim.m.sp.1.v4[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.sp.1.v4[i,]=ti.thresh(sim.m.sp.1.v4[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.sp.1.v4[i,]=bayesmooth(sim.m.sp.1.v4[i,],sigma=sigma.sp.1.v4,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.sp.1.v4[i,]=bayesmooth(sim.m.sp.1.v4[i,],sigma=sigma.sp.1.v4,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.sp.1.v4[i,]=ti.thresh(sim.m.sp.1.v4[i,],sigma=sigma.sp.1.v4,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.sp.1.v4[i,]=ti.thresh(sim.m.sp.1.v4[i,],sigma=sigma.sp.1.v4,filter.number=8,family="DaubLeAsymm")
  sigma.sp.3.v4.est=sig.est.func(sim.m.sp.3.v4[i,],n)
  mu.est.ash.haar.sp.3.v4[i,]=bayesmooth(sim.m.sp.3.v4[i,],sigma.sp.3.v4.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.sp.3.v4[i,]=bayesmooth(sim.m.sp.3.v4[i,],sigma.sp.3.v4.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.sp.3.v4[i,]=bayesmooth(sim.m.sp.3.v4[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.sp.3.v4[i,]=bayesmooth(sim.m.sp.3.v4[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.sp.3.v4[i,]=ti.thresh(sim.m.sp.3.v4[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.sp.3.v4[i,]=ti.thresh(sim.m.sp.3.v4[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.sp.3.v4[i,]=bayesmooth(sim.m.sp.3.v4[i,],sigma=sigma.sp.3.v4,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.sp.3.v4[i,]=bayesmooth(sim.m.sp.3.v4[i,],sigma=sigma.sp.3.v4,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.sp.3.v4[i,]=ti.thresh(sim.m.sp.3.v4[i,],sigma=sigma.sp.3.v4,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.sp.3.v4[i,]=ti.thresh(sim.m.sp.3.v4[i,],sigma=sigma.sp.3.v4,filter.number=8,family="DaubLeAsymm")

  sigma.bump.1.v4.est=sig.est.func(sim.m.bump.1.v4[i,],n)
  mu.est.ash.haar.bump.1.v4[i,]=bayesmooth(sim.m.bump.1.v4[i,],sigma.bump.1.v4.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.bump.1.v4[i,]=bayesmooth(sim.m.bump.1.v4[i,],sigma.bump.1.v4.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.bump.1.v4[i,]=bayesmooth(sim.m.bump.1.v4[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.bump.1.v4[i,]=bayesmooth(sim.m.bump.1.v4[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.bump.1.v4[i,]=ti.thresh(sim.m.bump.1.v4[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.bump.1.v4[i,]=ti.thresh(sim.m.bump.1.v4[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.bump.1.v4[i,]=bayesmooth(sim.m.bump.1.v4[i,],sigma=sigma.bump.1.v4,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.bump.1.v4[i,]=bayesmooth(sim.m.bump.1.v4[i,],sigma=sigma.bump.1.v4,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.bump.1.v4[i,]=ti.thresh(sim.m.bump.1.v4[i,],sigma=sigma.bump.1.v4,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.bump.1.v4[i,]=ti.thresh(sim.m.bump.1.v4[i,],sigma=sigma.bump.1.v4,filter.number=8,family="DaubLeAsymm")
  sigma.bump.3.v4.est=sig.est.func(sim.m.bump.3.v4[i,],n)
  mu.est.ash.haar.bump.3.v4[i,]=bayesmooth(sim.m.bump.3.v4[i,],sigma.bump.3.v4.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.bump.3.v4[i,]=bayesmooth(sim.m.bump.3.v4[i,],sigma.bump.3.v4.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.bump.3.v4[i,]=bayesmooth(sim.m.bump.3.v4[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.bump.3.v4[i,]=bayesmooth(sim.m.bump.3.v4[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.bump.3.v4[i,]=ti.thresh(sim.m.bump.3.v4[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.bump.3.v4[i,]=ti.thresh(sim.m.bump.3.v4[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.bump.3.v4[i,]=bayesmooth(sim.m.bump.3.v4[i,],sigma=sigma.bump.3.v4,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.bump.3.v4[i,]=bayesmooth(sim.m.bump.3.v4[i,],sigma=sigma.bump.3.v4,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.bump.3.v4[i,]=ti.thresh(sim.m.bump.3.v4[i,],sigma=sigma.bump.3.v4,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.bump.3.v4[i,]=ti.thresh(sim.m.bump.3.v4[i,],sigma=sigma.bump.3.v4,filter.number=8,family="DaubLeAsymm")
  
  sigma.blk.1.v4.est=sig.est.func(sim.m.blk.1.v4[i,],n)
  mu.est.ash.haar.blk.1.v4[i,]=bayesmooth(sim.m.blk.1.v4[i,],sigma.blk.1.v4.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.blk.1.v4[i,]=bayesmooth(sim.m.blk.1.v4[i,],sigma.blk.1.v4.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.blk.1.v4[i,]=bayesmooth(sim.m.blk.1.v4[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.blk.1.v4[i,]=bayesmooth(sim.m.blk.1.v4[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.blk.1.v4[i,]=ti.thresh(sim.m.blk.1.v4[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.blk.1.v4[i,]=ti.thresh(sim.m.blk.1.v4[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.blk.1.v4[i,]=bayesmooth(sim.m.blk.1.v4[i,],sigma=sigma.blk.1.v4,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.blk.1.v4[i,]=bayesmooth(sim.m.blk.1.v4[i,],sigma=sigma.blk.1.v4,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.blk.1.v4[i,]=ti.thresh(sim.m.blk.1.v4[i,],sigma=sigma.blk.1.v4,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.blk.1.v4[i,]=ti.thresh(sim.m.blk.1.v4[i,],sigma=sigma.blk.1.v4,filter.number=8,family="DaubLeAsymm")
  sigma.blk.3.v4.est=sig.est.func(sim.m.blk.3.v4[i,],n)
  mu.est.ash.haar.blk.3.v4[i,]=bayesmooth(sim.m.blk.3.v4[i,],sigma.blk.3.v4.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.blk.3.v4[i,]=bayesmooth(sim.m.blk.3.v4[i,],sigma.blk.3.v4.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.blk.3.v4[i,]=bayesmooth(sim.m.blk.3.v4[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.blk.3.v4[i,]=bayesmooth(sim.m.blk.3.v4[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.blk.3.v4[i,]=ti.thresh(sim.m.blk.3.v4[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.blk.3.v4[i,]=ti.thresh(sim.m.blk.3.v4[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.blk.3.v4[i,]=bayesmooth(sim.m.blk.3.v4[i,],sigma=sigma.blk.3.v4,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.blk.3.v4[i,]=bayesmooth(sim.m.blk.3.v4[i,],sigma=sigma.blk.3.v4,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.blk.3.v4[i,]=ti.thresh(sim.m.blk.3.v4[i,],sigma=sigma.blk.3.v4,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.blk.3.v4[i,]=ti.thresh(sim.m.blk.3.v4[i,],sigma=sigma.blk.3.v4,filter.number=8,family="DaubLeAsymm")

  sigma.ang.1.v4.est=sig.est.func(sim.m.ang.1.v4[i,],n)
  mu.est.ash.haar.ang.1.v4[i,]=bayesmooth(sim.m.ang.1.v4[i,],sigma.ang.1.v4.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.ang.1.v4[i,]=bayesmooth(sim.m.ang.1.v4[i,],sigma.ang.1.v4.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.ang.1.v4[i,]=bayesmooth(sim.m.ang.1.v4[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.ang.1.v4[i,]=bayesmooth(sim.m.ang.1.v4[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.ang.1.v4[i,]=ti.thresh(sim.m.ang.1.v4[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.ang.1.v4[i,]=ti.thresh(sim.m.ang.1.v4[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.ang.1.v4[i,]=bayesmooth(sim.m.ang.1.v4[i,],sigma=sigma.ang.1.v4,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.ang.1.v4[i,]=bayesmooth(sim.m.ang.1.v4[i,],sigma=sigma.ang.1.v4,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.ang.1.v4[i,]=ti.thresh(sim.m.ang.1.v4[i,],sigma=sigma.ang.1.v4,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.ang.1.v4[i,]=ti.thresh(sim.m.ang.1.v4[i,],sigma=sigma.ang.1.v4,filter.number=8,family="DaubLeAsymm")
  sigma.ang.3.v4.est=sig.est.func(sim.m.ang.3.v4[i,],n)
  mu.est.ash.haar.ang.3.v4[i,]=bayesmooth(sim.m.ang.3.v4[i,],sigma.ang.3.v4.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.ang.3.v4[i,]=bayesmooth(sim.m.ang.3.v4[i,],sigma.ang.3.v4.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.ang.3.v4[i,]=bayesmooth(sim.m.ang.3.v4[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.ang.3.v4[i,]=bayesmooth(sim.m.ang.3.v4[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.ang.3.v4[i,]=ti.thresh(sim.m.ang.3.v4[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.ang.3.v4[i,]=ti.thresh(sim.m.ang.3.v4[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.ang.3.v4[i,]=bayesmooth(sim.m.ang.3.v4[i,],sigma=sigma.ang.3.v4,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.ang.3.v4[i,]=bayesmooth(sim.m.ang.3.v4[i,],sigma=sigma.ang.3.v4,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.ang.3.v4[i,]=ti.thresh(sim.m.ang.3.v4[i,],sigma=sigma.ang.3.v4,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.ang.3.v4[i,]=ti.thresh(sim.m.ang.3.v4[i,],sigma=sigma.ang.3.v4,filter.number=8,family="DaubLeAsymm")

  sigma.dop.1.v4.est=sig.est.func(sim.m.dop.1.v4[i,],n)
  mu.est.ash.haar.dop.1.v4[i,]=bayesmooth(sim.m.dop.1.v4[i,],sigma.dop.1.v4.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.dop.1.v4[i,]=bayesmooth(sim.m.dop.1.v4[i,],sigma.dop.1.v4.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.dop.1.v4[i,]=bayesmooth(sim.m.dop.1.v4[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.dop.1.v4[i,]=bayesmooth(sim.m.dop.1.v4[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.dop.1.v4[i,]=ti.thresh(sim.m.dop.1.v4[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.dop.1.v4[i,]=ti.thresh(sim.m.dop.1.v4[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.dop.1.v4[i,]=bayesmooth(sim.m.dop.1.v4[i,],sigma=sigma.dop.1.v4,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.dop.1.v4[i,]=bayesmooth(sim.m.dop.1.v4[i,],sigma=sigma.dop.1.v4,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.dop.1.v4[i,]=ti.thresh(sim.m.dop.1.v4[i,],sigma=sigma.dop.1.v4,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.dop.1.v4[i,]=ti.thresh(sim.m.dop.1.v4[i,],sigma=sigma.dop.1.v4,filter.number=8,family="DaubLeAsymm")
  sigma.dop.3.v4.est=sig.est.func(sim.m.dop.3.v4[i,],n)
  mu.est.ash.haar.dop.3.v4[i,]=bayesmooth(sim.m.dop.3.v4[i,],sigma.dop.3.v4.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.dop.3.v4[i,]=bayesmooth(sim.m.dop.3.v4[i,],sigma.dop.3.v4.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.dop.3.v4[i,]=bayesmooth(sim.m.dop.3.v4[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.dop.3.v4[i,]=bayesmooth(sim.m.dop.3.v4[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.dop.3.v4[i,]=ti.thresh(sim.m.dop.3.v4[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.dop.3.v4[i,]=ti.thresh(sim.m.dop.3.v4[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.dop.3.v4[i,]=bayesmooth(sim.m.dop.3.v4[i,],sigma=sigma.dop.3.v4,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.dop.3.v4[i,]=bayesmooth(sim.m.dop.3.v4[i,],sigma=sigma.dop.3.v4,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.dop.3.v4[i,]=ti.thresh(sim.m.dop.3.v4[i,],sigma=sigma.dop.3.v4,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.dop.3.v4[i,]=ti.thresh(sim.m.dop.3.v4[i,],sigma=sigma.dop.3.v4,filter.number=8,family="DaubLeAsymm")

  sigma.blip.1.v4.est=sig.est.func(sim.m.blip.1.v4[i,],n)
  mu.est.ash.haar.blip.1.v4[i,]=bayesmooth(sim.m.blip.1.v4[i,],sigma.blip.1.v4.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.blip.1.v4[i,]=bayesmooth(sim.m.blip.1.v4[i,],sigma.blip.1.v4.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.blip.1.v4[i,]=bayesmooth(sim.m.blip.1.v4[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.blip.1.v4[i,]=bayesmooth(sim.m.blip.1.v4[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.blip.1.v4[i,]=ti.thresh(sim.m.blip.1.v4[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.blip.1.v4[i,]=ti.thresh(sim.m.blip.1.v4[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.blip.1.v4[i,]=bayesmooth(sim.m.blip.1.v4[i,],sigma=sigma.blip.1.v4,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.blip.1.v4[i,]=bayesmooth(sim.m.blip.1.v4[i,],sigma=sigma.blip.1.v4,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.blip.1.v4[i,]=ti.thresh(sim.m.blip.1.v4[i,],sigma=sigma.blip.1.v4,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.blip.1.v4[i,]=ti.thresh(sim.m.blip.1.v4[i,],sigma=sigma.blip.1.v4,filter.number=8,family="DaubLeAsymm")
  sigma.blip.3.v4.est=sig.est.func(sim.m.blip.3.v4[i,],n)
  mu.est.ash.haar.blip.3.v4[i,]=bayesmooth(sim.m.blip.3.v4[i,],sigma.blip.3.v4.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.blip.3.v4[i,]=bayesmooth(sim.m.blip.3.v4[i,],sigma.blip.3.v4.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.blip.3.v4[i,]=bayesmooth(sim.m.blip.3.v4[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.blip.3.v4[i,]=bayesmooth(sim.m.blip.3.v4[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.blip.3.v4[i,]=ti.thresh(sim.m.blip.3.v4[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.blip.3.v4[i,]=ti.thresh(sim.m.blip.3.v4[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.blip.3.v4[i,]=bayesmooth(sim.m.blip.3.v4[i,],sigma=sigma.blip.3.v4,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.blip.3.v4[i,]=bayesmooth(sim.m.blip.3.v4[i,],sigma=sigma.blip.3.v4,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.blip.3.v4[i,]=ti.thresh(sim.m.blip.3.v4[i,],sigma=sigma.blip.3.v4,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.blip.3.v4[i,]=ti.thresh(sim.m.blip.3.v4[i,],sigma=sigma.blip.3.v4,filter.number=8,family="DaubLeAsymm")

  sigma.cor.1.v4.est=sig.est.func(sim.m.cor.1.v4[i,],n)
  mu.est.ash.haar.cor.1.v4[i,]=bayesmooth(sim.m.cor.1.v4[i,],sigma.cor.1.v4.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.cor.1.v4[i,]=bayesmooth(sim.m.cor.1.v4[i,],sigma.cor.1.v4.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.cor.1.v4[i,]=bayesmooth(sim.m.cor.1.v4[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.cor.1.v4[i,]=bayesmooth(sim.m.cor.1.v4[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.cor.1.v4[i,]=ti.thresh(sim.m.cor.1.v4[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.cor.1.v4[i,]=ti.thresh(sim.m.cor.1.v4[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.cor.1.v4[i,]=bayesmooth(sim.m.cor.1.v4[i,],sigma=sigma.cor.1.v4,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.cor.1.v4[i,]=bayesmooth(sim.m.cor.1.v4[i,],sigma=sigma.cor.1.v4,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.cor.1.v4[i,]=ti.thresh(sim.m.cor.1.v4[i,],sigma=sigma.cor.1.v4,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.cor.1.v4[i,]=ti.thresh(sim.m.cor.1.v4[i,],sigma=sigma.cor.1.v4,filter.number=8,family="DaubLeAsymm")
  sigma.cor.3.v4.est=sig.est.func(sim.m.cor.3.v4[i,],n)
  mu.est.ash.haar.cor.3.v4[i,]=bayesmooth(sim.m.cor.3.v4[i,],sigma.cor.3.v4.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.cor.3.v4[i,]=bayesmooth(sim.m.cor.3.v4[i,],sigma.cor.3.v4.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.cor.3.v4[i,]=bayesmooth(sim.m.cor.3.v4[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.cor.3.v4[i,]=bayesmooth(sim.m.cor.3.v4[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.cor.3.v4[i,]=ti.thresh(sim.m.cor.3.v4[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.cor.3.v4[i,]=ti.thresh(sim.m.cor.3.v4[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.cor.3.v4[i,]=bayesmooth(sim.m.cor.3.v4[i,],sigma=sigma.cor.3.v4,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.cor.3.v4[i,]=bayesmooth(sim.m.cor.3.v4[i,],sigma=sigma.cor.3.v4,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.cor.3.v4[i,]=ti.thresh(sim.m.cor.3.v4[i,],sigma=sigma.cor.3.v4,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.cor.3.v4[i,]=ti.thresh(sim.m.cor.3.v4[i,],sigma=sigma.cor.3.v4,filter.number=8,family="DaubLeAsymm")


  sigma.sp.1.v5.est=sig.est.func(sim.m.sp.1.v5[i,],n)
  mu.est.ash.haar.sp.1.v5[i,]=bayesmooth(sim.m.sp.1.v5[i,],sigma.sp.1.v5.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.sp.1.v5[i,]=bayesmooth(sim.m.sp.1.v5[i,],sigma.sp.1.v5.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.sp.1.v5[i,]=bayesmooth(sim.m.sp.1.v5[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.sp.1.v5[i,]=bayesmooth(sim.m.sp.1.v5[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.sp.1.v5[i,]=ti.thresh(sim.m.sp.1.v5[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.sp.1.v5[i,]=ti.thresh(sim.m.sp.1.v5[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.sp.1.v5[i,]=bayesmooth(sim.m.sp.1.v5[i,],sigma=sigma.sp.1.v5,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.sp.1.v5[i,]=bayesmooth(sim.m.sp.1.v5[i,],sigma=sigma.sp.1.v5,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.sp.1.v5[i,]=ti.thresh(sim.m.sp.1.v5[i,],sigma=sigma.sp.1.v5,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.sp.1.v5[i,]=ti.thresh(sim.m.sp.1.v5[i,],sigma=sigma.sp.1.v5,filter.number=8,family="DaubLeAsymm")
  sigma.sp.3.v5.est=sig.est.func(sim.m.sp.3.v5[i,],n)
  mu.est.ash.haar.sp.3.v5[i,]=bayesmooth(sim.m.sp.3.v5[i,],sigma.sp.3.v5.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.sp.3.v5[i,]=bayesmooth(sim.m.sp.3.v5[i,],sigma.sp.3.v5.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.sp.3.v5[i,]=bayesmooth(sim.m.sp.3.v5[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.sp.3.v5[i,]=bayesmooth(sim.m.sp.3.v5[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.sp.3.v5[i,]=ti.thresh(sim.m.sp.3.v5[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.sp.3.v5[i,]=ti.thresh(sim.m.sp.3.v5[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.sp.3.v5[i,]=bayesmooth(sim.m.sp.3.v5[i,],sigma=sigma.sp.3.v5,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.sp.3.v5[i,]=bayesmooth(sim.m.sp.3.v5[i,],sigma=sigma.sp.3.v5,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.sp.3.v5[i,]=ti.thresh(sim.m.sp.3.v5[i,],sigma=sigma.sp.3.v5,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.sp.3.v5[i,]=ti.thresh(sim.m.sp.3.v5[i,],sigma=sigma.sp.3.v5,filter.number=8,family="DaubLeAsymm")

  sigma.bump.1.v5.est=sig.est.func(sim.m.bump.1.v5[i,],n)
  mu.est.ash.haar.bump.1.v5[i,]=bayesmooth(sim.m.bump.1.v5[i,],sigma.bump.1.v5.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.bump.1.v5[i,]=bayesmooth(sim.m.bump.1.v5[i,],sigma.bump.1.v5.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.bump.1.v5[i,]=bayesmooth(sim.m.bump.1.v5[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.bump.1.v5[i,]=bayesmooth(sim.m.bump.1.v5[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.bump.1.v5[i,]=ti.thresh(sim.m.bump.1.v5[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.bump.1.v5[i,]=ti.thresh(sim.m.bump.1.v5[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.bump.1.v5[i,]=bayesmooth(sim.m.bump.1.v5[i,],sigma=sigma.bump.1.v5,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.bump.1.v5[i,]=bayesmooth(sim.m.bump.1.v5[i,],sigma=sigma.bump.1.v5,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.bump.1.v5[i,]=ti.thresh(sim.m.bump.1.v5[i,],sigma=sigma.bump.1.v5,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.bump.1.v5[i,]=ti.thresh(sim.m.bump.1.v5[i,],sigma=sigma.bump.1.v5,filter.number=8,family="DaubLeAsymm")
  sigma.bump.3.v5.est=sig.est.func(sim.m.bump.3.v5[i,],n)
  mu.est.ash.haar.bump.3.v5[i,]=bayesmooth(sim.m.bump.3.v5[i,],sigma.bump.3.v5.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.bump.3.v5[i,]=bayesmooth(sim.m.bump.3.v5[i,],sigma.bump.3.v5.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.bump.3.v5[i,]=bayesmooth(sim.m.bump.3.v5[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.bump.3.v5[i,]=bayesmooth(sim.m.bump.3.v5[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.bump.3.v5[i,]=ti.thresh(sim.m.bump.3.v5[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.bump.3.v5[i,]=ti.thresh(sim.m.bump.3.v5[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.bump.3.v5[i,]=bayesmooth(sim.m.bump.3.v5[i,],sigma=sigma.bump.3.v5,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.bump.3.v5[i,]=bayesmooth(sim.m.bump.3.v5[i,],sigma=sigma.bump.3.v5,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.bump.3.v5[i,]=ti.thresh(sim.m.bump.3.v5[i,],sigma=sigma.bump.3.v5,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.bump.3.v5[i,]=ti.thresh(sim.m.bump.3.v5[i,],sigma=sigma.bump.3.v5,filter.number=8,family="DaubLeAsymm")
  
  sigma.blk.1.v5.est=sig.est.func(sim.m.blk.1.v5[i,],n)
  mu.est.ash.haar.blk.1.v5[i,]=bayesmooth(sim.m.blk.1.v5[i,],sigma.blk.1.v5.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.blk.1.v5[i,]=bayesmooth(sim.m.blk.1.v5[i,],sigma.blk.1.v5.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.blk.1.v5[i,]=bayesmooth(sim.m.blk.1.v5[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.blk.1.v5[i,]=bayesmooth(sim.m.blk.1.v5[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.blk.1.v5[i,]=ti.thresh(sim.m.blk.1.v5[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.blk.1.v5[i,]=ti.thresh(sim.m.blk.1.v5[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.blk.1.v5[i,]=bayesmooth(sim.m.blk.1.v5[i,],sigma=sigma.blk.1.v5,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.blk.1.v5[i,]=bayesmooth(sim.m.blk.1.v5[i,],sigma=sigma.blk.1.v5,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.blk.1.v5[i,]=ti.thresh(sim.m.blk.1.v5[i,],sigma=sigma.blk.1.v5,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.blk.1.v5[i,]=ti.thresh(sim.m.blk.1.v5[i,],sigma=sigma.blk.1.v5,filter.number=8,family="DaubLeAsymm")
  sigma.blk.3.v5.est=sig.est.func(sim.m.blk.3.v5[i,],n)
  mu.est.ash.haar.blk.3.v5[i,]=bayesmooth(sim.m.blk.3.v5[i,],sigma.blk.3.v5.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.blk.3.v5[i,]=bayesmooth(sim.m.blk.3.v5[i,],sigma.blk.3.v5.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.blk.3.v5[i,]=bayesmooth(sim.m.blk.3.v5[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.blk.3.v5[i,]=bayesmooth(sim.m.blk.3.v5[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.blk.3.v5[i,]=ti.thresh(sim.m.blk.3.v5[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.blk.3.v5[i,]=ti.thresh(sim.m.blk.3.v5[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.blk.3.v5[i,]=bayesmooth(sim.m.blk.3.v5[i,],sigma=sigma.blk.3.v5,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.blk.3.v5[i,]=bayesmooth(sim.m.blk.3.v5[i,],sigma=sigma.blk.3.v5,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.blk.3.v5[i,]=ti.thresh(sim.m.blk.3.v5[i,],sigma=sigma.blk.3.v5,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.blk.3.v5[i,]=ti.thresh(sim.m.blk.3.v5[i,],sigma=sigma.blk.3.v5,filter.number=8,family="DaubLeAsymm")

  sigma.ang.1.v5.est=sig.est.func(sim.m.ang.1.v5[i,],n)
  mu.est.ash.haar.ang.1.v5[i,]=bayesmooth(sim.m.ang.1.v5[i,],sigma.ang.1.v5.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.ang.1.v5[i,]=bayesmooth(sim.m.ang.1.v5[i,],sigma.ang.1.v5.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.ang.1.v5[i,]=bayesmooth(sim.m.ang.1.v5[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.ang.1.v5[i,]=bayesmooth(sim.m.ang.1.v5[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.ang.1.v5[i,]=ti.thresh(sim.m.ang.1.v5[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.ang.1.v5[i,]=ti.thresh(sim.m.ang.1.v5[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.ang.1.v5[i,]=bayesmooth(sim.m.ang.1.v5[i,],sigma=sigma.ang.1.v5,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.ang.1.v5[i,]=bayesmooth(sim.m.ang.1.v5[i,],sigma=sigma.ang.1.v5,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.ang.1.v5[i,]=ti.thresh(sim.m.ang.1.v5[i,],sigma=sigma.ang.1.v5,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.ang.1.v5[i,]=ti.thresh(sim.m.ang.1.v5[i,],sigma=sigma.ang.1.v5,filter.number=8,family="DaubLeAsymm")
  sigma.ang.3.v5.est=sig.est.func(sim.m.ang.3.v5[i,],n)
  mu.est.ash.haar.ang.3.v5[i,]=bayesmooth(sim.m.ang.3.v5[i,],sigma.ang.3.v5.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.ang.3.v5[i,]=bayesmooth(sim.m.ang.3.v5[i,],sigma.ang.3.v5.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.ang.3.v5[i,]=bayesmooth(sim.m.ang.3.v5[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.ang.3.v5[i,]=bayesmooth(sim.m.ang.3.v5[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.ang.3.v5[i,]=ti.thresh(sim.m.ang.3.v5[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.ang.3.v5[i,]=ti.thresh(sim.m.ang.3.v5[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.ang.3.v5[i,]=bayesmooth(sim.m.ang.3.v5[i,],sigma=sigma.ang.3.v5,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.ang.3.v5[i,]=bayesmooth(sim.m.ang.3.v5[i,],sigma=sigma.ang.3.v5,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.ang.3.v5[i,]=ti.thresh(sim.m.ang.3.v5[i,],sigma=sigma.ang.3.v5,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.ang.3.v5[i,]=ti.thresh(sim.m.ang.3.v5[i,],sigma=sigma.ang.3.v5,filter.number=8,family="DaubLeAsymm")

  sigma.dop.1.v5.est=sig.est.func(sim.m.dop.1.v5[i,],n)
  mu.est.ash.haar.dop.1.v5[i,]=bayesmooth(sim.m.dop.1.v5[i,],sigma.dop.1.v5.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.dop.1.v5[i,]=bayesmooth(sim.m.dop.1.v5[i,],sigma.dop.1.v5.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.dop.1.v5[i,]=bayesmooth(sim.m.dop.1.v5[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.dop.1.v5[i,]=bayesmooth(sim.m.dop.1.v5[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.dop.1.v5[i,]=ti.thresh(sim.m.dop.1.v5[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.dop.1.v5[i,]=ti.thresh(sim.m.dop.1.v5[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.dop.1.v5[i,]=bayesmooth(sim.m.dop.1.v5[i,],sigma=sigma.dop.1.v5,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.dop.1.v5[i,]=bayesmooth(sim.m.dop.1.v5[i,],sigma=sigma.dop.1.v5,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.dop.1.v5[i,]=ti.thresh(sim.m.dop.1.v5[i,],sigma=sigma.dop.1.v5,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.dop.1.v5[i,]=ti.thresh(sim.m.dop.1.v5[i,],sigma=sigma.dop.1.v5,filter.number=8,family="DaubLeAsymm")
  sigma.dop.3.v5.est=sig.est.func(sim.m.dop.3.v5[i,],n)
  mu.est.ash.haar.dop.3.v5[i,]=bayesmooth(sim.m.dop.3.v5[i,],sigma.dop.3.v5.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.dop.3.v5[i,]=bayesmooth(sim.m.dop.3.v5[i,],sigma.dop.3.v5.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.dop.3.v5[i,]=bayesmooth(sim.m.dop.3.v5[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.dop.3.v5[i,]=bayesmooth(sim.m.dop.3.v5[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.dop.3.v5[i,]=ti.thresh(sim.m.dop.3.v5[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.dop.3.v5[i,]=ti.thresh(sim.m.dop.3.v5[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.dop.3.v5[i,]=bayesmooth(sim.m.dop.3.v5[i,],sigma=sigma.dop.3.v5,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.dop.3.v5[i,]=bayesmooth(sim.m.dop.3.v5[i,],sigma=sigma.dop.3.v5,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.dop.3.v5[i,]=ti.thresh(sim.m.dop.3.v5[i,],sigma=sigma.dop.3.v5,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.dop.3.v5[i,]=ti.thresh(sim.m.dop.3.v5[i,],sigma=sigma.dop.3.v5,filter.number=8,family="DaubLeAsymm")

  sigma.blip.1.v5.est=sig.est.func(sim.m.blip.1.v5[i,],n)
  mu.est.ash.haar.blip.1.v5[i,]=bayesmooth(sim.m.blip.1.v5[i,],sigma.blip.1.v5.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.blip.1.v5[i,]=bayesmooth(sim.m.blip.1.v5[i,],sigma.blip.1.v5.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.blip.1.v5[i,]=bayesmooth(sim.m.blip.1.v5[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.blip.1.v5[i,]=bayesmooth(sim.m.blip.1.v5[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.blip.1.v5[i,]=ti.thresh(sim.m.blip.1.v5[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.blip.1.v5[i,]=ti.thresh(sim.m.blip.1.v5[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.blip.1.v5[i,]=bayesmooth(sim.m.blip.1.v5[i,],sigma=sigma.blip.1.v5,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.blip.1.v5[i,]=bayesmooth(sim.m.blip.1.v5[i,],sigma=sigma.blip.1.v5,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.blip.1.v5[i,]=ti.thresh(sim.m.blip.1.v5[i,],sigma=sigma.blip.1.v5,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.blip.1.v5[i,]=ti.thresh(sim.m.blip.1.v5[i,],sigma=sigma.blip.1.v5,filter.number=8,family="DaubLeAsymm")
  sigma.blip.3.v5.est=sig.est.func(sim.m.blip.3.v5[i,],n)
  mu.est.ash.haar.blip.3.v5[i,]=bayesmooth(sim.m.blip.3.v5[i,],sigma.blip.3.v5.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.blip.3.v5[i,]=bayesmooth(sim.m.blip.3.v5[i,],sigma.blip.3.v5.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.blip.3.v5[i,]=bayesmooth(sim.m.blip.3.v5[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.blip.3.v5[i,]=bayesmooth(sim.m.blip.3.v5[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.blip.3.v5[i,]=ti.thresh(sim.m.blip.3.v5[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.blip.3.v5[i,]=ti.thresh(sim.m.blip.3.v5[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.blip.3.v5[i,]=bayesmooth(sim.m.blip.3.v5[i,],sigma=sigma.blip.3.v5,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.blip.3.v5[i,]=bayesmooth(sim.m.blip.3.v5[i,],sigma=sigma.blip.3.v5,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.blip.3.v5[i,]=ti.thresh(sim.m.blip.3.v5[i,],sigma=sigma.blip.3.v5,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.blip.3.v5[i,]=ti.thresh(sim.m.blip.3.v5[i,],sigma=sigma.blip.3.v5,filter.number=8,family="DaubLeAsymm")

  sigma.cor.1.v5.est=sig.est.func(sim.m.cor.1.v5[i,],n)
  mu.est.ash.haar.cor.1.v5[i,]=bayesmooth(sim.m.cor.1.v5[i,],sigma.cor.1.v5.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.cor.1.v5[i,]=bayesmooth(sim.m.cor.1.v5[i,],sigma.cor.1.v5.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.cor.1.v5[i,]=bayesmooth(sim.m.cor.1.v5[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.cor.1.v5[i,]=bayesmooth(sim.m.cor.1.v5[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.cor.1.v5[i,]=ti.thresh(sim.m.cor.1.v5[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.cor.1.v5[i,]=ti.thresh(sim.m.cor.1.v5[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.cor.1.v5[i,]=bayesmooth(sim.m.cor.1.v5[i,],sigma=sigma.cor.1.v5,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.cor.1.v5[i,]=bayesmooth(sim.m.cor.1.v5[i,],sigma=sigma.cor.1.v5,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.cor.1.v5[i,]=ti.thresh(sim.m.cor.1.v5[i,],sigma=sigma.cor.1.v5,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.cor.1.v5[i,]=ti.thresh(sim.m.cor.1.v5[i,],sigma=sigma.cor.1.v5,filter.number=8,family="DaubLeAsymm")
  sigma.cor.3.v5.est=sig.est.func(sim.m.cor.3.v5[i,],n)
  mu.est.ash.haar.cor.3.v5[i,]=bayesmooth(sim.m.cor.3.v5[i,],sigma.cor.3.v5.est,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ash.s8.cor.3.v5[i,]=bayesmooth(sim.m.cor.3.v5[i,],sigma.cor.3.v5.est,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.ashe.haar.cor.3.v5[i,]=bayesmooth(sim.m.cor.3.v5[i,],filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.ashe.s8.cor.3.v5[i,]=bayesmooth(sim.m.cor.3.v5[i,],filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tie.haar.cor.3.v5[i,]=ti.thresh(sim.m.cor.3.v5[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.cor.3.v5[i,]=ti.thresh(sim.m.cor.3.v5[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.cor.3.v5[i,]=bayesmooth(sim.m.cor.3.v5[i,],sigma=sigma.cor.3.v5,filter.number=1,family="DaubExPhase",gridmult=0)
  mu.est.asht.s8.cor.3.v5[i,]=bayesmooth(sim.m.cor.3.v5[i,],sigma=sigma.cor.3.v5,filter.number=8,family="DaubLeAsymm",gridmult=0)
  mu.est.tit.haar.cor.3.v5[i,]=ti.thresh(sim.m.cor.3.v5[i,],sigma=sigma.cor.3.v5,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.cor.3.v5[i,]=ti.thresh(sim.m.cor.3.v5[i,],sigma=sigma.cor.3.v5,filter.number=8,family="DaubLeAsymm")


  mu.est.ebayes.sp.1.v1[i,]=waveti.ebayes(sim.m.sp.1.v1[i,],noise.level=sigma.sp.1.v1.est)
  mu.est.ebayes.sp.1.v2[i,]=waveti.ebayes(sim.m.sp.1.v2[i,],noise.level=sigma.sp.1.v2.est)
  mu.est.ebayes.sp.1.v3[i,]=waveti.ebayes(sim.m.sp.1.v3[i,],noise.level=sigma.sp.1.v3.est)
  mu.est.ebayes.sp.1.v4[i,]=waveti.ebayes(sim.m.sp.1.v4[i,],noise.level=sigma.sp.1.v4.est)
  mu.est.ebayes.sp.1.v5[i,]=waveti.ebayes(sim.m.sp.1.v5[i,],noise.level=sigma.sp.1.v5.est)
  mu.est.ebayes.sp.3.v1[i,]=waveti.ebayes(sim.m.sp.3.v1[i,],noise.level=sigma.sp.3.v1.est)
  mu.est.ebayes.sp.3.v2[i,]=waveti.ebayes(sim.m.sp.3.v2[i,],noise.level=sigma.sp.3.v2.est)
  mu.est.ebayes.sp.3.v3[i,]=waveti.ebayes(sim.m.sp.3.v3[i,],noise.level=sigma.sp.3.v3.est)
  mu.est.ebayes.sp.3.v4[i,]=waveti.ebayes(sim.m.sp.3.v4[i,],noise.level=sigma.sp.3.v4.est)
  mu.est.ebayes.sp.3.v5[i,]=waveti.ebayes(sim.m.sp.3.v5[i,],noise.level=sigma.sp.3.v5.est)
  mu.est.ebayes.bump.1.v1[i,]=waveti.ebayes(sim.m.bump.1.v1[i,],noise.level=sigma.bump.1.v1.est)
  mu.est.ebayes.bump.1.v2[i,]=waveti.ebayes(sim.m.bump.1.v2[i,],noise.level=sigma.bump.1.v2.est)
  mu.est.ebayes.bump.1.v3[i,]=waveti.ebayes(sim.m.bump.1.v3[i,],noise.level=sigma.bump.1.v3.est)
  mu.est.ebayes.bump.1.v4[i,]=waveti.ebayes(sim.m.bump.1.v4[i,],noise.level=sigma.bump.1.v4.est)
  mu.est.ebayes.bump.1.v5[i,]=waveti.ebayes(sim.m.bump.1.v5[i,],noise.level=sigma.bump.1.v5.est)
  mu.est.ebayes.bump.3.v1[i,]=waveti.ebayes(sim.m.bump.3.v1[i,],noise.level=sigma.bump.3.v1.est)
  mu.est.ebayes.bump.3.v2[i,]=waveti.ebayes(sim.m.bump.3.v2[i,],noise.level=sigma.bump.3.v2.est)
  mu.est.ebayes.bump.3.v3[i,]=waveti.ebayes(sim.m.bump.3.v3[i,],noise.level=sigma.bump.3.v3.est)
  mu.est.ebayes.bump.3.v4[i,]=waveti.ebayes(sim.m.bump.3.v4[i,],noise.level=sigma.bump.3.v4.est)
  mu.est.ebayes.bump.3.v5[i,]=waveti.ebayes(sim.m.bump.3.v5[i,],noise.level=sigma.bump.3.v5.est)
  mu.est.ebayes.blk.1.v1[i,]=waveti.ebayes(sim.m.blk.1.v1[i,],noise.level=sigma.blk.1.v1.est)
  mu.est.ebayes.blk.1.v2[i,]=waveti.ebayes(sim.m.blk.1.v2[i,],noise.level=sigma.blk.1.v2.est)
  mu.est.ebayes.blk.1.v3[i,]=waveti.ebayes(sim.m.blk.1.v3[i,],noise.level=sigma.blk.1.v3.est)
  mu.est.ebayes.blk.1.v4[i,]=waveti.ebayes(sim.m.blk.1.v4[i,],noise.level=sigma.blk.1.v4.est)
  mu.est.ebayes.blk.1.v5[i,]=waveti.ebayes(sim.m.blk.1.v5[i,],noise.level=sigma.blk.1.v5.est)
  mu.est.ebayes.blk.3.v1[i,]=waveti.ebayes(sim.m.blk.3.v1[i,],noise.level=sigma.blk.3.v1.est)
  mu.est.ebayes.blk.3.v2[i,]=waveti.ebayes(sim.m.blk.3.v2[i,],noise.level=sigma.blk.3.v2.est)
  mu.est.ebayes.blk.3.v3[i,]=waveti.ebayes(sim.m.blk.3.v3[i,],noise.level=sigma.blk.3.v3.est)
  mu.est.ebayes.blk.3.v4[i,]=waveti.ebayes(sim.m.blk.3.v4[i,],noise.level=sigma.blk.3.v4.est)
  mu.est.ebayes.blk.3.v5[i,]=waveti.ebayes(sim.m.blk.3.v5[i,],noise.level=sigma.blk.3.v5.est)
  mu.est.ebayes.ang.1.v1[i,]=waveti.ebayes(sim.m.ang.1.v1[i,],noise.level=sigma.ang.1.v1.est)
  mu.est.ebayes.ang.1.v2[i,]=waveti.ebayes(sim.m.ang.1.v2[i,],noise.level=sigma.ang.1.v2.est)
  mu.est.ebayes.ang.1.v3[i,]=waveti.ebayes(sim.m.ang.1.v3[i,],noise.level=sigma.ang.1.v3.est)
  mu.est.ebayes.ang.1.v4[i,]=waveti.ebayes(sim.m.ang.1.v4[i,],noise.level=sigma.ang.1.v4.est)
  mu.est.ebayes.ang.1.v5[i,]=waveti.ebayes(sim.m.ang.1.v5[i,],noise.level=sigma.ang.1.v5.est)
  mu.est.ebayes.ang.3.v1[i,]=waveti.ebayes(sim.m.ang.3.v1[i,],noise.level=sigma.ang.3.v1.est)
  mu.est.ebayes.ang.3.v2[i,]=waveti.ebayes(sim.m.ang.3.v2[i,],noise.level=sigma.ang.3.v2.est)
  mu.est.ebayes.ang.3.v3[i,]=waveti.ebayes(sim.m.ang.3.v3[i,],noise.level=sigma.ang.3.v3.est)
  mu.est.ebayes.ang.3.v4[i,]=waveti.ebayes(sim.m.ang.3.v4[i,],noise.level=sigma.ang.3.v4.est)
  mu.est.ebayes.ang.3.v5[i,]=waveti.ebayes(sim.m.ang.3.v5[i,],noise.level=sigma.ang.3.v5.est)
  mu.est.ebayes.dop.1.v1[i,]=waveti.ebayes(sim.m.dop.1.v1[i,],noise.level=sigma.dop.1.v1.est)
  mu.est.ebayes.dop.1.v2[i,]=waveti.ebayes(sim.m.dop.1.v2[i,],noise.level=sigma.dop.1.v2.est)
  mu.est.ebayes.dop.1.v3[i,]=waveti.ebayes(sim.m.dop.1.v3[i,],noise.level=sigma.dop.1.v3.est)
  mu.est.ebayes.dop.1.v4[i,]=waveti.ebayes(sim.m.dop.1.v4[i,],noise.level=sigma.dop.1.v4.est)
  mu.est.ebayes.dop.1.v5[i,]=waveti.ebayes(sim.m.dop.1.v5[i,],noise.level=sigma.dop.1.v5.est)
  mu.est.ebayes.dop.3.v1[i,]=waveti.ebayes(sim.m.dop.3.v1[i,],noise.level=sigma.dop.3.v1.est)
  mu.est.ebayes.dop.3.v2[i,]=waveti.ebayes(sim.m.dop.3.v2[i,],noise.level=sigma.dop.3.v2.est)
  mu.est.ebayes.dop.3.v3[i,]=waveti.ebayes(sim.m.dop.3.v3[i,],noise.level=sigma.dop.3.v3.est)
  mu.est.ebayes.dop.3.v4[i,]=waveti.ebayes(sim.m.dop.3.v4[i,],noise.level=sigma.dop.3.v4.est)
  mu.est.ebayes.dop.3.v5[i,]=waveti.ebayes(sim.m.dop.3.v5[i,],noise.level=sigma.dop.3.v5.est)
  mu.est.ebayes.blip.1.v1[i,]=waveti.ebayes(sim.m.blip.1.v1[i,],noise.level=sigma.blip.1.v1.est)
  mu.est.ebayes.blip.1.v2[i,]=waveti.ebayes(sim.m.blip.1.v2[i,],noise.level=sigma.blip.1.v2.est)
  mu.est.ebayes.blip.1.v3[i,]=waveti.ebayes(sim.m.blip.1.v3[i,],noise.level=sigma.blip.1.v3.est)
  mu.est.ebayes.blip.1.v4[i,]=waveti.ebayes(sim.m.blip.1.v4[i,],noise.level=sigma.blip.1.v4.est)
  mu.est.ebayes.blip.1.v5[i,]=waveti.ebayes(sim.m.blip.1.v5[i,],noise.level=sigma.blip.1.v5.est)
  mu.est.ebayes.blip.3.v1[i,]=waveti.ebayes(sim.m.blip.3.v1[i,],noise.level=sigma.blip.3.v1.est)
  mu.est.ebayes.blip.3.v2[i,]=waveti.ebayes(sim.m.blip.3.v2[i,],noise.level=sigma.blip.3.v2.est)
  mu.est.ebayes.blip.3.v3[i,]=waveti.ebayes(sim.m.blip.3.v3[i,],noise.level=sigma.blip.3.v3.est)
  mu.est.ebayes.blip.3.v4[i,]=waveti.ebayes(sim.m.blip.3.v4[i,],noise.level=sigma.blip.3.v4.est)
  mu.est.ebayes.blip.3.v5[i,]=waveti.ebayes(sim.m.blip.3.v5[i,],noise.level=sigma.blip.3.v5.est)
  mu.est.ebayes.cor.1.v1[i,]=waveti.ebayes(sim.m.cor.1.v1[i,],noise.level=sigma.cor.1.v1.est)
  mu.est.ebayes.cor.1.v2[i,]=waveti.ebayes(sim.m.cor.1.v2[i,],noise.level=sigma.cor.1.v2.est)
  mu.est.ebayes.cor.1.v3[i,]=waveti.ebayes(sim.m.cor.1.v3[i,],noise.level=sigma.cor.1.v3.est)
  mu.est.ebayes.cor.1.v4[i,]=waveti.ebayes(sim.m.cor.1.v4[i,],noise.level=sigma.cor.1.v4.est)
  mu.est.ebayes.cor.1.v5[i,]=waveti.ebayes(sim.m.cor.1.v5[i,],noise.level=sigma.cor.1.v5.est)
  mu.est.ebayes.cor.3.v1[i,]=waveti.ebayes(sim.m.cor.3.v1[i,],noise.level=sigma.cor.3.v1.est)
  mu.est.ebayes.cor.3.v2[i,]=waveti.ebayes(sim.m.cor.3.v2[i,],noise.level=sigma.cor.3.v2.est)
  mu.est.ebayes.cor.3.v3[i,]=waveti.ebayes(sim.m.cor.3.v3[i,],noise.level=sigma.cor.3.v3.est)
  mu.est.ebayes.cor.3.v4[i,]=waveti.ebayes(sim.m.cor.3.v4[i,],noise.level=sigma.cor.3.v4.est)
  mu.est.ebayes.cor.3.v5[i,]=waveti.ebayes(sim.m.cor.3.v5[i,],noise.level=sigma.cor.3.v5.est)

  mu.est.tieb.haar.sp.1.v1[i,]=ti.thresh(sim.m.sp.1.v1[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.sp.1.v1[i,]=ti.thresh(sim.m.sp.1.v1[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.sp.1.v2[i,]=ti.thresh(sim.m.sp.1.v2[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.sp.1.v2[i,]=ti.thresh(sim.m.sp.1.v2[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.sp.1.v3[i,]=ti.thresh(sim.m.sp.1.v3[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.sp.1.v3[i,]=ti.thresh(sim.m.sp.1.v3[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.sp.1.v4[i,]=ti.thresh(sim.m.sp.1.v4[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.sp.1.v4[i,]=ti.thresh(sim.m.sp.1.v4[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.sp.1.v5[i,]=ti.thresh(sim.m.sp.1.v5[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.sp.1.v5[i,]=ti.thresh(sim.m.sp.1.v5[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.sp.3.v1[i,]=ti.thresh(sim.m.sp.3.v1[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.sp.3.v1[i,]=ti.thresh(sim.m.sp.3.v1[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.sp.3.v2[i,]=ti.thresh(sim.m.sp.3.v2[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.sp.3.v2[i,]=ti.thresh(sim.m.sp.3.v2[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.sp.3.v3[i,]=ti.thresh(sim.m.sp.3.v3[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.sp.3.v3[i,]=ti.thresh(sim.m.sp.3.v3[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.sp.3.v4[i,]=ti.thresh(sim.m.sp.3.v4[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.sp.3.v4[i,]=ti.thresh(sim.m.sp.3.v4[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.sp.3.v5[i,]=ti.thresh(sim.m.sp.3.v5[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.sp.3.v5[i,]=ti.thresh(sim.m.sp.3.v5[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")


  mu.est.tieb.haar.bump.1.v1[i,]=ti.thresh(sim.m.bump.1.v1[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.bump.1.v1[i,]=ti.thresh(sim.m.bump.1.v1[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.bump.1.v2[i,]=ti.thresh(sim.m.bump.1.v2[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.bump.1.v2[i,]=ti.thresh(sim.m.bump.1.v2[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.bump.1.v3[i,]=ti.thresh(sim.m.bump.1.v3[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.bump.1.v3[i,]=ti.thresh(sim.m.bump.1.v3[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.bump.1.v4[i,]=ti.thresh(sim.m.bump.1.v4[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.bump.1.v4[i,]=ti.thresh(sim.m.bump.1.v4[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.bump.1.v5[i,]=ti.thresh(sim.m.bump.1.v5[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.bump.1.v5[i,]=ti.thresh(sim.m.bump.1.v5[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.bump.3.v1[i,]=ti.thresh(sim.m.bump.3.v1[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.bump.3.v1[i,]=ti.thresh(sim.m.bump.3.v1[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.bump.3.v2[i,]=ti.thresh(sim.m.bump.3.v2[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.bump.3.v2[i,]=ti.thresh(sim.m.bump.3.v2[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.bump.3.v3[i,]=ti.thresh(sim.m.bump.3.v3[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.bump.3.v3[i,]=ti.thresh(sim.m.bump.3.v3[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.bump.3.v4[i,]=ti.thresh(sim.m.bump.3.v4[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.bump.3.v4[i,]=ti.thresh(sim.m.bump.3.v4[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.bump.3.v5[i,]=ti.thresh(sim.m.bump.3.v5[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.bump.3.v5[i,]=ti.thresh(sim.m.bump.3.v5[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.blk.1.v1[i,]=ti.thresh(sim.m.blk.1.v1[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.blk.1.v1[i,]=ti.thresh(sim.m.blk.1.v1[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.blk.1.v2[i,]=ti.thresh(sim.m.blk.1.v2[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.blk.1.v2[i,]=ti.thresh(sim.m.blk.1.v2[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.blk.1.v3[i,]=ti.thresh(sim.m.blk.1.v3[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.blk.1.v3[i,]=ti.thresh(sim.m.blk.1.v3[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.blk.1.v4[i,]=ti.thresh(sim.m.blk.1.v4[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.blk.1.v4[i,]=ti.thresh(sim.m.blk.1.v4[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.blk.1.v5[i,]=ti.thresh(sim.m.blk.1.v5[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.blk.1.v5[i,]=ti.thresh(sim.m.blk.1.v5[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.blk.3.v1[i,]=ti.thresh(sim.m.blk.3.v1[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.blk.3.v1[i,]=ti.thresh(sim.m.blk.3.v1[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.blk.3.v2[i,]=ti.thresh(sim.m.blk.3.v2[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.blk.3.v2[i,]=ti.thresh(sim.m.blk.3.v2[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.blk.3.v3[i,]=ti.thresh(sim.m.blk.3.v3[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.blk.3.v3[i,]=ti.thresh(sim.m.blk.3.v3[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.blk.3.v4[i,]=ti.thresh(sim.m.blk.3.v4[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.blk.3.v4[i,]=ti.thresh(sim.m.blk.3.v4[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.blk.3.v5[i,]=ti.thresh(sim.m.blk.3.v5[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.blk.3.v5[i,]=ti.thresh(sim.m.blk.3.v5[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")


  mu.est.tieb.haar.ang.1.v1[i,]=ti.thresh(sim.m.ang.1.v1[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.ang.1.v1[i,]=ti.thresh(sim.m.ang.1.v1[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.ang.1.v2[i,]=ti.thresh(sim.m.ang.1.v2[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.ang.1.v2[i,]=ti.thresh(sim.m.ang.1.v2[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.ang.1.v3[i,]=ti.thresh(sim.m.ang.1.v3[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.ang.1.v3[i,]=ti.thresh(sim.m.ang.1.v3[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.ang.1.v4[i,]=ti.thresh(sim.m.ang.1.v4[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.ang.1.v4[i,]=ti.thresh(sim.m.ang.1.v4[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.ang.1.v5[i,]=ti.thresh(sim.m.ang.1.v5[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.ang.1.v5[i,]=ti.thresh(sim.m.ang.1.v5[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.ang.3.v1[i,]=ti.thresh(sim.m.ang.3.v1[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.ang.3.v1[i,]=ti.thresh(sim.m.ang.3.v1[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.ang.3.v2[i,]=ti.thresh(sim.m.ang.3.v2[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.ang.3.v2[i,]=ti.thresh(sim.m.ang.3.v2[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.ang.3.v3[i,]=ti.thresh(sim.m.ang.3.v3[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.ang.3.v3[i,]=ti.thresh(sim.m.ang.3.v3[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.ang.3.v4[i,]=ti.thresh(sim.m.ang.3.v4[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.ang.3.v4[i,]=ti.thresh(sim.m.ang.3.v4[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.ang.3.v5[i,]=ti.thresh(sim.m.ang.3.v5[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.ang.3.v5[i,]=ti.thresh(sim.m.ang.3.v5[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")


  mu.est.tieb.haar.dop.1.v1[i,]=ti.thresh(sim.m.dop.1.v1[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.dop.1.v1[i,]=ti.thresh(sim.m.dop.1.v1[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.dop.1.v2[i,]=ti.thresh(sim.m.dop.1.v2[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.dop.1.v2[i,]=ti.thresh(sim.m.dop.1.v2[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.dop.1.v3[i,]=ti.thresh(sim.m.dop.1.v3[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.dop.1.v3[i,]=ti.thresh(sim.m.dop.1.v3[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.dop.1.v4[i,]=ti.thresh(sim.m.dop.1.v4[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.dop.1.v4[i,]=ti.thresh(sim.m.dop.1.v4[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.dop.1.v5[i,]=ti.thresh(sim.m.dop.1.v5[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.dop.1.v5[i,]=ti.thresh(sim.m.dop.1.v5[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.dop.3.v1[i,]=ti.thresh(sim.m.dop.3.v1[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.dop.3.v1[i,]=ti.thresh(sim.m.dop.3.v1[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.dop.3.v2[i,]=ti.thresh(sim.m.dop.3.v2[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.dop.3.v2[i,]=ti.thresh(sim.m.dop.3.v2[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.dop.3.v3[i,]=ti.thresh(sim.m.dop.3.v3[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.dop.3.v3[i,]=ti.thresh(sim.m.dop.3.v3[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.dop.3.v4[i,]=ti.thresh(sim.m.dop.3.v4[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.dop.3.v4[i,]=ti.thresh(sim.m.dop.3.v4[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.dop.3.v5[i,]=ti.thresh(sim.m.dop.3.v5[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.dop.3.v5[i,]=ti.thresh(sim.m.dop.3.v5[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")


  mu.est.tieb.haar.blip.1.v1[i,]=ti.thresh(sim.m.blip.1.v1[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.blip.1.v1[i,]=ti.thresh(sim.m.blip.1.v1[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.blip.1.v2[i,]=ti.thresh(sim.m.blip.1.v2[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.blip.1.v2[i,]=ti.thresh(sim.m.blip.1.v2[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.blip.1.v3[i,]=ti.thresh(sim.m.blip.1.v3[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.blip.1.v3[i,]=ti.thresh(sim.m.blip.1.v3[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.blip.1.v4[i,]=ti.thresh(sim.m.blip.1.v4[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.blip.1.v4[i,]=ti.thresh(sim.m.blip.1.v4[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.blip.1.v5[i,]=ti.thresh(sim.m.blip.1.v5[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.blip.1.v5[i,]=ti.thresh(sim.m.blip.1.v5[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.blip.3.v1[i,]=ti.thresh(sim.m.blip.3.v1[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.blip.3.v1[i,]=ti.thresh(sim.m.blip.3.v1[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.blip.3.v2[i,]=ti.thresh(sim.m.blip.3.v2[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.blip.3.v2[i,]=ti.thresh(sim.m.blip.3.v2[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.blip.3.v3[i,]=ti.thresh(sim.m.blip.3.v3[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.blip.3.v3[i,]=ti.thresh(sim.m.blip.3.v3[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.blip.3.v4[i,]=ti.thresh(sim.m.blip.3.v4[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.blip.3.v4[i,]=ti.thresh(sim.m.blip.3.v4[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.blip.3.v5[i,]=ti.thresh(sim.m.blip.3.v5[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.blip.3.v5[i,]=ti.thresh(sim.m.blip.3.v5[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")



  mu.est.tieb.haar.cor.1.v1[i,]=ti.thresh(sim.m.cor.1.v1[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.cor.1.v1[i,]=ti.thresh(sim.m.cor.1.v1[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.cor.1.v2[i,]=ti.thresh(sim.m.cor.1.v2[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.cor.1.v2[i,]=ti.thresh(sim.m.cor.1.v2[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.cor.1.v3[i,]=ti.thresh(sim.m.cor.1.v3[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.cor.1.v3[i,]=ti.thresh(sim.m.cor.1.v3[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.cor.1.v4[i,]=ti.thresh(sim.m.cor.1.v4[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.cor.1.v4[i,]=ti.thresh(sim.m.cor.1.v4[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.cor.1.v5[i,]=ti.thresh(sim.m.cor.1.v5[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.cor.1.v5[i,]=ti.thresh(sim.m.cor.1.v5[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.cor.3.v1[i,]=ti.thresh(sim.m.cor.3.v1[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.cor.3.v1[i,]=ti.thresh(sim.m.cor.3.v1[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.cor.3.v2[i,]=ti.thresh(sim.m.cor.3.v2[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.cor.3.v2[i,]=ti.thresh(sim.m.cor.3.v2[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.cor.3.v3[i,]=ti.thresh(sim.m.cor.3.v3[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.cor.3.v3[i,]=ti.thresh(sim.m.cor.3.v3[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.cor.3.v4[i,]=ti.thresh(sim.m.cor.3.v4[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.cor.3.v4[i,]=ti.thresh(sim.m.cor.3.v4[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.cor.3.v5[i,]=ti.thresh(sim.m.cor.3.v5[i,],method="bayesm",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.cor.3.v5[i,]=ti.thresh(sim.m.cor.3.v5[i,],method="bayesm",filter.number=8,family="DaubLeAsymm")




  mu.est.ashe.haar.j.sp.1.v1[i,]=bayesmooth(sim.m.sp.1.v1[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.sp.1.v1[i,]=bayesmooth(sim.m.sp.1.v1[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.sp.1.v2[i,]=bayesmooth(sim.m.sp.1.v2[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.sp.1.v2[i,]=bayesmooth(sim.m.sp.1.v2[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.sp.1.v3[i,]=bayesmooth(sim.m.sp.1.v3[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.sp.1.v3[i,]=bayesmooth(sim.m.sp.1.v3[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.sp.1.v4[i,]=bayesmooth(sim.m.sp.1.v4[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.sp.1.v4[i,]=bayesmooth(sim.m.sp.1.v4[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.sp.1.v5[i,]=bayesmooth(sim.m.sp.1.v5[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.sp.1.v5[i,]=bayesmooth(sim.m.sp.1.v5[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.sp.3.v1[i,]=bayesmooth(sim.m.sp.3.v1[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.sp.3.v1[i,]=bayesmooth(sim.m.sp.3.v1[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.sp.3.v2[i,]=bayesmooth(sim.m.sp.3.v2[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.sp.3.v2[i,]=bayesmooth(sim.m.sp.3.v2[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.sp.3.v3[i,]=bayesmooth(sim.m.sp.3.v3[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.sp.3.v3[i,]=bayesmooth(sim.m.sp.3.v3[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.sp.3.v4[i,]=bayesmooth(sim.m.sp.3.v4[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.sp.3.v4[i,]=bayesmooth(sim.m.sp.3.v4[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.sp.3.v5[i,]=bayesmooth(sim.m.sp.3.v5[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.sp.3.v5[i,]=bayesmooth(sim.m.sp.3.v5[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)


  mu.est.ashe.haar.j.bump.1.v1[i,]=bayesmooth(sim.m.bump.1.v1[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.bump.1.v1[i,]=bayesmooth(sim.m.bump.1.v1[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.bump.1.v2[i,]=bayesmooth(sim.m.bump.1.v2[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.bump.1.v2[i,]=bayesmooth(sim.m.bump.1.v2[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.bump.1.v3[i,]=bayesmooth(sim.m.bump.1.v3[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.bump.1.v3[i,]=bayesmooth(sim.m.bump.1.v3[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.bump.1.v4[i,]=bayesmooth(sim.m.bump.1.v4[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.bump.1.v4[i,]=bayesmooth(sim.m.bump.1.v4[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.bump.1.v5[i,]=bayesmooth(sim.m.bump.1.v5[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.bump.1.v5[i,]=bayesmooth(sim.m.bump.1.v5[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.bump.3.v1[i,]=bayesmooth(sim.m.bump.3.v1[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.bump.3.v1[i,]=bayesmooth(sim.m.bump.3.v1[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.bump.3.v2[i,]=bayesmooth(sim.m.bump.3.v2[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.bump.3.v2[i,]=bayesmooth(sim.m.bump.3.v2[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.bump.3.v3[i,]=bayesmooth(sim.m.bump.3.v3[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.bump.3.v3[i,]=bayesmooth(sim.m.bump.3.v3[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.bump.3.v4[i,]=bayesmooth(sim.m.bump.3.v4[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.bump.3.v4[i,]=bayesmooth(sim.m.bump.3.v4[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.bump.3.v5[i,]=bayesmooth(sim.m.bump.3.v5[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.bump.3.v5[i,]=bayesmooth(sim.m.bump.3.v5[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.blk.1.v1[i,]=bayesmooth(sim.m.blk.1.v1[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.blk.1.v1[i,]=bayesmooth(sim.m.blk.1.v1[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.blk.1.v2[i,]=bayesmooth(sim.m.blk.1.v2[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.blk.1.v2[i,]=bayesmooth(sim.m.blk.1.v2[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.blk.1.v3[i,]=bayesmooth(sim.m.blk.1.v3[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.blk.1.v3[i,]=bayesmooth(sim.m.blk.1.v3[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.blk.1.v4[i,]=bayesmooth(sim.m.blk.1.v4[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.blk.1.v4[i,]=bayesmooth(sim.m.blk.1.v4[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.blk.1.v5[i,]=bayesmooth(sim.m.blk.1.v5[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.blk.1.v5[i,]=bayesmooth(sim.m.blk.1.v5[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.blk.3.v1[i,]=bayesmooth(sim.m.blk.3.v1[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.blk.3.v1[i,]=bayesmooth(sim.m.blk.3.v1[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.blk.3.v2[i,]=bayesmooth(sim.m.blk.3.v2[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.blk.3.v2[i,]=bayesmooth(sim.m.blk.3.v2[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.blk.3.v3[i,]=bayesmooth(sim.m.blk.3.v3[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.blk.3.v3[i,]=bayesmooth(sim.m.blk.3.v3[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.blk.3.v4[i,]=bayesmooth(sim.m.blk.3.v4[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.blk.3.v4[i,]=bayesmooth(sim.m.blk.3.v4[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.blk.3.v5[i,]=bayesmooth(sim.m.blk.3.v5[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.blk.3.v5[i,]=bayesmooth(sim.m.blk.3.v5[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)


  mu.est.ashe.haar.j.ang.1.v1[i,]=bayesmooth(sim.m.ang.1.v1[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.ang.1.v1[i,]=bayesmooth(sim.m.ang.1.v1[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.ang.1.v2[i,]=bayesmooth(sim.m.ang.1.v2[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.ang.1.v2[i,]=bayesmooth(sim.m.ang.1.v2[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.ang.1.v3[i,]=bayesmooth(sim.m.ang.1.v3[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.ang.1.v3[i,]=bayesmooth(sim.m.ang.1.v3[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.ang.1.v4[i,]=bayesmooth(sim.m.ang.1.v4[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.ang.1.v4[i,]=bayesmooth(sim.m.ang.1.v4[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.ang.1.v5[i,]=bayesmooth(sim.m.ang.1.v5[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.ang.1.v5[i,]=bayesmooth(sim.m.ang.1.v5[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.ang.3.v1[i,]=bayesmooth(sim.m.ang.3.v1[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.ang.3.v1[i,]=bayesmooth(sim.m.ang.3.v1[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.ang.3.v2[i,]=bayesmooth(sim.m.ang.3.v2[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.ang.3.v2[i,]=bayesmooth(sim.m.ang.3.v2[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.ang.3.v3[i,]=bayesmooth(sim.m.ang.3.v3[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.ang.3.v3[i,]=bayesmooth(sim.m.ang.3.v3[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.ang.3.v4[i,]=bayesmooth(sim.m.ang.3.v4[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.ang.3.v4[i,]=bayesmooth(sim.m.ang.3.v4[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.ang.3.v5[i,]=bayesmooth(sim.m.ang.3.v5[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.ang.3.v5[i,]=bayesmooth(sim.m.ang.3.v5[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)


  mu.est.ashe.haar.j.dop.1.v1[i,]=bayesmooth(sim.m.dop.1.v1[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.dop.1.v1[i,]=bayesmooth(sim.m.dop.1.v1[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.dop.1.v2[i,]=bayesmooth(sim.m.dop.1.v2[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.dop.1.v2[i,]=bayesmooth(sim.m.dop.1.v2[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.dop.1.v3[i,]=bayesmooth(sim.m.dop.1.v3[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.dop.1.v3[i,]=bayesmooth(sim.m.dop.1.v3[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.dop.1.v4[i,]=bayesmooth(sim.m.dop.1.v4[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.dop.1.v4[i,]=bayesmooth(sim.m.dop.1.v4[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.dop.1.v5[i,]=bayesmooth(sim.m.dop.1.v5[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.dop.1.v5[i,]=bayesmooth(sim.m.dop.1.v5[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.dop.3.v1[i,]=bayesmooth(sim.m.dop.3.v1[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.dop.3.v1[i,]=bayesmooth(sim.m.dop.3.v1[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.dop.3.v2[i,]=bayesmooth(sim.m.dop.3.v2[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.dop.3.v2[i,]=bayesmooth(sim.m.dop.3.v2[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.dop.3.v3[i,]=bayesmooth(sim.m.dop.3.v3[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.dop.3.v3[i,]=bayesmooth(sim.m.dop.3.v3[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.dop.3.v4[i,]=bayesmooth(sim.m.dop.3.v4[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.dop.3.v4[i,]=bayesmooth(sim.m.dop.3.v4[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.dop.3.v5[i,]=bayesmooth(sim.m.dop.3.v5[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.dop.3.v5[i,]=bayesmooth(sim.m.dop.3.v5[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)


  mu.est.ashe.haar.j.blip.1.v1[i,]=bayesmooth(sim.m.blip.1.v1[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.blip.1.v1[i,]=bayesmooth(sim.m.blip.1.v1[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.blip.1.v2[i,]=bayesmooth(sim.m.blip.1.v2[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.blip.1.v2[i,]=bayesmooth(sim.m.blip.1.v2[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.blip.1.v3[i,]=bayesmooth(sim.m.blip.1.v3[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.blip.1.v3[i,]=bayesmooth(sim.m.blip.1.v3[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.blip.1.v4[i,]=bayesmooth(sim.m.blip.1.v4[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.blip.1.v4[i,]=bayesmooth(sim.m.blip.1.v4[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.blip.1.v5[i,]=bayesmooth(sim.m.blip.1.v5[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.blip.1.v5[i,]=bayesmooth(sim.m.blip.1.v5[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.blip.3.v1[i,]=bayesmooth(sim.m.blip.3.v1[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.blip.3.v1[i,]=bayesmooth(sim.m.blip.3.v1[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.blip.3.v2[i,]=bayesmooth(sim.m.blip.3.v2[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.blip.3.v2[i,]=bayesmooth(sim.m.blip.3.v2[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.blip.3.v3[i,]=bayesmooth(sim.m.blip.3.v3[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.blip.3.v3[i,]=bayesmooth(sim.m.blip.3.v3[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.blip.3.v4[i,]=bayesmooth(sim.m.blip.3.v4[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.blip.3.v4[i,]=bayesmooth(sim.m.blip.3.v4[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.blip.3.v5[i,]=bayesmooth(sim.m.blip.3.v5[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.blip.3.v5[i,]=bayesmooth(sim.m.blip.3.v5[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)



  mu.est.ashe.haar.j.cor.1.v1[i,]=bayesmooth(sim.m.cor.1.v1[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.cor.1.v1[i,]=bayesmooth(sim.m.cor.1.v1[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.cor.1.v2[i,]=bayesmooth(sim.m.cor.1.v2[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.cor.1.v2[i,]=bayesmooth(sim.m.cor.1.v2[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.cor.1.v3[i,]=bayesmooth(sim.m.cor.1.v3[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.cor.1.v3[i,]=bayesmooth(sim.m.cor.1.v3[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.cor.1.v4[i,]=bayesmooth(sim.m.cor.1.v4[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.cor.1.v4[i,]=bayesmooth(sim.m.cor.1.v4[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.cor.1.v5[i,]=bayesmooth(sim.m.cor.1.v5[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.cor.1.v5[i,]=bayesmooth(sim.m.cor.1.v5[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.cor.3.v1[i,]=bayesmooth(sim.m.cor.3.v1[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.cor.3.v1[i,]=bayesmooth(sim.m.cor.3.v1[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.cor.3.v2[i,]=bayesmooth(sim.m.cor.3.v2[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.cor.3.v2[i,]=bayesmooth(sim.m.cor.3.v2[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.cor.3.v3[i,]=bayesmooth(sim.m.cor.3.v3[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.cor.3.v3[i,]=bayesmooth(sim.m.cor.3.v3[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.cor.3.v4[i,]=bayesmooth(sim.m.cor.3.v4[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.cor.3.v4[i,]=bayesmooth(sim.m.cor.3.v4[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.haar.j.cor.3.v5[i,]=bayesmooth(sim.m.cor.3.v5[i,],filter.number=1,family="DaubExPhase",jash=TRUE)
  mu.est.ashe.s8.j.h.cor.3.v5[i,]=bayesmooth(sim.m.cor.3.v5[i,],filter.number=8,family="DaubLeAsymm",jash=TRUE)




  mu.est.ashe.s8.j.sp.1.v1[i,]=bayesmooth(sim.m.sp.1.v1[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.sp.1.v2[i,]=bayesmooth(sim.m.sp.1.v2[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.sp.1.v3[i,]=bayesmooth(sim.m.sp.1.v3[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.sp.1.v4[i,]=bayesmooth(sim.m.sp.1.v4[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.sp.1.v5[i,]=bayesmooth(sim.m.sp.1.v5[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.sp.3.v1[i,]=bayesmooth(sim.m.sp.3.v1[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.sp.3.v2[i,]=bayesmooth(sim.m.sp.3.v2[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.sp.3.v3[i,]=bayesmooth(sim.m.sp.3.v3[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.sp.3.v4[i,]=bayesmooth(sim.m.sp.3.v4[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.sp.3.v5[i,]=bayesmooth(sim.m.sp.3.v5[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)


  mu.est.ashe.s8.j.bump.1.v1[i,]=bayesmooth(sim.m.bump.1.v1[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.bump.1.v2[i,]=bayesmooth(sim.m.bump.1.v2[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.bump.1.v3[i,]=bayesmooth(sim.m.bump.1.v3[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.bump.1.v4[i,]=bayesmooth(sim.m.bump.1.v4[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.bump.1.v5[i,]=bayesmooth(sim.m.bump.1.v5[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.bump.3.v1[i,]=bayesmooth(sim.m.bump.3.v1[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.bump.3.v2[i,]=bayesmooth(sim.m.bump.3.v2[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.bump.3.v3[i,]=bayesmooth(sim.m.bump.3.v3[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.bump.3.v4[i,]=bayesmooth(sim.m.bump.3.v4[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.bump.3.v5[i,]=bayesmooth(sim.m.bump.3.v5[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.blk.1.v1[i,]=bayesmooth(sim.m.blk.1.v1[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.blk.1.v2[i,]=bayesmooth(sim.m.blk.1.v2[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.blk.1.v3[i,]=bayesmooth(sim.m.blk.1.v3[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.blk.1.v4[i,]=bayesmooth(sim.m.blk.1.v4[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.blk.1.v5[i,]=bayesmooth(sim.m.blk.1.v5[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.blk.3.v1[i,]=bayesmooth(sim.m.blk.3.v1[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.blk.3.v2[i,]=bayesmooth(sim.m.blk.3.v2[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.blk.3.v3[i,]=bayesmooth(sim.m.blk.3.v3[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.blk.3.v4[i,]=bayesmooth(sim.m.blk.3.v4[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.blk.3.v5[i,]=bayesmooth(sim.m.blk.3.v5[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)


  mu.est.ashe.s8.j.ang.1.v1[i,]=bayesmooth(sim.m.ang.1.v1[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.ang.1.v2[i,]=bayesmooth(sim.m.ang.1.v2[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.ang.1.v3[i,]=bayesmooth(sim.m.ang.1.v3[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.ang.1.v4[i,]=bayesmooth(sim.m.ang.1.v4[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.ang.1.v5[i,]=bayesmooth(sim.m.ang.1.v5[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.ang.3.v1[i,]=bayesmooth(sim.m.ang.3.v1[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.ang.3.v2[i,]=bayesmooth(sim.m.ang.3.v2[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.ang.3.v3[i,]=bayesmooth(sim.m.ang.3.v3[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.ang.3.v4[i,]=bayesmooth(sim.m.ang.3.v4[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.ang.3.v5[i,]=bayesmooth(sim.m.ang.3.v5[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)


  mu.est.ashe.s8.j.dop.1.v1[i,]=bayesmooth(sim.m.dop.1.v1[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.dop.1.v2[i,]=bayesmooth(sim.m.dop.1.v2[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.dop.1.v3[i,]=bayesmooth(sim.m.dop.1.v3[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.dop.1.v4[i,]=bayesmooth(sim.m.dop.1.v4[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.dop.1.v5[i,]=bayesmooth(sim.m.dop.1.v5[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.dop.3.v1[i,]=bayesmooth(sim.m.dop.3.v1[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.dop.3.v2[i,]=bayesmooth(sim.m.dop.3.v2[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.dop.3.v3[i,]=bayesmooth(sim.m.dop.3.v3[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.dop.3.v4[i,]=bayesmooth(sim.m.dop.3.v4[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.dop.3.v5[i,]=bayesmooth(sim.m.dop.3.v5[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)


  mu.est.ashe.s8.j.blip.1.v1[i,]=bayesmooth(sim.m.blip.1.v1[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.blip.1.v2[i,]=bayesmooth(sim.m.blip.1.v2[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.blip.1.v3[i,]=bayesmooth(sim.m.blip.1.v3[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.blip.1.v4[i,]=bayesmooth(sim.m.blip.1.v4[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.blip.1.v5[i,]=bayesmooth(sim.m.blip.1.v5[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.blip.3.v1[i,]=bayesmooth(sim.m.blip.3.v1[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.blip.3.v2[i,]=bayesmooth(sim.m.blip.3.v2[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.blip.3.v3[i,]=bayesmooth(sim.m.blip.3.v3[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.blip.3.v4[i,]=bayesmooth(sim.m.blip.3.v4[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.blip.3.v5[i,]=bayesmooth(sim.m.blip.3.v5[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)



  mu.est.ashe.s8.j.cor.1.v1[i,]=bayesmooth(sim.m.cor.1.v1[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.cor.1.v2[i,]=bayesmooth(sim.m.cor.1.v2[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.cor.1.v3[i,]=bayesmooth(sim.m.cor.1.v3[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.cor.1.v4[i,]=bayesmooth(sim.m.cor.1.v4[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.cor.1.v5[i,]=bayesmooth(sim.m.cor.1.v5[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.cor.3.v1[i,]=bayesmooth(sim.m.cor.3.v1[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.cor.3.v2[i,]=bayesmooth(sim.m.cor.3.v2[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.cor.3.v3[i,]=bayesmooth(sim.m.cor.3.v3[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.cor.3.v4[i,]=bayesmooth(sim.m.cor.3.v4[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.cor.3.v5[i,]=bayesmooth(sim.m.cor.3.v5[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashef.s8.sp.1.v1[i,]=bayesmooth(sim.m.sp.1.v1[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.sp.1.v2[i,]=bayesmooth(sim.m.sp.1.v2[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.sp.1.v3[i,]=bayesmooth(sim.m.sp.1.v3[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.sp.1.v4[i,]=bayesmooth(sim.m.sp.1.v4[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.sp.1.v5[i,]=bayesmooth(sim.m.sp.1.v5[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.sp.3.v1[i,]=bayesmooth(sim.m.sp.3.v1[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.sp.3.v2[i,]=bayesmooth(sim.m.sp.3.v2[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.sp.3.v3[i,]=bayesmooth(sim.m.sp.3.v3[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.sp.3.v4[i,]=bayesmooth(sim.m.sp.3.v4[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.sp.3.v5[i,]=bayesmooth(sim.m.sp.3.v5[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)


  mu.est.ashef.s8.bump.1.v1[i,]=bayesmooth(sim.m.bump.1.v1[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.bump.1.v2[i,]=bayesmooth(sim.m.bump.1.v2[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.bump.1.v3[i,]=bayesmooth(sim.m.bump.1.v3[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.bump.1.v4[i,]=bayesmooth(sim.m.bump.1.v4[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.bump.1.v5[i,]=bayesmooth(sim.m.bump.1.v5[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.bump.3.v1[i,]=bayesmooth(sim.m.bump.3.v1[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.bump.3.v2[i,]=bayesmooth(sim.m.bump.3.v2[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.bump.3.v3[i,]=bayesmooth(sim.m.bump.3.v3[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.bump.3.v4[i,]=bayesmooth(sim.m.bump.3.v4[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.bump.3.v5[i,]=bayesmooth(sim.m.bump.3.v5[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.blk.1.v1[i,]=bayesmooth(sim.m.blk.1.v1[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.blk.1.v2[i,]=bayesmooth(sim.m.blk.1.v2[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.blk.1.v3[i,]=bayesmooth(sim.m.blk.1.v3[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.blk.1.v4[i,]=bayesmooth(sim.m.blk.1.v4[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.blk.1.v5[i,]=bayesmooth(sim.m.blk.1.v5[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.blk.3.v1[i,]=bayesmooth(sim.m.blk.3.v1[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.blk.3.v2[i,]=bayesmooth(sim.m.blk.3.v2[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.blk.3.v3[i,]=bayesmooth(sim.m.blk.3.v3[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.blk.3.v4[i,]=bayesmooth(sim.m.blk.3.v4[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.blk.3.v5[i,]=bayesmooth(sim.m.blk.3.v5[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)


  mu.est.ashef.s8.ang.1.v1[i,]=bayesmooth(sim.m.ang.1.v1[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.ang.1.v2[i,]=bayesmooth(sim.m.ang.1.v2[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.ang.1.v3[i,]=bayesmooth(sim.m.ang.1.v3[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.ang.1.v4[i,]=bayesmooth(sim.m.ang.1.v4[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.ang.1.v5[i,]=bayesmooth(sim.m.ang.1.v5[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.ang.3.v1[i,]=bayesmooth(sim.m.ang.3.v1[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.ang.3.v2[i,]=bayesmooth(sim.m.ang.3.v2[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.ang.3.v3[i,]=bayesmooth(sim.m.ang.3.v3[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.ang.3.v4[i,]=bayesmooth(sim.m.ang.3.v4[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.ang.3.v5[i,]=bayesmooth(sim.m.ang.3.v5[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)


  mu.est.ashef.s8.dop.1.v1[i,]=bayesmooth(sim.m.dop.1.v1[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.dop.1.v2[i,]=bayesmooth(sim.m.dop.1.v2[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.dop.1.v3[i,]=bayesmooth(sim.m.dop.1.v3[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.dop.1.v4[i,]=bayesmooth(sim.m.dop.1.v4[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.dop.1.v5[i,]=bayesmooth(sim.m.dop.1.v5[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.dop.3.v1[i,]=bayesmooth(sim.m.dop.3.v1[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.dop.3.v2[i,]=bayesmooth(sim.m.dop.3.v2[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.dop.3.v3[i,]=bayesmooth(sim.m.dop.3.v3[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.dop.3.v4[i,]=bayesmooth(sim.m.dop.3.v4[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.dop.3.v5[i,]=bayesmooth(sim.m.dop.3.v5[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)


  mu.est.ashef.s8.blip.1.v1[i,]=bayesmooth(sim.m.blip.1.v1[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.blip.1.v2[i,]=bayesmooth(sim.m.blip.1.v2[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.blip.1.v3[i,]=bayesmooth(sim.m.blip.1.v3[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.blip.1.v4[i,]=bayesmooth(sim.m.blip.1.v4[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.blip.1.v5[i,]=bayesmooth(sim.m.blip.1.v5[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.blip.3.v1[i,]=bayesmooth(sim.m.blip.3.v1[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.blip.3.v2[i,]=bayesmooth(sim.m.blip.3.v2[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.blip.3.v3[i,]=bayesmooth(sim.m.blip.3.v3[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.blip.3.v4[i,]=bayesmooth(sim.m.blip.3.v4[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.blip.3.v5[i,]=bayesmooth(sim.m.blip.3.v5[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)



  mu.est.ashef.s8.cor.1.v1[i,]=bayesmooth(sim.m.cor.1.v1[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.cor.1.v2[i,]=bayesmooth(sim.m.cor.1.v2[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.cor.1.v3[i,]=bayesmooth(sim.m.cor.1.v3[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.cor.1.v4[i,]=bayesmooth(sim.m.cor.1.v4[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.cor.1.v5[i,]=bayesmooth(sim.m.cor.1.v5[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.cor.3.v1[i,]=bayesmooth(sim.m.cor.3.v1[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.cor.3.v2[i,]=bayesmooth(sim.m.cor.3.v2[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.cor.3.v3[i,]=bayesmooth(sim.m.cor.3.v3[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.cor.3.v4[i,]=bayesmooth(sim.m.cor.3.v4[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)

  mu.est.ashef.s8.cor.3.v5[i,]=bayesmooth(sim.m.cor.3.v5[i,],filter.number=8,family="DaubLeAsymm",gridmult=2)


  mu.est.ashef.haar.sp.1.v1[i,]=bayesmooth(sim.m.sp.1.v1[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.sp.1.v2[i,]=bayesmooth(sim.m.sp.1.v2[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.sp.1.v3[i,]=bayesmooth(sim.m.sp.1.v3[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.sp.1.v4[i,]=bayesmooth(sim.m.sp.1.v4[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.sp.1.v5[i,]=bayesmooth(sim.m.sp.1.v5[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.sp.3.v1[i,]=bayesmooth(sim.m.sp.3.v1[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.sp.3.v2[i,]=bayesmooth(sim.m.sp.3.v2[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.sp.3.v3[i,]=bayesmooth(sim.m.sp.3.v3[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.sp.3.v4[i,]=bayesmooth(sim.m.sp.3.v4[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.sp.3.v5[i,]=bayesmooth(sim.m.sp.3.v5[i,],filter.number=1,family="DaubExPhase",gridmult=2)


  mu.est.ashef.haar.bump.1.v1[i,]=bayesmooth(sim.m.bump.1.v1[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.bump.1.v2[i,]=bayesmooth(sim.m.bump.1.v2[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.bump.1.v3[i,]=bayesmooth(sim.m.bump.1.v3[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.bump.1.v4[i,]=bayesmooth(sim.m.bump.1.v4[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.bump.1.v5[i,]=bayesmooth(sim.m.bump.1.v5[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.bump.3.v1[i,]=bayesmooth(sim.m.bump.3.v1[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.bump.3.v2[i,]=bayesmooth(sim.m.bump.3.v2[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.bump.3.v3[i,]=bayesmooth(sim.m.bump.3.v3[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.bump.3.v4[i,]=bayesmooth(sim.m.bump.3.v4[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.bump.3.v5[i,]=bayesmooth(sim.m.bump.3.v5[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.blk.1.v1[i,]=bayesmooth(sim.m.blk.1.v1[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.blk.1.v2[i,]=bayesmooth(sim.m.blk.1.v2[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.blk.1.v3[i,]=bayesmooth(sim.m.blk.1.v3[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.blk.1.v4[i,]=bayesmooth(sim.m.blk.1.v4[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.blk.1.v5[i,]=bayesmooth(sim.m.blk.1.v5[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.blk.3.v1[i,]=bayesmooth(sim.m.blk.3.v1[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.blk.3.v2[i,]=bayesmooth(sim.m.blk.3.v2[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.blk.3.v3[i,]=bayesmooth(sim.m.blk.3.v3[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.blk.3.v4[i,]=bayesmooth(sim.m.blk.3.v4[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.blk.3.v5[i,]=bayesmooth(sim.m.blk.3.v5[i,],filter.number=1,family="DaubExPhase",gridmult=2)


  mu.est.ashef.haar.ang.1.v1[i,]=bayesmooth(sim.m.ang.1.v1[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.ang.1.v2[i,]=bayesmooth(sim.m.ang.1.v2[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.ang.1.v3[i,]=bayesmooth(sim.m.ang.1.v3[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.ang.1.v4[i,]=bayesmooth(sim.m.ang.1.v4[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.ang.1.v5[i,]=bayesmooth(sim.m.ang.1.v5[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.ang.3.v1[i,]=bayesmooth(sim.m.ang.3.v1[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.ang.3.v2[i,]=bayesmooth(sim.m.ang.3.v2[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.ang.3.v3[i,]=bayesmooth(sim.m.ang.3.v3[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.ang.3.v4[i,]=bayesmooth(sim.m.ang.3.v4[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.ang.3.v5[i,]=bayesmooth(sim.m.ang.3.v5[i,],filter.number=1,family="DaubExPhase",gridmult=2)


  mu.est.ashef.haar.dop.1.v1[i,]=bayesmooth(sim.m.dop.1.v1[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.dop.1.v2[i,]=bayesmooth(sim.m.dop.1.v2[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.dop.1.v3[i,]=bayesmooth(sim.m.dop.1.v3[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.dop.1.v4[i,]=bayesmooth(sim.m.dop.1.v4[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.dop.1.v5[i,]=bayesmooth(sim.m.dop.1.v5[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.dop.3.v1[i,]=bayesmooth(sim.m.dop.3.v1[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.dop.3.v2[i,]=bayesmooth(sim.m.dop.3.v2[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.dop.3.v3[i,]=bayesmooth(sim.m.dop.3.v3[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.dop.3.v4[i,]=bayesmooth(sim.m.dop.3.v4[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.dop.3.v5[i,]=bayesmooth(sim.m.dop.3.v5[i,],filter.number=1,family="DaubExPhase",gridmult=2)


  mu.est.ashef.haar.blip.1.v1[i,]=bayesmooth(sim.m.blip.1.v1[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.blip.1.v2[i,]=bayesmooth(sim.m.blip.1.v2[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.blip.1.v3[i,]=bayesmooth(sim.m.blip.1.v3[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.blip.1.v4[i,]=bayesmooth(sim.m.blip.1.v4[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.blip.1.v5[i,]=bayesmooth(sim.m.blip.1.v5[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.blip.3.v1[i,]=bayesmooth(sim.m.blip.3.v1[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.blip.3.v2[i,]=bayesmooth(sim.m.blip.3.v2[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.blip.3.v3[i,]=bayesmooth(sim.m.blip.3.v3[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.blip.3.v4[i,]=bayesmooth(sim.m.blip.3.v4[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.blip.3.v5[i,]=bayesmooth(sim.m.blip.3.v5[i,],filter.number=1,family="DaubExPhase",gridmult=2)



  mu.est.ashef.haar.cor.1.v1[i,]=bayesmooth(sim.m.cor.1.v1[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.cor.1.v2[i,]=bayesmooth(sim.m.cor.1.v2[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.cor.1.v3[i,]=bayesmooth(sim.m.cor.1.v3[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.cor.1.v4[i,]=bayesmooth(sim.m.cor.1.v4[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.cor.1.v5[i,]=bayesmooth(sim.m.cor.1.v5[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.cor.3.v1[i,]=bayesmooth(sim.m.cor.3.v1[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.cor.3.v2[i,]=bayesmooth(sim.m.cor.3.v2[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.cor.3.v3[i,]=bayesmooth(sim.m.cor.3.v3[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.cor.3.v4[i,]=bayesmooth(sim.m.cor.3.v4[i,],filter.number=1,family="DaubExPhase",gridmult=2)

  mu.est.ashef.haar.cor.3.v5[i,]=bayesmooth(sim.m.cor.3.v5[i,],filter.number=1,family="DaubExPhase",gridmult=2)


  print(i)
}


save.image("res_sim.RData")