source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/ashsmooth.R"))
library(wavethresh)
library(EbayesThresh)


n=1024
t=1:n/n


spike.f=function(x) (0.75*exp(-500*(x-0.23)^2)+1.5*exp(-2000*(x-0.33)^2)+3*exp(-8000*(x-0.47)^2)+2.25*exp(-16000*(x-0.69)^2)+0.5*exp(-32000*(x-0.83)^2))
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




#save data to import into matlab
write(t(sim.m.sp.1.v1),"sim_sp_1_v1.txt",ncolumns=n)
write(t(sim.m.sp.1.v2),"sim_sp_1_v2.txt",ncolumns=n)
write(t(sim.m.sp.1.v3),"sim_sp_1_v3.txt",ncolumns=n)
write(t(sim.m.sp.1.v4),"sim_sp_1_v4.txt",ncolumns=n)
write(t(sim.m.sp.1.v5),"sim_sp_1_v5.txt",ncolumns=n)

write(t(sim.m.sp.3.v1),"sim_sp_3_v1.txt",ncolumns=n)
write(t(sim.m.sp.3.v2),"sim_sp_3_v2.txt",ncolumns=n)
write(t(sim.m.sp.3.v3),"sim_sp_3_v3.txt",ncolumns=n)
write(t(sim.m.sp.3.v4),"sim_sp_3_v4.txt",ncolumns=n)
write(t(sim.m.sp.3.v5),"sim_sp_3_v5.txt",ncolumns=n)


write(t(sim.m.bump.1.v1),"sim_bump_1_v1.txt",ncolumns=n)
write(t(sim.m.bump.1.v2),"sim_bump_1_v2.txt",ncolumns=n)
write(t(sim.m.bump.1.v3),"sim_bump_1_v3.txt",ncolumns=n)
write(t(sim.m.bump.1.v4),"sim_bump_1_v4.txt",ncolumns=n)
write(t(sim.m.bump.1.v5),"sim_bump_1_v5.txt",ncolumns=n)

write(t(sim.m.bump.3.v1),"sim_bump_3_v1.txt",ncolumns=n)
write(t(sim.m.bump.3.v2),"sim_bump_3_v2.txt",ncolumns=n)
write(t(sim.m.bump.3.v3),"sim_bump_3_v3.txt",ncolumns=n)
write(t(sim.m.bump.3.v4),"sim_bump_3_v4.txt",ncolumns=n)
write(t(sim.m.bump.3.v5),"sim_bump_3_v5.txt",ncolumns=n)


write(t(sim.m.blk.1.v1),"sim_blk_1_v1.txt",ncolumns=n)
write(t(sim.m.blk.1.v2),"sim_blk_1_v2.txt",ncolumns=n)
write(t(sim.m.blk.1.v3),"sim_blk_1_v3.txt",ncolumns=n)
write(t(sim.m.blk.1.v4),"sim_blk_1_v4.txt",ncolumns=n)
write(t(sim.m.blk.1.v5),"sim_blk_1_v5.txt",ncolumns=n)

write(t(sim.m.blk.3.v1),"sim_blk_3_v1.txt",ncolumns=n)
write(t(sim.m.blk.3.v2),"sim_blk_3_v2.txt",ncolumns=n)
write(t(sim.m.blk.3.v3),"sim_blk_3_v3.txt",ncolumns=n)
write(t(sim.m.blk.3.v4),"sim_blk_3_v4.txt",ncolumns=n)
write(t(sim.m.blk.3.v5),"sim_blk_3_v5.txt",ncolumns=n)



write(t(sim.m.ang.1.v1),"sim_ang_1_v1.txt",ncolumns=n)
write(t(sim.m.ang.1.v2),"sim_ang_1_v2.txt",ncolumns=n)
write(t(sim.m.ang.1.v3),"sim_ang_1_v3.txt",ncolumns=n)
write(t(sim.m.ang.1.v4),"sim_ang_1_v4.txt",ncolumns=n)
write(t(sim.m.ang.1.v5),"sim_ang_1_v5.txt",ncolumns=n)

write(t(sim.m.ang.3.v1),"sim_ang_3_v1.txt",ncolumns=n)
write(t(sim.m.ang.3.v2),"sim_ang_3_v2.txt",ncolumns=n)
write(t(sim.m.ang.3.v3),"sim_ang_3_v3.txt",ncolumns=n)
write(t(sim.m.ang.3.v4),"sim_ang_3_v4.txt",ncolumns=n)
write(t(sim.m.ang.3.v5),"sim_ang_3_v5.txt",ncolumns=n)



write(t(sim.m.dop.1.v1),"sim_dop_1_v1.txt",ncolumns=n)
write(t(sim.m.dop.1.v2),"sim_dop_1_v2.txt",ncolumns=n)
write(t(sim.m.dop.1.v3),"sim_dop_1_v3.txt",ncolumns=n)
write(t(sim.m.dop.1.v4),"sim_dop_1_v4.txt",ncolumns=n)
write(t(sim.m.dop.1.v5),"sim_dop_1_v5.txt",ncolumns=n)

write(t(sim.m.dop.3.v1),"sim_dop_3_v1.txt",ncolumns=n)
write(t(sim.m.dop.3.v2),"sim_dop_3_v2.txt",ncolumns=n)
write(t(sim.m.dop.3.v3),"sim_dop_3_v3.txt",ncolumns=n)
write(t(sim.m.dop.3.v4),"sim_dop_3_v4.txt",ncolumns=n)
write(t(sim.m.dop.3.v5),"sim_dop_3_v5.txt",ncolumns=n)



write(t(sim.m.blip.1.v1),"sim_blip_1_v1.txt",ncolumns=n)
write(t(sim.m.blip.1.v2),"sim_blip_1_v2.txt",ncolumns=n)
write(t(sim.m.blip.1.v3),"sim_blip_1_v3.txt",ncolumns=n)
write(t(sim.m.blip.1.v4),"sim_blip_1_v4.txt",ncolumns=n)
write(t(sim.m.blip.1.v5),"sim_blip_1_v5.txt",ncolumns=n)

write(t(sim.m.blip.3.v1),"sim_blip_3_v1.txt",ncolumns=n)
write(t(sim.m.blip.3.v2),"sim_blip_3_v2.txt",ncolumns=n)
write(t(sim.m.blip.3.v3),"sim_blip_3_v3.txt",ncolumns=n)
write(t(sim.m.blip.3.v4),"sim_blip_3_v4.txt",ncolumns=n)
write(t(sim.m.blip.3.v5),"sim_blip_3_v5.txt",ncolumns=n)



write(t(sim.m.cor.1.v1),"sim_cor_1_v1.txt",ncolumns=n)
write(t(sim.m.cor.1.v2),"sim_cor_1_v2.txt",ncolumns=n)
write(t(sim.m.cor.1.v3),"sim_cor_1_v3.txt",ncolumns=n)
write(t(sim.m.cor.1.v4),"sim_cor_1_v4.txt",ncolumns=n)
write(t(sim.m.cor.1.v5),"sim_cor_1_v5.txt",ncolumns=n)

write(t(sim.m.cor.3.v1),"sim_cor_3_v1.txt",ncolumns=n)
write(t(sim.m.cor.3.v2),"sim_cor_3_v2.txt",ncolumns=n)
write(t(sim.m.cor.3.v3),"sim_cor_3_v3.txt",ncolumns=n)
write(t(sim.m.cor.3.v4),"sim_cor_3_v4.txt",ncolumns=n)
write(t(sim.m.cor.3.v5),"sim_cor_3_v5.txt",ncolumns=n)




pdf("signal_1d_gaus.pdf",height=14,width=14)
par(mfrow=c(4,2))
plot(mu.sp,ylab="Intesnity",main="Spikes",type='l')
plot(mu.bump,ylab="Intesnity",main="Bumps",type='l')
plot(mu.blk,ylab="Intesnity",main="Blocks",type='l')
plot(mu.ang,ylab="Intesnity",main="Angles",type='l')
plot(mu.dop,ylab="Intesnity",main="Doppler",type='l')
plot(mu.blip,ylab="Intesnity",main="Blip",type='l')
plot(mu.cor,ylab="Intesnity",main="Corner",type='l')
dev.off()


pdf("var_1d_gaus.pdf",height=14,width=14)
par(mfrow=c(3,2))
plot(sigma.sp.3.v1^2,ylab="Intesnity",main="V1",type='l')
plot(sigma.sp.3.v2^2,ylab="Intesnity",main="V2",type='l')
plot(sigma.sp.3.v3^2,ylab="Intesnity",main="v3",type='l')
plot(sigma.sp.3.v4^2,ylab="Intesnity",main="V4",type='l')
plot(sigma.sp.3.v5^2,ylab="Intesnity",main="V5",type='l')
dev.off()




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



for(i in 1:100){
  sigma.sp.1.v1.est=sig.est.func(sim.m.sp.1.v1[i,],n)
  mu.est.ash.haar.sp.1.v1[i,]=ashsmooth.gaus(sim.m.sp.1.v1[i,],sigma.sp.1.v1.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.sp.1.v1[i,]=ashsmooth.gaus(sim.m.sp.1.v1[i,],sigma.sp.1.v1.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.sp.1.v1[i,]=ashsmooth.gaus(sim.m.sp.1.v1[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.sp.1.v1[i,]=ashsmooth.gaus(sim.m.sp.1.v1[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.sp.1.v1[i,]=ti.thresh(sim.m.sp.1.v1[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.sp.1.v1[i,]=ti.thresh(sim.m.sp.1.v1[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.sp.1.v1[i,]=ashsmooth.gaus(sim.m.sp.1.v1[i,],sigma=sigma.sp.1.v1,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.sp.1.v1[i,]=ashsmooth.gaus(sim.m.sp.1.v1[i,],sigma=sigma.sp.1.v1,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.sp.1.v1[i,]=ti.thresh(sim.m.sp.1.v1[i,],sigma=sigma.sp.1.v1,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.sp.1.v1[i,]=ti.thresh(sim.m.sp.1.v1[i,],sigma=sigma.sp.1.v1,filter.number=8,family="DaubLeAsymm")
  sigma.sp.3.v1.est=sig.est.func(sim.m.sp.3.v1[i,],n)
  mu.est.ash.haar.sp.3.v1[i,]=ashsmooth.gaus(sim.m.sp.3.v1[i,],sigma.sp.3.v1.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.sp.3.v1[i,]=ashsmooth.gaus(sim.m.sp.3.v1[i,],sigma.sp.3.v1.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.sp.3.v1[i,]=ashsmooth.gaus(sim.m.sp.3.v1[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.sp.3.v1[i,]=ashsmooth.gaus(sim.m.sp.3.v1[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.sp.3.v1[i,]=ti.thresh(sim.m.sp.3.v1[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.sp.3.v1[i,]=ti.thresh(sim.m.sp.3.v1[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.sp.3.v1[i,]=ashsmooth.gaus(sim.m.sp.3.v1[i,],sigma=sigma.sp.3.v1,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.sp.3.v1[i,]=ashsmooth.gaus(sim.m.sp.3.v1[i,],sigma=sigma.sp.3.v1,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.sp.3.v1[i,]=ti.thresh(sim.m.sp.3.v1[i,],sigma=sigma.sp.3.v1,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.sp.3.v1[i,]=ti.thresh(sim.m.sp.3.v1[i,],sigma=sigma.sp.3.v1,filter.number=8,family="DaubLeAsymm")

  sigma.bump.1.v1.est=sig.est.func(sim.m.bump.1.v1[i,],n)
  mu.est.ash.haar.bump.1.v1[i,]=ashsmooth.gaus(sim.m.bump.1.v1[i,],sigma.bump.1.v1.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.bump.1.v1[i,]=ashsmooth.gaus(sim.m.bump.1.v1[i,],sigma.bump.1.v1.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.bump.1.v1[i,]=ashsmooth.gaus(sim.m.bump.1.v1[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.bump.1.v1[i,]=ashsmooth.gaus(sim.m.bump.1.v1[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.bump.1.v1[i,]=ti.thresh(sim.m.bump.1.v1[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.bump.1.v1[i,]=ti.thresh(sim.m.bump.1.v1[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.bump.1.v1[i,]=ashsmooth.gaus(sim.m.bump.1.v1[i,],sigma=sigma.bump.1.v1,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.bump.1.v1[i,]=ashsmooth.gaus(sim.m.bump.1.v1[i,],sigma=sigma.bump.1.v1,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.bump.1.v1[i,]=ti.thresh(sim.m.bump.1.v1[i,],sigma=sigma.bump.1.v1,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.bump.1.v1[i,]=ti.thresh(sim.m.bump.1.v1[i,],sigma=sigma.bump.1.v1,filter.number=8,family="DaubLeAsymm")
  sigma.bump.3.v1.est=sig.est.func(sim.m.bump.3.v1[i,],n)
  mu.est.ash.haar.bump.3.v1[i,]=ashsmooth.gaus(sim.m.bump.3.v1[i,],sigma.bump.3.v1.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.bump.3.v1[i,]=ashsmooth.gaus(sim.m.bump.3.v1[i,],sigma.bump.3.v1.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.bump.3.v1[i,]=ashsmooth.gaus(sim.m.bump.3.v1[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.bump.3.v1[i,]=ashsmooth.gaus(sim.m.bump.3.v1[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.bump.3.v1[i,]=ti.thresh(sim.m.bump.3.v1[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.bump.3.v1[i,]=ti.thresh(sim.m.bump.3.v1[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.bump.3.v1[i,]=ashsmooth.gaus(sim.m.bump.3.v1[i,],sigma=sigma.bump.3.v1,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.bump.3.v1[i,]=ashsmooth.gaus(sim.m.bump.3.v1[i,],sigma=sigma.bump.3.v1,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.bump.3.v1[i,]=ti.thresh(sim.m.bump.3.v1[i,],sigma=sigma.bump.3.v1,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.bump.3.v1[i,]=ti.thresh(sim.m.bump.3.v1[i,],sigma=sigma.bump.3.v1,filter.number=8,family="DaubLeAsymm")
  
  sigma.blk.1.v1.est=sig.est.func(sim.m.blk.1.v1[i,],n)
  mu.est.ash.haar.blk.1.v1[i,]=ashsmooth.gaus(sim.m.blk.1.v1[i,],sigma.blk.1.v1.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.blk.1.v1[i,]=ashsmooth.gaus(sim.m.blk.1.v1[i,],sigma.blk.1.v1.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.blk.1.v1[i,]=ashsmooth.gaus(sim.m.blk.1.v1[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.blk.1.v1[i,]=ashsmooth.gaus(sim.m.blk.1.v1[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.blk.1.v1[i,]=ti.thresh(sim.m.blk.1.v1[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.blk.1.v1[i,]=ti.thresh(sim.m.blk.1.v1[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.blk.1.v1[i,]=ashsmooth.gaus(sim.m.blk.1.v1[i,],sigma=sigma.blk.1.v1,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.blk.1.v1[i,]=ashsmooth.gaus(sim.m.blk.1.v1[i,],sigma=sigma.blk.1.v1,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.blk.1.v1[i,]=ti.thresh(sim.m.blk.1.v1[i,],sigma=sigma.blk.1.v1,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.blk.1.v1[i,]=ti.thresh(sim.m.blk.1.v1[i,],sigma=sigma.blk.1.v1,filter.number=8,family="DaubLeAsymm")
  sigma.blk.3.v1.est=sig.est.func(sim.m.blk.3.v1[i,],n)
  mu.est.ash.haar.blk.3.v1[i,]=ashsmooth.gaus(sim.m.blk.3.v1[i,],sigma.blk.3.v1.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.blk.3.v1[i,]=ashsmooth.gaus(sim.m.blk.3.v1[i,],sigma.blk.3.v1.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.blk.3.v1[i,]=ashsmooth.gaus(sim.m.blk.3.v1[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.blk.3.v1[i,]=ashsmooth.gaus(sim.m.blk.3.v1[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.blk.3.v1[i,]=ti.thresh(sim.m.blk.3.v1[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.blk.3.v1[i,]=ti.thresh(sim.m.blk.3.v1[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.blk.3.v1[i,]=ashsmooth.gaus(sim.m.blk.3.v1[i,],sigma=sigma.blk.3.v1,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.blk.3.v1[i,]=ashsmooth.gaus(sim.m.blk.3.v1[i,],sigma=sigma.blk.3.v1,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.blk.3.v1[i,]=ti.thresh(sim.m.blk.3.v1[i,],sigma=sigma.blk.3.v1,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.blk.3.v1[i,]=ti.thresh(sim.m.blk.3.v1[i,],sigma=sigma.blk.3.v1,filter.number=8,family="DaubLeAsymm")

  sigma.ang.1.v1.est=sig.est.func(sim.m.ang.1.v1[i,],n)
  mu.est.ash.haar.ang.1.v1[i,]=ashsmooth.gaus(sim.m.ang.1.v1[i,],sigma.ang.1.v1.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.ang.1.v1[i,]=ashsmooth.gaus(sim.m.ang.1.v1[i,],sigma.ang.1.v1.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.ang.1.v1[i,]=ashsmooth.gaus(sim.m.ang.1.v1[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.ang.1.v1[i,]=ashsmooth.gaus(sim.m.ang.1.v1[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.ang.1.v1[i,]=ti.thresh(sim.m.ang.1.v1[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.ang.1.v1[i,]=ti.thresh(sim.m.ang.1.v1[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.ang.1.v1[i,]=ashsmooth.gaus(sim.m.ang.1.v1[i,],sigma=sigma.ang.1.v1,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.ang.1.v1[i,]=ashsmooth.gaus(sim.m.ang.1.v1[i,],sigma=sigma.ang.1.v1,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.ang.1.v1[i,]=ti.thresh(sim.m.ang.1.v1[i,],sigma=sigma.ang.1.v1,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.ang.1.v1[i,]=ti.thresh(sim.m.ang.1.v1[i,],sigma=sigma.ang.1.v1,filter.number=8,family="DaubLeAsymm")
  sigma.ang.3.v1.est=sig.est.func(sim.m.ang.3.v1[i,],n)
  mu.est.ash.haar.ang.3.v1[i,]=ashsmooth.gaus(sim.m.ang.3.v1[i,],sigma.ang.3.v1.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.ang.3.v1[i,]=ashsmooth.gaus(sim.m.ang.3.v1[i,],sigma.ang.3.v1.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.ang.3.v1[i,]=ashsmooth.gaus(sim.m.ang.3.v1[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.ang.3.v1[i,]=ashsmooth.gaus(sim.m.ang.3.v1[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.ang.3.v1[i,]=ti.thresh(sim.m.ang.3.v1[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.ang.3.v1[i,]=ti.thresh(sim.m.ang.3.v1[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.ang.3.v1[i,]=ashsmooth.gaus(sim.m.ang.3.v1[i,],sigma=sigma.ang.3.v1,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.ang.3.v1[i,]=ashsmooth.gaus(sim.m.ang.3.v1[i,],sigma=sigma.ang.3.v1,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.ang.3.v1[i,]=ti.thresh(sim.m.ang.3.v1[i,],sigma=sigma.ang.3.v1,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.ang.3.v1[i,]=ti.thresh(sim.m.ang.3.v1[i,],sigma=sigma.ang.3.v1,filter.number=8,family="DaubLeAsymm")

  sigma.dop.1.v1.est=sig.est.func(sim.m.dop.1.v1[i,],n)
  mu.est.ash.haar.dop.1.v1[i,]=ashsmooth.gaus(sim.m.dop.1.v1[i,],sigma.dop.1.v1.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.dop.1.v1[i,]=ashsmooth.gaus(sim.m.dop.1.v1[i,],sigma.dop.1.v1.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.dop.1.v1[i,]=ashsmooth.gaus(sim.m.dop.1.v1[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.dop.1.v1[i,]=ashsmooth.gaus(sim.m.dop.1.v1[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.dop.1.v1[i,]=ti.thresh(sim.m.dop.1.v1[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.dop.1.v1[i,]=ti.thresh(sim.m.dop.1.v1[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.dop.1.v1[i,]=ashsmooth.gaus(sim.m.dop.1.v1[i,],sigma=sigma.dop.1.v1,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.dop.1.v1[i,]=ashsmooth.gaus(sim.m.dop.1.v1[i,],sigma=sigma.dop.1.v1,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.dop.1.v1[i,]=ti.thresh(sim.m.dop.1.v1[i,],sigma=sigma.dop.1.v1,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.dop.1.v1[i,]=ti.thresh(sim.m.dop.1.v1[i,],sigma=sigma.dop.1.v1,filter.number=8,family="DaubLeAsymm")
  sigma.dop.3.v1.est=sig.est.func(sim.m.dop.3.v1[i,],n)
  mu.est.ash.haar.dop.3.v1[i,]=ashsmooth.gaus(sim.m.dop.3.v1[i,],sigma.dop.3.v1.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.dop.3.v1[i,]=ashsmooth.gaus(sim.m.dop.3.v1[i,],sigma.dop.3.v1.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.dop.3.v1[i,]=ashsmooth.gaus(sim.m.dop.3.v1[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.dop.3.v1[i,]=ashsmooth.gaus(sim.m.dop.3.v1[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.dop.3.v1[i,]=ti.thresh(sim.m.dop.3.v1[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.dop.3.v1[i,]=ti.thresh(sim.m.dop.3.v1[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.dop.3.v1[i,]=ashsmooth.gaus(sim.m.dop.3.v1[i,],sigma=sigma.dop.3.v1,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.dop.3.v1[i,]=ashsmooth.gaus(sim.m.dop.3.v1[i,],sigma=sigma.dop.3.v1,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.dop.3.v1[i,]=ti.thresh(sim.m.dop.3.v1[i,],sigma=sigma.dop.3.v1,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.dop.3.v1[i,]=ti.thresh(sim.m.dop.3.v1[i,],sigma=sigma.dop.3.v1,filter.number=8,family="DaubLeAsymm")

  sigma.blip.1.v1.est=sig.est.func(sim.m.blip.1.v1[i,],n)
  mu.est.ash.haar.blip.1.v1[i,]=ashsmooth.gaus(sim.m.blip.1.v1[i,],sigma.blip.1.v1.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.blip.1.v1[i,]=ashsmooth.gaus(sim.m.blip.1.v1[i,],sigma.blip.1.v1.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.blip.1.v1[i,]=ashsmooth.gaus(sim.m.blip.1.v1[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.blip.1.v1[i,]=ashsmooth.gaus(sim.m.blip.1.v1[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.blip.1.v1[i,]=ti.thresh(sim.m.blip.1.v1[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.blip.1.v1[i,]=ti.thresh(sim.m.blip.1.v1[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.blip.1.v1[i,]=ashsmooth.gaus(sim.m.blip.1.v1[i,],sigma=sigma.blip.1.v1,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.blip.1.v1[i,]=ashsmooth.gaus(sim.m.blip.1.v1[i,],sigma=sigma.blip.1.v1,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.blip.1.v1[i,]=ti.thresh(sim.m.blip.1.v1[i,],sigma=sigma.blip.1.v1,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.blip.1.v1[i,]=ti.thresh(sim.m.blip.1.v1[i,],sigma=sigma.blip.1.v1,filter.number=8,family="DaubLeAsymm")
  sigma.blip.3.v1.est=sig.est.func(sim.m.blip.3.v1[i,],n)
  mu.est.ash.haar.blip.3.v1[i,]=ashsmooth.gaus(sim.m.blip.3.v1[i,],sigma.blip.3.v1.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.blip.3.v1[i,]=ashsmooth.gaus(sim.m.blip.3.v1[i,],sigma.blip.3.v1.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.blip.3.v1[i,]=ashsmooth.gaus(sim.m.blip.3.v1[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.blip.3.v1[i,]=ashsmooth.gaus(sim.m.blip.3.v1[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.blip.3.v1[i,]=ti.thresh(sim.m.blip.3.v1[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.blip.3.v1[i,]=ti.thresh(sim.m.blip.3.v1[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.blip.3.v1[i,]=ashsmooth.gaus(sim.m.blip.3.v1[i,],sigma=sigma.blip.3.v1,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.blip.3.v1[i,]=ashsmooth.gaus(sim.m.blip.3.v1[i,],sigma=sigma.blip.3.v1,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.blip.3.v1[i,]=ti.thresh(sim.m.blip.3.v1[i,],sigma=sigma.blip.3.v1,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.blip.3.v1[i,]=ti.thresh(sim.m.blip.3.v1[i,],sigma=sigma.blip.3.v1,filter.number=8,family="DaubLeAsymm")

  sigma.cor.1.v1.est=sig.est.func(sim.m.cor.1.v1[i,],n)
  mu.est.ash.haar.cor.1.v1[i,]=ashsmooth.gaus(sim.m.cor.1.v1[i,],sigma.cor.1.v1.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.cor.1.v1[i,]=ashsmooth.gaus(sim.m.cor.1.v1[i,],sigma.cor.1.v1.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.cor.1.v1[i,]=ashsmooth.gaus(sim.m.cor.1.v1[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.cor.1.v1[i,]=ashsmooth.gaus(sim.m.cor.1.v1[i,],filter.number=8,family="DaubLeAsymm",)
  mu.est.tie.haar.cor.1.v1[i,]=ti.thresh(sim.m.cor.1.v1[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.cor.1.v1[i,]=ti.thresh(sim.m.cor.1.v1[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.cor.1.v1[i,]=ashsmooth.gaus(sim.m.cor.1.v1[i,],sigma=sigma.cor.1.v1,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.cor.1.v1[i,]=ashsmooth.gaus(sim.m.cor.1.v1[i,],sigma=sigma.cor.1.v1,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.cor.1.v1[i,]=ti.thresh(sim.m.cor.1.v1[i,],sigma=sigma.cor.1.v1,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.cor.1.v1[i,]=ti.thresh(sim.m.cor.1.v1[i,],sigma=sigma.cor.1.v1,filter.number=8,family="DaubLeAsymm")
  sigma.cor.3.v1.est=sig.est.func(sim.m.cor.3.v1[i,],n)
  mu.est.ash.haar.cor.3.v1[i,]=ashsmooth.gaus(sim.m.cor.3.v1[i,],sigma.cor.3.v1.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.cor.3.v1[i,]=ashsmooth.gaus(sim.m.cor.3.v1[i,],sigma.cor.3.v1.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.cor.3.v1[i,]=ashsmooth.gaus(sim.m.cor.3.v1[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.cor.3.v1[i,]=ashsmooth.gaus(sim.m.cor.3.v1[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.cor.3.v1[i,]=ti.thresh(sim.m.cor.3.v1[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.cor.3.v1[i,]=ti.thresh(sim.m.cor.3.v1[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.cor.3.v1[i,]=ashsmooth.gaus(sim.m.cor.3.v1[i,],sigma=sigma.cor.3.v1,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.cor.3.v1[i,]=ashsmooth.gaus(sim.m.cor.3.v1[i,],sigma=sigma.cor.3.v1,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.cor.3.v1[i,]=ti.thresh(sim.m.cor.3.v1[i,],sigma=sigma.cor.3.v1,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.cor.3.v1[i,]=ti.thresh(sim.m.cor.3.v1[i,],sigma=sigma.cor.3.v1,filter.number=8,family="DaubLeAsymm")



  sigma.sp.1.v2.est=sig.est.func(sim.m.sp.1.v2[i,],n)
  mu.est.ash.haar.sp.1.v2[i,]=ashsmooth.gaus(sim.m.sp.1.v2[i,],sigma.sp.1.v2.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.sp.1.v2[i,]=ashsmooth.gaus(sim.m.sp.1.v2[i,],sigma.sp.1.v2.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.sp.1.v2[i,]=ashsmooth.gaus(sim.m.sp.1.v2[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.sp.1.v2[i,]=ashsmooth.gaus(sim.m.sp.1.v2[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.sp.1.v2[i,]=ti.thresh(sim.m.sp.1.v2[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.sp.1.v2[i,]=ti.thresh(sim.m.sp.1.v2[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.sp.1.v2[i,]=ashsmooth.gaus(sim.m.sp.1.v2[i,],sigma=sigma.sp.1.v2,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.sp.1.v2[i,]=ashsmooth.gaus(sim.m.sp.1.v2[i,],sigma=sigma.sp.1.v2,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.sp.1.v2[i,]=ti.thresh(sim.m.sp.1.v2[i,],sigma=sigma.sp.1.v2,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.sp.1.v2[i,]=ti.thresh(sim.m.sp.1.v2[i,],sigma=sigma.sp.1.v2,filter.number=8,family="DaubLeAsymm")
  sigma.sp.3.v2.est=sig.est.func(sim.m.sp.3.v2[i,],n)
  mu.est.ash.haar.sp.3.v2[i,]=ashsmooth.gaus(sim.m.sp.3.v2[i,],sigma.sp.3.v2.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.sp.3.v2[i,]=ashsmooth.gaus(sim.m.sp.3.v2[i,],sigma.sp.3.v2.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.sp.3.v2[i,]=ashsmooth.gaus(sim.m.sp.3.v2[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.sp.3.v2[i,]=ashsmooth.gaus(sim.m.sp.3.v2[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.sp.3.v2[i,]=ti.thresh(sim.m.sp.3.v2[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.sp.3.v2[i,]=ti.thresh(sim.m.sp.3.v2[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.sp.3.v2[i,]=ashsmooth.gaus(sim.m.sp.3.v2[i,],sigma=sigma.sp.3.v2,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.sp.3.v2[i,]=ashsmooth.gaus(sim.m.sp.3.v2[i,],sigma=sigma.sp.3.v2,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.sp.3.v2[i,]=ti.thresh(sim.m.sp.3.v2[i,],sigma=sigma.sp.3.v2,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.sp.3.v2[i,]=ti.thresh(sim.m.sp.3.v2[i,],sigma=sigma.sp.3.v2,filter.number=8,family="DaubLeAsymm")

  sigma.bump.1.v2.est=sig.est.func(sim.m.bump.1.v2[i,],n)
  mu.est.ash.haar.bump.1.v2[i,]=ashsmooth.gaus(sim.m.bump.1.v2[i,],sigma.bump.1.v2.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.bump.1.v2[i,]=ashsmooth.gaus(sim.m.bump.1.v2[i,],sigma.bump.1.v2.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.bump.1.v2[i,]=ashsmooth.gaus(sim.m.bump.1.v2[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.bump.1.v2[i,]=ashsmooth.gaus(sim.m.bump.1.v2[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.bump.1.v2[i,]=ti.thresh(sim.m.bump.1.v2[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.bump.1.v2[i,]=ti.thresh(sim.m.bump.1.v2[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.bump.1.v2[i,]=ashsmooth.gaus(sim.m.bump.1.v2[i,],sigma=sigma.bump.1.v2,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.bump.1.v2[i,]=ashsmooth.gaus(sim.m.bump.1.v2[i,],sigma=sigma.bump.1.v2,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.bump.1.v2[i,]=ti.thresh(sim.m.bump.1.v2[i,],sigma=sigma.bump.1.v2,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.bump.1.v2[i,]=ti.thresh(sim.m.bump.1.v2[i,],sigma=sigma.bump.1.v2,filter.number=8,family="DaubLeAsymm")
  sigma.bump.3.v2.est=sig.est.func(sim.m.bump.3.v2[i,],n)
  mu.est.ash.haar.bump.3.v2[i,]=ashsmooth.gaus(sim.m.bump.3.v2[i,],sigma.bump.3.v2.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.bump.3.v2[i,]=ashsmooth.gaus(sim.m.bump.3.v2[i,],sigma.bump.3.v2.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.bump.3.v2[i,]=ashsmooth.gaus(sim.m.bump.3.v2[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.bump.3.v2[i,]=ashsmooth.gaus(sim.m.bump.3.v2[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.bump.3.v2[i,]=ti.thresh(sim.m.bump.3.v2[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.bump.3.v2[i,]=ti.thresh(sim.m.bump.3.v2[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.bump.3.v2[i,]=ashsmooth.gaus(sim.m.bump.3.v2[i,],sigma=sigma.bump.3.v2,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.bump.3.v2[i,]=ashsmooth.gaus(sim.m.bump.3.v2[i,],sigma=sigma.bump.3.v2,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.bump.3.v2[i,]=ti.thresh(sim.m.bump.3.v2[i,],sigma=sigma.bump.3.v2,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.bump.3.v2[i,]=ti.thresh(sim.m.bump.3.v2[i,],sigma=sigma.bump.3.v2,filter.number=8,family="DaubLeAsymm")
  
  sigma.blk.1.v2.est=sig.est.func(sim.m.blk.1.v2[i,],n)
  mu.est.ash.haar.blk.1.v2[i,]=ashsmooth.gaus(sim.m.blk.1.v2[i,],sigma.blk.1.v2.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.blk.1.v2[i,]=ashsmooth.gaus(sim.m.blk.1.v2[i,],sigma.blk.1.v2.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.blk.1.v2[i,]=ashsmooth.gaus(sim.m.blk.1.v2[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.blk.1.v2[i,]=ashsmooth.gaus(sim.m.blk.1.v2[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.blk.1.v2[i,]=ti.thresh(sim.m.blk.1.v2[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.blk.1.v2[i,]=ti.thresh(sim.m.blk.1.v2[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.blk.1.v2[i,]=ashsmooth.gaus(sim.m.blk.1.v2[i,],sigma=sigma.blk.1.v2,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.blk.1.v2[i,]=ashsmooth.gaus(sim.m.blk.1.v2[i,],sigma=sigma.blk.1.v2,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.blk.1.v2[i,]=ti.thresh(sim.m.blk.1.v2[i,],sigma=sigma.blk.1.v2,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.blk.1.v2[i,]=ti.thresh(sim.m.blk.1.v2[i,],sigma=sigma.blk.1.v2,filter.number=8,family="DaubLeAsymm")
  sigma.blk.3.v2.est=sig.est.func(sim.m.blk.3.v2[i,],n)
  mu.est.ash.haar.blk.3.v2[i,]=ashsmooth.gaus(sim.m.blk.3.v2[i,],sigma.blk.3.v2.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.blk.3.v2[i,]=ashsmooth.gaus(sim.m.blk.3.v2[i,],sigma.blk.3.v2.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.blk.3.v2[i,]=ashsmooth.gaus(sim.m.blk.3.v2[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.blk.3.v2[i,]=ashsmooth.gaus(sim.m.blk.3.v2[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.blk.3.v2[i,]=ti.thresh(sim.m.blk.3.v2[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.blk.3.v2[i,]=ti.thresh(sim.m.blk.3.v2[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.blk.3.v2[i,]=ashsmooth.gaus(sim.m.blk.3.v2[i,],sigma=sigma.blk.3.v2,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.blk.3.v2[i,]=ashsmooth.gaus(sim.m.blk.3.v2[i,],sigma=sigma.blk.3.v2,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.blk.3.v2[i,]=ti.thresh(sim.m.blk.3.v2[i,],sigma=sigma.blk.3.v2,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.blk.3.v2[i,]=ti.thresh(sim.m.blk.3.v2[i,],sigma=sigma.blk.3.v2,filter.number=8,family="DaubLeAsymm")

  sigma.ang.1.v2.est=sig.est.func(sim.m.ang.1.v2[i,],n)
  mu.est.ash.haar.ang.1.v2[i,]=ashsmooth.gaus(sim.m.ang.1.v2[i,],sigma.ang.1.v2.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.ang.1.v2[i,]=ashsmooth.gaus(sim.m.ang.1.v2[i,],sigma.ang.1.v2.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.ang.1.v2[i,]=ashsmooth.gaus(sim.m.ang.1.v2[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.ang.1.v2[i,]=ashsmooth.gaus(sim.m.ang.1.v2[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.ang.1.v2[i,]=ti.thresh(sim.m.ang.1.v2[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.ang.1.v2[i,]=ti.thresh(sim.m.ang.1.v2[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.ang.1.v2[i,]=ashsmooth.gaus(sim.m.ang.1.v2[i,],sigma=sigma.ang.1.v2,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.ang.1.v2[i,]=ashsmooth.gaus(sim.m.ang.1.v2[i,],sigma=sigma.ang.1.v2,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.ang.1.v2[i,]=ti.thresh(sim.m.ang.1.v2[i,],sigma=sigma.ang.1.v2,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.ang.1.v2[i,]=ti.thresh(sim.m.ang.1.v2[i,],sigma=sigma.ang.1.v2,filter.number=8,family="DaubLeAsymm")
  sigma.ang.3.v2.est=sig.est.func(sim.m.ang.3.v2[i,],n)
  mu.est.ash.haar.ang.3.v2[i,]=ashsmooth.gaus(sim.m.ang.3.v2[i,],sigma.ang.3.v2.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.ang.3.v2[i,]=ashsmooth.gaus(sim.m.ang.3.v2[i,],sigma.ang.3.v2.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.ang.3.v2[i,]=ashsmooth.gaus(sim.m.ang.3.v2[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.ang.3.v2[i,]=ashsmooth.gaus(sim.m.ang.3.v2[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.ang.3.v2[i,]=ti.thresh(sim.m.ang.3.v2[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.ang.3.v2[i,]=ti.thresh(sim.m.ang.3.v2[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.ang.3.v2[i,]=ashsmooth.gaus(sim.m.ang.3.v2[i,],sigma=sigma.ang.3.v2,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.ang.3.v2[i,]=ashsmooth.gaus(sim.m.ang.3.v2[i,],sigma=sigma.ang.3.v2,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.ang.3.v2[i,]=ti.thresh(sim.m.ang.3.v2[i,],sigma=sigma.ang.3.v2,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.ang.3.v2[i,]=ti.thresh(sim.m.ang.3.v2[i,],sigma=sigma.ang.3.v2,filter.number=8,family="DaubLeAsymm")

  sigma.dop.1.v2.est=sig.est.func(sim.m.dop.1.v2[i,],n)
  mu.est.ash.haar.dop.1.v2[i,]=ashsmooth.gaus(sim.m.dop.1.v2[i,],sigma.dop.1.v2.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.dop.1.v2[i,]=ashsmooth.gaus(sim.m.dop.1.v2[i,],sigma.dop.1.v2.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.dop.1.v2[i,]=ashsmooth.gaus(sim.m.dop.1.v2[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.dop.1.v2[i,]=ashsmooth.gaus(sim.m.dop.1.v2[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.dop.1.v2[i,]=ti.thresh(sim.m.dop.1.v2[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.dop.1.v2[i,]=ti.thresh(sim.m.dop.1.v2[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.dop.1.v2[i,]=ashsmooth.gaus(sim.m.dop.1.v2[i,],sigma=sigma.dop.1.v2,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.dop.1.v2[i,]=ashsmooth.gaus(sim.m.dop.1.v2[i,],sigma=sigma.dop.1.v2,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.dop.1.v2[i,]=ti.thresh(sim.m.dop.1.v2[i,],sigma=sigma.dop.1.v2,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.dop.1.v2[i,]=ti.thresh(sim.m.dop.1.v2[i,],sigma=sigma.dop.1.v2,filter.number=8,family="DaubLeAsymm")
  sigma.dop.3.v2.est=sig.est.func(sim.m.dop.3.v2[i,],n)
  mu.est.ash.haar.dop.3.v2[i,]=ashsmooth.gaus(sim.m.dop.3.v2[i,],sigma.dop.3.v2.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.dop.3.v2[i,]=ashsmooth.gaus(sim.m.dop.3.v2[i,],sigma.dop.3.v2.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.dop.3.v2[i,]=ashsmooth.gaus(sim.m.dop.3.v2[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.dop.3.v2[i,]=ashsmooth.gaus(sim.m.dop.3.v2[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.dop.3.v2[i,]=ti.thresh(sim.m.dop.3.v2[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.dop.3.v2[i,]=ti.thresh(sim.m.dop.3.v2[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.dop.3.v2[i,]=ashsmooth.gaus(sim.m.dop.3.v2[i,],sigma=sigma.dop.3.v2,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.dop.3.v2[i,]=ashsmooth.gaus(sim.m.dop.3.v2[i,],sigma=sigma.dop.3.v2,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.dop.3.v2[i,]=ti.thresh(sim.m.dop.3.v2[i,],sigma=sigma.dop.3.v2,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.dop.3.v2[i,]=ti.thresh(sim.m.dop.3.v2[i,],sigma=sigma.dop.3.v2,filter.number=8,family="DaubLeAsymm")

  sigma.blip.1.v2.est=sig.est.func(sim.m.blip.1.v2[i,],n)
  mu.est.ash.haar.blip.1.v2[i,]=ashsmooth.gaus(sim.m.blip.1.v2[i,],sigma.blip.1.v2.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.blip.1.v2[i,]=ashsmooth.gaus(sim.m.blip.1.v2[i,],sigma.blip.1.v2.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.blip.1.v2[i,]=ashsmooth.gaus(sim.m.blip.1.v2[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.blip.1.v2[i,]=ashsmooth.gaus(sim.m.blip.1.v2[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.blip.1.v2[i,]=ti.thresh(sim.m.blip.1.v2[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.blip.1.v2[i,]=ti.thresh(sim.m.blip.1.v2[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.blip.1.v2[i,]=ashsmooth.gaus(sim.m.blip.1.v2[i,],sigma=sigma.blip.1.v2,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.blip.1.v2[i,]=ashsmooth.gaus(sim.m.blip.1.v2[i,],sigma=sigma.blip.1.v2,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.blip.1.v2[i,]=ti.thresh(sim.m.blip.1.v2[i,],sigma=sigma.blip.1.v2,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.blip.1.v2[i,]=ti.thresh(sim.m.blip.1.v2[i,],sigma=sigma.blip.1.v2,filter.number=8,family="DaubLeAsymm")
  sigma.blip.3.v2.est=sig.est.func(sim.m.blip.3.v2[i,],n)
  mu.est.ash.haar.blip.3.v2[i,]=ashsmooth.gaus(sim.m.blip.3.v2[i,],sigma.blip.3.v2.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.blip.3.v2[i,]=ashsmooth.gaus(sim.m.blip.3.v2[i,],sigma.blip.3.v2.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.blip.3.v2[i,]=ashsmooth.gaus(sim.m.blip.3.v2[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.blip.3.v2[i,]=ashsmooth.gaus(sim.m.blip.3.v2[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.blip.3.v2[i,]=ti.thresh(sim.m.blip.3.v2[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.blip.3.v2[i,]=ti.thresh(sim.m.blip.3.v2[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.blip.3.v2[i,]=ashsmooth.gaus(sim.m.blip.3.v2[i,],sigma=sigma.blip.3.v2,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.blip.3.v2[i,]=ashsmooth.gaus(sim.m.blip.3.v2[i,],sigma=sigma.blip.3.v2,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.blip.3.v2[i,]=ti.thresh(sim.m.blip.3.v2[i,],sigma=sigma.blip.3.v2,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.blip.3.v2[i,]=ti.thresh(sim.m.blip.3.v2[i,],sigma=sigma.blip.3.v2,filter.number=8,family="DaubLeAsymm")

  sigma.cor.1.v2.est=sig.est.func(sim.m.cor.1.v2[i,],n)
  mu.est.ash.haar.cor.1.v2[i,]=ashsmooth.gaus(sim.m.cor.1.v2[i,],sigma.cor.1.v2.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.cor.1.v2[i,]=ashsmooth.gaus(sim.m.cor.1.v2[i,],sigma.cor.1.v2.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.cor.1.v2[i,]=ashsmooth.gaus(sim.m.cor.1.v2[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.cor.1.v2[i,]=ashsmooth.gaus(sim.m.cor.1.v2[i,],filter.number=8,family="DaubLeAsymm",)
  mu.est.tie.haar.cor.1.v2[i,]=ti.thresh(sim.m.cor.1.v2[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.cor.1.v2[i,]=ti.thresh(sim.m.cor.1.v2[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.cor.1.v2[i,]=ashsmooth.gaus(sim.m.cor.1.v2[i,],sigma=sigma.cor.1.v2,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.cor.1.v2[i,]=ashsmooth.gaus(sim.m.cor.1.v2[i,],sigma=sigma.cor.1.v2,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.cor.1.v2[i,]=ti.thresh(sim.m.cor.1.v2[i,],sigma=sigma.cor.1.v2,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.cor.1.v2[i,]=ti.thresh(sim.m.cor.1.v2[i,],sigma=sigma.cor.1.v2,filter.number=8,family="DaubLeAsymm")
  sigma.cor.3.v2.est=sig.est.func(sim.m.cor.3.v2[i,],n)
  mu.est.ash.haar.cor.3.v2[i,]=ashsmooth.gaus(sim.m.cor.3.v2[i,],sigma.cor.3.v2.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.cor.3.v2[i,]=ashsmooth.gaus(sim.m.cor.3.v2[i,],sigma.cor.3.v2.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.cor.3.v2[i,]=ashsmooth.gaus(sim.m.cor.3.v2[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.cor.3.v2[i,]=ashsmooth.gaus(sim.m.cor.3.v2[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.cor.3.v2[i,]=ti.thresh(sim.m.cor.3.v2[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.cor.3.v2[i,]=ti.thresh(sim.m.cor.3.v2[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.cor.3.v2[i,]=ashsmooth.gaus(sim.m.cor.3.v2[i,],sigma=sigma.cor.3.v2,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.cor.3.v2[i,]=ashsmooth.gaus(sim.m.cor.3.v2[i,],sigma=sigma.cor.3.v2,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.cor.3.v2[i,]=ti.thresh(sim.m.cor.3.v2[i,],sigma=sigma.cor.3.v2,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.cor.3.v2[i,]=ti.thresh(sim.m.cor.3.v2[i,],sigma=sigma.cor.3.v2,filter.number=8,family="DaubLeAsymm")





  sigma.sp.1.v3.est=sig.est.func(sim.m.sp.1.v3[i,],n)
  mu.est.ash.haar.sp.1.v3[i,]=ashsmooth.gaus(sim.m.sp.1.v3[i,],sigma.sp.1.v3.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.sp.1.v3[i,]=ashsmooth.gaus(sim.m.sp.1.v3[i,],sigma.sp.1.v3.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.sp.1.v3[i,]=ashsmooth.gaus(sim.m.sp.1.v3[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.sp.1.v3[i,]=ashsmooth.gaus(sim.m.sp.1.v3[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.sp.1.v3[i,]=ti.thresh(sim.m.sp.1.v3[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.sp.1.v3[i,]=ti.thresh(sim.m.sp.1.v3[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.sp.1.v3[i,]=ashsmooth.gaus(sim.m.sp.1.v3[i,],sigma=sigma.sp.1.v3,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.sp.1.v3[i,]=ashsmooth.gaus(sim.m.sp.1.v3[i,],sigma=sigma.sp.1.v3,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.sp.1.v3[i,]=ti.thresh(sim.m.sp.1.v3[i,],sigma=sigma.sp.1.v3,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.sp.1.v3[i,]=ti.thresh(sim.m.sp.1.v3[i,],sigma=sigma.sp.1.v3,filter.number=8,family="DaubLeAsymm")
  sigma.sp.3.v3.est=sig.est.func(sim.m.sp.3.v3[i,],n)
  mu.est.ash.haar.sp.3.v3[i,]=ashsmooth.gaus(sim.m.sp.3.v3[i,],sigma.sp.3.v3.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.sp.3.v3[i,]=ashsmooth.gaus(sim.m.sp.3.v3[i,],sigma.sp.3.v3.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.sp.3.v3[i,]=ashsmooth.gaus(sim.m.sp.3.v3[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.sp.3.v3[i,]=ashsmooth.gaus(sim.m.sp.3.v3[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.sp.3.v3[i,]=ti.thresh(sim.m.sp.3.v3[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.sp.3.v3[i,]=ti.thresh(sim.m.sp.3.v3[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.sp.3.v3[i,]=ashsmooth.gaus(sim.m.sp.3.v3[i,],sigma=sigma.sp.3.v3,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.sp.3.v3[i,]=ashsmooth.gaus(sim.m.sp.3.v3[i,],sigma=sigma.sp.3.v3,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.sp.3.v3[i,]=ti.thresh(sim.m.sp.3.v3[i,],sigma=sigma.sp.3.v3,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.sp.3.v3[i,]=ti.thresh(sim.m.sp.3.v3[i,],sigma=sigma.sp.3.v3,filter.number=8,family="DaubLeAsymm")

  sigma.bump.1.v3.est=sig.est.func(sim.m.bump.1.v3[i,],n)
  mu.est.ash.haar.bump.1.v3[i,]=ashsmooth.gaus(sim.m.bump.1.v3[i,],sigma.bump.1.v3.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.bump.1.v3[i,]=ashsmooth.gaus(sim.m.bump.1.v3[i,],sigma.bump.1.v3.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.bump.1.v3[i,]=ashsmooth.gaus(sim.m.bump.1.v3[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.bump.1.v3[i,]=ashsmooth.gaus(sim.m.bump.1.v3[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.bump.1.v3[i,]=ti.thresh(sim.m.bump.1.v3[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.bump.1.v3[i,]=ti.thresh(sim.m.bump.1.v3[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.bump.1.v3[i,]=ashsmooth.gaus(sim.m.bump.1.v3[i,],sigma=sigma.bump.1.v3,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.bump.1.v3[i,]=ashsmooth.gaus(sim.m.bump.1.v3[i,],sigma=sigma.bump.1.v3,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.bump.1.v3[i,]=ti.thresh(sim.m.bump.1.v3[i,],sigma=sigma.bump.1.v3,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.bump.1.v3[i,]=ti.thresh(sim.m.bump.1.v3[i,],sigma=sigma.bump.1.v3,filter.number=8,family="DaubLeAsymm")
  sigma.bump.3.v3.est=sig.est.func(sim.m.bump.3.v3[i,],n)
  mu.est.ash.haar.bump.3.v3[i,]=ashsmooth.gaus(sim.m.bump.3.v3[i,],sigma.bump.3.v3.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.bump.3.v3[i,]=ashsmooth.gaus(sim.m.bump.3.v3[i,],sigma.bump.3.v3.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.bump.3.v3[i,]=ashsmooth.gaus(sim.m.bump.3.v3[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.bump.3.v3[i,]=ashsmooth.gaus(sim.m.bump.3.v3[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.bump.3.v3[i,]=ti.thresh(sim.m.bump.3.v3[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.bump.3.v3[i,]=ti.thresh(sim.m.bump.3.v3[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.bump.3.v3[i,]=ashsmooth.gaus(sim.m.bump.3.v3[i,],sigma=sigma.bump.3.v3,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.bump.3.v3[i,]=ashsmooth.gaus(sim.m.bump.3.v3[i,],sigma=sigma.bump.3.v3,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.bump.3.v3[i,]=ti.thresh(sim.m.bump.3.v3[i,],sigma=sigma.bump.3.v3,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.bump.3.v3[i,]=ti.thresh(sim.m.bump.3.v3[i,],sigma=sigma.bump.3.v3,filter.number=8,family="DaubLeAsymm")
  
  sigma.blk.1.v3.est=sig.est.func(sim.m.blk.1.v3[i,],n)
  mu.est.ash.haar.blk.1.v3[i,]=ashsmooth.gaus(sim.m.blk.1.v3[i,],sigma.blk.1.v3.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.blk.1.v3[i,]=ashsmooth.gaus(sim.m.blk.1.v3[i,],sigma.blk.1.v3.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.blk.1.v3[i,]=ashsmooth.gaus(sim.m.blk.1.v3[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.blk.1.v3[i,]=ashsmooth.gaus(sim.m.blk.1.v3[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.blk.1.v3[i,]=ti.thresh(sim.m.blk.1.v3[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.blk.1.v3[i,]=ti.thresh(sim.m.blk.1.v3[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.blk.1.v3[i,]=ashsmooth.gaus(sim.m.blk.1.v3[i,],sigma=sigma.blk.1.v3,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.blk.1.v3[i,]=ashsmooth.gaus(sim.m.blk.1.v3[i,],sigma=sigma.blk.1.v3,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.blk.1.v3[i,]=ti.thresh(sim.m.blk.1.v3[i,],sigma=sigma.blk.1.v3,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.blk.1.v3[i,]=ti.thresh(sim.m.blk.1.v3[i,],sigma=sigma.blk.1.v3,filter.number=8,family="DaubLeAsymm")
  sigma.blk.3.v3.est=sig.est.func(sim.m.blk.3.v3[i,],n)
  mu.est.ash.haar.blk.3.v3[i,]=ashsmooth.gaus(sim.m.blk.3.v3[i,],sigma.blk.3.v3.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.blk.3.v3[i,]=ashsmooth.gaus(sim.m.blk.3.v3[i,],sigma.blk.3.v3.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.blk.3.v3[i,]=ashsmooth.gaus(sim.m.blk.3.v3[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.blk.3.v3[i,]=ashsmooth.gaus(sim.m.blk.3.v3[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.blk.3.v3[i,]=ti.thresh(sim.m.blk.3.v3[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.blk.3.v3[i,]=ti.thresh(sim.m.blk.3.v3[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.blk.3.v3[i,]=ashsmooth.gaus(sim.m.blk.3.v3[i,],sigma=sigma.blk.3.v3,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.blk.3.v3[i,]=ashsmooth.gaus(sim.m.blk.3.v3[i,],sigma=sigma.blk.3.v3,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.blk.3.v3[i,]=ti.thresh(sim.m.blk.3.v3[i,],sigma=sigma.blk.3.v3,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.blk.3.v3[i,]=ti.thresh(sim.m.blk.3.v3[i,],sigma=sigma.blk.3.v3,filter.number=8,family="DaubLeAsymm")

  sigma.ang.1.v3.est=sig.est.func(sim.m.ang.1.v3[i,],n)
  mu.est.ash.haar.ang.1.v3[i,]=ashsmooth.gaus(sim.m.ang.1.v3[i,],sigma.ang.1.v3.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.ang.1.v3[i,]=ashsmooth.gaus(sim.m.ang.1.v3[i,],sigma.ang.1.v3.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.ang.1.v3[i,]=ashsmooth.gaus(sim.m.ang.1.v3[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.ang.1.v3[i,]=ashsmooth.gaus(sim.m.ang.1.v3[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.ang.1.v3[i,]=ti.thresh(sim.m.ang.1.v3[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.ang.1.v3[i,]=ti.thresh(sim.m.ang.1.v3[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.ang.1.v3[i,]=ashsmooth.gaus(sim.m.ang.1.v3[i,],sigma=sigma.ang.1.v3,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.ang.1.v3[i,]=ashsmooth.gaus(sim.m.ang.1.v3[i,],sigma=sigma.ang.1.v3,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.ang.1.v3[i,]=ti.thresh(sim.m.ang.1.v3[i,],sigma=sigma.ang.1.v3,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.ang.1.v3[i,]=ti.thresh(sim.m.ang.1.v3[i,],sigma=sigma.ang.1.v3,filter.number=8,family="DaubLeAsymm")
  sigma.ang.3.v3.est=sig.est.func(sim.m.ang.3.v3[i,],n)
  mu.est.ash.haar.ang.3.v3[i,]=ashsmooth.gaus(sim.m.ang.3.v3[i,],sigma.ang.3.v3.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.ang.3.v3[i,]=ashsmooth.gaus(sim.m.ang.3.v3[i,],sigma.ang.3.v3.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.ang.3.v3[i,]=ashsmooth.gaus(sim.m.ang.3.v3[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.ang.3.v3[i,]=ashsmooth.gaus(sim.m.ang.3.v3[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.ang.3.v3[i,]=ti.thresh(sim.m.ang.3.v3[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.ang.3.v3[i,]=ti.thresh(sim.m.ang.3.v3[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.ang.3.v3[i,]=ashsmooth.gaus(sim.m.ang.3.v3[i,],sigma=sigma.ang.3.v3,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.ang.3.v3[i,]=ashsmooth.gaus(sim.m.ang.3.v3[i,],sigma=sigma.ang.3.v3,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.ang.3.v3[i,]=ti.thresh(sim.m.ang.3.v3[i,],sigma=sigma.ang.3.v3,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.ang.3.v3[i,]=ti.thresh(sim.m.ang.3.v3[i,],sigma=sigma.ang.3.v3,filter.number=8,family="DaubLeAsymm")

  sigma.dop.1.v3.est=sig.est.func(sim.m.dop.1.v3[i,],n)
  mu.est.ash.haar.dop.1.v3[i,]=ashsmooth.gaus(sim.m.dop.1.v3[i,],sigma.dop.1.v3.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.dop.1.v3[i,]=ashsmooth.gaus(sim.m.dop.1.v3[i,],sigma.dop.1.v3.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.dop.1.v3[i,]=ashsmooth.gaus(sim.m.dop.1.v3[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.dop.1.v3[i,]=ashsmooth.gaus(sim.m.dop.1.v3[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.dop.1.v3[i,]=ti.thresh(sim.m.dop.1.v3[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.dop.1.v3[i,]=ti.thresh(sim.m.dop.1.v3[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.dop.1.v3[i,]=ashsmooth.gaus(sim.m.dop.1.v3[i,],sigma=sigma.dop.1.v3,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.dop.1.v3[i,]=ashsmooth.gaus(sim.m.dop.1.v3[i,],sigma=sigma.dop.1.v3,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.dop.1.v3[i,]=ti.thresh(sim.m.dop.1.v3[i,],sigma=sigma.dop.1.v3,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.dop.1.v3[i,]=ti.thresh(sim.m.dop.1.v3[i,],sigma=sigma.dop.1.v3,filter.number=8,family="DaubLeAsymm")
  sigma.dop.3.v3.est=sig.est.func(sim.m.dop.3.v3[i,],n)
  mu.est.ash.haar.dop.3.v3[i,]=ashsmooth.gaus(sim.m.dop.3.v3[i,],sigma.dop.3.v3.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.dop.3.v3[i,]=ashsmooth.gaus(sim.m.dop.3.v3[i,],sigma.dop.3.v3.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.dop.3.v3[i,]=ashsmooth.gaus(sim.m.dop.3.v3[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.dop.3.v3[i,]=ashsmooth.gaus(sim.m.dop.3.v3[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.dop.3.v3[i,]=ti.thresh(sim.m.dop.3.v3[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.dop.3.v3[i,]=ti.thresh(sim.m.dop.3.v3[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.dop.3.v3[i,]=ashsmooth.gaus(sim.m.dop.3.v3[i,],sigma=sigma.dop.3.v3,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.dop.3.v3[i,]=ashsmooth.gaus(sim.m.dop.3.v3[i,],sigma=sigma.dop.3.v3,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.dop.3.v3[i,]=ti.thresh(sim.m.dop.3.v3[i,],sigma=sigma.dop.3.v3,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.dop.3.v3[i,]=ti.thresh(sim.m.dop.3.v3[i,],sigma=sigma.dop.3.v3,filter.number=8,family="DaubLeAsymm")

  sigma.blip.1.v3.est=sig.est.func(sim.m.blip.1.v3[i,],n)
  mu.est.ash.haar.blip.1.v3[i,]=ashsmooth.gaus(sim.m.blip.1.v3[i,],sigma.blip.1.v3.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.blip.1.v3[i,]=ashsmooth.gaus(sim.m.blip.1.v3[i,],sigma.blip.1.v3.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.blip.1.v3[i,]=ashsmooth.gaus(sim.m.blip.1.v3[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.blip.1.v3[i,]=ashsmooth.gaus(sim.m.blip.1.v3[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.blip.1.v3[i,]=ti.thresh(sim.m.blip.1.v3[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.blip.1.v3[i,]=ti.thresh(sim.m.blip.1.v3[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.blip.1.v3[i,]=ashsmooth.gaus(sim.m.blip.1.v3[i,],sigma=sigma.blip.1.v3,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.blip.1.v3[i,]=ashsmooth.gaus(sim.m.blip.1.v3[i,],sigma=sigma.blip.1.v3,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.blip.1.v3[i,]=ti.thresh(sim.m.blip.1.v3[i,],sigma=sigma.blip.1.v3,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.blip.1.v3[i,]=ti.thresh(sim.m.blip.1.v3[i,],sigma=sigma.blip.1.v3,filter.number=8,family="DaubLeAsymm")
  sigma.blip.3.v3.est=sig.est.func(sim.m.blip.3.v3[i,],n)
  mu.est.ash.haar.blip.3.v3[i,]=ashsmooth.gaus(sim.m.blip.3.v3[i,],sigma.blip.3.v3.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.blip.3.v3[i,]=ashsmooth.gaus(sim.m.blip.3.v3[i,],sigma.blip.3.v3.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.blip.3.v3[i,]=ashsmooth.gaus(sim.m.blip.3.v3[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.blip.3.v3[i,]=ashsmooth.gaus(sim.m.blip.3.v3[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.blip.3.v3[i,]=ti.thresh(sim.m.blip.3.v3[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.blip.3.v3[i,]=ti.thresh(sim.m.blip.3.v3[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.blip.3.v3[i,]=ashsmooth.gaus(sim.m.blip.3.v3[i,],sigma=sigma.blip.3.v3,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.blip.3.v3[i,]=ashsmooth.gaus(sim.m.blip.3.v3[i,],sigma=sigma.blip.3.v3,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.blip.3.v3[i,]=ti.thresh(sim.m.blip.3.v3[i,],sigma=sigma.blip.3.v3,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.blip.3.v3[i,]=ti.thresh(sim.m.blip.3.v3[i,],sigma=sigma.blip.3.v3,filter.number=8,family="DaubLeAsymm")

  sigma.cor.1.v3.est=sig.est.func(sim.m.cor.1.v3[i,],n)
  mu.est.ash.haar.cor.1.v3[i,]=ashsmooth.gaus(sim.m.cor.1.v3[i,],sigma.cor.1.v3.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.cor.1.v3[i,]=ashsmooth.gaus(sim.m.cor.1.v3[i,],sigma.cor.1.v3.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.cor.1.v3[i,]=ashsmooth.gaus(sim.m.cor.1.v3[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.cor.1.v3[i,]=ashsmooth.gaus(sim.m.cor.1.v3[i,],filter.number=8,family="DaubLeAsymm",)
  mu.est.tie.haar.cor.1.v3[i,]=ti.thresh(sim.m.cor.1.v3[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.cor.1.v3[i,]=ti.thresh(sim.m.cor.1.v3[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.cor.1.v3[i,]=ashsmooth.gaus(sim.m.cor.1.v3[i,],sigma=sigma.cor.1.v3,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.cor.1.v3[i,]=ashsmooth.gaus(sim.m.cor.1.v3[i,],sigma=sigma.cor.1.v3,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.cor.1.v3[i,]=ti.thresh(sim.m.cor.1.v3[i,],sigma=sigma.cor.1.v3,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.cor.1.v3[i,]=ti.thresh(sim.m.cor.1.v3[i,],sigma=sigma.cor.1.v3,filter.number=8,family="DaubLeAsymm")
  sigma.cor.3.v3.est=sig.est.func(sim.m.cor.3.v3[i,],n)
  mu.est.ash.haar.cor.3.v3[i,]=ashsmooth.gaus(sim.m.cor.3.v3[i,],sigma.cor.3.v3.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.cor.3.v3[i,]=ashsmooth.gaus(sim.m.cor.3.v3[i,],sigma.cor.3.v3.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.cor.3.v3[i,]=ashsmooth.gaus(sim.m.cor.3.v3[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.cor.3.v3[i,]=ashsmooth.gaus(sim.m.cor.3.v3[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.cor.3.v3[i,]=ti.thresh(sim.m.cor.3.v3[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.cor.3.v3[i,]=ti.thresh(sim.m.cor.3.v3[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.cor.3.v3[i,]=ashsmooth.gaus(sim.m.cor.3.v3[i,],sigma=sigma.cor.3.v3,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.cor.3.v3[i,]=ashsmooth.gaus(sim.m.cor.3.v3[i,],sigma=sigma.cor.3.v3,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.cor.3.v3[i,]=ti.thresh(sim.m.cor.3.v3[i,],sigma=sigma.cor.3.v3,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.cor.3.v3[i,]=ti.thresh(sim.m.cor.3.v3[i,],sigma=sigma.cor.3.v3,filter.number=8,family="DaubLeAsymm")



  sigma.sp.1.v4.est=sig.est.func(sim.m.sp.1.v4[i,],n)
  mu.est.ash.haar.sp.1.v4[i,]=ashsmooth.gaus(sim.m.sp.1.v4[i,],sigma.sp.1.v4.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.sp.1.v4[i,]=ashsmooth.gaus(sim.m.sp.1.v4[i,],sigma.sp.1.v4.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.sp.1.v4[i,]=ashsmooth.gaus(sim.m.sp.1.v4[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.sp.1.v4[i,]=ashsmooth.gaus(sim.m.sp.1.v4[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.sp.1.v4[i,]=ti.thresh(sim.m.sp.1.v4[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.sp.1.v4[i,]=ti.thresh(sim.m.sp.1.v4[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.sp.1.v4[i,]=ashsmooth.gaus(sim.m.sp.1.v4[i,],sigma=sigma.sp.1.v4,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.sp.1.v4[i,]=ashsmooth.gaus(sim.m.sp.1.v4[i,],sigma=sigma.sp.1.v4,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.sp.1.v4[i,]=ti.thresh(sim.m.sp.1.v4[i,],sigma=sigma.sp.1.v4,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.sp.1.v4[i,]=ti.thresh(sim.m.sp.1.v4[i,],sigma=sigma.sp.1.v4,filter.number=8,family="DaubLeAsymm")
  sigma.sp.3.v4.est=sig.est.func(sim.m.sp.3.v4[i,],n)
  mu.est.ash.haar.sp.3.v4[i,]=ashsmooth.gaus(sim.m.sp.3.v4[i,],sigma.sp.3.v4.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.sp.3.v4[i,]=ashsmooth.gaus(sim.m.sp.3.v4[i,],sigma.sp.3.v4.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.sp.3.v4[i,]=ashsmooth.gaus(sim.m.sp.3.v4[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.sp.3.v4[i,]=ashsmooth.gaus(sim.m.sp.3.v4[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.sp.3.v4[i,]=ti.thresh(sim.m.sp.3.v4[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.sp.3.v4[i,]=ti.thresh(sim.m.sp.3.v4[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.sp.3.v4[i,]=ashsmooth.gaus(sim.m.sp.3.v4[i,],sigma=sigma.sp.3.v4,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.sp.3.v4[i,]=ashsmooth.gaus(sim.m.sp.3.v4[i,],sigma=sigma.sp.3.v4,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.sp.3.v4[i,]=ti.thresh(sim.m.sp.3.v4[i,],sigma=sigma.sp.3.v4,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.sp.3.v4[i,]=ti.thresh(sim.m.sp.3.v4[i,],sigma=sigma.sp.3.v4,filter.number=8,family="DaubLeAsymm")

  sigma.bump.1.v4.est=sig.est.func(sim.m.bump.1.v4[i,],n)
  mu.est.ash.haar.bump.1.v4[i,]=ashsmooth.gaus(sim.m.bump.1.v4[i,],sigma.bump.1.v4.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.bump.1.v4[i,]=ashsmooth.gaus(sim.m.bump.1.v4[i,],sigma.bump.1.v4.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.bump.1.v4[i,]=ashsmooth.gaus(sim.m.bump.1.v4[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.bump.1.v4[i,]=ashsmooth.gaus(sim.m.bump.1.v4[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.bump.1.v4[i,]=ti.thresh(sim.m.bump.1.v4[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.bump.1.v4[i,]=ti.thresh(sim.m.bump.1.v4[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.bump.1.v4[i,]=ashsmooth.gaus(sim.m.bump.1.v4[i,],sigma=sigma.bump.1.v4,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.bump.1.v4[i,]=ashsmooth.gaus(sim.m.bump.1.v4[i,],sigma=sigma.bump.1.v4,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.bump.1.v4[i,]=ti.thresh(sim.m.bump.1.v4[i,],sigma=sigma.bump.1.v4,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.bump.1.v4[i,]=ti.thresh(sim.m.bump.1.v4[i,],sigma=sigma.bump.1.v4,filter.number=8,family="DaubLeAsymm")
  sigma.bump.3.v4.est=sig.est.func(sim.m.bump.3.v4[i,],n)
  mu.est.ash.haar.bump.3.v4[i,]=ashsmooth.gaus(sim.m.bump.3.v4[i,],sigma.bump.3.v4.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.bump.3.v4[i,]=ashsmooth.gaus(sim.m.bump.3.v4[i,],sigma.bump.3.v4.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.bump.3.v4[i,]=ashsmooth.gaus(sim.m.bump.3.v4[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.bump.3.v4[i,]=ashsmooth.gaus(sim.m.bump.3.v4[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.bump.3.v4[i,]=ti.thresh(sim.m.bump.3.v4[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.bump.3.v4[i,]=ti.thresh(sim.m.bump.3.v4[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.bump.3.v4[i,]=ashsmooth.gaus(sim.m.bump.3.v4[i,],sigma=sigma.bump.3.v4,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.bump.3.v4[i,]=ashsmooth.gaus(sim.m.bump.3.v4[i,],sigma=sigma.bump.3.v4,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.bump.3.v4[i,]=ti.thresh(sim.m.bump.3.v4[i,],sigma=sigma.bump.3.v4,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.bump.3.v4[i,]=ti.thresh(sim.m.bump.3.v4[i,],sigma=sigma.bump.3.v4,filter.number=8,family="DaubLeAsymm")
  
  sigma.blk.1.v4.est=sig.est.func(sim.m.blk.1.v4[i,],n)
  mu.est.ash.haar.blk.1.v4[i,]=ashsmooth.gaus(sim.m.blk.1.v4[i,],sigma.blk.1.v4.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.blk.1.v4[i,]=ashsmooth.gaus(sim.m.blk.1.v4[i,],sigma.blk.1.v4.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.blk.1.v4[i,]=ashsmooth.gaus(sim.m.blk.1.v4[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.blk.1.v4[i,]=ashsmooth.gaus(sim.m.blk.1.v4[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.blk.1.v4[i,]=ti.thresh(sim.m.blk.1.v4[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.blk.1.v4[i,]=ti.thresh(sim.m.blk.1.v4[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.blk.1.v4[i,]=ashsmooth.gaus(sim.m.blk.1.v4[i,],sigma=sigma.blk.1.v4,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.blk.1.v4[i,]=ashsmooth.gaus(sim.m.blk.1.v4[i,],sigma=sigma.blk.1.v4,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.blk.1.v4[i,]=ti.thresh(sim.m.blk.1.v4[i,],sigma=sigma.blk.1.v4,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.blk.1.v4[i,]=ti.thresh(sim.m.blk.1.v4[i,],sigma=sigma.blk.1.v4,filter.number=8,family="DaubLeAsymm")
  sigma.blk.3.v4.est=sig.est.func(sim.m.blk.3.v4[i,],n)
  mu.est.ash.haar.blk.3.v4[i,]=ashsmooth.gaus(sim.m.blk.3.v4[i,],sigma.blk.3.v4.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.blk.3.v4[i,]=ashsmooth.gaus(sim.m.blk.3.v4[i,],sigma.blk.3.v4.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.blk.3.v4[i,]=ashsmooth.gaus(sim.m.blk.3.v4[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.blk.3.v4[i,]=ashsmooth.gaus(sim.m.blk.3.v4[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.blk.3.v4[i,]=ti.thresh(sim.m.blk.3.v4[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.blk.3.v4[i,]=ti.thresh(sim.m.blk.3.v4[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.blk.3.v4[i,]=ashsmooth.gaus(sim.m.blk.3.v4[i,],sigma=sigma.blk.3.v4,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.blk.3.v4[i,]=ashsmooth.gaus(sim.m.blk.3.v4[i,],sigma=sigma.blk.3.v4,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.blk.3.v4[i,]=ti.thresh(sim.m.blk.3.v4[i,],sigma=sigma.blk.3.v4,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.blk.3.v4[i,]=ti.thresh(sim.m.blk.3.v4[i,],sigma=sigma.blk.3.v4,filter.number=8,family="DaubLeAsymm")

  sigma.ang.1.v4.est=sig.est.func(sim.m.ang.1.v4[i,],n)
  mu.est.ash.haar.ang.1.v4[i,]=ashsmooth.gaus(sim.m.ang.1.v4[i,],sigma.ang.1.v4.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.ang.1.v4[i,]=ashsmooth.gaus(sim.m.ang.1.v4[i,],sigma.ang.1.v4.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.ang.1.v4[i,]=ashsmooth.gaus(sim.m.ang.1.v4[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.ang.1.v4[i,]=ashsmooth.gaus(sim.m.ang.1.v4[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.ang.1.v4[i,]=ti.thresh(sim.m.ang.1.v4[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.ang.1.v4[i,]=ti.thresh(sim.m.ang.1.v4[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.ang.1.v4[i,]=ashsmooth.gaus(sim.m.ang.1.v4[i,],sigma=sigma.ang.1.v4,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.ang.1.v4[i,]=ashsmooth.gaus(sim.m.ang.1.v4[i,],sigma=sigma.ang.1.v4,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.ang.1.v4[i,]=ti.thresh(sim.m.ang.1.v4[i,],sigma=sigma.ang.1.v4,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.ang.1.v4[i,]=ti.thresh(sim.m.ang.1.v4[i,],sigma=sigma.ang.1.v4,filter.number=8,family="DaubLeAsymm")
  sigma.ang.3.v4.est=sig.est.func(sim.m.ang.3.v4[i,],n)
  mu.est.ash.haar.ang.3.v4[i,]=ashsmooth.gaus(sim.m.ang.3.v4[i,],sigma.ang.3.v4.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.ang.3.v4[i,]=ashsmooth.gaus(sim.m.ang.3.v4[i,],sigma.ang.3.v4.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.ang.3.v4[i,]=ashsmooth.gaus(sim.m.ang.3.v4[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.ang.3.v4[i,]=ashsmooth.gaus(sim.m.ang.3.v4[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.ang.3.v4[i,]=ti.thresh(sim.m.ang.3.v4[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.ang.3.v4[i,]=ti.thresh(sim.m.ang.3.v4[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.ang.3.v4[i,]=ashsmooth.gaus(sim.m.ang.3.v4[i,],sigma=sigma.ang.3.v4,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.ang.3.v4[i,]=ashsmooth.gaus(sim.m.ang.3.v4[i,],sigma=sigma.ang.3.v4,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.ang.3.v4[i,]=ti.thresh(sim.m.ang.3.v4[i,],sigma=sigma.ang.3.v4,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.ang.3.v4[i,]=ti.thresh(sim.m.ang.3.v4[i,],sigma=sigma.ang.3.v4,filter.number=8,family="DaubLeAsymm")

  sigma.dop.1.v4.est=sig.est.func(sim.m.dop.1.v4[i,],n)
  mu.est.ash.haar.dop.1.v4[i,]=ashsmooth.gaus(sim.m.dop.1.v4[i,],sigma.dop.1.v4.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.dop.1.v4[i,]=ashsmooth.gaus(sim.m.dop.1.v4[i,],sigma.dop.1.v4.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.dop.1.v4[i,]=ashsmooth.gaus(sim.m.dop.1.v4[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.dop.1.v4[i,]=ashsmooth.gaus(sim.m.dop.1.v4[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.dop.1.v4[i,]=ti.thresh(sim.m.dop.1.v4[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.dop.1.v4[i,]=ti.thresh(sim.m.dop.1.v4[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.dop.1.v4[i,]=ashsmooth.gaus(sim.m.dop.1.v4[i,],sigma=sigma.dop.1.v4,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.dop.1.v4[i,]=ashsmooth.gaus(sim.m.dop.1.v4[i,],sigma=sigma.dop.1.v4,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.dop.1.v4[i,]=ti.thresh(sim.m.dop.1.v4[i,],sigma=sigma.dop.1.v4,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.dop.1.v4[i,]=ti.thresh(sim.m.dop.1.v4[i,],sigma=sigma.dop.1.v4,filter.number=8,family="DaubLeAsymm")
  sigma.dop.3.v4.est=sig.est.func(sim.m.dop.3.v4[i,],n)
  mu.est.ash.haar.dop.3.v4[i,]=ashsmooth.gaus(sim.m.dop.3.v4[i,],sigma.dop.3.v4.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.dop.3.v4[i,]=ashsmooth.gaus(sim.m.dop.3.v4[i,],sigma.dop.3.v4.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.dop.3.v4[i,]=ashsmooth.gaus(sim.m.dop.3.v4[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.dop.3.v4[i,]=ashsmooth.gaus(sim.m.dop.3.v4[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.dop.3.v4[i,]=ti.thresh(sim.m.dop.3.v4[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.dop.3.v4[i,]=ti.thresh(sim.m.dop.3.v4[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.dop.3.v4[i,]=ashsmooth.gaus(sim.m.dop.3.v4[i,],sigma=sigma.dop.3.v4,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.dop.3.v4[i,]=ashsmooth.gaus(sim.m.dop.3.v4[i,],sigma=sigma.dop.3.v4,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.dop.3.v4[i,]=ti.thresh(sim.m.dop.3.v4[i,],sigma=sigma.dop.3.v4,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.dop.3.v4[i,]=ti.thresh(sim.m.dop.3.v4[i,],sigma=sigma.dop.3.v4,filter.number=8,family="DaubLeAsymm")

  sigma.blip.1.v4.est=sig.est.func(sim.m.blip.1.v4[i,],n)
  mu.est.ash.haar.blip.1.v4[i,]=ashsmooth.gaus(sim.m.blip.1.v4[i,],sigma.blip.1.v4.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.blip.1.v4[i,]=ashsmooth.gaus(sim.m.blip.1.v4[i,],sigma.blip.1.v4.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.blip.1.v4[i,]=ashsmooth.gaus(sim.m.blip.1.v4[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.blip.1.v4[i,]=ashsmooth.gaus(sim.m.blip.1.v4[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.blip.1.v4[i,]=ti.thresh(sim.m.blip.1.v4[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.blip.1.v4[i,]=ti.thresh(sim.m.blip.1.v4[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.blip.1.v4[i,]=ashsmooth.gaus(sim.m.blip.1.v4[i,],sigma=sigma.blip.1.v4,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.blip.1.v4[i,]=ashsmooth.gaus(sim.m.blip.1.v4[i,],sigma=sigma.blip.1.v4,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.blip.1.v4[i,]=ti.thresh(sim.m.blip.1.v4[i,],sigma=sigma.blip.1.v4,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.blip.1.v4[i,]=ti.thresh(sim.m.blip.1.v4[i,],sigma=sigma.blip.1.v4,filter.number=8,family="DaubLeAsymm")
  sigma.blip.3.v4.est=sig.est.func(sim.m.blip.3.v4[i,],n)
  mu.est.ash.haar.blip.3.v4[i,]=ashsmooth.gaus(sim.m.blip.3.v4[i,],sigma.blip.3.v4.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.blip.3.v4[i,]=ashsmooth.gaus(sim.m.blip.3.v4[i,],sigma.blip.3.v4.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.blip.3.v4[i,]=ashsmooth.gaus(sim.m.blip.3.v4[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.blip.3.v4[i,]=ashsmooth.gaus(sim.m.blip.3.v4[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.blip.3.v4[i,]=ti.thresh(sim.m.blip.3.v4[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.blip.3.v4[i,]=ti.thresh(sim.m.blip.3.v4[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.blip.3.v4[i,]=ashsmooth.gaus(sim.m.blip.3.v4[i,],sigma=sigma.blip.3.v4,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.blip.3.v4[i,]=ashsmooth.gaus(sim.m.blip.3.v4[i,],sigma=sigma.blip.3.v4,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.blip.3.v4[i,]=ti.thresh(sim.m.blip.3.v4[i,],sigma=sigma.blip.3.v4,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.blip.3.v4[i,]=ti.thresh(sim.m.blip.3.v4[i,],sigma=sigma.blip.3.v4,filter.number=8,family="DaubLeAsymm")

  sigma.cor.1.v4.est=sig.est.func(sim.m.cor.1.v4[i,],n)
  mu.est.ash.haar.cor.1.v4[i,]=ashsmooth.gaus(sim.m.cor.1.v4[i,],sigma.cor.1.v4.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.cor.1.v4[i,]=ashsmooth.gaus(sim.m.cor.1.v4[i,],sigma.cor.1.v4.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.cor.1.v4[i,]=ashsmooth.gaus(sim.m.cor.1.v4[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.cor.1.v4[i,]=ashsmooth.gaus(sim.m.cor.1.v4[i,],filter.number=8,family="DaubLeAsymm",)
  mu.est.tie.haar.cor.1.v4[i,]=ti.thresh(sim.m.cor.1.v4[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.cor.1.v4[i,]=ti.thresh(sim.m.cor.1.v4[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.cor.1.v4[i,]=ashsmooth.gaus(sim.m.cor.1.v4[i,],sigma=sigma.cor.1.v4,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.cor.1.v4[i,]=ashsmooth.gaus(sim.m.cor.1.v4[i,],sigma=sigma.cor.1.v4,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.cor.1.v4[i,]=ti.thresh(sim.m.cor.1.v4[i,],sigma=sigma.cor.1.v4,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.cor.1.v4[i,]=ti.thresh(sim.m.cor.1.v4[i,],sigma=sigma.cor.1.v4,filter.number=8,family="DaubLeAsymm")
  sigma.cor.3.v4.est=sig.est.func(sim.m.cor.3.v4[i,],n)
  mu.est.ash.haar.cor.3.v4[i,]=ashsmooth.gaus(sim.m.cor.3.v4[i,],sigma.cor.3.v4.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.cor.3.v4[i,]=ashsmooth.gaus(sim.m.cor.3.v4[i,],sigma.cor.3.v4.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.cor.3.v4[i,]=ashsmooth.gaus(sim.m.cor.3.v4[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.cor.3.v4[i,]=ashsmooth.gaus(sim.m.cor.3.v4[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.cor.3.v4[i,]=ti.thresh(sim.m.cor.3.v4[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.cor.3.v4[i,]=ti.thresh(sim.m.cor.3.v4[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.cor.3.v4[i,]=ashsmooth.gaus(sim.m.cor.3.v4[i,],sigma=sigma.cor.3.v4,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.cor.3.v4[i,]=ashsmooth.gaus(sim.m.cor.3.v4[i,],sigma=sigma.cor.3.v4,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.cor.3.v4[i,]=ti.thresh(sim.m.cor.3.v4[i,],sigma=sigma.cor.3.v4,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.cor.3.v4[i,]=ti.thresh(sim.m.cor.3.v4[i,],sigma=sigma.cor.3.v4,filter.number=8,family="DaubLeAsymm")



  sigma.sp.1.v5.est=sig.est.func(sim.m.sp.1.v5[i,],n)
  mu.est.ash.haar.sp.1.v5[i,]=ashsmooth.gaus(sim.m.sp.1.v5[i,],sigma.sp.1.v5.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.sp.1.v5[i,]=ashsmooth.gaus(sim.m.sp.1.v5[i,],sigma.sp.1.v5.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.sp.1.v5[i,]=ashsmooth.gaus(sim.m.sp.1.v5[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.sp.1.v5[i,]=ashsmooth.gaus(sim.m.sp.1.v5[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.sp.1.v5[i,]=ti.thresh(sim.m.sp.1.v5[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.sp.1.v5[i,]=ti.thresh(sim.m.sp.1.v5[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.sp.1.v5[i,]=ashsmooth.gaus(sim.m.sp.1.v5[i,],sigma=sigma.sp.1.v5,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.sp.1.v5[i,]=ashsmooth.gaus(sim.m.sp.1.v5[i,],sigma=sigma.sp.1.v5,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.sp.1.v5[i,]=ti.thresh(sim.m.sp.1.v5[i,],sigma=sigma.sp.1.v5,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.sp.1.v5[i,]=ti.thresh(sim.m.sp.1.v5[i,],sigma=sigma.sp.1.v5,filter.number=8,family="DaubLeAsymm")
  sigma.sp.3.v5.est=sig.est.func(sim.m.sp.3.v5[i,],n)
  mu.est.ash.haar.sp.3.v5[i,]=ashsmooth.gaus(sim.m.sp.3.v5[i,],sigma.sp.3.v5.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.sp.3.v5[i,]=ashsmooth.gaus(sim.m.sp.3.v5[i,],sigma.sp.3.v5.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.sp.3.v5[i,]=ashsmooth.gaus(sim.m.sp.3.v5[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.sp.3.v5[i,]=ashsmooth.gaus(sim.m.sp.3.v5[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.sp.3.v5[i,]=ti.thresh(sim.m.sp.3.v5[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.sp.3.v5[i,]=ti.thresh(sim.m.sp.3.v5[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.sp.3.v5[i,]=ashsmooth.gaus(sim.m.sp.3.v5[i,],sigma=sigma.sp.3.v5,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.sp.3.v5[i,]=ashsmooth.gaus(sim.m.sp.3.v5[i,],sigma=sigma.sp.3.v5,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.sp.3.v5[i,]=ti.thresh(sim.m.sp.3.v5[i,],sigma=sigma.sp.3.v5,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.sp.3.v5[i,]=ti.thresh(sim.m.sp.3.v5[i,],sigma=sigma.sp.3.v5,filter.number=8,family="DaubLeAsymm")

  sigma.bump.1.v5.est=sig.est.func(sim.m.bump.1.v5[i,],n)
  mu.est.ash.haar.bump.1.v5[i,]=ashsmooth.gaus(sim.m.bump.1.v5[i,],sigma.bump.1.v5.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.bump.1.v5[i,]=ashsmooth.gaus(sim.m.bump.1.v5[i,],sigma.bump.1.v5.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.bump.1.v5[i,]=ashsmooth.gaus(sim.m.bump.1.v5[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.bump.1.v5[i,]=ashsmooth.gaus(sim.m.bump.1.v5[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.bump.1.v5[i,]=ti.thresh(sim.m.bump.1.v5[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.bump.1.v5[i,]=ti.thresh(sim.m.bump.1.v5[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.bump.1.v5[i,]=ashsmooth.gaus(sim.m.bump.1.v5[i,],sigma=sigma.bump.1.v5,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.bump.1.v5[i,]=ashsmooth.gaus(sim.m.bump.1.v5[i,],sigma=sigma.bump.1.v5,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.bump.1.v5[i,]=ti.thresh(sim.m.bump.1.v5[i,],sigma=sigma.bump.1.v5,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.bump.1.v5[i,]=ti.thresh(sim.m.bump.1.v5[i,],sigma=sigma.bump.1.v5,filter.number=8,family="DaubLeAsymm")
  sigma.bump.3.v5.est=sig.est.func(sim.m.bump.3.v5[i,],n)
  mu.est.ash.haar.bump.3.v5[i,]=ashsmooth.gaus(sim.m.bump.3.v5[i,],sigma.bump.3.v5.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.bump.3.v5[i,]=ashsmooth.gaus(sim.m.bump.3.v5[i,],sigma.bump.3.v5.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.bump.3.v5[i,]=ashsmooth.gaus(sim.m.bump.3.v5[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.bump.3.v5[i,]=ashsmooth.gaus(sim.m.bump.3.v5[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.bump.3.v5[i,]=ti.thresh(sim.m.bump.3.v5[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.bump.3.v5[i,]=ti.thresh(sim.m.bump.3.v5[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.bump.3.v5[i,]=ashsmooth.gaus(sim.m.bump.3.v5[i,],sigma=sigma.bump.3.v5,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.bump.3.v5[i,]=ashsmooth.gaus(sim.m.bump.3.v5[i,],sigma=sigma.bump.3.v5,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.bump.3.v5[i,]=ti.thresh(sim.m.bump.3.v5[i,],sigma=sigma.bump.3.v5,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.bump.3.v5[i,]=ti.thresh(sim.m.bump.3.v5[i,],sigma=sigma.bump.3.v5,filter.number=8,family="DaubLeAsymm")
  
  sigma.blk.1.v5.est=sig.est.func(sim.m.blk.1.v5[i,],n)
  mu.est.ash.haar.blk.1.v5[i,]=ashsmooth.gaus(sim.m.blk.1.v5[i,],sigma.blk.1.v5.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.blk.1.v5[i,]=ashsmooth.gaus(sim.m.blk.1.v5[i,],sigma.blk.1.v5.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.blk.1.v5[i,]=ashsmooth.gaus(sim.m.blk.1.v5[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.blk.1.v5[i,]=ashsmooth.gaus(sim.m.blk.1.v5[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.blk.1.v5[i,]=ti.thresh(sim.m.blk.1.v5[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.blk.1.v5[i,]=ti.thresh(sim.m.blk.1.v5[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.blk.1.v5[i,]=ashsmooth.gaus(sim.m.blk.1.v5[i,],sigma=sigma.blk.1.v5,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.blk.1.v5[i,]=ashsmooth.gaus(sim.m.blk.1.v5[i,],sigma=sigma.blk.1.v5,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.blk.1.v5[i,]=ti.thresh(sim.m.blk.1.v5[i,],sigma=sigma.blk.1.v5,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.blk.1.v5[i,]=ti.thresh(sim.m.blk.1.v5[i,],sigma=sigma.blk.1.v5,filter.number=8,family="DaubLeAsymm")
  sigma.blk.3.v5.est=sig.est.func(sim.m.blk.3.v5[i,],n)
  mu.est.ash.haar.blk.3.v5[i,]=ashsmooth.gaus(sim.m.blk.3.v5[i,],sigma.blk.3.v5.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.blk.3.v5[i,]=ashsmooth.gaus(sim.m.blk.3.v5[i,],sigma.blk.3.v5.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.blk.3.v5[i,]=ashsmooth.gaus(sim.m.blk.3.v5[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.blk.3.v5[i,]=ashsmooth.gaus(sim.m.blk.3.v5[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.blk.3.v5[i,]=ti.thresh(sim.m.blk.3.v5[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.blk.3.v5[i,]=ti.thresh(sim.m.blk.3.v5[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.blk.3.v5[i,]=ashsmooth.gaus(sim.m.blk.3.v5[i,],sigma=sigma.blk.3.v5,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.blk.3.v5[i,]=ashsmooth.gaus(sim.m.blk.3.v5[i,],sigma=sigma.blk.3.v5,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.blk.3.v5[i,]=ti.thresh(sim.m.blk.3.v5[i,],sigma=sigma.blk.3.v5,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.blk.3.v5[i,]=ti.thresh(sim.m.blk.3.v5[i,],sigma=sigma.blk.3.v5,filter.number=8,family="DaubLeAsymm")

  sigma.ang.1.v5.est=sig.est.func(sim.m.ang.1.v5[i,],n)
  mu.est.ash.haar.ang.1.v5[i,]=ashsmooth.gaus(sim.m.ang.1.v5[i,],sigma.ang.1.v5.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.ang.1.v5[i,]=ashsmooth.gaus(sim.m.ang.1.v5[i,],sigma.ang.1.v5.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.ang.1.v5[i,]=ashsmooth.gaus(sim.m.ang.1.v5[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.ang.1.v5[i,]=ashsmooth.gaus(sim.m.ang.1.v5[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.ang.1.v5[i,]=ti.thresh(sim.m.ang.1.v5[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.ang.1.v5[i,]=ti.thresh(sim.m.ang.1.v5[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.ang.1.v5[i,]=ashsmooth.gaus(sim.m.ang.1.v5[i,],sigma=sigma.ang.1.v5,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.ang.1.v5[i,]=ashsmooth.gaus(sim.m.ang.1.v5[i,],sigma=sigma.ang.1.v5,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.ang.1.v5[i,]=ti.thresh(sim.m.ang.1.v5[i,],sigma=sigma.ang.1.v5,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.ang.1.v5[i,]=ti.thresh(sim.m.ang.1.v5[i,],sigma=sigma.ang.1.v5,filter.number=8,family="DaubLeAsymm")
  sigma.ang.3.v5.est=sig.est.func(sim.m.ang.3.v5[i,],n)
  mu.est.ash.haar.ang.3.v5[i,]=ashsmooth.gaus(sim.m.ang.3.v5[i,],sigma.ang.3.v5.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.ang.3.v5[i,]=ashsmooth.gaus(sim.m.ang.3.v5[i,],sigma.ang.3.v5.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.ang.3.v5[i,]=ashsmooth.gaus(sim.m.ang.3.v5[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.ang.3.v5[i,]=ashsmooth.gaus(sim.m.ang.3.v5[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.ang.3.v5[i,]=ti.thresh(sim.m.ang.3.v5[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.ang.3.v5[i,]=ti.thresh(sim.m.ang.3.v5[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.ang.3.v5[i,]=ashsmooth.gaus(sim.m.ang.3.v5[i,],sigma=sigma.ang.3.v5,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.ang.3.v5[i,]=ashsmooth.gaus(sim.m.ang.3.v5[i,],sigma=sigma.ang.3.v5,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.ang.3.v5[i,]=ti.thresh(sim.m.ang.3.v5[i,],sigma=sigma.ang.3.v5,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.ang.3.v5[i,]=ti.thresh(sim.m.ang.3.v5[i,],sigma=sigma.ang.3.v5,filter.number=8,family="DaubLeAsymm")

  sigma.dop.1.v5.est=sig.est.func(sim.m.dop.1.v5[i,],n)
  mu.est.ash.haar.dop.1.v5[i,]=ashsmooth.gaus(sim.m.dop.1.v5[i,],sigma.dop.1.v5.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.dop.1.v5[i,]=ashsmooth.gaus(sim.m.dop.1.v5[i,],sigma.dop.1.v5.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.dop.1.v5[i,]=ashsmooth.gaus(sim.m.dop.1.v5[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.dop.1.v5[i,]=ashsmooth.gaus(sim.m.dop.1.v5[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.dop.1.v5[i,]=ti.thresh(sim.m.dop.1.v5[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.dop.1.v5[i,]=ti.thresh(sim.m.dop.1.v5[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.dop.1.v5[i,]=ashsmooth.gaus(sim.m.dop.1.v5[i,],sigma=sigma.dop.1.v5,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.dop.1.v5[i,]=ashsmooth.gaus(sim.m.dop.1.v5[i,],sigma=sigma.dop.1.v5,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.dop.1.v5[i,]=ti.thresh(sim.m.dop.1.v5[i,],sigma=sigma.dop.1.v5,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.dop.1.v5[i,]=ti.thresh(sim.m.dop.1.v5[i,],sigma=sigma.dop.1.v5,filter.number=8,family="DaubLeAsymm")
  sigma.dop.3.v5.est=sig.est.func(sim.m.dop.3.v5[i,],n)
  mu.est.ash.haar.dop.3.v5[i,]=ashsmooth.gaus(sim.m.dop.3.v5[i,],sigma.dop.3.v5.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.dop.3.v5[i,]=ashsmooth.gaus(sim.m.dop.3.v5[i,],sigma.dop.3.v5.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.dop.3.v5[i,]=ashsmooth.gaus(sim.m.dop.3.v5[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.dop.3.v5[i,]=ashsmooth.gaus(sim.m.dop.3.v5[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.dop.3.v5[i,]=ti.thresh(sim.m.dop.3.v5[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.dop.3.v5[i,]=ti.thresh(sim.m.dop.3.v5[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.dop.3.v5[i,]=ashsmooth.gaus(sim.m.dop.3.v5[i,],sigma=sigma.dop.3.v5,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.dop.3.v5[i,]=ashsmooth.gaus(sim.m.dop.3.v5[i,],sigma=sigma.dop.3.v5,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.dop.3.v5[i,]=ti.thresh(sim.m.dop.3.v5[i,],sigma=sigma.dop.3.v5,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.dop.3.v5[i,]=ti.thresh(sim.m.dop.3.v5[i,],sigma=sigma.dop.3.v5,filter.number=8,family="DaubLeAsymm")

  sigma.blip.1.v5.est=sig.est.func(sim.m.blip.1.v5[i,],n)
  mu.est.ash.haar.blip.1.v5[i,]=ashsmooth.gaus(sim.m.blip.1.v5[i,],sigma.blip.1.v5.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.blip.1.v5[i,]=ashsmooth.gaus(sim.m.blip.1.v5[i,],sigma.blip.1.v5.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.blip.1.v5[i,]=ashsmooth.gaus(sim.m.blip.1.v5[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.blip.1.v5[i,]=ashsmooth.gaus(sim.m.blip.1.v5[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.blip.1.v5[i,]=ti.thresh(sim.m.blip.1.v5[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.blip.1.v5[i,]=ti.thresh(sim.m.blip.1.v5[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.blip.1.v5[i,]=ashsmooth.gaus(sim.m.blip.1.v5[i,],sigma=sigma.blip.1.v5,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.blip.1.v5[i,]=ashsmooth.gaus(sim.m.blip.1.v5[i,],sigma=sigma.blip.1.v5,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.blip.1.v5[i,]=ti.thresh(sim.m.blip.1.v5[i,],sigma=sigma.blip.1.v5,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.blip.1.v5[i,]=ti.thresh(sim.m.blip.1.v5[i,],sigma=sigma.blip.1.v5,filter.number=8,family="DaubLeAsymm")
  sigma.blip.3.v5.est=sig.est.func(sim.m.blip.3.v5[i,],n)
  mu.est.ash.haar.blip.3.v5[i,]=ashsmooth.gaus(sim.m.blip.3.v5[i,],sigma.blip.3.v5.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.blip.3.v5[i,]=ashsmooth.gaus(sim.m.blip.3.v5[i,],sigma.blip.3.v5.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.blip.3.v5[i,]=ashsmooth.gaus(sim.m.blip.3.v5[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.blip.3.v5[i,]=ashsmooth.gaus(sim.m.blip.3.v5[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.blip.3.v5[i,]=ti.thresh(sim.m.blip.3.v5[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.blip.3.v5[i,]=ti.thresh(sim.m.blip.3.v5[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.blip.3.v5[i,]=ashsmooth.gaus(sim.m.blip.3.v5[i,],sigma=sigma.blip.3.v5,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.blip.3.v5[i,]=ashsmooth.gaus(sim.m.blip.3.v5[i,],sigma=sigma.blip.3.v5,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.blip.3.v5[i,]=ti.thresh(sim.m.blip.3.v5[i,],sigma=sigma.blip.3.v5,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.blip.3.v5[i,]=ti.thresh(sim.m.blip.3.v5[i,],sigma=sigma.blip.3.v5,filter.number=8,family="DaubLeAsymm")

  sigma.cor.1.v5.est=sig.est.func(sim.m.cor.1.v5[i,],n)
  mu.est.ash.haar.cor.1.v5[i,]=ashsmooth.gaus(sim.m.cor.1.v5[i,],sigma.cor.1.v5.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.cor.1.v5[i,]=ashsmooth.gaus(sim.m.cor.1.v5[i,],sigma.cor.1.v5.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.cor.1.v5[i,]=ashsmooth.gaus(sim.m.cor.1.v5[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.cor.1.v5[i,]=ashsmooth.gaus(sim.m.cor.1.v5[i,],filter.number=8,family="DaubLeAsymm",)
  mu.est.tie.haar.cor.1.v5[i,]=ti.thresh(sim.m.cor.1.v5[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.cor.1.v5[i,]=ti.thresh(sim.m.cor.1.v5[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.cor.1.v5[i,]=ashsmooth.gaus(sim.m.cor.1.v5[i,],sigma=sigma.cor.1.v5,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.cor.1.v5[i,]=ashsmooth.gaus(sim.m.cor.1.v5[i,],sigma=sigma.cor.1.v5,filter.number=8,family="DaubLeAsymm")
  mu.est.tit.haar.cor.1.v5[i,]=ti.thresh(sim.m.cor.1.v5[i,],sigma=sigma.cor.1.v5,filter.number=1,family="DaubExPhase")
  mu.est.tit.s8.cor.1.v5[i,]=ti.thresh(sim.m.cor.1.v5[i,],sigma=sigma.cor.1.v5,filter.number=8,family="DaubLeAsymm")
  sigma.cor.3.v5.est=sig.est.func(sim.m.cor.3.v5[i,],n)
  mu.est.ash.haar.cor.3.v5[i,]=ashsmooth.gaus(sim.m.cor.3.v5[i,],sigma.cor.3.v5.est,filter.number=1,family="DaubExPhase")
  mu.est.ash.s8.cor.3.v5[i,]=ashsmooth.gaus(sim.m.cor.3.v5[i,],sigma.cor.3.v5.est,filter.number=8,family="DaubLeAsymm")
  mu.est.ashe.haar.cor.3.v5[i,]=ashsmooth.gaus(sim.m.cor.3.v5[i,],filter.number=1,family="DaubExPhase")
  mu.est.ashe.s8.cor.3.v5[i,]=ashsmooth.gaus(sim.m.cor.3.v5[i,],filter.number=8,family="DaubLeAsymm")
  mu.est.tie.haar.cor.3.v5[i,]=ti.thresh(sim.m.cor.3.v5[i,],method="rmad",filter.number=1,family="DaubExPhase")
  mu.est.tie.s8.cor.3.v5[i,]=ti.thresh(sim.m.cor.3.v5[i,],method="rmad",filter.number=8,family="DaubLeAsymm")
  mu.est.asht.haar.cor.3.v5[i,]=ashsmooth.gaus(sim.m.cor.3.v5[i,],sigma=sigma.cor.3.v5,filter.number=1,family="DaubExPhase")
  mu.est.asht.s8.cor.3.v5[i,]=ashsmooth.gaus(sim.m.cor.3.v5[i,],sigma=sigma.cor.3.v5,filter.number=8,family="DaubLeAsymm")
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

  mu.est.tieb.haar.sp.1.v1[i,]=ti.thresh(sim.m.sp.1.v1[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.sp.1.v1[i,]=ti.thresh(sim.m.sp.1.v1[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.sp.1.v2[i,]=ti.thresh(sim.m.sp.1.v2[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.sp.1.v2[i,]=ti.thresh(sim.m.sp.1.v2[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.sp.1.v3[i,]=ti.thresh(sim.m.sp.1.v3[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.sp.1.v3[i,]=ti.thresh(sim.m.sp.1.v3[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.sp.1.v4[i,]=ti.thresh(sim.m.sp.1.v4[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.sp.1.v4[i,]=ti.thresh(sim.m.sp.1.v4[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.sp.1.v5[i,]=ti.thresh(sim.m.sp.1.v5[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.sp.1.v5[i,]=ti.thresh(sim.m.sp.1.v5[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.sp.3.v1[i,]=ti.thresh(sim.m.sp.3.v1[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.sp.3.v1[i,]=ti.thresh(sim.m.sp.3.v1[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.sp.3.v2[i,]=ti.thresh(sim.m.sp.3.v2[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.sp.3.v2[i,]=ti.thresh(sim.m.sp.3.v2[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.sp.3.v3[i,]=ti.thresh(sim.m.sp.3.v3[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.sp.3.v3[i,]=ti.thresh(sim.m.sp.3.v3[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.sp.3.v4[i,]=ti.thresh(sim.m.sp.3.v4[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.sp.3.v4[i,]=ti.thresh(sim.m.sp.3.v4[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.sp.3.v5[i,]=ti.thresh(sim.m.sp.3.v5[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.sp.3.v5[i,]=ti.thresh(sim.m.sp.3.v5[i,],method="smash",filter.number=8,family="DaubLeAsymm")


  mu.est.tieb.haar.bump.1.v1[i,]=ti.thresh(sim.m.bump.1.v1[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.bump.1.v1[i,]=ti.thresh(sim.m.bump.1.v1[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.bump.1.v2[i,]=ti.thresh(sim.m.bump.1.v2[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.bump.1.v2[i,]=ti.thresh(sim.m.bump.1.v2[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.bump.1.v3[i,]=ti.thresh(sim.m.bump.1.v3[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.bump.1.v3[i,]=ti.thresh(sim.m.bump.1.v3[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.bump.1.v4[i,]=ti.thresh(sim.m.bump.1.v4[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.bump.1.v4[i,]=ti.thresh(sim.m.bump.1.v4[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.bump.1.v5[i,]=ti.thresh(sim.m.bump.1.v5[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.bump.1.v5[i,]=ti.thresh(sim.m.bump.1.v5[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.bump.3.v1[i,]=ti.thresh(sim.m.bump.3.v1[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.bump.3.v1[i,]=ti.thresh(sim.m.bump.3.v1[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.bump.3.v2[i,]=ti.thresh(sim.m.bump.3.v2[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.bump.3.v2[i,]=ti.thresh(sim.m.bump.3.v2[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.bump.3.v3[i,]=ti.thresh(sim.m.bump.3.v3[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.bump.3.v3[i,]=ti.thresh(sim.m.bump.3.v3[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.bump.3.v4[i,]=ti.thresh(sim.m.bump.3.v4[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.bump.3.v4[i,]=ti.thresh(sim.m.bump.3.v4[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.bump.3.v5[i,]=ti.thresh(sim.m.bump.3.v5[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.bump.3.v5[i,]=ti.thresh(sim.m.bump.3.v5[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.blk.1.v1[i,]=ti.thresh(sim.m.blk.1.v1[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.blk.1.v1[i,]=ti.thresh(sim.m.blk.1.v1[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.blk.1.v2[i,]=ti.thresh(sim.m.blk.1.v2[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.blk.1.v2[i,]=ti.thresh(sim.m.blk.1.v2[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.blk.1.v3[i,]=ti.thresh(sim.m.blk.1.v3[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.blk.1.v3[i,]=ti.thresh(sim.m.blk.1.v3[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.blk.1.v4[i,]=ti.thresh(sim.m.blk.1.v4[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.blk.1.v4[i,]=ti.thresh(sim.m.blk.1.v4[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.blk.1.v5[i,]=ti.thresh(sim.m.blk.1.v5[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.blk.1.v5[i,]=ti.thresh(sim.m.blk.1.v5[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.blk.3.v1[i,]=ti.thresh(sim.m.blk.3.v1[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.blk.3.v1[i,]=ti.thresh(sim.m.blk.3.v1[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.blk.3.v2[i,]=ti.thresh(sim.m.blk.3.v2[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.blk.3.v2[i,]=ti.thresh(sim.m.blk.3.v2[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.blk.3.v3[i,]=ti.thresh(sim.m.blk.3.v3[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.blk.3.v3[i,]=ti.thresh(sim.m.blk.3.v3[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.blk.3.v4[i,]=ti.thresh(sim.m.blk.3.v4[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.blk.3.v4[i,]=ti.thresh(sim.m.blk.3.v4[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.blk.3.v5[i,]=ti.thresh(sim.m.blk.3.v5[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.blk.3.v5[i,]=ti.thresh(sim.m.blk.3.v5[i,],method="smash",filter.number=8,family="DaubLeAsymm")


  mu.est.tieb.haar.ang.1.v1[i,]=ti.thresh(sim.m.ang.1.v1[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.ang.1.v1[i,]=ti.thresh(sim.m.ang.1.v1[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.ang.1.v2[i,]=ti.thresh(sim.m.ang.1.v2[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.ang.1.v2[i,]=ti.thresh(sim.m.ang.1.v2[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.ang.1.v3[i,]=ti.thresh(sim.m.ang.1.v3[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.ang.1.v3[i,]=ti.thresh(sim.m.ang.1.v3[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.ang.1.v4[i,]=ti.thresh(sim.m.ang.1.v4[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.ang.1.v4[i,]=ti.thresh(sim.m.ang.1.v4[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.ang.1.v5[i,]=ti.thresh(sim.m.ang.1.v5[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.ang.1.v5[i,]=ti.thresh(sim.m.ang.1.v5[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.ang.3.v1[i,]=ti.thresh(sim.m.ang.3.v1[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.ang.3.v1[i,]=ti.thresh(sim.m.ang.3.v1[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.ang.3.v2[i,]=ti.thresh(sim.m.ang.3.v2[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.ang.3.v2[i,]=ti.thresh(sim.m.ang.3.v2[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.ang.3.v3[i,]=ti.thresh(sim.m.ang.3.v3[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.ang.3.v3[i,]=ti.thresh(sim.m.ang.3.v3[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.ang.3.v4[i,]=ti.thresh(sim.m.ang.3.v4[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.ang.3.v4[i,]=ti.thresh(sim.m.ang.3.v4[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.ang.3.v5[i,]=ti.thresh(sim.m.ang.3.v5[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.ang.3.v5[i,]=ti.thresh(sim.m.ang.3.v5[i,],method="smash",filter.number=8,family="DaubLeAsymm")


  mu.est.tieb.haar.dop.1.v1[i,]=ti.thresh(sim.m.dop.1.v1[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.dop.1.v1[i,]=ti.thresh(sim.m.dop.1.v1[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.dop.1.v2[i,]=ti.thresh(sim.m.dop.1.v2[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.dop.1.v2[i,]=ti.thresh(sim.m.dop.1.v2[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.dop.1.v3[i,]=ti.thresh(sim.m.dop.1.v3[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.dop.1.v3[i,]=ti.thresh(sim.m.dop.1.v3[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.dop.1.v4[i,]=ti.thresh(sim.m.dop.1.v4[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.dop.1.v4[i,]=ti.thresh(sim.m.dop.1.v4[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.dop.1.v5[i,]=ti.thresh(sim.m.dop.1.v5[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.dop.1.v5[i,]=ti.thresh(sim.m.dop.1.v5[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.dop.3.v1[i,]=ti.thresh(sim.m.dop.3.v1[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.dop.3.v1[i,]=ti.thresh(sim.m.dop.3.v1[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.dop.3.v2[i,]=ti.thresh(sim.m.dop.3.v2[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.dop.3.v2[i,]=ti.thresh(sim.m.dop.3.v2[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.dop.3.v3[i,]=ti.thresh(sim.m.dop.3.v3[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.dop.3.v3[i,]=ti.thresh(sim.m.dop.3.v3[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.dop.3.v4[i,]=ti.thresh(sim.m.dop.3.v4[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.dop.3.v4[i,]=ti.thresh(sim.m.dop.3.v4[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.dop.3.v5[i,]=ti.thresh(sim.m.dop.3.v5[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.dop.3.v5[i,]=ti.thresh(sim.m.dop.3.v5[i,],method="smash",filter.number=8,family="DaubLeAsymm")


  mu.est.tieb.haar.blip.1.v1[i,]=ti.thresh(sim.m.blip.1.v1[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.blip.1.v1[i,]=ti.thresh(sim.m.blip.1.v1[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.blip.1.v2[i,]=ti.thresh(sim.m.blip.1.v2[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.blip.1.v2[i,]=ti.thresh(sim.m.blip.1.v2[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.blip.1.v3[i,]=ti.thresh(sim.m.blip.1.v3[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.blip.1.v3[i,]=ti.thresh(sim.m.blip.1.v3[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.blip.1.v4[i,]=ti.thresh(sim.m.blip.1.v4[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.blip.1.v4[i,]=ti.thresh(sim.m.blip.1.v4[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.blip.1.v5[i,]=ti.thresh(sim.m.blip.1.v5[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.blip.1.v5[i,]=ti.thresh(sim.m.blip.1.v5[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.blip.3.v1[i,]=ti.thresh(sim.m.blip.3.v1[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.blip.3.v1[i,]=ti.thresh(sim.m.blip.3.v1[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.blip.3.v2[i,]=ti.thresh(sim.m.blip.3.v2[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.blip.3.v2[i,]=ti.thresh(sim.m.blip.3.v2[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.blip.3.v3[i,]=ti.thresh(sim.m.blip.3.v3[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.blip.3.v3[i,]=ti.thresh(sim.m.blip.3.v3[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.blip.3.v4[i,]=ti.thresh(sim.m.blip.3.v4[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.blip.3.v4[i,]=ti.thresh(sim.m.blip.3.v4[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.blip.3.v5[i,]=ti.thresh(sim.m.blip.3.v5[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.blip.3.v5[i,]=ti.thresh(sim.m.blip.3.v5[i,],method="smash",filter.number=8,family="DaubLeAsymm")



  mu.est.tieb.haar.cor.1.v1[i,]=ti.thresh(sim.m.cor.1.v1[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.cor.1.v1[i,]=ti.thresh(sim.m.cor.1.v1[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.cor.1.v2[i,]=ti.thresh(sim.m.cor.1.v2[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.cor.1.v2[i,]=ti.thresh(sim.m.cor.1.v2[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.cor.1.v3[i,]=ti.thresh(sim.m.cor.1.v3[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.cor.1.v3[i,]=ti.thresh(sim.m.cor.1.v3[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.cor.1.v4[i,]=ti.thresh(sim.m.cor.1.v4[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.cor.1.v4[i,]=ti.thresh(sim.m.cor.1.v4[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.cor.1.v5[i,]=ti.thresh(sim.m.cor.1.v5[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.cor.1.v5[i,]=ti.thresh(sim.m.cor.1.v5[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.cor.3.v1[i,]=ti.thresh(sim.m.cor.3.v1[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.cor.3.v1[i,]=ti.thresh(sim.m.cor.3.v1[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.cor.3.v2[i,]=ti.thresh(sim.m.cor.3.v2[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.cor.3.v2[i,]=ti.thresh(sim.m.cor.3.v2[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.cor.3.v3[i,]=ti.thresh(sim.m.cor.3.v3[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.cor.3.v3[i,]=ti.thresh(sim.m.cor.3.v3[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.cor.3.v4[i,]=ti.thresh(sim.m.cor.3.v4[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.cor.3.v4[i,]=ti.thresh(sim.m.cor.3.v4[i,],method="smash",filter.number=8,family="DaubLeAsymm")

  mu.est.tieb.haar.cor.3.v5[i,]=ti.thresh(sim.m.cor.3.v5[i,],method="smash",filter.number=1,family="DaubExPhase")
  mu.est.tieb.s8.cor.3.v5[i,]=ti.thresh(sim.m.cor.3.v5[i,],method="smash",filter.number=8,family="DaubLeAsymm")



  mu.est.ashe.s8.j.sp.1.v1[i,]=ashsmooth.gaus(sim.m.sp.1.v1[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.sp.1.v2[i,]=ashsmooth.gaus(sim.m.sp.1.v2[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.sp.1.v3[i,]=ashsmooth.gaus(sim.m.sp.1.v3[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.sp.1.v4[i,]=ashsmooth.gaus(sim.m.sp.1.v4[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.sp.1.v5[i,]=ashsmooth.gaus(sim.m.sp.1.v5[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.sp.3.v1[i,]=ashsmooth.gaus(sim.m.sp.3.v1[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.sp.3.v2[i,]=ashsmooth.gaus(sim.m.sp.3.v2[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.sp.3.v3[i,]=ashsmooth.gaus(sim.m.sp.3.v3[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.sp.3.v4[i,]=ashsmooth.gaus(sim.m.sp.3.v4[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.sp.3.v5[i,]=ashsmooth.gaus(sim.m.sp.3.v5[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)


  mu.est.ashe.s8.j.bump.1.v1[i,]=ashsmooth.gaus(sim.m.bump.1.v1[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.bump.1.v2[i,]=ashsmooth.gaus(sim.m.bump.1.v2[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.bump.1.v3[i,]=ashsmooth.gaus(sim.m.bump.1.v3[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.bump.1.v4[i,]=ashsmooth.gaus(sim.m.bump.1.v4[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.bump.1.v5[i,]=ashsmooth.gaus(sim.m.bump.1.v5[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.bump.3.v1[i,]=ashsmooth.gaus(sim.m.bump.3.v1[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.bump.3.v2[i,]=ashsmooth.gaus(sim.m.bump.3.v2[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.bump.3.v3[i,]=ashsmooth.gaus(sim.m.bump.3.v3[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.bump.3.v4[i,]=ashsmooth.gaus(sim.m.bump.3.v4[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.bump.3.v5[i,]=ashsmooth.gaus(sim.m.bump.3.v5[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.blk.1.v1[i,]=ashsmooth.gaus(sim.m.blk.1.v1[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.blk.1.v2[i,]=ashsmooth.gaus(sim.m.blk.1.v2[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.blk.1.v3[i,]=ashsmooth.gaus(sim.m.blk.1.v3[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.blk.1.v4[i,]=ashsmooth.gaus(sim.m.blk.1.v4[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.blk.1.v5[i,]=ashsmooth.gaus(sim.m.blk.1.v5[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.blk.3.v1[i,]=ashsmooth.gaus(sim.m.blk.3.v1[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.blk.3.v2[i,]=ashsmooth.gaus(sim.m.blk.3.v2[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.blk.3.v3[i,]=ashsmooth.gaus(sim.m.blk.3.v3[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.blk.3.v4[i,]=ashsmooth.gaus(sim.m.blk.3.v4[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.blk.3.v5[i,]=ashsmooth.gaus(sim.m.blk.3.v5[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)


  mu.est.ashe.s8.j.ang.1.v1[i,]=ashsmooth.gaus(sim.m.ang.1.v1[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.ang.1.v2[i,]=ashsmooth.gaus(sim.m.ang.1.v2[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.ang.1.v3[i,]=ashsmooth.gaus(sim.m.ang.1.v3[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.ang.1.v4[i,]=ashsmooth.gaus(sim.m.ang.1.v4[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.ang.1.v5[i,]=ashsmooth.gaus(sim.m.ang.1.v5[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.ang.3.v1[i,]=ashsmooth.gaus(sim.m.ang.3.v1[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.ang.3.v2[i,]=ashsmooth.gaus(sim.m.ang.3.v2[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.ang.3.v3[i,]=ashsmooth.gaus(sim.m.ang.3.v3[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.ang.3.v4[i,]=ashsmooth.gaus(sim.m.ang.3.v4[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.ang.3.v5[i,]=ashsmooth.gaus(sim.m.ang.3.v5[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)


  mu.est.ashe.s8.j.dop.1.v1[i,]=ashsmooth.gaus(sim.m.dop.1.v1[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.dop.1.v2[i,]=ashsmooth.gaus(sim.m.dop.1.v2[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.dop.1.v3[i,]=ashsmooth.gaus(sim.m.dop.1.v3[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.dop.1.v4[i,]=ashsmooth.gaus(sim.m.dop.1.v4[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.dop.1.v5[i,]=ashsmooth.gaus(sim.m.dop.1.v5[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.dop.3.v1[i,]=ashsmooth.gaus(sim.m.dop.3.v1[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.dop.3.v2[i,]=ashsmooth.gaus(sim.m.dop.3.v2[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.dop.3.v3[i,]=ashsmooth.gaus(sim.m.dop.3.v3[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.dop.3.v4[i,]=ashsmooth.gaus(sim.m.dop.3.v4[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.dop.3.v5[i,]=ashsmooth.gaus(sim.m.dop.3.v5[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)


  mu.est.ashe.s8.j.blip.1.v1[i,]=ashsmooth.gaus(sim.m.blip.1.v1[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.blip.1.v2[i,]=ashsmooth.gaus(sim.m.blip.1.v2[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.blip.1.v3[i,]=ashsmooth.gaus(sim.m.blip.1.v3[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.blip.1.v4[i,]=ashsmooth.gaus(sim.m.blip.1.v4[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.blip.1.v5[i,]=ashsmooth.gaus(sim.m.blip.1.v5[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.blip.3.v1[i,]=ashsmooth.gaus(sim.m.blip.3.v1[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.blip.3.v2[i,]=ashsmooth.gaus(sim.m.blip.3.v2[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.blip.3.v3[i,]=ashsmooth.gaus(sim.m.blip.3.v3[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.blip.3.v4[i,]=ashsmooth.gaus(sim.m.blip.3.v4[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.blip.3.v5[i,]=ashsmooth.gaus(sim.m.blip.3.v5[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)



  mu.est.ashe.s8.j.cor.1.v1[i,]=ashsmooth.gaus(sim.m.cor.1.v1[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.cor.1.v2[i,]=ashsmooth.gaus(sim.m.cor.1.v2[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.cor.1.v3[i,]=ashsmooth.gaus(sim.m.cor.1.v3[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.cor.1.v4[i,]=ashsmooth.gaus(sim.m.cor.1.v4[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.cor.1.v5[i,]=ashsmooth.gaus(sim.m.cor.1.v5[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.cor.3.v1[i,]=ashsmooth.gaus(sim.m.cor.3.v1[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.cor.3.v2[i,]=ashsmooth.gaus(sim.m.cor.3.v2[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.cor.3.v3[i,]=ashsmooth.gaus(sim.m.cor.3.v3[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.cor.3.v4[i,]=ashsmooth.gaus(sim.m.cor.3.v4[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  mu.est.ashe.s8.j.cor.3.v5[i,]=ashsmooth.gaus(sim.m.cor.3.v5[i,],v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE)

  print(i)
}



setwd("D:/Grad School/Spring 2013/multiscale_ash/simulation_1d_g/n_1024")

mu.est.ti.sp.1.v1=data.matrix(read.csv("est_ti_sp_1_v1.csv",header=FALSE))
mu.est.ti.sp.1.v2=data.matrix(read.csv("est_ti_sp_1_v2.csv",header=FALSE))
mu.est.ti.sp.1.v3=data.matrix(read.csv("est_ti_sp_1_v3.csv",header=FALSE))
mu.est.ti.sp.1.v4=data.matrix(read.csv("est_ti_sp_1_v4.csv",header=FALSE))
mu.est.ti.sp.1.v5=data.matrix(read.csv("est_ti_sp_1_v5.csv",header=FALSE))
mu.est.ti.sp.3.v1=data.matrix(read.csv("est_ti_sp_3_v1.csv",header=FALSE))
mu.est.ti.sp.3.v2=data.matrix(read.csv("est_ti_sp_3_v2.csv",header=FALSE))
mu.est.ti.sp.3.v3=data.matrix(read.csv("est_ti_sp_3_v3.csv",header=FALSE))
mu.est.ti.sp.3.v4=data.matrix(read.csv("est_ti_sp_3_v4.csv",header=FALSE))
mu.est.ti.sp.3.v5=data.matrix(read.csv("est_ti_sp_3_v5.csv",header=FALSE))

mu.est.ti.bump.1.v1=data.matrix(read.csv("est_ti_bump_1_v1.csv",header=FALSE))
mu.est.ti.bump.1.v2=data.matrix(read.csv("est_ti_bump_1_v2.csv",header=FALSE))
mu.est.ti.bump.1.v3=data.matrix(read.csv("est_ti_bump_1_v3.csv",header=FALSE))
mu.est.ti.bump.1.v4=data.matrix(read.csv("est_ti_bump_1_v4.csv",header=FALSE))
mu.est.ti.bump.1.v5=data.matrix(read.csv("est_ti_bump_1_v5.csv",header=FALSE))
mu.est.ti.bump.3.v1=data.matrix(read.csv("est_ti_bump_3_v1.csv",header=FALSE))
mu.est.ti.bump.3.v2=data.matrix(read.csv("est_ti_bump_3_v2.csv",header=FALSE))
mu.est.ti.bump.3.v3=data.matrix(read.csv("est_ti_bump_3_v3.csv",header=FALSE))
mu.est.ti.bump.3.v4=data.matrix(read.csv("est_ti_bump_3_v4.csv",header=FALSE))
mu.est.ti.bump.3.v5=data.matrix(read.csv("est_ti_bump_3_v5.csv",header=FALSE))

mu.est.ti.blk.1.v1=data.matrix(read.csv("est_ti_blk_1_v1.csv",header=FALSE))
mu.est.ti.blk.1.v2=data.matrix(read.csv("est_ti_blk_1_v2.csv",header=FALSE))
mu.est.ti.blk.1.v3=data.matrix(read.csv("est_ti_blk_1_v3.csv",header=FALSE))
mu.est.ti.blk.1.v4=data.matrix(read.csv("est_ti_blk_1_v4.csv",header=FALSE))
mu.est.ti.blk.1.v5=data.matrix(read.csv("est_ti_blk_1_v5.csv",header=FALSE))
mu.est.ti.blk.3.v1=data.matrix(read.csv("est_ti_blk_3_v1.csv",header=FALSE))
mu.est.ti.blk.3.v2=data.matrix(read.csv("est_ti_blk_3_v2.csv",header=FALSE))
mu.est.ti.blk.3.v3=data.matrix(read.csv("est_ti_blk_3_v3.csv",header=FALSE))
mu.est.ti.blk.3.v4=data.matrix(read.csv("est_ti_blk_3_v4.csv",header=FALSE))
mu.est.ti.blk.3.v5=data.matrix(read.csv("est_ti_blk_3_v5.csv",header=FALSE))

mu.est.ti.ang.1.v1=data.matrix(read.csv("est_ti_ang_1_v1.csv",header=FALSE))
mu.est.ti.ang.1.v2=data.matrix(read.csv("est_ti_ang_1_v2.csv",header=FALSE))
mu.est.ti.ang.1.v3=data.matrix(read.csv("est_ti_ang_1_v3.csv",header=FALSE))
mu.est.ti.ang.1.v4=data.matrix(read.csv("est_ti_ang_1_v4.csv",header=FALSE))
mu.est.ti.ang.1.v5=data.matrix(read.csv("est_ti_ang_1_v5.csv",header=FALSE))
mu.est.ti.ang.3.v1=data.matrix(read.csv("est_ti_ang_3_v1.csv",header=FALSE))
mu.est.ti.ang.3.v2=data.matrix(read.csv("est_ti_ang_3_v2.csv",header=FALSE))
mu.est.ti.ang.3.v3=data.matrix(read.csv("est_ti_ang_3_v3.csv",header=FALSE))
mu.est.ti.ang.3.v4=data.matrix(read.csv("est_ti_ang_3_v4.csv",header=FALSE))
mu.est.ti.ang.3.v5=data.matrix(read.csv("est_ti_ang_3_v5.csv",header=FALSE))

mu.est.ti.dop.1.v1=data.matrix(read.csv("est_ti_dop_1_v1.csv",header=FALSE))
mu.est.ti.dop.1.v2=data.matrix(read.csv("est_ti_dop_1_v2.csv",header=FALSE))
mu.est.ti.dop.1.v3=data.matrix(read.csv("est_ti_dop_1_v3.csv",header=FALSE))
mu.est.ti.dop.1.v4=data.matrix(read.csv("est_ti_dop_1_v4.csv",header=FALSE))
mu.est.ti.dop.1.v5=data.matrix(read.csv("est_ti_dop_1_v5.csv",header=FALSE))
mu.est.ti.dop.3.v1=data.matrix(read.csv("est_ti_dop_3_v1.csv",header=FALSE))
mu.est.ti.dop.3.v2=data.matrix(read.csv("est_ti_dop_3_v2.csv",header=FALSE))
mu.est.ti.dop.3.v3=data.matrix(read.csv("est_ti_dop_3_v3.csv",header=FALSE))
mu.est.ti.dop.3.v4=data.matrix(read.csv("est_ti_dop_3_v4.csv",header=FALSE))
mu.est.ti.dop.3.v5=data.matrix(read.csv("est_ti_dop_3_v5.csv",header=FALSE))

mu.est.ti.blip.1.v1=data.matrix(read.csv("est_ti_blip_1_v1.csv",header=FALSE))
mu.est.ti.blip.1.v2=data.matrix(read.csv("est_ti_blip_1_v2.csv",header=FALSE))
mu.est.ti.blip.1.v3=data.matrix(read.csv("est_ti_blip_1_v3.csv",header=FALSE))
mu.est.ti.blip.1.v4=data.matrix(read.csv("est_ti_blip_1_v4.csv",header=FALSE))
mu.est.ti.blip.1.v5=data.matrix(read.csv("est_ti_blip_1_v5.csv",header=FALSE))
mu.est.ti.blip.3.v1=data.matrix(read.csv("est_ti_blip_3_v1.csv",header=FALSE))
mu.est.ti.blip.3.v2=data.matrix(read.csv("est_ti_blip_3_v2.csv",header=FALSE))
mu.est.ti.blip.3.v3=data.matrix(read.csv("est_ti_blip_3_v3.csv",header=FALSE))
mu.est.ti.blip.3.v4=data.matrix(read.csv("est_ti_blip_3_v4.csv",header=FALSE))
mu.est.ti.blip.3.v5=data.matrix(read.csv("est_ti_blip_3_v5.csv",header=FALSE))

mu.est.ti.cor.1.v1=data.matrix(read.csv("est_ti_cor_1_v1.csv",header=FALSE))
mu.est.ti.cor.1.v2=data.matrix(read.csv("est_ti_cor_1_v2.csv",header=FALSE))
mu.est.ti.cor.1.v3=data.matrix(read.csv("est_ti_cor_1_v3.csv",header=FALSE))
mu.est.ti.cor.1.v4=data.matrix(read.csv("est_ti_cor_1_v4.csv",header=FALSE))
mu.est.ti.cor.1.v5=data.matrix(read.csv("est_ti_cor_1_v5.csv",header=FALSE))
mu.est.ti.cor.3.v1=data.matrix(read.csv("est_ti_cor_3_v1.csv",header=FALSE))
mu.est.ti.cor.3.v2=data.matrix(read.csv("est_ti_cor_3_v2.csv",header=FALSE))
mu.est.ti.cor.3.v3=data.matrix(read.csv("est_ti_cor_3_v3.csv",header=FALSE))
mu.est.ti.cor.3.v4=data.matrix(read.csv("est_ti_cor_3_v4.csv",header=FALSE))
mu.est.ti.cor.3.v5=data.matrix(read.csv("est_ti_cor_3_v5.csv",header=FALSE))



mu.est.js.sp.1.v1=data.matrix(read.csv("est_js_sp_1_v1.csv",header=FALSE))
mu.est.js.sp.1.v2=data.matrix(read.csv("est_js_sp_1_v2.csv",header=FALSE))
mu.est.js.sp.1.v3=data.matrix(read.csv("est_js_sp_1_v3.csv",header=FALSE))
mu.est.js.sp.1.v4=data.matrix(read.csv("est_js_sp_1_v4.csv",header=FALSE))
mu.est.js.sp.1.v5=data.matrix(read.csv("est_js_sp_1_v5.csv",header=FALSE))
mu.est.js.sp.3.v1=data.matrix(read.csv("est_js_sp_3_v1.csv",header=FALSE))
mu.est.js.sp.3.v2=data.matrix(read.csv("est_js_sp_3_v2.csv",header=FALSE))
mu.est.js.sp.3.v3=data.matrix(read.csv("est_js_sp_3_v3.csv",header=FALSE))
mu.est.js.sp.3.v4=data.matrix(read.csv("est_js_sp_3_v4.csv",header=FALSE))
mu.est.js.sp.3.v5=data.matrix(read.csv("est_js_sp_3_v5.csv",header=FALSE))

mu.est.js.bump.1.v1=data.matrix(read.csv("est_js_bump_1_v1.csv",header=FALSE))
mu.est.js.bump.1.v2=data.matrix(read.csv("est_js_bump_1_v2.csv",header=FALSE))
mu.est.js.bump.1.v3=data.matrix(read.csv("est_js_bump_1_v3.csv",header=FALSE))
mu.est.js.bump.1.v4=data.matrix(read.csv("est_js_bump_1_v4.csv",header=FALSE))
mu.est.js.bump.1.v5=data.matrix(read.csv("est_js_bump_1_v5.csv",header=FALSE))
mu.est.js.bump.3.v1=data.matrix(read.csv("est_js_bump_3_v1.csv",header=FALSE))
mu.est.js.bump.3.v2=data.matrix(read.csv("est_js_bump_3_v2.csv",header=FALSE))
mu.est.js.bump.3.v3=data.matrix(read.csv("est_js_bump_3_v3.csv",header=FALSE))
mu.est.js.bump.3.v4=data.matrix(read.csv("est_js_bump_3_v4.csv",header=FALSE))
mu.est.js.bump.3.v5=data.matrix(read.csv("est_js_bump_3_v5.csv",header=FALSE))

mu.est.js.blk.1.v1=data.matrix(read.csv("est_js_blk_1_v1.csv",header=FALSE))
mu.est.js.blk.1.v2=data.matrix(read.csv("est_js_blk_1_v2.csv",header=FALSE))
mu.est.js.blk.1.v3=data.matrix(read.csv("est_js_blk_1_v3.csv",header=FALSE))
mu.est.js.blk.1.v4=data.matrix(read.csv("est_js_blk_1_v4.csv",header=FALSE))
mu.est.js.blk.1.v5=data.matrix(read.csv("est_js_blk_1_v5.csv",header=FALSE))
mu.est.js.blk.3.v1=data.matrix(read.csv("est_js_blk_3_v1.csv",header=FALSE))
mu.est.js.blk.3.v2=data.matrix(read.csv("est_js_blk_3_v2.csv",header=FALSE))
mu.est.js.blk.3.v3=data.matrix(read.csv("est_js_blk_3_v3.csv",header=FALSE))
mu.est.js.blk.3.v4=data.matrix(read.csv("est_js_blk_3_v4.csv",header=FALSE))
mu.est.js.blk.3.v5=data.matrix(read.csv("est_js_blk_3_v5.csv",header=FALSE))

mu.est.js.ang.1.v1=data.matrix(read.csv("est_js_ang_1_v1.csv",header=FALSE))
mu.est.js.ang.1.v2=data.matrix(read.csv("est_js_ang_1_v2.csv",header=FALSE))
mu.est.js.ang.1.v3=data.matrix(read.csv("est_js_ang_1_v3.csv",header=FALSE))
mu.est.js.ang.1.v4=data.matrix(read.csv("est_js_ang_1_v4.csv",header=FALSE))
mu.est.js.ang.1.v5=data.matrix(read.csv("est_js_ang_1_v5.csv",header=FALSE))
mu.est.js.ang.3.v1=data.matrix(read.csv("est_js_ang_3_v1.csv",header=FALSE))
mu.est.js.ang.3.v2=data.matrix(read.csv("est_js_ang_3_v2.csv",header=FALSE))
mu.est.js.ang.3.v3=data.matrix(read.csv("est_js_ang_3_v3.csv",header=FALSE))
mu.est.js.ang.3.v4=data.matrix(read.csv("est_js_ang_3_v4.csv",header=FALSE))
mu.est.js.ang.3.v5=data.matrix(read.csv("est_js_ang_3_v5.csv",header=FALSE))

mu.est.js.dop.1.v1=data.matrix(read.csv("est_js_dop_1_v1.csv",header=FALSE))
mu.est.js.dop.1.v2=data.matrix(read.csv("est_js_dop_1_v2.csv",header=FALSE))
mu.est.js.dop.1.v3=data.matrix(read.csv("est_js_dop_1_v3.csv",header=FALSE))
mu.est.js.dop.1.v4=data.matrix(read.csv("est_js_dop_1_v4.csv",header=FALSE))
mu.est.js.dop.1.v5=data.matrix(read.csv("est_js_dop_1_v5.csv",header=FALSE))
mu.est.js.dop.3.v1=data.matrix(read.csv("est_js_dop_3_v1.csv",header=FALSE))
mu.est.js.dop.3.v2=data.matrix(read.csv("est_js_dop_3_v2.csv",header=FALSE))
mu.est.js.dop.3.v3=data.matrix(read.csv("est_js_dop_3_v3.csv",header=FALSE))
mu.est.js.dop.3.v4=data.matrix(read.csv("est_js_dop_3_v4.csv",header=FALSE))
mu.est.js.dop.3.v5=data.matrix(read.csv("est_js_dop_3_v5.csv",header=FALSE))

mu.est.js.blip.1.v1=data.matrix(read.csv("est_js_blip_1_v1.csv",header=FALSE))
mu.est.js.blip.1.v2=data.matrix(read.csv("est_js_blip_1_v2.csv",header=FALSE))
mu.est.js.blip.1.v3=data.matrix(read.csv("est_js_blip_1_v3.csv",header=FALSE))
mu.est.js.blip.1.v4=data.matrix(read.csv("est_js_blip_1_v4.csv",header=FALSE))
mu.est.js.blip.1.v5=data.matrix(read.csv("est_js_blip_1_v5.csv",header=FALSE))
mu.est.js.blip.3.v1=data.matrix(read.csv("est_js_blip_3_v1.csv",header=FALSE))
mu.est.js.blip.3.v2=data.matrix(read.csv("est_js_blip_3_v2.csv",header=FALSE))
mu.est.js.blip.3.v3=data.matrix(read.csv("est_js_blip_3_v3.csv",header=FALSE))
mu.est.js.blip.3.v4=data.matrix(read.csv("est_js_blip_3_v4.csv",header=FALSE))
mu.est.js.blip.3.v5=data.matrix(read.csv("est_js_blip_3_v5.csv",header=FALSE))

mu.est.js.cor.1.v1=data.matrix(read.csv("est_js_cor_1_v1.csv",header=FALSE))
mu.est.js.cor.1.v2=data.matrix(read.csv("est_js_cor_1_v2.csv",header=FALSE))
mu.est.js.cor.1.v3=data.matrix(read.csv("est_js_cor_1_v3.csv",header=FALSE))
mu.est.js.cor.1.v4=data.matrix(read.csv("est_js_cor_1_v4.csv",header=FALSE))
mu.est.js.cor.1.v5=data.matrix(read.csv("est_js_cor_1_v5.csv",header=FALSE))
mu.est.js.cor.3.v1=data.matrix(read.csv("est_js_cor_3_v1.csv",header=FALSE))
mu.est.js.cor.3.v2=data.matrix(read.csv("est_js_cor_3_v2.csv",header=FALSE))
mu.est.js.cor.3.v3=data.matrix(read.csv("est_js_cor_3_v3.csv",header=FALSE))
mu.est.js.cor.3.v4=data.matrix(read.csv("est_js_cor_3_v4.csv",header=FALSE))
mu.est.js.cor.3.v5=data.matrix(read.csv("est_js_cor_3_v5.csv",header=FALSE))



mu.est.bams.sp.1.v1=data.matrix(read.csv("est_bams_sp_1_v1.csv",header=FALSE))
mu.est.bams.sp.1.v2=data.matrix(read.csv("est_bams_sp_1_v2.csv",header=FALSE))
mu.est.bams.sp.1.v3=data.matrix(read.csv("est_bams_sp_1_v3.csv",header=FALSE))
mu.est.bams.sp.1.v4=data.matrix(read.csv("est_bams_sp_1_v4.csv",header=FALSE))
mu.est.bams.sp.1.v5=data.matrix(read.csv("est_bams_sp_1_v5.csv",header=FALSE))
mu.est.bams.sp.3.v1=data.matrix(read.csv("est_bams_sp_3_v1.csv",header=FALSE))
mu.est.bams.sp.3.v2=data.matrix(read.csv("est_bams_sp_3_v2.csv",header=FALSE))
mu.est.bams.sp.3.v3=data.matrix(read.csv("est_bams_sp_3_v3.csv",header=FALSE))
mu.est.bams.sp.3.v4=data.matrix(read.csv("est_bams_sp_3_v4.csv",header=FALSE))
mu.est.bams.sp.3.v5=data.matrix(read.csv("est_bams_sp_3_v5.csv",header=FALSE))

mu.est.bams.bump.1.v1=data.matrix(read.csv("est_bams_bump_1_v1.csv",header=FALSE))
mu.est.bams.bump.1.v2=data.matrix(read.csv("est_bams_bump_1_v2.csv",header=FALSE))
mu.est.bams.bump.1.v3=data.matrix(read.csv("est_bams_bump_1_v3.csv",header=FALSE))
mu.est.bams.bump.1.v4=data.matrix(read.csv("est_bams_bump_1_v4.csv",header=FALSE))
mu.est.bams.bump.1.v5=data.matrix(read.csv("est_bams_bump_1_v5.csv",header=FALSE))
mu.est.bams.bump.3.v1=data.matrix(read.csv("est_bams_bump_3_v1.csv",header=FALSE))
mu.est.bams.bump.3.v2=data.matrix(read.csv("est_bams_bump_3_v2.csv",header=FALSE))
mu.est.bams.bump.3.v3=data.matrix(read.csv("est_bams_bump_3_v3.csv",header=FALSE))
mu.est.bams.bump.3.v4=data.matrix(read.csv("est_bams_bump_3_v4.csv",header=FALSE))
mu.est.bams.bump.3.v5=data.matrix(read.csv("est_bams_bump_3_v5.csv",header=FALSE))

mu.est.bams.blk.1.v1=data.matrix(read.csv("est_bams_blk_1_v1.csv",header=FALSE))
mu.est.bams.blk.1.v2=data.matrix(read.csv("est_bams_blk_1_v2.csv",header=FALSE))
mu.est.bams.blk.1.v3=data.matrix(read.csv("est_bams_blk_1_v3.csv",header=FALSE))
mu.est.bams.blk.1.v4=data.matrix(read.csv("est_bams_blk_1_v4.csv",header=FALSE))
mu.est.bams.blk.1.v5=data.matrix(read.csv("est_bams_blk_1_v5.csv",header=FALSE))
mu.est.bams.blk.3.v1=data.matrix(read.csv("est_bams_blk_3_v1.csv",header=FALSE))
mu.est.bams.blk.3.v2=data.matrix(read.csv("est_bams_blk_3_v2.csv",header=FALSE))
mu.est.bams.blk.3.v3=data.matrix(read.csv("est_bams_blk_3_v3.csv",header=FALSE))
mu.est.bams.blk.3.v4=data.matrix(read.csv("est_bams_blk_3_v4.csv",header=FALSE))
mu.est.bams.blk.3.v5=data.matrix(read.csv("est_bams_blk_3_v5.csv",header=FALSE))

mu.est.bams.ang.1.v1=data.matrix(read.csv("est_bams_ang_1_v1.csv",header=FALSE))
mu.est.bams.ang.1.v2=data.matrix(read.csv("est_bams_ang_1_v2.csv",header=FALSE))
mu.est.bams.ang.1.v3=data.matrix(read.csv("est_bams_ang_1_v3.csv",header=FALSE))
mu.est.bams.ang.1.v4=data.matrix(read.csv("est_bams_ang_1_v4.csv",header=FALSE))
mu.est.bams.ang.1.v5=data.matrix(read.csv("est_bams_ang_1_v5.csv",header=FALSE))
mu.est.bams.ang.3.v1=data.matrix(read.csv("est_bams_ang_3_v1.csv",header=FALSE))
mu.est.bams.ang.3.v2=data.matrix(read.csv("est_bams_ang_3_v2.csv",header=FALSE))
mu.est.bams.ang.3.v3=data.matrix(read.csv("est_bams_ang_3_v3.csv",header=FALSE))
mu.est.bams.ang.3.v4=data.matrix(read.csv("est_bams_ang_3_v4.csv",header=FALSE))
mu.est.bams.ang.3.v5=data.matrix(read.csv("est_bams_ang_3_v5.csv",header=FALSE))

mu.est.bams.dop.1.v1=data.matrix(read.csv("est_bams_dop_1_v1.csv",header=FALSE))
mu.est.bams.dop.1.v2=data.matrix(read.csv("est_bams_dop_1_v2.csv",header=FALSE))
mu.est.bams.dop.1.v3=data.matrix(read.csv("est_bams_dop_1_v3.csv",header=FALSE))
mu.est.bams.dop.1.v4=data.matrix(read.csv("est_bams_dop_1_v4.csv",header=FALSE))
mu.est.bams.dop.1.v5=data.matrix(read.csv("est_bams_dop_1_v5.csv",header=FALSE))
mu.est.bams.dop.3.v1=data.matrix(read.csv("est_bams_dop_3_v1.csv",header=FALSE))
mu.est.bams.dop.3.v2=data.matrix(read.csv("est_bams_dop_3_v2.csv",header=FALSE))
mu.est.bams.dop.3.v3=data.matrix(read.csv("est_bams_dop_3_v3.csv",header=FALSE))
mu.est.bams.dop.3.v4=data.matrix(read.csv("est_bams_dop_3_v4.csv",header=FALSE))
mu.est.bams.dop.3.v5=data.matrix(read.csv("est_bams_dop_3_v5.csv",header=FALSE))

mu.est.bams.blip.1.v1=data.matrix(read.csv("est_bams_blip_1_v1.csv",header=FALSE))
mu.est.bams.blip.1.v2=data.matrix(read.csv("est_bams_blip_1_v2.csv",header=FALSE))
mu.est.bams.blip.1.v3=data.matrix(read.csv("est_bams_blip_1_v3.csv",header=FALSE))
mu.est.bams.blip.1.v4=data.matrix(read.csv("est_bams_blip_1_v4.csv",header=FALSE))
mu.est.bams.blip.1.v5=data.matrix(read.csv("est_bams_blip_1_v5.csv",header=FALSE))
mu.est.bams.blip.3.v1=data.matrix(read.csv("est_bams_blip_3_v1.csv",header=FALSE))
mu.est.bams.blip.3.v2=data.matrix(read.csv("est_bams_blip_3_v2.csv",header=FALSE))
mu.est.bams.blip.3.v3=data.matrix(read.csv("est_bams_blip_3_v3.csv",header=FALSE))
mu.est.bams.blip.3.v4=data.matrix(read.csv("est_bams_blip_3_v4.csv",header=FALSE))
mu.est.bams.blip.3.v5=data.matrix(read.csv("est_bams_blip_3_v5.csv",header=FALSE))

mu.est.bams.cor.1.v1=data.matrix(read.csv("est_bams_cor_1_v1.csv",header=FALSE))
mu.est.bams.cor.1.v2=data.matrix(read.csv("est_bams_cor_1_v2.csv",header=FALSE))
mu.est.bams.cor.1.v3=data.matrix(read.csv("est_bams_cor_1_v3.csv",header=FALSE))
mu.est.bams.cor.1.v4=data.matrix(read.csv("est_bams_cor_1_v4.csv",header=FALSE))
mu.est.bams.cor.1.v5=data.matrix(read.csv("est_bams_cor_1_v5.csv",header=FALSE))
mu.est.bams.cor.3.v1=data.matrix(read.csv("est_bams_cor_3_v1.csv",header=FALSE))
mu.est.bams.cor.3.v2=data.matrix(read.csv("est_bams_cor_3_v2.csv",header=FALSE))
mu.est.bams.cor.3.v3=data.matrix(read.csv("est_bams_cor_3_v3.csv",header=FALSE))
mu.est.bams.cor.3.v4=data.matrix(read.csv("est_bams_cor_3_v4.csv",header=FALSE))
mu.est.bams.cor.3.v5=data.matrix(read.csv("est_bams_cor_3_v5.csv",header=FALSE))


mu.est.nblk.sp.1.v1=data.matrix(read.csv("est_nblk_sp_1_v1.csv",header=FALSE))
mu.est.nblk.sp.1.v2=data.matrix(read.csv("est_nblk_sp_1_v2.csv",header=FALSE))
mu.est.nblk.sp.1.v3=data.matrix(read.csv("est_nblk_sp_1_v3.csv",header=FALSE))
mu.est.nblk.sp.1.v4=data.matrix(read.csv("est_nblk_sp_1_v4.csv",header=FALSE))
mu.est.nblk.sp.1.v5=data.matrix(read.csv("est_nblk_sp_1_v5.csv",header=FALSE))
mu.est.nblk.sp.3.v1=data.matrix(read.csv("est_nblk_sp_3_v1.csv",header=FALSE))
mu.est.nblk.sp.3.v2=data.matrix(read.csv("est_nblk_sp_3_v2.csv",header=FALSE))
mu.est.nblk.sp.3.v3=data.matrix(read.csv("est_nblk_sp_3_v3.csv",header=FALSE))
mu.est.nblk.sp.3.v4=data.matrix(read.csv("est_nblk_sp_3_v4.csv",header=FALSE))
mu.est.nblk.sp.3.v5=data.matrix(read.csv("est_nblk_sp_3_v5.csv",header=FALSE))

mu.est.nblk.bump.1.v1=data.matrix(read.csv("est_nblk_bump_1_v1.csv",header=FALSE))
mu.est.nblk.bump.1.v2=data.matrix(read.csv("est_nblk_bump_1_v2.csv",header=FALSE))
mu.est.nblk.bump.1.v3=data.matrix(read.csv("est_nblk_bump_1_v3.csv",header=FALSE))
mu.est.nblk.bump.1.v4=data.matrix(read.csv("est_nblk_bump_1_v4.csv",header=FALSE))
mu.est.nblk.bump.1.v5=data.matrix(read.csv("est_nblk_bump_1_v5.csv",header=FALSE))
mu.est.nblk.bump.3.v1=data.matrix(read.csv("est_nblk_bump_3_v1.csv",header=FALSE))
mu.est.nblk.bump.3.v2=data.matrix(read.csv("est_nblk_bump_3_v2.csv",header=FALSE))
mu.est.nblk.bump.3.v3=data.matrix(read.csv("est_nblk_bump_3_v3.csv",header=FALSE))
mu.est.nblk.bump.3.v4=data.matrix(read.csv("est_nblk_bump_3_v4.csv",header=FALSE))
mu.est.nblk.bump.3.v5=data.matrix(read.csv("est_nblk_bump_3_v5.csv",header=FALSE))

mu.est.nblk.blk.1.v1=data.matrix(read.csv("est_nblk_blk_1_v1.csv",header=FALSE))
mu.est.nblk.blk.1.v2=data.matrix(read.csv("est_nblk_blk_1_v2.csv",header=FALSE))
mu.est.nblk.blk.1.v3=data.matrix(read.csv("est_nblk_blk_1_v3.csv",header=FALSE))
mu.est.nblk.blk.1.v4=data.matrix(read.csv("est_nblk_blk_1_v4.csv",header=FALSE))
mu.est.nblk.blk.1.v5=data.matrix(read.csv("est_nblk_blk_1_v5.csv",header=FALSE))
mu.est.nblk.blk.3.v1=data.matrix(read.csv("est_nblk_blk_3_v1.csv",header=FALSE))
mu.est.nblk.blk.3.v2=data.matrix(read.csv("est_nblk_blk_3_v2.csv",header=FALSE))
mu.est.nblk.blk.3.v3=data.matrix(read.csv("est_nblk_blk_3_v3.csv",header=FALSE))
mu.est.nblk.blk.3.v4=data.matrix(read.csv("est_nblk_blk_3_v4.csv",header=FALSE))
mu.est.nblk.blk.3.v5=data.matrix(read.csv("est_nblk_blk_3_v5.csv",header=FALSE))

mu.est.nblk.ang.1.v1=data.matrix(read.csv("est_nblk_ang_1_v1.csv",header=FALSE))
mu.est.nblk.ang.1.v2=data.matrix(read.csv("est_nblk_ang_1_v2.csv",header=FALSE))
mu.est.nblk.ang.1.v3=data.matrix(read.csv("est_nblk_ang_1_v3.csv",header=FALSE))
mu.est.nblk.ang.1.v4=data.matrix(read.csv("est_nblk_ang_1_v4.csv",header=FALSE))
mu.est.nblk.ang.1.v5=data.matrix(read.csv("est_nblk_ang_1_v5.csv",header=FALSE))
mu.est.nblk.ang.3.v1=data.matrix(read.csv("est_nblk_ang_3_v1.csv",header=FALSE))
mu.est.nblk.ang.3.v2=data.matrix(read.csv("est_nblk_ang_3_v2.csv",header=FALSE))
mu.est.nblk.ang.3.v3=data.matrix(read.csv("est_nblk_ang_3_v3.csv",header=FALSE))
mu.est.nblk.ang.3.v4=data.matrix(read.csv("est_nblk_ang_3_v4.csv",header=FALSE))
mu.est.nblk.ang.3.v5=data.matrix(read.csv("est_nblk_ang_3_v5.csv",header=FALSE))

mu.est.nblk.dop.1.v1=data.matrix(read.csv("est_nblk_dop_1_v1.csv",header=FALSE))
mu.est.nblk.dop.1.v2=data.matrix(read.csv("est_nblk_dop_1_v2.csv",header=FALSE))
mu.est.nblk.dop.1.v3=data.matrix(read.csv("est_nblk_dop_1_v3.csv",header=FALSE))
mu.est.nblk.dop.1.v4=data.matrix(read.csv("est_nblk_dop_1_v4.csv",header=FALSE))
mu.est.nblk.dop.1.v5=data.matrix(read.csv("est_nblk_dop_1_v5.csv",header=FALSE))
mu.est.nblk.dop.3.v1=data.matrix(read.csv("est_nblk_dop_3_v1.csv",header=FALSE))
mu.est.nblk.dop.3.v2=data.matrix(read.csv("est_nblk_dop_3_v2.csv",header=FALSE))
mu.est.nblk.dop.3.v3=data.matrix(read.csv("est_nblk_dop_3_v3.csv",header=FALSE))
mu.est.nblk.dop.3.v4=data.matrix(read.csv("est_nblk_dop_3_v4.csv",header=FALSE))
mu.est.nblk.dop.3.v5=data.matrix(read.csv("est_nblk_dop_3_v5.csv",header=FALSE))

mu.est.nblk.blip.1.v1=data.matrix(read.csv("est_nblk_blip_1_v1.csv",header=FALSE))
mu.est.nblk.blip.1.v2=data.matrix(read.csv("est_nblk_blip_1_v2.csv",header=FALSE))
mu.est.nblk.blip.1.v3=data.matrix(read.csv("est_nblk_blip_1_v3.csv",header=FALSE))
mu.est.nblk.blip.1.v4=data.matrix(read.csv("est_nblk_blip_1_v4.csv",header=FALSE))
mu.est.nblk.blip.1.v5=data.matrix(read.csv("est_nblk_blip_1_v5.csv",header=FALSE))
mu.est.nblk.blip.3.v1=data.matrix(read.csv("est_nblk_blip_3_v1.csv",header=FALSE))
mu.est.nblk.blip.3.v2=data.matrix(read.csv("est_nblk_blip_3_v2.csv",header=FALSE))
mu.est.nblk.blip.3.v3=data.matrix(read.csv("est_nblk_blip_3_v3.csv",header=FALSE))
mu.est.nblk.blip.3.v4=data.matrix(read.csv("est_nblk_blip_3_v4.csv",header=FALSE))
mu.est.nblk.blip.3.v5=data.matrix(read.csv("est_nblk_blip_3_v5.csv",header=FALSE))

mu.est.nblk.cor.1.v1=data.matrix(read.csv("est_nblk_cor_1_v1.csv",header=FALSE))
mu.est.nblk.cor.1.v2=data.matrix(read.csv("est_nblk_cor_1_v2.csv",header=FALSE))
mu.est.nblk.cor.1.v3=data.matrix(read.csv("est_nblk_cor_1_v3.csv",header=FALSE))
mu.est.nblk.cor.1.v4=data.matrix(read.csv("est_nblk_cor_1_v4.csv",header=FALSE))
mu.est.nblk.cor.1.v5=data.matrix(read.csv("est_nblk_cor_1_v5.csv",header=FALSE))
mu.est.nblk.cor.3.v1=data.matrix(read.csv("est_nblk_cor_3_v1.csv",header=FALSE))
mu.est.nblk.cor.3.v2=data.matrix(read.csv("est_nblk_cor_3_v2.csv",header=FALSE))
mu.est.nblk.cor.3.v3=data.matrix(read.csv("est_nblk_cor_3_v3.csv",header=FALSE))
mu.est.nblk.cor.3.v4=data.matrix(read.csv("est_nblk_cor_3_v4.csv",header=FALSE))
mu.est.nblk.cor.3.v5=data.matrix(read.csv("est_nblk_cor_3_v5.csv",header=FALSE))



mu.est.sure.sp.1.v1=data.matrix(read.csv("est_sure_sp_1_v1.csv",header=FALSE))
mu.est.sure.sp.1.v2=data.matrix(read.csv("est_sure_sp_1_v2.csv",header=FALSE))
mu.est.sure.sp.1.v3=data.matrix(read.csv("est_sure_sp_1_v3.csv",header=FALSE))
mu.est.sure.sp.1.v4=data.matrix(read.csv("est_sure_sp_1_v4.csv",header=FALSE))
mu.est.sure.sp.1.v5=data.matrix(read.csv("est_sure_sp_1_v5.csv",header=FALSE))
mu.est.sure.sp.3.v1=data.matrix(read.csv("est_sure_sp_3_v1.csv",header=FALSE))
mu.est.sure.sp.3.v2=data.matrix(read.csv("est_sure_sp_3_v2.csv",header=FALSE))
mu.est.sure.sp.3.v3=data.matrix(read.csv("est_sure_sp_3_v3.csv",header=FALSE))
mu.est.sure.sp.3.v4=data.matrix(read.csv("est_sure_sp_3_v4.csv",header=FALSE))
mu.est.sure.sp.3.v5=data.matrix(read.csv("est_sure_sp_3_v5.csv",header=FALSE))

mu.est.sure.bump.1.v1=data.matrix(read.csv("est_sure_bump_1_v1.csv",header=FALSE))
mu.est.sure.bump.1.v2=data.matrix(read.csv("est_sure_bump_1_v2.csv",header=FALSE))
mu.est.sure.bump.1.v3=data.matrix(read.csv("est_sure_bump_1_v3.csv",header=FALSE))
mu.est.sure.bump.1.v4=data.matrix(read.csv("est_sure_bump_1_v4.csv",header=FALSE))
mu.est.sure.bump.1.v5=data.matrix(read.csv("est_sure_bump_1_v5.csv",header=FALSE))
mu.est.sure.bump.3.v1=data.matrix(read.csv("est_sure_bump_3_v1.csv",header=FALSE))
mu.est.sure.bump.3.v2=data.matrix(read.csv("est_sure_bump_3_v2.csv",header=FALSE))
mu.est.sure.bump.3.v3=data.matrix(read.csv("est_sure_bump_3_v3.csv",header=FALSE))
mu.est.sure.bump.3.v4=data.matrix(read.csv("est_sure_bump_3_v4.csv",header=FALSE))
mu.est.sure.bump.3.v5=data.matrix(read.csv("est_sure_bump_3_v5.csv",header=FALSE))

mu.est.sure.blk.1.v1=data.matrix(read.csv("est_sure_blk_1_v1.csv",header=FALSE))
mu.est.sure.blk.1.v2=data.matrix(read.csv("est_sure_blk_1_v2.csv",header=FALSE))
mu.est.sure.blk.1.v3=data.matrix(read.csv("est_sure_blk_1_v3.csv",header=FALSE))
mu.est.sure.blk.1.v4=data.matrix(read.csv("est_sure_blk_1_v4.csv",header=FALSE))
mu.est.sure.blk.1.v5=data.matrix(read.csv("est_sure_blk_1_v5.csv",header=FALSE))
mu.est.sure.blk.3.v1=data.matrix(read.csv("est_sure_blk_3_v1.csv",header=FALSE))
mu.est.sure.blk.3.v2=data.matrix(read.csv("est_sure_blk_3_v2.csv",header=FALSE))
mu.est.sure.blk.3.v3=data.matrix(read.csv("est_sure_blk_3_v3.csv",header=FALSE))
mu.est.sure.blk.3.v4=data.matrix(read.csv("est_sure_blk_3_v4.csv",header=FALSE))
mu.est.sure.blk.3.v5=data.matrix(read.csv("est_sure_blk_3_v5.csv",header=FALSE))

mu.est.sure.ang.1.v1=data.matrix(read.csv("est_sure_ang_1_v1.csv",header=FALSE))
mu.est.sure.ang.1.v2=data.matrix(read.csv("est_sure_ang_1_v2.csv",header=FALSE))
mu.est.sure.ang.1.v3=data.matrix(read.csv("est_sure_ang_1_v3.csv",header=FALSE))
mu.est.sure.ang.1.v4=data.matrix(read.csv("est_sure_ang_1_v4.csv",header=FALSE))
mu.est.sure.ang.1.v5=data.matrix(read.csv("est_sure_ang_1_v5.csv",header=FALSE))
mu.est.sure.ang.3.v1=data.matrix(read.csv("est_sure_ang_3_v1.csv",header=FALSE))
mu.est.sure.ang.3.v2=data.matrix(read.csv("est_sure_ang_3_v2.csv",header=FALSE))
mu.est.sure.ang.3.v3=data.matrix(read.csv("est_sure_ang_3_v3.csv",header=FALSE))
mu.est.sure.ang.3.v4=data.matrix(read.csv("est_sure_ang_3_v4.csv",header=FALSE))
mu.est.sure.ang.3.v5=data.matrix(read.csv("est_sure_ang_3_v5.csv",header=FALSE))

mu.est.sure.dop.1.v1=data.matrix(read.csv("est_sure_dop_1_v1.csv",header=FALSE))
mu.est.sure.dop.1.v2=data.matrix(read.csv("est_sure_dop_1_v2.csv",header=FALSE))
mu.est.sure.dop.1.v3=data.matrix(read.csv("est_sure_dop_1_v3.csv",header=FALSE))
mu.est.sure.dop.1.v4=data.matrix(read.csv("est_sure_dop_1_v4.csv",header=FALSE))
mu.est.sure.dop.1.v5=data.matrix(read.csv("est_sure_dop_1_v5.csv",header=FALSE))
mu.est.sure.dop.3.v1=data.matrix(read.csv("est_sure_dop_3_v1.csv",header=FALSE))
mu.est.sure.dop.3.v2=data.matrix(read.csv("est_sure_dop_3_v2.csv",header=FALSE))
mu.est.sure.dop.3.v3=data.matrix(read.csv("est_sure_dop_3_v3.csv",header=FALSE))
mu.est.sure.dop.3.v4=data.matrix(read.csv("est_sure_dop_3_v4.csv",header=FALSE))
mu.est.sure.dop.3.v5=data.matrix(read.csv("est_sure_dop_3_v5.csv",header=FALSE))

mu.est.sure.blip.1.v1=data.matrix(read.csv("est_sure_blip_1_v1.csv",header=FALSE))
mu.est.sure.blip.1.v2=data.matrix(read.csv("est_sure_blip_1_v2.csv",header=FALSE))
mu.est.sure.blip.1.v3=data.matrix(read.csv("est_sure_blip_1_v3.csv",header=FALSE))
mu.est.sure.blip.1.v4=data.matrix(read.csv("est_sure_blip_1_v4.csv",header=FALSE))
mu.est.sure.blip.1.v5=data.matrix(read.csv("est_sure_blip_1_v5.csv",header=FALSE))
mu.est.sure.blip.3.v1=data.matrix(read.csv("est_sure_blip_3_v1.csv",header=FALSE))
mu.est.sure.blip.3.v2=data.matrix(read.csv("est_sure_blip_3_v2.csv",header=FALSE))
mu.est.sure.blip.3.v3=data.matrix(read.csv("est_sure_blip_3_v3.csv",header=FALSE))
mu.est.sure.blip.3.v4=data.matrix(read.csv("est_sure_blip_3_v4.csv",header=FALSE))
mu.est.sure.blip.3.v5=data.matrix(read.csv("est_sure_blip_3_v5.csv",header=FALSE))

mu.est.sure.cor.1.v1=data.matrix(read.csv("est_sure_cor_1_v1.csv",header=FALSE))
mu.est.sure.cor.1.v2=data.matrix(read.csv("est_sure_cor_1_v2.csv",header=FALSE))
mu.est.sure.cor.1.v3=data.matrix(read.csv("est_sure_cor_1_v3.csv",header=FALSE))
mu.est.sure.cor.1.v4=data.matrix(read.csv("est_sure_cor_1_v4.csv",header=FALSE))
mu.est.sure.cor.1.v5=data.matrix(read.csv("est_sure_cor_1_v5.csv",header=FALSE))
mu.est.sure.cor.3.v1=data.matrix(read.csv("est_sure_cor_3_v1.csv",header=FALSE))
mu.est.sure.cor.3.v2=data.matrix(read.csv("est_sure_cor_3_v2.csv",header=FALSE))
mu.est.sure.cor.3.v3=data.matrix(read.csv("est_sure_cor_3_v3.csv",header=FALSE))
mu.est.sure.cor.3.v4=data.matrix(read.csv("est_sure_cor_3_v4.csv",header=FALSE))
mu.est.sure.cor.3.v5=data.matrix(read.csv("est_sure_cor_3_v5.csv",header=FALSE))





mu.est.postmean.sp.1.v1=data.matrix(read.csv("est_postmean_sp_1_v1.csv",header=FALSE))
mu.est.postmean.sp.1.v2=data.matrix(read.csv("est_postmean_sp_1_v2.csv",header=FALSE))
mu.est.postmean.sp.1.v3=data.matrix(read.csv("est_postmean_sp_1_v3.csv",header=FALSE))
mu.est.postmean.sp.1.v4=data.matrix(read.csv("est_postmean_sp_1_v4.csv",header=FALSE))
mu.est.postmean.sp.1.v5=data.matrix(read.csv("est_postmean_sp_1_v5.csv",header=FALSE))
mu.est.postmean.sp.3.v1=data.matrix(read.csv("est_postmean_sp_3_v1.csv",header=FALSE))
mu.est.postmean.sp.3.v2=data.matrix(read.csv("est_postmean_sp_3_v2.csv",header=FALSE))
mu.est.postmean.sp.3.v3=data.matrix(read.csv("est_postmean_sp_3_v3.csv",header=FALSE))
mu.est.postmean.sp.3.v4=data.matrix(read.csv("est_postmean_sp_3_v4.csv",header=FALSE))
mu.est.postmean.sp.3.v5=data.matrix(read.csv("est_postmean_sp_3_v5.csv",header=FALSE))

mu.est.postmean.bump.1.v1=data.matrix(read.csv("est_postmean_bump_1_v1.csv",header=FALSE))
mu.est.postmean.bump.1.v2=data.matrix(read.csv("est_postmean_bump_1_v2.csv",header=FALSE))
mu.est.postmean.bump.1.v3=data.matrix(read.csv("est_postmean_bump_1_v3.csv",header=FALSE))
mu.est.postmean.bump.1.v4=data.matrix(read.csv("est_postmean_bump_1_v4.csv",header=FALSE))
mu.est.postmean.bump.1.v5=data.matrix(read.csv("est_postmean_bump_1_v5.csv",header=FALSE))
mu.est.postmean.bump.3.v1=data.matrix(read.csv("est_postmean_bump_3_v1.csv",header=FALSE))
mu.est.postmean.bump.3.v2=data.matrix(read.csv("est_postmean_bump_3_v2.csv",header=FALSE))
mu.est.postmean.bump.3.v3=data.matrix(read.csv("est_postmean_bump_3_v3.csv",header=FALSE))
mu.est.postmean.bump.3.v4=data.matrix(read.csv("est_postmean_bump_3_v4.csv",header=FALSE))
mu.est.postmean.bump.3.v5=data.matrix(read.csv("est_postmean_bump_3_v5.csv",header=FALSE))

mu.est.postmean.blk.1.v1=data.matrix(read.csv("est_postmean_blk_1_v1.csv",header=FALSE))
mu.est.postmean.blk.1.v2=data.matrix(read.csv("est_postmean_blk_1_v2.csv",header=FALSE))
mu.est.postmean.blk.1.v3=data.matrix(read.csv("est_postmean_blk_1_v3.csv",header=FALSE))
mu.est.postmean.blk.1.v4=data.matrix(read.csv("est_postmean_blk_1_v4.csv",header=FALSE))
mu.est.postmean.blk.1.v5=data.matrix(read.csv("est_postmean_blk_1_v5.csv",header=FALSE))
mu.est.postmean.blk.3.v1=data.matrix(read.csv("est_postmean_blk_3_v1.csv",header=FALSE))
mu.est.postmean.blk.3.v2=data.matrix(read.csv("est_postmean_blk_3_v2.csv",header=FALSE))
mu.est.postmean.blk.3.v3=data.matrix(read.csv("est_postmean_blk_3_v3.csv",header=FALSE))
mu.est.postmean.blk.3.v4=data.matrix(read.csv("est_postmean_blk_3_v4.csv",header=FALSE))
mu.est.postmean.blk.3.v5=data.matrix(read.csv("est_postmean_blk_3_v5.csv",header=FALSE))

mu.est.postmean.ang.1.v1=data.matrix(read.csv("est_postmean_ang_1_v1.csv",header=FALSE))
mu.est.postmean.ang.1.v2=data.matrix(read.csv("est_postmean_ang_1_v2.csv",header=FALSE))
mu.est.postmean.ang.1.v3=data.matrix(read.csv("est_postmean_ang_1_v3.csv",header=FALSE))
mu.est.postmean.ang.1.v4=data.matrix(read.csv("est_postmean_ang_1_v4.csv",header=FALSE))
mu.est.postmean.ang.1.v5=data.matrix(read.csv("est_postmean_ang_1_v5.csv",header=FALSE))
mu.est.postmean.ang.3.v1=data.matrix(read.csv("est_postmean_ang_3_v1.csv",header=FALSE))
mu.est.postmean.ang.3.v2=data.matrix(read.csv("est_postmean_ang_3_v2.csv",header=FALSE))
mu.est.postmean.ang.3.v3=data.matrix(read.csv("est_postmean_ang_3_v3.csv",header=FALSE))
mu.est.postmean.ang.3.v4=data.matrix(read.csv("est_postmean_ang_3_v4.csv",header=FALSE))
mu.est.postmean.ang.3.v5=data.matrix(read.csv("est_postmean_ang_3_v5.csv",header=FALSE))

mu.est.postmean.dop.1.v1=data.matrix(read.csv("est_postmean_dop_1_v1.csv",header=FALSE))
mu.est.postmean.dop.1.v2=data.matrix(read.csv("est_postmean_dop_1_v2.csv",header=FALSE))
mu.est.postmean.dop.1.v3=data.matrix(read.csv("est_postmean_dop_1_v3.csv",header=FALSE))
mu.est.postmean.dop.1.v4=data.matrix(read.csv("est_postmean_dop_1_v4.csv",header=FALSE))
mu.est.postmean.dop.1.v5=data.matrix(read.csv("est_postmean_dop_1_v5.csv",header=FALSE))
mu.est.postmean.dop.3.v1=data.matrix(read.csv("est_postmean_dop_3_v1.csv",header=FALSE))
mu.est.postmean.dop.3.v2=data.matrix(read.csv("est_postmean_dop_3_v2.csv",header=FALSE))
mu.est.postmean.dop.3.v3=data.matrix(read.csv("est_postmean_dop_3_v3.csv",header=FALSE))
mu.est.postmean.dop.3.v4=data.matrix(read.csv("est_postmean_dop_3_v4.csv",header=FALSE))
mu.est.postmean.dop.3.v5=data.matrix(read.csv("est_postmean_dop_3_v5.csv",header=FALSE))

mu.est.postmean.blip.1.v1=data.matrix(read.csv("est_postmean_blip_1_v1.csv",header=FALSE))
mu.est.postmean.blip.1.v2=data.matrix(read.csv("est_postmean_blip_1_v2.csv",header=FALSE))
mu.est.postmean.blip.1.v3=data.matrix(read.csv("est_postmean_blip_1_v3.csv",header=FALSE))
mu.est.postmean.blip.1.v4=data.matrix(read.csv("est_postmean_blip_1_v4.csv",header=FALSE))
mu.est.postmean.blip.1.v5=data.matrix(read.csv("est_postmean_blip_1_v5.csv",header=FALSE))
mu.est.postmean.blip.3.v1=data.matrix(read.csv("est_postmean_blip_3_v1.csv",header=FALSE))
mu.est.postmean.blip.3.v2=data.matrix(read.csv("est_postmean_blip_3_v2.csv",header=FALSE))
mu.est.postmean.blip.3.v3=data.matrix(read.csv("est_postmean_blip_3_v3.csv",header=FALSE))
mu.est.postmean.blip.3.v4=data.matrix(read.csv("est_postmean_blip_3_v4.csv",header=FALSE))
mu.est.postmean.blip.3.v5=data.matrix(read.csv("est_postmean_blip_3_v5.csv",header=FALSE))

mu.est.postmean.cor.1.v1=data.matrix(read.csv("est_postmean_cor_1_v1.csv",header=FALSE))
mu.est.postmean.cor.1.v2=data.matrix(read.csv("est_postmean_cor_1_v2.csv",header=FALSE))
mu.est.postmean.cor.1.v3=data.matrix(read.csv("est_postmean_cor_1_v3.csv",header=FALSE))
mu.est.postmean.cor.1.v4=data.matrix(read.csv("est_postmean_cor_1_v4.csv",header=FALSE))
mu.est.postmean.cor.1.v5=data.matrix(read.csv("est_postmean_cor_1_v5.csv",header=FALSE))
mu.est.postmean.cor.3.v1=data.matrix(read.csv("est_postmean_cor_3_v1.csv",header=FALSE))
mu.est.postmean.cor.3.v2=data.matrix(read.csv("est_postmean_cor_3_v2.csv",header=FALSE))
mu.est.postmean.cor.3.v3=data.matrix(read.csv("est_postmean_cor_3_v3.csv",header=FALSE))
mu.est.postmean.cor.3.v4=data.matrix(read.csv("est_postmean_cor_3_v4.csv",header=FALSE))
mu.est.postmean.cor.3.v5=data.matrix(read.csv("est_postmean_cor_3_v5.csv",header=FALSE))




mise.ti.sp.1.v1=mise(mu.est.ti.sp.1.v1,mu.sp)
mise.ti.sp.1.v2=mise(mu.est.ti.sp.1.v2,mu.sp)
mise.ti.sp.1.v3=mise(mu.est.ti.sp.1.v3,mu.sp)
mise.ti.sp.1.v4=mise(mu.est.ti.sp.1.v4,mu.sp)
mise.ti.sp.1.v5=mise(mu.est.ti.sp.1.v5,mu.sp)
mise.ti.sp.3.v1=mise(mu.est.ti.sp.3.v1,mu.sp)
mise.ti.sp.3.v2=mise(mu.est.ti.sp.3.v2,mu.sp)
mise.ti.sp.3.v3=mise(mu.est.ti.sp.3.v3,mu.sp)
mise.ti.sp.3.v4=mise(mu.est.ti.sp.3.v4,mu.sp)
mise.ti.sp.3.v5=mise(mu.est.ti.sp.3.v5,mu.sp)

mise.ti.bump.1.v1=mise(mu.est.ti.bump.1.v1,mu.bump)
mise.ti.bump.1.v2=mise(mu.est.ti.bump.1.v2,mu.bump)
mise.ti.bump.1.v3=mise(mu.est.ti.bump.1.v3,mu.bump)
mise.ti.bump.1.v4=mise(mu.est.ti.bump.1.v4,mu.bump)
mise.ti.bump.1.v5=mise(mu.est.ti.bump.1.v5,mu.bump)
mise.ti.bump.3.v1=mise(mu.est.ti.bump.3.v1,mu.bump)
mise.ti.bump.3.v2=mise(mu.est.ti.bump.3.v2,mu.bump)
mise.ti.bump.3.v3=mise(mu.est.ti.bump.3.v3,mu.bump)
mise.ti.bump.3.v4=mise(mu.est.ti.bump.3.v4,mu.bump)
mise.ti.bump.3.v5=mise(mu.est.ti.bump.3.v5,mu.bump)

mise.ti.blk.1.v1=mise(mu.est.ti.blk.1.v1,mu.blk)
mise.ti.blk.1.v2=mise(mu.est.ti.blk.1.v2,mu.blk)
mise.ti.blk.1.v3=mise(mu.est.ti.blk.1.v3,mu.blk)
mise.ti.blk.1.v4=mise(mu.est.ti.blk.1.v4,mu.blk)
mise.ti.blk.1.v5=mise(mu.est.ti.blk.1.v5,mu.blk)
mise.ti.blk.3.v1=mise(mu.est.ti.blk.3.v1,mu.blk)
mise.ti.blk.3.v2=mise(mu.est.ti.blk.3.v2,mu.blk)
mise.ti.blk.3.v3=mise(mu.est.ti.blk.3.v3,mu.blk)
mise.ti.blk.3.v4=mise(mu.est.ti.blk.3.v4,mu.blk)
mise.ti.blk.3.v5=mise(mu.est.ti.blk.3.v5,mu.blk)

mise.ti.ang.1.v1=mise(mu.est.ti.ang.1.v1,mu.ang)
mise.ti.ang.1.v2=mise(mu.est.ti.ang.1.v2,mu.ang)
mise.ti.ang.1.v3=mise(mu.est.ti.ang.1.v3,mu.ang)
mise.ti.ang.1.v4=mise(mu.est.ti.ang.1.v4,mu.ang)
mise.ti.ang.1.v5=mise(mu.est.ti.ang.1.v5,mu.ang)
mise.ti.ang.3.v1=mise(mu.est.ti.ang.3.v1,mu.ang)
mise.ti.ang.3.v2=mise(mu.est.ti.ang.3.v2,mu.ang)
mise.ti.ang.3.v3=mise(mu.est.ti.ang.3.v3,mu.ang)
mise.ti.ang.3.v4=mise(mu.est.ti.ang.3.v4,mu.ang)
mise.ti.ang.3.v5=mise(mu.est.ti.ang.3.v5,mu.ang)

mise.ti.dop.1.v1=mise(mu.est.ti.dop.1.v1,mu.dop)
mise.ti.dop.1.v2=mise(mu.est.ti.dop.1.v2,mu.dop)
mise.ti.dop.1.v3=mise(mu.est.ti.dop.1.v3,mu.dop)
mise.ti.dop.1.v4=mise(mu.est.ti.dop.1.v4,mu.dop)
mise.ti.dop.1.v5=mise(mu.est.ti.dop.1.v5,mu.dop)
mise.ti.dop.3.v1=mise(mu.est.ti.dop.3.v1,mu.dop)
mise.ti.dop.3.v2=mise(mu.est.ti.dop.3.v2,mu.dop)
mise.ti.dop.3.v3=mise(mu.est.ti.dop.3.v3,mu.dop)
mise.ti.dop.3.v4=mise(mu.est.ti.dop.3.v4,mu.dop)
mise.ti.dop.3.v5=mise(mu.est.ti.dop.3.v5,mu.dop)

mise.ti.blip.1.v1=mise(mu.est.ti.blip.1.v1,mu.blip)
mise.ti.blip.1.v2=mise(mu.est.ti.blip.1.v2,mu.blip)
mise.ti.blip.1.v3=mise(mu.est.ti.blip.1.v3,mu.blip)
mise.ti.blip.1.v4=mise(mu.est.ti.blip.1.v4,mu.blip)
mise.ti.blip.1.v5=mise(mu.est.ti.blip.1.v5,mu.blip)
mise.ti.blip.3.v1=mise(mu.est.ti.blip.3.v1,mu.blip)
mise.ti.blip.3.v2=mise(mu.est.ti.blip.3.v2,mu.blip)
mise.ti.blip.3.v3=mise(mu.est.ti.blip.3.v3,mu.blip)
mise.ti.blip.3.v4=mise(mu.est.ti.blip.3.v4,mu.blip)
mise.ti.blip.3.v5=mise(mu.est.ti.blip.3.v5,mu.blip)

mise.ti.cor.1.v1=mise(mu.est.ti.cor.1.v1,mu.cor)
mise.ti.cor.1.v2=mise(mu.est.ti.cor.1.v2,mu.cor)
mise.ti.cor.1.v3=mise(mu.est.ti.cor.1.v3,mu.cor)
mise.ti.cor.1.v4=mise(mu.est.ti.cor.1.v4,mu.cor)
mise.ti.cor.1.v5=mise(mu.est.ti.cor.1.v5,mu.cor)
mise.ti.cor.3.v1=mise(mu.est.ti.cor.3.v1,mu.cor)
mise.ti.cor.3.v2=mise(mu.est.ti.cor.3.v2,mu.cor)
mise.ti.cor.3.v3=mise(mu.est.ti.cor.3.v3,mu.cor)
mise.ti.cor.3.v4=mise(mu.est.ti.cor.3.v4,mu.cor)
mise.ti.cor.3.v5=mise(mu.est.ti.cor.3.v5,mu.cor)



mise.js.sp.1.v1=mise(mu.est.js.sp.1.v1,mu.sp)
mise.js.sp.1.v2=mise(mu.est.js.sp.1.v2,mu.sp)
mise.js.sp.1.v3=mise(mu.est.js.sp.1.v3,mu.sp)
mise.js.sp.1.v4=mise(mu.est.js.sp.1.v4,mu.sp)
mise.js.sp.1.v5=mise(mu.est.js.sp.1.v5,mu.sp)
mise.js.sp.3.v1=mise(mu.est.js.sp.3.v1,mu.sp)
mise.js.sp.3.v2=mise(mu.est.js.sp.3.v2,mu.sp)
mise.js.sp.3.v3=mise(mu.est.js.sp.3.v3,mu.sp)
mise.js.sp.3.v4=mise(mu.est.js.sp.3.v4,mu.sp)
mise.js.sp.3.v5=mise(mu.est.js.sp.3.v5,mu.sp)

mise.js.bump.1.v1=mise(mu.est.js.bump.1.v1,mu.bump)
mise.js.bump.1.v2=mise(mu.est.js.bump.1.v2,mu.bump)
mise.js.bump.1.v3=mise(mu.est.js.bump.1.v3,mu.bump)
mise.js.bump.1.v4=mise(mu.est.js.bump.1.v4,mu.bump)
mise.js.bump.1.v5=mise(mu.est.js.bump.1.v5,mu.bump)
mise.js.bump.3.v1=mise(mu.est.js.bump.3.v1,mu.bump)
mise.js.bump.3.v2=mise(mu.est.js.bump.3.v2,mu.bump)
mise.js.bump.3.v3=mise(mu.est.js.bump.3.v3,mu.bump)
mise.js.bump.3.v4=mise(mu.est.js.bump.3.v4,mu.bump)
mise.js.bump.3.v5=mise(mu.est.js.bump.3.v5,mu.bump)

mise.js.blk.1.v1=mise(mu.est.js.blk.1.v1,mu.blk)
mise.js.blk.1.v2=mise(mu.est.js.blk.1.v2,mu.blk)
mise.js.blk.1.v3=mise(mu.est.js.blk.1.v3,mu.blk)
mise.js.blk.1.v4=mise(mu.est.js.blk.1.v4,mu.blk)
mise.js.blk.1.v5=mise(mu.est.js.blk.1.v5,mu.blk)
mise.js.blk.3.v1=mise(mu.est.js.blk.3.v1,mu.blk)
mise.js.blk.3.v2=mise(mu.est.js.blk.3.v2,mu.blk)
mise.js.blk.3.v3=mise(mu.est.js.blk.3.v3,mu.blk)
mise.js.blk.3.v4=mise(mu.est.js.blk.3.v4,mu.blk)
mise.js.blk.3.v5=mise(mu.est.js.blk.3.v5,mu.blk)

mise.js.ang.1.v1=mise(mu.est.js.ang.1.v1,mu.ang)
mise.js.ang.1.v2=mise(mu.est.js.ang.1.v2,mu.ang)
mise.js.ang.1.v3=mise(mu.est.js.ang.1.v3,mu.ang)
mise.js.ang.1.v4=mise(mu.est.js.ang.1.v4,mu.ang)
mise.js.ang.1.v5=mise(mu.est.js.ang.1.v5,mu.ang)
mise.js.ang.3.v1=mise(mu.est.js.ang.3.v1,mu.ang)
mise.js.ang.3.v2=mise(mu.est.js.ang.3.v2,mu.ang)
mise.js.ang.3.v3=mise(mu.est.js.ang.3.v3,mu.ang)
mise.js.ang.3.v4=mise(mu.est.js.ang.3.v4,mu.ang)
mise.js.ang.3.v5=mise(mu.est.js.ang.3.v5,mu.ang)

mise.js.dop.1.v1=mise(mu.est.js.dop.1.v1,mu.dop)
mise.js.dop.1.v2=mise(mu.est.js.dop.1.v2,mu.dop)
mise.js.dop.1.v3=mise(mu.est.js.dop.1.v3,mu.dop)
mise.js.dop.1.v4=mise(mu.est.js.dop.1.v4,mu.dop)
mise.js.dop.1.v5=mise(mu.est.js.dop.1.v5,mu.dop)
mise.js.dop.3.v1=mise(mu.est.js.dop.3.v1,mu.dop)
mise.js.dop.3.v2=mise(mu.est.js.dop.3.v2,mu.dop)
mise.js.dop.3.v3=mise(mu.est.js.dop.3.v3,mu.dop)
mise.js.dop.3.v4=mise(mu.est.js.dop.3.v4,mu.dop)
mise.js.dop.3.v5=mise(mu.est.js.dop.3.v5,mu.dop)

mise.js.blip.1.v1=mise(mu.est.js.blip.1.v1,mu.blip)
mise.js.blip.1.v2=mise(mu.est.js.blip.1.v2,mu.blip)
mise.js.blip.1.v3=mise(mu.est.js.blip.1.v3,mu.blip)
mise.js.blip.1.v4=mise(mu.est.js.blip.1.v4,mu.blip)
mise.js.blip.1.v5=mise(mu.est.js.blip.1.v5,mu.blip)
mise.js.blip.3.v1=mise(mu.est.js.blip.3.v1,mu.blip)
mise.js.blip.3.v2=mise(mu.est.js.blip.3.v2,mu.blip)
mise.js.blip.3.v3=mise(mu.est.js.blip.3.v3,mu.blip)
mise.js.blip.3.v4=mise(mu.est.js.blip.3.v4,mu.blip)
mise.js.blip.3.v5=mise(mu.est.js.blip.3.v5,mu.blip)

mise.js.cor.1.v1=mise(mu.est.js.cor.1.v1,mu.cor)
mise.js.cor.1.v2=mise(mu.est.js.cor.1.v2,mu.cor)
mise.js.cor.1.v3=mise(mu.est.js.cor.1.v3,mu.cor)
mise.js.cor.1.v4=mise(mu.est.js.cor.1.v4,mu.cor)
mise.js.cor.1.v5=mise(mu.est.js.cor.1.v5,mu.cor)
mise.js.cor.3.v1=mise(mu.est.js.cor.3.v1,mu.cor)
mise.js.cor.3.v2=mise(mu.est.js.cor.3.v2,mu.cor)
mise.js.cor.3.v3=mise(mu.est.js.cor.3.v3,mu.cor)
mise.js.cor.3.v4=mise(mu.est.js.cor.3.v4,mu.cor)
mise.js.cor.3.v5=mise(mu.est.js.cor.3.v5,mu.cor)



mise.bams.sp.1.v1=mise(mu.est.bams.sp.1.v1,mu.sp)
mise.bams.sp.1.v2=mise(mu.est.bams.sp.1.v2,mu.sp)
mise.bams.sp.1.v3=mise(mu.est.bams.sp.1.v3,mu.sp)
mise.bams.sp.1.v4=mise(mu.est.bams.sp.1.v4,mu.sp)
mise.bams.sp.1.v5=mise(mu.est.bams.sp.1.v5,mu.sp)
mise.bams.sp.3.v1=mise(mu.est.bams.sp.3.v1,mu.sp)
mise.bams.sp.3.v2=mise(mu.est.bams.sp.3.v2,mu.sp)
mise.bams.sp.3.v3=mise(mu.est.bams.sp.3.v3,mu.sp)
mise.bams.sp.3.v4=mise(mu.est.bams.sp.3.v4,mu.sp)
mise.bams.sp.3.v5=mise(mu.est.bams.sp.3.v5,mu.sp)

mise.bams.bump.1.v1=mise(mu.est.bams.bump.1.v1,mu.bump)
mise.bams.bump.1.v2=mise(mu.est.bams.bump.1.v2,mu.bump)
mise.bams.bump.1.v3=mise(mu.est.bams.bump.1.v3,mu.bump)
mise.bams.bump.1.v4=mise(mu.est.bams.bump.1.v4,mu.bump)
mise.bams.bump.1.v5=mise(mu.est.bams.bump.1.v5,mu.bump)
mise.bams.bump.3.v1=mise(mu.est.bams.bump.3.v1,mu.bump)
mise.bams.bump.3.v2=mise(mu.est.bams.bump.3.v2,mu.bump)
mise.bams.bump.3.v3=mise(mu.est.bams.bump.3.v3,mu.bump)
mise.bams.bump.3.v4=mise(mu.est.bams.bump.3.v4,mu.bump)
mise.bams.bump.3.v5=mise(mu.est.bams.bump.3.v5,mu.bump)

mise.bams.blk.1.v1=mise(mu.est.bams.blk.1.v1,mu.blk)
mise.bams.blk.1.v2=mise(mu.est.bams.blk.1.v2,mu.blk)
mise.bams.blk.1.v3=mise(mu.est.bams.blk.1.v3,mu.blk)
mise.bams.blk.1.v4=mise(mu.est.bams.blk.1.v4,mu.blk)
mise.bams.blk.1.v5=mise(mu.est.bams.blk.1.v5,mu.blk)
mise.bams.blk.3.v1=mise(mu.est.bams.blk.3.v1,mu.blk)
mise.bams.blk.3.v2=mise(mu.est.bams.blk.3.v2,mu.blk)
mise.bams.blk.3.v3=mise(mu.est.bams.blk.3.v3,mu.blk)
mise.bams.blk.3.v4=mise(mu.est.bams.blk.3.v4,mu.blk)
mise.bams.blk.3.v5=mise(mu.est.bams.blk.3.v5,mu.blk)

mise.bams.ang.1.v1=mise(mu.est.bams.ang.1.v1,mu.ang)
mise.bams.ang.1.v2=mise(mu.est.bams.ang.1.v2,mu.ang)
mise.bams.ang.1.v3=mise(mu.est.bams.ang.1.v3,mu.ang)
mise.bams.ang.1.v4=mise(mu.est.bams.ang.1.v4,mu.ang)
mise.bams.ang.1.v5=mise(mu.est.bams.ang.1.v5,mu.ang)
mise.bams.ang.3.v1=mise(mu.est.bams.ang.3.v1,mu.ang)
mise.bams.ang.3.v2=mise(mu.est.bams.ang.3.v2,mu.ang)
mise.bams.ang.3.v3=mise(mu.est.bams.ang.3.v3,mu.ang)
mise.bams.ang.3.v4=mise(mu.est.bams.ang.3.v4,mu.ang)
mise.bams.ang.3.v5=mise(mu.est.bams.ang.3.v5,mu.ang)

mise.bams.dop.1.v1=mise(mu.est.bams.dop.1.v1,mu.dop)
mise.bams.dop.1.v2=mise(mu.est.bams.dop.1.v2,mu.dop)
mise.bams.dop.1.v3=mise(mu.est.bams.dop.1.v3,mu.dop)
mise.bams.dop.1.v4=mise(mu.est.bams.dop.1.v4,mu.dop)
mise.bams.dop.1.v5=mise(mu.est.bams.dop.1.v5,mu.dop)
mise.bams.dop.3.v1=mise(mu.est.bams.dop.3.v1,mu.dop)
mise.bams.dop.3.v2=mise(mu.est.bams.dop.3.v2,mu.dop)
mise.bams.dop.3.v3=mise(mu.est.bams.dop.3.v3,mu.dop)
mise.bams.dop.3.v4=mise(mu.est.bams.dop.3.v4,mu.dop)
mise.bams.dop.3.v5=mise(mu.est.bams.dop.3.v5,mu.dop)

mise.bams.blip.1.v1=mise(mu.est.bams.blip.1.v1,mu.blip)
mise.bams.blip.1.v2=mise(mu.est.bams.blip.1.v2,mu.blip)
mise.bams.blip.1.v3=mise(mu.est.bams.blip.1.v3,mu.blip)
mise.bams.blip.1.v4=mise(mu.est.bams.blip.1.v4,mu.blip)
mise.bams.blip.1.v5=mise(mu.est.bams.blip.1.v5,mu.blip)
mise.bams.blip.3.v1=mise(mu.est.bams.blip.3.v1,mu.blip)
mise.bams.blip.3.v2=mise(mu.est.bams.blip.3.v2,mu.blip)
mise.bams.blip.3.v3=mise(mu.est.bams.blip.3.v3,mu.blip)
mise.bams.blip.3.v4=mise(mu.est.bams.blip.3.v4,mu.blip)
mise.bams.blip.3.v5=mise(mu.est.bams.blip.3.v5,mu.blip)

mise.bams.cor.1.v1=mise(mu.est.bams.cor.1.v1,mu.cor)
mise.bams.cor.1.v2=mise(mu.est.bams.cor.1.v2,mu.cor)
mise.bams.cor.1.v3=mise(mu.est.bams.cor.1.v3,mu.cor)
mise.bams.cor.1.v4=mise(mu.est.bams.cor.1.v4,mu.cor)
mise.bams.cor.1.v5=mise(mu.est.bams.cor.1.v5,mu.cor)
mise.bams.cor.3.v1=mise(mu.est.bams.cor.3.v1,mu.cor)
mise.bams.cor.3.v2=mise(mu.est.bams.cor.3.v2,mu.cor)
mise.bams.cor.3.v3=mise(mu.est.bams.cor.3.v3,mu.cor)
mise.bams.cor.3.v4=mise(mu.est.bams.cor.3.v4,mu.cor)
mise.bams.cor.3.v5=mise(mu.est.bams.cor.3.v5,mu.cor)



mise.nblk.sp.1.v1=mise(mu.est.nblk.sp.1.v1,mu.sp)
mise.nblk.sp.1.v2=mise(mu.est.nblk.sp.1.v2,mu.sp)
mise.nblk.sp.1.v3=mise(mu.est.nblk.sp.1.v3,mu.sp)
mise.nblk.sp.1.v4=mise(mu.est.nblk.sp.1.v4,mu.sp)
mise.nblk.sp.1.v5=mise(mu.est.nblk.sp.1.v5,mu.sp)
mise.nblk.sp.3.v1=mise(mu.est.nblk.sp.3.v1,mu.sp)
mise.nblk.sp.3.v2=mise(mu.est.nblk.sp.3.v2,mu.sp)
mise.nblk.sp.3.v3=mise(mu.est.nblk.sp.3.v3,mu.sp)
mise.nblk.sp.3.v4=mise(mu.est.nblk.sp.3.v4,mu.sp)
mise.nblk.sp.3.v5=mise(mu.est.nblk.sp.3.v5,mu.sp)

mise.nblk.bump.1.v1=mise(mu.est.nblk.bump.1.v1,mu.bump)
mise.nblk.bump.1.v2=mise(mu.est.nblk.bump.1.v2,mu.bump)
mise.nblk.bump.1.v3=mise(mu.est.nblk.bump.1.v3,mu.bump)
mise.nblk.bump.1.v4=mise(mu.est.nblk.bump.1.v4,mu.bump)
mise.nblk.bump.1.v5=mise(mu.est.nblk.bump.1.v5,mu.bump)
mise.nblk.bump.3.v1=mise(mu.est.nblk.bump.3.v1,mu.bump)
mise.nblk.bump.3.v2=mise(mu.est.nblk.bump.3.v2,mu.bump)
mise.nblk.bump.3.v3=mise(mu.est.nblk.bump.3.v3,mu.bump)
mise.nblk.bump.3.v4=mise(mu.est.nblk.bump.3.v4,mu.bump)
mise.nblk.bump.3.v5=mise(mu.est.nblk.bump.3.v5,mu.bump)

mise.nblk.blk.1.v1=mise(mu.est.nblk.blk.1.v1,mu.blk)
mise.nblk.blk.1.v2=mise(mu.est.nblk.blk.1.v2,mu.blk)
mise.nblk.blk.1.v3=mise(mu.est.nblk.blk.1.v3,mu.blk)
mise.nblk.blk.1.v4=mise(mu.est.nblk.blk.1.v4,mu.blk)
mise.nblk.blk.1.v5=mise(mu.est.nblk.blk.1.v5,mu.blk)
mise.nblk.blk.3.v1=mise(mu.est.nblk.blk.3.v1,mu.blk)
mise.nblk.blk.3.v2=mise(mu.est.nblk.blk.3.v2,mu.blk)
mise.nblk.blk.3.v3=mise(mu.est.nblk.blk.3.v3,mu.blk)
mise.nblk.blk.3.v4=mise(mu.est.nblk.blk.3.v4,mu.blk)
mise.nblk.blk.3.v5=mise(mu.est.nblk.blk.3.v5,mu.blk)

mise.nblk.ang.1.v1=mise(mu.est.nblk.ang.1.v1,mu.ang)
mise.nblk.ang.1.v2=mise(mu.est.nblk.ang.1.v2,mu.ang)
mise.nblk.ang.1.v3=mise(mu.est.nblk.ang.1.v3,mu.ang)
mise.nblk.ang.1.v4=mise(mu.est.nblk.ang.1.v4,mu.ang)
mise.nblk.ang.1.v5=mise(mu.est.nblk.ang.1.v5,mu.ang)
mise.nblk.ang.3.v1=mise(mu.est.nblk.ang.3.v1,mu.ang)
mise.nblk.ang.3.v2=mise(mu.est.nblk.ang.3.v2,mu.ang)
mise.nblk.ang.3.v3=mise(mu.est.nblk.ang.3.v3,mu.ang)
mise.nblk.ang.3.v4=mise(mu.est.nblk.ang.3.v4,mu.ang)
mise.nblk.ang.3.v5=mise(mu.est.nblk.ang.3.v5,mu.ang)

mise.nblk.dop.1.v1=mise(mu.est.nblk.dop.1.v1,mu.dop)
mise.nblk.dop.1.v2=mise(mu.est.nblk.dop.1.v2,mu.dop)
mise.nblk.dop.1.v3=mise(mu.est.nblk.dop.1.v3,mu.dop)
mise.nblk.dop.1.v4=mise(mu.est.nblk.dop.1.v4,mu.dop)
mise.nblk.dop.1.v5=mise(mu.est.nblk.dop.1.v5,mu.dop)
mise.nblk.dop.3.v1=mise(mu.est.nblk.dop.3.v1,mu.dop)
mise.nblk.dop.3.v2=mise(mu.est.nblk.dop.3.v2,mu.dop)
mise.nblk.dop.3.v3=mise(mu.est.nblk.dop.3.v3,mu.dop)
mise.nblk.dop.3.v4=mise(mu.est.nblk.dop.3.v4,mu.dop)
mise.nblk.dop.3.v5=mise(mu.est.nblk.dop.3.v5,mu.dop)

mise.nblk.blip.1.v1=mise(mu.est.nblk.blip.1.v1,mu.blip)
mise.nblk.blip.1.v2=mise(mu.est.nblk.blip.1.v2,mu.blip)
mise.nblk.blip.1.v3=mise(mu.est.nblk.blip.1.v3,mu.blip)
mise.nblk.blip.1.v4=mise(mu.est.nblk.blip.1.v4,mu.blip)
mise.nblk.blip.1.v5=mise(mu.est.nblk.blip.1.v5,mu.blip)
mise.nblk.blip.3.v1=mise(mu.est.nblk.blip.3.v1,mu.blip)
mise.nblk.blip.3.v2=mise(mu.est.nblk.blip.3.v2,mu.blip)
mise.nblk.blip.3.v3=mise(mu.est.nblk.blip.3.v3,mu.blip)
mise.nblk.blip.3.v4=mise(mu.est.nblk.blip.3.v4,mu.blip)
mise.nblk.blip.3.v5=mise(mu.est.nblk.blip.3.v5,mu.blip)

mise.nblk.cor.1.v1=mise(mu.est.nblk.cor.1.v1,mu.cor)
mise.nblk.cor.1.v2=mise(mu.est.nblk.cor.1.v2,mu.cor)
mise.nblk.cor.1.v3=mise(mu.est.nblk.cor.1.v3,mu.cor)
mise.nblk.cor.1.v4=mise(mu.est.nblk.cor.1.v4,mu.cor)
mise.nblk.cor.1.v5=mise(mu.est.nblk.cor.1.v5,mu.cor)
mise.nblk.cor.3.v1=mise(mu.est.nblk.cor.3.v1,mu.cor)
mise.nblk.cor.3.v2=mise(mu.est.nblk.cor.3.v2,mu.cor)
mise.nblk.cor.3.v3=mise(mu.est.nblk.cor.3.v3,mu.cor)
mise.nblk.cor.3.v4=mise(mu.est.nblk.cor.3.v4,mu.cor)
mise.nblk.cor.3.v5=mise(mu.est.nblk.cor.3.v5,mu.cor)


mise.sure.sp.1.v1=mise(mu.est.sure.sp.1.v1,mu.sp)
mise.sure.sp.1.v2=mise(mu.est.sure.sp.1.v2,mu.sp)
mise.sure.sp.1.v3=mise(mu.est.sure.sp.1.v3,mu.sp)
mise.sure.sp.1.v4=mise(mu.est.sure.sp.1.v4,mu.sp)
mise.sure.sp.1.v5=mise(mu.est.sure.sp.1.v5,mu.sp)
mise.sure.sp.3.v1=mise(mu.est.sure.sp.3.v1,mu.sp)
mise.sure.sp.3.v2=mise(mu.est.sure.sp.3.v2,mu.sp)
mise.sure.sp.3.v3=mise(mu.est.sure.sp.3.v3,mu.sp)
mise.sure.sp.3.v4=mise(mu.est.sure.sp.3.v4,mu.sp)
mise.sure.sp.3.v5=mise(mu.est.sure.sp.3.v5,mu.sp)

mise.sure.bump.1.v1=mise(mu.est.sure.bump.1.v1,mu.bump)
mise.sure.bump.1.v2=mise(mu.est.sure.bump.1.v2,mu.bump)
mise.sure.bump.1.v3=mise(mu.est.sure.bump.1.v3,mu.bump)
mise.sure.bump.1.v4=mise(mu.est.sure.bump.1.v4,mu.bump)
mise.sure.bump.1.v5=mise(mu.est.sure.bump.1.v5,mu.bump)
mise.sure.bump.3.v1=mise(mu.est.sure.bump.3.v1,mu.bump)
mise.sure.bump.3.v2=mise(mu.est.sure.bump.3.v2,mu.bump)
mise.sure.bump.3.v3=mise(mu.est.sure.bump.3.v3,mu.bump)
mise.sure.bump.3.v4=mise(mu.est.sure.bump.3.v4,mu.bump)
mise.sure.bump.3.v5=mise(mu.est.sure.bump.3.v5,mu.bump)

mise.sure.blk.1.v1=mise(mu.est.sure.blk.1.v1,mu.blk)
mise.sure.blk.1.v2=mise(mu.est.sure.blk.1.v2,mu.blk)
mise.sure.blk.1.v3=mise(mu.est.sure.blk.1.v3,mu.blk)
mise.sure.blk.1.v4=mise(mu.est.sure.blk.1.v4,mu.blk)
mise.sure.blk.1.v5=mise(mu.est.sure.blk.1.v5,mu.blk)
mise.sure.blk.3.v1=mise(mu.est.sure.blk.3.v1,mu.blk)
mise.sure.blk.3.v2=mise(mu.est.sure.blk.3.v2,mu.blk)
mise.sure.blk.3.v3=mise(mu.est.sure.blk.3.v3,mu.blk)
mise.sure.blk.3.v4=mise(mu.est.sure.blk.3.v4,mu.blk)
mise.sure.blk.3.v5=mise(mu.est.sure.blk.3.v5,mu.blk)

mise.sure.ang.1.v1=mise(mu.est.sure.ang.1.v1,mu.ang)
mise.sure.ang.1.v2=mise(mu.est.sure.ang.1.v2,mu.ang)
mise.sure.ang.1.v3=mise(mu.est.sure.ang.1.v3,mu.ang)
mise.sure.ang.1.v4=mise(mu.est.sure.ang.1.v4,mu.ang)
mise.sure.ang.1.v5=mise(mu.est.sure.ang.1.v5,mu.ang)
mise.sure.ang.3.v1=mise(mu.est.sure.ang.3.v1,mu.ang)
mise.sure.ang.3.v2=mise(mu.est.sure.ang.3.v2,mu.ang)
mise.sure.ang.3.v3=mise(mu.est.sure.ang.3.v3,mu.ang)
mise.sure.ang.3.v4=mise(mu.est.sure.ang.3.v4,mu.ang)
mise.sure.ang.3.v5=mise(mu.est.sure.ang.3.v5,mu.ang)

mise.sure.dop.1.v1=mise(mu.est.sure.dop.1.v1,mu.dop)
mise.sure.dop.1.v2=mise(mu.est.sure.dop.1.v2,mu.dop)
mise.sure.dop.1.v3=mise(mu.est.sure.dop.1.v3,mu.dop)
mise.sure.dop.1.v4=mise(mu.est.sure.dop.1.v4,mu.dop)
mise.sure.dop.1.v5=mise(mu.est.sure.dop.1.v5,mu.dop)
mise.sure.dop.3.v1=mise(mu.est.sure.dop.3.v1,mu.dop)
mise.sure.dop.3.v2=mise(mu.est.sure.dop.3.v2,mu.dop)
mise.sure.dop.3.v3=mise(mu.est.sure.dop.3.v3,mu.dop)
mise.sure.dop.3.v4=mise(mu.est.sure.dop.3.v4,mu.dop)
mise.sure.dop.3.v5=mise(mu.est.sure.dop.3.v5,mu.dop)

mise.sure.blip.1.v1=mise(mu.est.sure.blip.1.v1,mu.blip)
mise.sure.blip.1.v2=mise(mu.est.sure.blip.1.v2,mu.blip)
mise.sure.blip.1.v3=mise(mu.est.sure.blip.1.v3,mu.blip)
mise.sure.blip.1.v4=mise(mu.est.sure.blip.1.v4,mu.blip)
mise.sure.blip.1.v5=mise(mu.est.sure.blip.1.v5,mu.blip)
mise.sure.blip.3.v1=mise(mu.est.sure.blip.3.v1,mu.blip)
mise.sure.blip.3.v2=mise(mu.est.sure.blip.3.v2,mu.blip)
mise.sure.blip.3.v3=mise(mu.est.sure.blip.3.v3,mu.blip)
mise.sure.blip.3.v4=mise(mu.est.sure.blip.3.v4,mu.blip)
mise.sure.blip.3.v5=mise(mu.est.sure.blip.3.v5,mu.blip)

mise.sure.cor.1.v1=mise(mu.est.sure.cor.1.v1,mu.cor)
mise.sure.cor.1.v2=mise(mu.est.sure.cor.1.v2,mu.cor)
mise.sure.cor.1.v3=mise(mu.est.sure.cor.1.v3,mu.cor)
mise.sure.cor.1.v4=mise(mu.est.sure.cor.1.v4,mu.cor)
mise.sure.cor.1.v5=mise(mu.est.sure.cor.1.v5,mu.cor)
mise.sure.cor.3.v1=mise(mu.est.sure.cor.3.v1,mu.cor)
mise.sure.cor.3.v2=mise(mu.est.sure.cor.3.v2,mu.cor)
mise.sure.cor.3.v3=mise(mu.est.sure.cor.3.v3,mu.cor)
mise.sure.cor.3.v4=mise(mu.est.sure.cor.3.v4,mu.cor)
mise.sure.cor.3.v5=mise(mu.est.sure.cor.3.v5,mu.cor)



mise.postmean.sp.1.v1=mise(mu.est.postmean.sp.1.v1,mu.sp)
mise.postmean.sp.1.v2=mise(mu.est.postmean.sp.1.v2,mu.sp)
mise.postmean.sp.1.v3=mise(mu.est.postmean.sp.1.v3,mu.sp)
mise.postmean.sp.1.v4=mise(mu.est.postmean.sp.1.v4,mu.sp)
mise.postmean.sp.1.v5=mise(mu.est.postmean.sp.1.v5,mu.sp)
mise.postmean.sp.3.v1=mise(mu.est.postmean.sp.3.v1,mu.sp)
mise.postmean.sp.3.v2=mise(mu.est.postmean.sp.3.v2,mu.sp)
mise.postmean.sp.3.v3=mise(mu.est.postmean.sp.3.v3,mu.sp)
mise.postmean.sp.3.v4=mise(mu.est.postmean.sp.3.v4,mu.sp)
mise.postmean.sp.3.v5=mise(mu.est.postmean.sp.3.v5,mu.sp)

mise.postmean.bump.1.v1=mise(mu.est.postmean.bump.1.v1,mu.bump)
mise.postmean.bump.1.v2=mise(mu.est.postmean.bump.1.v2,mu.bump)
mise.postmean.bump.1.v3=mise(mu.est.postmean.bump.1.v3,mu.bump)
mise.postmean.bump.1.v4=mise(mu.est.postmean.bump.1.v4,mu.bump)
mise.postmean.bump.1.v5=mise(mu.est.postmean.bump.1.v5,mu.bump)
mise.postmean.bump.3.v1=mise(mu.est.postmean.bump.3.v1,mu.bump)
mise.postmean.bump.3.v2=mise(mu.est.postmean.bump.3.v2,mu.bump)
mise.postmean.bump.3.v3=mise(mu.est.postmean.bump.3.v3,mu.bump)
mise.postmean.bump.3.v4=mise(mu.est.postmean.bump.3.v4,mu.bump)
mise.postmean.bump.3.v5=mise(mu.est.postmean.bump.3.v5,mu.bump)

mise.postmean.blk.1.v1=mise(mu.est.postmean.blk.1.v1,mu.blk)
mise.postmean.blk.1.v2=mise(mu.est.postmean.blk.1.v2,mu.blk)
mise.postmean.blk.1.v3=mise(mu.est.postmean.blk.1.v3,mu.blk)
mise.postmean.blk.1.v4=mise(mu.est.postmean.blk.1.v4,mu.blk)
mise.postmean.blk.1.v5=mise(mu.est.postmean.blk.1.v5,mu.blk)
mise.postmean.blk.3.v1=mise(mu.est.postmean.blk.3.v1,mu.blk)
mise.postmean.blk.3.v2=mise(mu.est.postmean.blk.3.v2,mu.blk)
mise.postmean.blk.3.v3=mise(mu.est.postmean.blk.3.v3,mu.blk)
mise.postmean.blk.3.v4=mise(mu.est.postmean.blk.3.v4,mu.blk)
mise.postmean.blk.3.v5=mise(mu.est.postmean.blk.3.v5,mu.blk)

mise.postmean.ang.1.v1=mise(mu.est.postmean.ang.1.v1,mu.ang)
mise.postmean.ang.1.v2=mise(mu.est.postmean.ang.1.v2,mu.ang)
mise.postmean.ang.1.v3=mise(mu.est.postmean.ang.1.v3,mu.ang)
mise.postmean.ang.1.v4=mise(mu.est.postmean.ang.1.v4,mu.ang)
mise.postmean.ang.1.v5=mise(mu.est.postmean.ang.1.v5,mu.ang)
mise.postmean.ang.3.v1=mise(mu.est.postmean.ang.3.v1,mu.ang)
mise.postmean.ang.3.v2=mise(mu.est.postmean.ang.3.v2,mu.ang)
mise.postmean.ang.3.v3=mise(mu.est.postmean.ang.3.v3,mu.ang)
mise.postmean.ang.3.v4=mise(mu.est.postmean.ang.3.v4,mu.ang)
mise.postmean.ang.3.v5=mise(mu.est.postmean.ang.3.v5,mu.ang)

mise.postmean.dop.1.v1=mise(mu.est.postmean.dop.1.v1,mu.dop)
mise.postmean.dop.1.v2=mise(mu.est.postmean.dop.1.v2,mu.dop)
mise.postmean.dop.1.v3=mise(mu.est.postmean.dop.1.v3,mu.dop)
mise.postmean.dop.1.v4=mise(mu.est.postmean.dop.1.v4,mu.dop)
mise.postmean.dop.1.v5=mise(mu.est.postmean.dop.1.v5,mu.dop)
mise.postmean.dop.3.v1=mise(mu.est.postmean.dop.3.v1,mu.dop)
mise.postmean.dop.3.v2=mise(mu.est.postmean.dop.3.v2,mu.dop)
mise.postmean.dop.3.v3=mise(mu.est.postmean.dop.3.v3,mu.dop)
mise.postmean.dop.3.v4=mise(mu.est.postmean.dop.3.v4,mu.dop)
mise.postmean.dop.3.v5=mise(mu.est.postmean.dop.3.v5,mu.dop)

mise.postmean.blip.1.v1=mise(mu.est.postmean.blip.1.v1,mu.blip)
mise.postmean.blip.1.v2=mise(mu.est.postmean.blip.1.v2,mu.blip)
mise.postmean.blip.1.v3=mise(mu.est.postmean.blip.1.v3,mu.blip)
mise.postmean.blip.1.v4=mise(mu.est.postmean.blip.1.v4,mu.blip)
mise.postmean.blip.1.v5=mise(mu.est.postmean.blip.1.v5,mu.blip)
mise.postmean.blip.3.v1=mise(mu.est.postmean.blip.3.v1,mu.blip)
mise.postmean.blip.3.v2=mise(mu.est.postmean.blip.3.v2,mu.blip)
mise.postmean.blip.3.v3=mise(mu.est.postmean.blip.3.v3,mu.blip)
mise.postmean.blip.3.v4=mise(mu.est.postmean.blip.3.v4,mu.blip)
mise.postmean.blip.3.v5=mise(mu.est.postmean.blip.3.v5,mu.blip)

mise.postmean.cor.1.v1=mise(mu.est.postmean.cor.1.v1,mu.cor)
mise.postmean.cor.1.v2=mise(mu.est.postmean.cor.1.v2,mu.cor)
mise.postmean.cor.1.v3=mise(mu.est.postmean.cor.1.v3,mu.cor)
mise.postmean.cor.1.v4=mise(mu.est.postmean.cor.1.v4,mu.cor)
mise.postmean.cor.1.v5=mise(mu.est.postmean.cor.1.v5,mu.cor)
mise.postmean.cor.3.v1=mise(mu.est.postmean.cor.3.v1,mu.cor)
mise.postmean.cor.3.v2=mise(mu.est.postmean.cor.3.v2,mu.cor)
mise.postmean.cor.3.v3=mise(mu.est.postmean.cor.3.v3,mu.cor)
mise.postmean.cor.3.v4=mise(mu.est.postmean.cor.3.v4,mu.cor)
mise.postmean.cor.3.v5=mise(mu.est.postmean.cor.3.v5,mu.cor)




mise.ash.haar.sp.1.v1=mise(mu.est.ash.haar.sp.1.v1,mu.sp)
mise.ash.haar.sp.1.v2=mise(mu.est.ash.haar.sp.1.v2,mu.sp)
mise.ash.haar.sp.1.v3=mise(mu.est.ash.haar.sp.1.v3,mu.sp)
mise.ash.haar.sp.1.v4=mise(mu.est.ash.haar.sp.1.v4,mu.sp)
mise.ash.haar.sp.1.v5=mise(mu.est.ash.haar.sp.1.v5,mu.sp)
mise.ash.haar.sp.3.v1=mise(mu.est.ash.haar.sp.3.v1,mu.sp)
mise.ash.haar.sp.3.v2=mise(mu.est.ash.haar.sp.3.v2,mu.sp)
mise.ash.haar.sp.3.v3=mise(mu.est.ash.haar.sp.3.v3,mu.sp)
mise.ash.haar.sp.3.v4=mise(mu.est.ash.haar.sp.3.v4,mu.sp)
mise.ash.haar.sp.3.v5=mise(mu.est.ash.haar.sp.3.v5,mu.sp)

mise.ash.haar.bump.1.v1=mise(mu.est.ash.haar.bump.1.v1,mu.bump)
mise.ash.haar.bump.1.v2=mise(mu.est.ash.haar.bump.1.v2,mu.bump)
mise.ash.haar.bump.1.v3=mise(mu.est.ash.haar.bump.1.v3,mu.bump)
mise.ash.haar.bump.1.v4=mise(mu.est.ash.haar.bump.1.v4,mu.bump)
mise.ash.haar.bump.1.v5=mise(mu.est.ash.haar.bump.1.v5,mu.bump)
mise.ash.haar.bump.3.v1=mise(mu.est.ash.haar.bump.3.v1,mu.bump)
mise.ash.haar.bump.3.v2=mise(mu.est.ash.haar.bump.3.v2,mu.bump)
mise.ash.haar.bump.3.v3=mise(mu.est.ash.haar.bump.3.v3,mu.bump)
mise.ash.haar.bump.3.v4=mise(mu.est.ash.haar.bump.3.v4,mu.bump)
mise.ash.haar.bump.3.v5=mise(mu.est.ash.haar.bump.3.v5,mu.bump)

mise.ash.haar.blk.1.v1=mise(mu.est.ash.haar.blk.1.v1,mu.blk)
mise.ash.haar.blk.1.v2=mise(mu.est.ash.haar.blk.1.v2,mu.blk)
mise.ash.haar.blk.1.v3=mise(mu.est.ash.haar.blk.1.v3,mu.blk)
mise.ash.haar.blk.1.v4=mise(mu.est.ash.haar.blk.1.v4,mu.blk)
mise.ash.haar.blk.1.v5=mise(mu.est.ash.haar.blk.1.v5,mu.blk)
mise.ash.haar.blk.3.v1=mise(mu.est.ash.haar.blk.3.v1,mu.blk)
mise.ash.haar.blk.3.v2=mise(mu.est.ash.haar.blk.3.v2,mu.blk)
mise.ash.haar.blk.3.v3=mise(mu.est.ash.haar.blk.3.v3,mu.blk)
mise.ash.haar.blk.3.v4=mise(mu.est.ash.haar.blk.3.v4,mu.blk)
mise.ash.haar.blk.3.v5=mise(mu.est.ash.haar.blk.3.v5,mu.blk)

mise.ash.haar.ang.1.v1=mise(mu.est.ash.haar.ang.1.v1,mu.ang)
mise.ash.haar.ang.1.v2=mise(mu.est.ash.haar.ang.1.v2,mu.ang)
mise.ash.haar.ang.1.v3=mise(mu.est.ash.haar.ang.1.v3,mu.ang)
mise.ash.haar.ang.1.v4=mise(mu.est.ash.haar.ang.1.v4,mu.ang)
mise.ash.haar.ang.1.v5=mise(mu.est.ash.haar.ang.1.v5,mu.ang)
mise.ash.haar.ang.3.v1=mise(mu.est.ash.haar.ang.3.v1,mu.ang)
mise.ash.haar.ang.3.v2=mise(mu.est.ash.haar.ang.3.v2,mu.ang)
mise.ash.haar.ang.3.v3=mise(mu.est.ash.haar.ang.3.v3,mu.ang)
mise.ash.haar.ang.3.v4=mise(mu.est.ash.haar.ang.3.v4,mu.ang)
mise.ash.haar.ang.3.v5=mise(mu.est.ash.haar.ang.3.v5,mu.ang)

mise.ash.haar.dop.1.v1=mise(mu.est.ash.haar.dop.1.v1,mu.dop)
mise.ash.haar.dop.1.v2=mise(mu.est.ash.haar.dop.1.v2,mu.dop)
mise.ash.haar.dop.1.v3=mise(mu.est.ash.haar.dop.1.v3,mu.dop)
mise.ash.haar.dop.1.v4=mise(mu.est.ash.haar.dop.1.v4,mu.dop)
mise.ash.haar.dop.1.v5=mise(mu.est.ash.haar.dop.1.v5,mu.dop)
mise.ash.haar.dop.3.v1=mise(mu.est.ash.haar.dop.3.v1,mu.dop)
mise.ash.haar.dop.3.v2=mise(mu.est.ash.haar.dop.3.v2,mu.dop)
mise.ash.haar.dop.3.v3=mise(mu.est.ash.haar.dop.3.v3,mu.dop)
mise.ash.haar.dop.3.v4=mise(mu.est.ash.haar.dop.3.v4,mu.dop)
mise.ash.haar.dop.3.v5=mise(mu.est.ash.haar.dop.3.v5,mu.dop)

mise.ash.haar.blip.1.v1=mise(mu.est.ash.haar.blip.1.v1,mu.blip)
mise.ash.haar.blip.1.v2=mise(mu.est.ash.haar.blip.1.v2,mu.blip)
mise.ash.haar.blip.1.v3=mise(mu.est.ash.haar.blip.1.v3,mu.blip)
mise.ash.haar.blip.1.v4=mise(mu.est.ash.haar.blip.1.v4,mu.blip)
mise.ash.haar.blip.1.v5=mise(mu.est.ash.haar.blip.1.v5,mu.blip)
mise.ash.haar.blip.3.v1=mise(mu.est.ash.haar.blip.3.v1,mu.blip)
mise.ash.haar.blip.3.v2=mise(mu.est.ash.haar.blip.3.v2,mu.blip)
mise.ash.haar.blip.3.v3=mise(mu.est.ash.haar.blip.3.v3,mu.blip)
mise.ash.haar.blip.3.v4=mise(mu.est.ash.haar.blip.3.v4,mu.blip)
mise.ash.haar.blip.3.v5=mise(mu.est.ash.haar.blip.3.v5,mu.blip)

mise.ash.haar.cor.1.v1=mise(mu.est.ash.haar.cor.1.v1,mu.cor)
mise.ash.haar.cor.1.v2=mise(mu.est.ash.haar.cor.1.v2,mu.cor)
mise.ash.haar.cor.1.v3=mise(mu.est.ash.haar.cor.1.v3,mu.cor)
mise.ash.haar.cor.1.v4=mise(mu.est.ash.haar.cor.1.v4,mu.cor)
mise.ash.haar.cor.1.v5=mise(mu.est.ash.haar.cor.1.v5,mu.cor)
mise.ash.haar.cor.3.v1=mise(mu.est.ash.haar.cor.3.v1,mu.cor)
mise.ash.haar.cor.3.v2=mise(mu.est.ash.haar.cor.3.v2,mu.cor)
mise.ash.haar.cor.3.v3=mise(mu.est.ash.haar.cor.3.v3,mu.cor)
mise.ash.haar.cor.3.v4=mise(mu.est.ash.haar.cor.3.v4,mu.cor)
mise.ash.haar.cor.3.v5=mise(mu.est.ash.haar.cor.3.v5,mu.cor)


mise.ash.s8.sp.1.v1=mise(mu.est.ash.s8.sp.1.v1,mu.sp)
mise.ash.s8.sp.1.v2=mise(mu.est.ash.s8.sp.1.v2,mu.sp)
mise.ash.s8.sp.1.v3=mise(mu.est.ash.s8.sp.1.v3,mu.sp)
mise.ash.s8.sp.1.v4=mise(mu.est.ash.s8.sp.1.v4,mu.sp)
mise.ash.s8.sp.1.v5=mise(mu.est.ash.s8.sp.1.v5,mu.sp)
mise.ash.s8.sp.3.v1=mise(mu.est.ash.s8.sp.3.v1,mu.sp)
mise.ash.s8.sp.3.v2=mise(mu.est.ash.s8.sp.3.v2,mu.sp)
mise.ash.s8.sp.3.v3=mise(mu.est.ash.s8.sp.3.v3,mu.sp)
mise.ash.s8.sp.3.v4=mise(mu.est.ash.s8.sp.3.v4,mu.sp)
mise.ash.s8.sp.3.v5=mise(mu.est.ash.s8.sp.3.v5,mu.sp)

mise.ash.s8.bump.1.v1=mise(mu.est.ash.s8.bump.1.v1,mu.bump)
mise.ash.s8.bump.1.v2=mise(mu.est.ash.s8.bump.1.v2,mu.bump)
mise.ash.s8.bump.1.v3=mise(mu.est.ash.s8.bump.1.v3,mu.bump)
mise.ash.s8.bump.1.v4=mise(mu.est.ash.s8.bump.1.v4,mu.bump)
mise.ash.s8.bump.1.v5=mise(mu.est.ash.s8.bump.1.v5,mu.bump)
mise.ash.s8.bump.3.v1=mise(mu.est.ash.s8.bump.3.v1,mu.bump)
mise.ash.s8.bump.3.v2=mise(mu.est.ash.s8.bump.3.v2,mu.bump)
mise.ash.s8.bump.3.v3=mise(mu.est.ash.s8.bump.3.v3,mu.bump)
mise.ash.s8.bump.3.v4=mise(mu.est.ash.s8.bump.3.v4,mu.bump)
mise.ash.s8.bump.3.v5=mise(mu.est.ash.s8.bump.3.v5,mu.bump)

mise.ash.s8.blk.1.v1=mise(mu.est.ash.s8.blk.1.v1,mu.blk)
mise.ash.s8.blk.1.v2=mise(mu.est.ash.s8.blk.1.v2,mu.blk)
mise.ash.s8.blk.1.v3=mise(mu.est.ash.s8.blk.1.v3,mu.blk)
mise.ash.s8.blk.1.v4=mise(mu.est.ash.s8.blk.1.v4,mu.blk)
mise.ash.s8.blk.1.v5=mise(mu.est.ash.s8.blk.1.v5,mu.blk)
mise.ash.s8.blk.3.v1=mise(mu.est.ash.s8.blk.3.v1,mu.blk)
mise.ash.s8.blk.3.v2=mise(mu.est.ash.s8.blk.3.v2,mu.blk)
mise.ash.s8.blk.3.v3=mise(mu.est.ash.s8.blk.3.v3,mu.blk)
mise.ash.s8.blk.3.v4=mise(mu.est.ash.s8.blk.3.v4,mu.blk)
mise.ash.s8.blk.3.v5=mise(mu.est.ash.s8.blk.3.v5,mu.blk)

mise.ash.s8.ang.1.v1=mise(mu.est.ash.s8.ang.1.v1,mu.ang)
mise.ash.s8.ang.1.v2=mise(mu.est.ash.s8.ang.1.v2,mu.ang)
mise.ash.s8.ang.1.v3=mise(mu.est.ash.s8.ang.1.v3,mu.ang)
mise.ash.s8.ang.1.v4=mise(mu.est.ash.s8.ang.1.v4,mu.ang)
mise.ash.s8.ang.1.v5=mise(mu.est.ash.s8.ang.1.v5,mu.ang)
mise.ash.s8.ang.3.v1=mise(mu.est.ash.s8.ang.3.v1,mu.ang)
mise.ash.s8.ang.3.v2=mise(mu.est.ash.s8.ang.3.v2,mu.ang)
mise.ash.s8.ang.3.v3=mise(mu.est.ash.s8.ang.3.v3,mu.ang)
mise.ash.s8.ang.3.v4=mise(mu.est.ash.s8.ang.3.v4,mu.ang)
mise.ash.s8.ang.3.v5=mise(mu.est.ash.s8.ang.3.v5,mu.ang)

mise.ash.s8.dop.1.v1=mise(mu.est.ash.s8.dop.1.v1,mu.dop)
mise.ash.s8.dop.1.v2=mise(mu.est.ash.s8.dop.1.v2,mu.dop)
mise.ash.s8.dop.1.v3=mise(mu.est.ash.s8.dop.1.v3,mu.dop)
mise.ash.s8.dop.1.v4=mise(mu.est.ash.s8.dop.1.v4,mu.dop)
mise.ash.s8.dop.1.v5=mise(mu.est.ash.s8.dop.1.v5,mu.dop)
mise.ash.s8.dop.3.v1=mise(mu.est.ash.s8.dop.3.v1,mu.dop)
mise.ash.s8.dop.3.v2=mise(mu.est.ash.s8.dop.3.v2,mu.dop)
mise.ash.s8.dop.3.v3=mise(mu.est.ash.s8.dop.3.v3,mu.dop)
mise.ash.s8.dop.3.v4=mise(mu.est.ash.s8.dop.3.v4,mu.dop)
mise.ash.s8.dop.3.v5=mise(mu.est.ash.s8.dop.3.v5,mu.dop)

mise.ash.s8.blip.1.v1=mise(mu.est.ash.s8.blip.1.v1,mu.blip)
mise.ash.s8.blip.1.v2=mise(mu.est.ash.s8.blip.1.v2,mu.blip)
mise.ash.s8.blip.1.v3=mise(mu.est.ash.s8.blip.1.v3,mu.blip)
mise.ash.s8.blip.1.v4=mise(mu.est.ash.s8.blip.1.v4,mu.blip)
mise.ash.s8.blip.1.v5=mise(mu.est.ash.s8.blip.1.v5,mu.blip)
mise.ash.s8.blip.3.v1=mise(mu.est.ash.s8.blip.3.v1,mu.blip)
mise.ash.s8.blip.3.v2=mise(mu.est.ash.s8.blip.3.v2,mu.blip)
mise.ash.s8.blip.3.v3=mise(mu.est.ash.s8.blip.3.v3,mu.blip)
mise.ash.s8.blip.3.v4=mise(mu.est.ash.s8.blip.3.v4,mu.blip)
mise.ash.s8.blip.3.v5=mise(mu.est.ash.s8.blip.3.v5,mu.blip)

mise.ash.s8.cor.1.v1=mise(mu.est.ash.s8.cor.1.v1,mu.cor)
mise.ash.s8.cor.1.v2=mise(mu.est.ash.s8.cor.1.v2,mu.cor)
mise.ash.s8.cor.1.v3=mise(mu.est.ash.s8.cor.1.v3,mu.cor)
mise.ash.s8.cor.1.v4=mise(mu.est.ash.s8.cor.1.v4,mu.cor)
mise.ash.s8.cor.1.v5=mise(mu.est.ash.s8.cor.1.v5,mu.cor)
mise.ash.s8.cor.3.v1=mise(mu.est.ash.s8.cor.3.v1,mu.cor)
mise.ash.s8.cor.3.v2=mise(mu.est.ash.s8.cor.3.v2,mu.cor)
mise.ash.s8.cor.3.v3=mise(mu.est.ash.s8.cor.3.v3,mu.cor)
mise.ash.s8.cor.3.v4=mise(mu.est.ash.s8.cor.3.v4,mu.cor)
mise.ash.s8.cor.3.v5=mise(mu.est.ash.s8.cor.3.v5,mu.cor)



mise.ashe.haar.sp.1.v1=mise(mu.est.ashe.haar.sp.1.v1,mu.sp)
mise.ashe.haar.sp.1.v2=mise(mu.est.ashe.haar.sp.1.v2,mu.sp)
mise.ashe.haar.sp.1.v3=mise(mu.est.ashe.haar.sp.1.v3,mu.sp)
mise.ashe.haar.sp.1.v4=mise(mu.est.ashe.haar.sp.1.v4,mu.sp)
mise.ashe.haar.sp.1.v5=mise(mu.est.ashe.haar.sp.1.v5,mu.sp)
mise.ashe.haar.sp.3.v1=mise(mu.est.ashe.haar.sp.3.v1,mu.sp)
mise.ashe.haar.sp.3.v2=mise(mu.est.ashe.haar.sp.3.v2,mu.sp)
mise.ashe.haar.sp.3.v3=mise(mu.est.ashe.haar.sp.3.v3,mu.sp)
mise.ashe.haar.sp.3.v4=mise(mu.est.ashe.haar.sp.3.v4,mu.sp)
mise.ashe.haar.sp.3.v5=mise(mu.est.ashe.haar.sp.3.v5,mu.sp)

mise.ashe.haar.bump.1.v1=mise(mu.est.ashe.haar.bump.1.v1,mu.bump)
mise.ashe.haar.bump.1.v2=mise(mu.est.ashe.haar.bump.1.v2,mu.bump)
mise.ashe.haar.bump.1.v3=mise(mu.est.ashe.haar.bump.1.v3,mu.bump)
mise.ashe.haar.bump.1.v4=mise(mu.est.ashe.haar.bump.1.v4,mu.bump)
mise.ashe.haar.bump.1.v5=mise(mu.est.ashe.haar.bump.1.v5,mu.bump)
mise.ashe.haar.bump.3.v1=mise(mu.est.ashe.haar.bump.3.v1,mu.bump)
mise.ashe.haar.bump.3.v2=mise(mu.est.ashe.haar.bump.3.v2,mu.bump)
mise.ashe.haar.bump.3.v3=mise(mu.est.ashe.haar.bump.3.v3,mu.bump)
mise.ashe.haar.bump.3.v4=mise(mu.est.ashe.haar.bump.3.v4,mu.bump)
mise.ashe.haar.bump.3.v5=mise(mu.est.ashe.haar.bump.3.v5,mu.bump)

mise.ashe.haar.blk.1.v1=mise(mu.est.ashe.haar.blk.1.v1,mu.blk)
mise.ashe.haar.blk.1.v2=mise(mu.est.ashe.haar.blk.1.v2,mu.blk)
mise.ashe.haar.blk.1.v3=mise(mu.est.ashe.haar.blk.1.v3,mu.blk)
mise.ashe.haar.blk.1.v4=mise(mu.est.ashe.haar.blk.1.v4,mu.blk)
mise.ashe.haar.blk.1.v5=mise(mu.est.ashe.haar.blk.1.v5,mu.blk)
mise.ashe.haar.blk.3.v1=mise(mu.est.ashe.haar.blk.3.v1,mu.blk)
mise.ashe.haar.blk.3.v2=mise(mu.est.ashe.haar.blk.3.v2,mu.blk)
mise.ashe.haar.blk.3.v3=mise(mu.est.ashe.haar.blk.3.v3,mu.blk)
mise.ashe.haar.blk.3.v4=mise(mu.est.ashe.haar.blk.3.v4,mu.blk)
mise.ashe.haar.blk.3.v5=mise(mu.est.ashe.haar.blk.3.v5,mu.blk)

mise.ashe.haar.ang.1.v1=mise(mu.est.ashe.haar.ang.1.v1,mu.ang)
mise.ashe.haar.ang.1.v2=mise(mu.est.ashe.haar.ang.1.v2,mu.ang)
mise.ashe.haar.ang.1.v3=mise(mu.est.ashe.haar.ang.1.v3,mu.ang)
mise.ashe.haar.ang.1.v4=mise(mu.est.ashe.haar.ang.1.v4,mu.ang)
mise.ashe.haar.ang.1.v5=mise(mu.est.ashe.haar.ang.1.v5,mu.ang)
mise.ashe.haar.ang.3.v1=mise(mu.est.ashe.haar.ang.3.v1,mu.ang)
mise.ashe.haar.ang.3.v2=mise(mu.est.ashe.haar.ang.3.v2,mu.ang)
mise.ashe.haar.ang.3.v3=mise(mu.est.ashe.haar.ang.3.v3,mu.ang)
mise.ashe.haar.ang.3.v4=mise(mu.est.ashe.haar.ang.3.v4,mu.ang)
mise.ashe.haar.ang.3.v5=mise(mu.est.ashe.haar.ang.3.v5,mu.ang)

mise.ashe.haar.dop.1.v1=mise(mu.est.ashe.haar.dop.1.v1,mu.dop)
mise.ashe.haar.dop.1.v2=mise(mu.est.ashe.haar.dop.1.v2,mu.dop)
mise.ashe.haar.dop.1.v3=mise(mu.est.ashe.haar.dop.1.v3,mu.dop)
mise.ashe.haar.dop.1.v4=mise(mu.est.ashe.haar.dop.1.v4,mu.dop)
mise.ashe.haar.dop.1.v5=mise(mu.est.ashe.haar.dop.1.v5,mu.dop)
mise.ashe.haar.dop.3.v1=mise(mu.est.ashe.haar.dop.3.v1,mu.dop)
mise.ashe.haar.dop.3.v2=mise(mu.est.ashe.haar.dop.3.v2,mu.dop)
mise.ashe.haar.dop.3.v3=mise(mu.est.ashe.haar.dop.3.v3,mu.dop)
mise.ashe.haar.dop.3.v4=mise(mu.est.ashe.haar.dop.3.v4,mu.dop)
mise.ashe.haar.dop.3.v5=mise(mu.est.ashe.haar.dop.3.v5,mu.dop)

mise.ashe.haar.blip.1.v1=mise(mu.est.ashe.haar.blip.1.v1,mu.blip)
mise.ashe.haar.blip.1.v2=mise(mu.est.ashe.haar.blip.1.v2,mu.blip)
mise.ashe.haar.blip.1.v3=mise(mu.est.ashe.haar.blip.1.v3,mu.blip)
mise.ashe.haar.blip.1.v4=mise(mu.est.ashe.haar.blip.1.v4,mu.blip)
mise.ashe.haar.blip.1.v5=mise(mu.est.ashe.haar.blip.1.v5,mu.blip)
mise.ashe.haar.blip.3.v1=mise(mu.est.ashe.haar.blip.3.v1,mu.blip)
mise.ashe.haar.blip.3.v2=mise(mu.est.ashe.haar.blip.3.v2,mu.blip)
mise.ashe.haar.blip.3.v3=mise(mu.est.ashe.haar.blip.3.v3,mu.blip)
mise.ashe.haar.blip.3.v4=mise(mu.est.ashe.haar.blip.3.v4,mu.blip)
mise.ashe.haar.blip.3.v5=mise(mu.est.ashe.haar.blip.3.v5,mu.blip)

mise.ashe.haar.cor.1.v1=mise(mu.est.ashe.haar.cor.1.v1,mu.cor)
mise.ashe.haar.cor.1.v2=mise(mu.est.ashe.haar.cor.1.v2,mu.cor)
mise.ashe.haar.cor.1.v3=mise(mu.est.ashe.haar.cor.1.v3,mu.cor)
mise.ashe.haar.cor.1.v4=mise(mu.est.ashe.haar.cor.1.v4,mu.cor)
mise.ashe.haar.cor.1.v5=mise(mu.est.ashe.haar.cor.1.v5,mu.cor)
mise.ashe.haar.cor.3.v1=mise(mu.est.ashe.haar.cor.3.v1,mu.cor)
mise.ashe.haar.cor.3.v2=mise(mu.est.ashe.haar.cor.3.v2,mu.cor)
mise.ashe.haar.cor.3.v3=mise(mu.est.ashe.haar.cor.3.v3,mu.cor)
mise.ashe.haar.cor.3.v4=mise(mu.est.ashe.haar.cor.3.v4,mu.cor)
mise.ashe.haar.cor.3.v5=mise(mu.est.ashe.haar.cor.3.v5,mu.cor)


mise.ashe.s8.sp.1.v1=mise(mu.est.ashe.s8.sp.1.v1,mu.sp)
mise.ashe.s8.sp.1.v2=mise(mu.est.ashe.s8.sp.1.v2,mu.sp)
mise.ashe.s8.sp.1.v3=mise(mu.est.ashe.s8.sp.1.v3,mu.sp)
mise.ashe.s8.sp.1.v4=mise(mu.est.ashe.s8.sp.1.v4,mu.sp)
mise.ashe.s8.sp.1.v5=mise(mu.est.ashe.s8.sp.1.v5,mu.sp)
mise.ashe.s8.sp.3.v1=mise(mu.est.ashe.s8.sp.3.v1,mu.sp)
mise.ashe.s8.sp.3.v2=mise(mu.est.ashe.s8.sp.3.v2,mu.sp)
mise.ashe.s8.sp.3.v3=mise(mu.est.ashe.s8.sp.3.v3,mu.sp)
mise.ashe.s8.sp.3.v4=mise(mu.est.ashe.s8.sp.3.v4,mu.sp)
mise.ashe.s8.sp.3.v5=mise(mu.est.ashe.s8.sp.3.v5,mu.sp)

mise.ashe.s8.bump.1.v1=mise(mu.est.ashe.s8.bump.1.v1,mu.bump)
mise.ashe.s8.bump.1.v2=mise(mu.est.ashe.s8.bump.1.v2,mu.bump)
mise.ashe.s8.bump.1.v3=mise(mu.est.ashe.s8.bump.1.v3,mu.bump)
mise.ashe.s8.bump.1.v4=mise(mu.est.ashe.s8.bump.1.v4,mu.bump)
mise.ashe.s8.bump.1.v5=mise(mu.est.ashe.s8.bump.1.v5,mu.bump)
mise.ashe.s8.bump.3.v1=mise(mu.est.ashe.s8.bump.3.v1,mu.bump)
mise.ashe.s8.bump.3.v2=mise(mu.est.ashe.s8.bump.3.v2,mu.bump)
mise.ashe.s8.bump.3.v3=mise(mu.est.ashe.s8.bump.3.v3,mu.bump)
mise.ashe.s8.bump.3.v4=mise(mu.est.ashe.s8.bump.3.v4,mu.bump)
mise.ashe.s8.bump.3.v5=mise(mu.est.ashe.s8.bump.3.v5,mu.bump)

mise.ashe.s8.blk.1.v1=mise(mu.est.ashe.s8.blk.1.v1,mu.blk)
mise.ashe.s8.blk.1.v2=mise(mu.est.ashe.s8.blk.1.v2,mu.blk)
mise.ashe.s8.blk.1.v3=mise(mu.est.ashe.s8.blk.1.v3,mu.blk)
mise.ashe.s8.blk.1.v4=mise(mu.est.ashe.s8.blk.1.v4,mu.blk)
mise.ashe.s8.blk.1.v5=mise(mu.est.ashe.s8.blk.1.v5,mu.blk)
mise.ashe.s8.blk.3.v1=mise(mu.est.ashe.s8.blk.3.v1,mu.blk)
mise.ashe.s8.blk.3.v2=mise(mu.est.ashe.s8.blk.3.v2,mu.blk)
mise.ashe.s8.blk.3.v3=mise(mu.est.ashe.s8.blk.3.v3,mu.blk)
mise.ashe.s8.blk.3.v4=mise(mu.est.ashe.s8.blk.3.v4,mu.blk)
mise.ashe.s8.blk.3.v5=mise(mu.est.ashe.s8.blk.3.v5,mu.blk)

mise.ashe.s8.ang.1.v1=mise(mu.est.ashe.s8.ang.1.v1,mu.ang)
mise.ashe.s8.ang.1.v2=mise(mu.est.ashe.s8.ang.1.v2,mu.ang)
mise.ashe.s8.ang.1.v3=mise(mu.est.ashe.s8.ang.1.v3,mu.ang)
mise.ashe.s8.ang.1.v4=mise(mu.est.ashe.s8.ang.1.v4,mu.ang)
mise.ashe.s8.ang.1.v5=mise(mu.est.ashe.s8.ang.1.v5,mu.ang)
mise.ashe.s8.ang.3.v1=mise(mu.est.ashe.s8.ang.3.v1,mu.ang)
mise.ashe.s8.ang.3.v2=mise(mu.est.ashe.s8.ang.3.v2,mu.ang)
mise.ashe.s8.ang.3.v3=mise(mu.est.ashe.s8.ang.3.v3,mu.ang)
mise.ashe.s8.ang.3.v4=mise(mu.est.ashe.s8.ang.3.v4,mu.ang)
mise.ashe.s8.ang.3.v5=mise(mu.est.ashe.s8.ang.3.v5,mu.ang)

mise.ashe.s8.dop.1.v1=mise(mu.est.ashe.s8.dop.1.v1,mu.dop)
mise.ashe.s8.dop.1.v2=mise(mu.est.ashe.s8.dop.1.v2,mu.dop)
mise.ashe.s8.dop.1.v3=mise(mu.est.ashe.s8.dop.1.v3,mu.dop)
mise.ashe.s8.dop.1.v4=mise(mu.est.ashe.s8.dop.1.v4,mu.dop)
mise.ashe.s8.dop.1.v5=mise(mu.est.ashe.s8.dop.1.v5,mu.dop)
mise.ashe.s8.dop.3.v1=mise(mu.est.ashe.s8.dop.3.v1,mu.dop)
mise.ashe.s8.dop.3.v2=mise(mu.est.ashe.s8.dop.3.v2,mu.dop)
mise.ashe.s8.dop.3.v3=mise(mu.est.ashe.s8.dop.3.v3,mu.dop)
mise.ashe.s8.dop.3.v4=mise(mu.est.ashe.s8.dop.3.v4,mu.dop)
mise.ashe.s8.dop.3.v5=mise(mu.est.ashe.s8.dop.3.v5,mu.dop)

mise.ashe.s8.blip.1.v1=mise(mu.est.ashe.s8.blip.1.v1,mu.blip)
mise.ashe.s8.blip.1.v2=mise(mu.est.ashe.s8.blip.1.v2,mu.blip)
mise.ashe.s8.blip.1.v3=mise(mu.est.ashe.s8.blip.1.v3,mu.blip)
mise.ashe.s8.blip.1.v4=mise(mu.est.ashe.s8.blip.1.v4,mu.blip)
mise.ashe.s8.blip.1.v5=mise(mu.est.ashe.s8.blip.1.v5,mu.blip)
mise.ashe.s8.blip.3.v1=mise(mu.est.ashe.s8.blip.3.v1,mu.blip)
mise.ashe.s8.blip.3.v2=mise(mu.est.ashe.s8.blip.3.v2,mu.blip)
mise.ashe.s8.blip.3.v3=mise(mu.est.ashe.s8.blip.3.v3,mu.blip)
mise.ashe.s8.blip.3.v4=mise(mu.est.ashe.s8.blip.3.v4,mu.blip)
mise.ashe.s8.blip.3.v5=mise(mu.est.ashe.s8.blip.3.v5,mu.blip)

mise.ashe.s8.cor.1.v1=mise(mu.est.ashe.s8.cor.1.v1,mu.cor)
mise.ashe.s8.cor.1.v2=mise(mu.est.ashe.s8.cor.1.v2,mu.cor)
mise.ashe.s8.cor.1.v3=mise(mu.est.ashe.s8.cor.1.v3,mu.cor)
mise.ashe.s8.cor.1.v4=mise(mu.est.ashe.s8.cor.1.v4,mu.cor)
mise.ashe.s8.cor.1.v5=mise(mu.est.ashe.s8.cor.1.v5,mu.cor)
mise.ashe.s8.cor.3.v1=mise(mu.est.ashe.s8.cor.3.v1,mu.cor)
mise.ashe.s8.cor.3.v2=mise(mu.est.ashe.s8.cor.3.v2,mu.cor)
mise.ashe.s8.cor.3.v3=mise(mu.est.ashe.s8.cor.3.v3,mu.cor)
mise.ashe.s8.cor.3.v4=mise(mu.est.ashe.s8.cor.3.v4,mu.cor)
mise.ashe.s8.cor.3.v5=mise(mu.est.ashe.s8.cor.3.v5,mu.cor)



mise.tie.haar.sp.1.v1=mise(mu.est.tie.haar.sp.1.v1,mu.sp)
mise.tie.haar.sp.1.v2=mise(mu.est.tie.haar.sp.1.v2,mu.sp)
mise.tie.haar.sp.1.v3=mise(mu.est.tie.haar.sp.1.v3,mu.sp)
mise.tie.haar.sp.1.v4=mise(mu.est.tie.haar.sp.1.v4,mu.sp)
mise.tie.haar.sp.1.v5=mise(mu.est.tie.haar.sp.1.v5,mu.sp)
mise.tie.haar.sp.3.v1=mise(mu.est.tie.haar.sp.3.v1,mu.sp)
mise.tie.haar.sp.3.v2=mise(mu.est.tie.haar.sp.3.v2,mu.sp)
mise.tie.haar.sp.3.v3=mise(mu.est.tie.haar.sp.3.v3,mu.sp)
mise.tie.haar.sp.3.v4=mise(mu.est.tie.haar.sp.3.v4,mu.sp)
mise.tie.haar.sp.3.v5=mise(mu.est.tie.haar.sp.3.v5,mu.sp)

mise.tie.haar.bump.1.v1=mise(mu.est.tie.haar.bump.1.v1,mu.bump)
mise.tie.haar.bump.1.v2=mise(mu.est.tie.haar.bump.1.v2,mu.bump)
mise.tie.haar.bump.1.v3=mise(mu.est.tie.haar.bump.1.v3,mu.bump)
mise.tie.haar.bump.1.v4=mise(mu.est.tie.haar.bump.1.v4,mu.bump)
mise.tie.haar.bump.1.v5=mise(mu.est.tie.haar.bump.1.v5,mu.bump)
mise.tie.haar.bump.3.v1=mise(mu.est.tie.haar.bump.3.v1,mu.bump)
mise.tie.haar.bump.3.v2=mise(mu.est.tie.haar.bump.3.v2,mu.bump)
mise.tie.haar.bump.3.v3=mise(mu.est.tie.haar.bump.3.v3,mu.bump)
mise.tie.haar.bump.3.v4=mise(mu.est.tie.haar.bump.3.v4,mu.bump)
mise.tie.haar.bump.3.v5=mise(mu.est.tie.haar.bump.3.v5,mu.bump)

mise.tie.haar.blk.1.v1=mise(mu.est.tie.haar.blk.1.v1,mu.blk)
mise.tie.haar.blk.1.v2=mise(mu.est.tie.haar.blk.1.v2,mu.blk)
mise.tie.haar.blk.1.v3=mise(mu.est.tie.haar.blk.1.v3,mu.blk)
mise.tie.haar.blk.1.v4=mise(mu.est.tie.haar.blk.1.v4,mu.blk)
mise.tie.haar.blk.1.v5=mise(mu.est.tie.haar.blk.1.v5,mu.blk)
mise.tie.haar.blk.3.v1=mise(mu.est.tie.haar.blk.3.v1,mu.blk)
mise.tie.haar.blk.3.v2=mise(mu.est.tie.haar.blk.3.v2,mu.blk)
mise.tie.haar.blk.3.v3=mise(mu.est.tie.haar.blk.3.v3,mu.blk)
mise.tie.haar.blk.3.v4=mise(mu.est.tie.haar.blk.3.v4,mu.blk)
mise.tie.haar.blk.3.v5=mise(mu.est.tie.haar.blk.3.v5,mu.blk)

mise.tie.haar.ang.1.v1=mise(mu.est.tie.haar.ang.1.v1,mu.ang)
mise.tie.haar.ang.1.v2=mise(mu.est.tie.haar.ang.1.v2,mu.ang)
mise.tie.haar.ang.1.v3=mise(mu.est.tie.haar.ang.1.v3,mu.ang)
mise.tie.haar.ang.1.v4=mise(mu.est.tie.haar.ang.1.v4,mu.ang)
mise.tie.haar.ang.1.v5=mise(mu.est.tie.haar.ang.1.v5,mu.ang)
mise.tie.haar.ang.3.v1=mise(mu.est.tie.haar.ang.3.v1,mu.ang)
mise.tie.haar.ang.3.v2=mise(mu.est.tie.haar.ang.3.v2,mu.ang)
mise.tie.haar.ang.3.v3=mise(mu.est.tie.haar.ang.3.v3,mu.ang)
mise.tie.haar.ang.3.v4=mise(mu.est.tie.haar.ang.3.v4,mu.ang)
mise.tie.haar.ang.3.v5=mise(mu.est.tie.haar.ang.3.v5,mu.ang)

mise.tie.haar.dop.1.v1=mise(mu.est.tie.haar.dop.1.v1,mu.dop)
mise.tie.haar.dop.1.v2=mise(mu.est.tie.haar.dop.1.v2,mu.dop)
mise.tie.haar.dop.1.v3=mise(mu.est.tie.haar.dop.1.v3,mu.dop)
mise.tie.haar.dop.1.v4=mise(mu.est.tie.haar.dop.1.v4,mu.dop)
mise.tie.haar.dop.1.v5=mise(mu.est.tie.haar.dop.1.v5,mu.dop)
mise.tie.haar.dop.3.v1=mise(mu.est.tie.haar.dop.3.v1,mu.dop)
mise.tie.haar.dop.3.v2=mise(mu.est.tie.haar.dop.3.v2,mu.dop)
mise.tie.haar.dop.3.v3=mise(mu.est.tie.haar.dop.3.v3,mu.dop)
mise.tie.haar.dop.3.v4=mise(mu.est.tie.haar.dop.3.v4,mu.dop)
mise.tie.haar.dop.3.v5=mise(mu.est.tie.haar.dop.3.v5,mu.dop)

mise.tie.haar.blip.1.v1=mise(mu.est.tie.haar.blip.1.v1,mu.blip)
mise.tie.haar.blip.1.v2=mise(mu.est.tie.haar.blip.1.v2,mu.blip)
mise.tie.haar.blip.1.v3=mise(mu.est.tie.haar.blip.1.v3,mu.blip)
mise.tie.haar.blip.1.v4=mise(mu.est.tie.haar.blip.1.v4,mu.blip)
mise.tie.haar.blip.1.v5=mise(mu.est.tie.haar.blip.1.v5,mu.blip)
mise.tie.haar.blip.3.v1=mise(mu.est.tie.haar.blip.3.v1,mu.blip)
mise.tie.haar.blip.3.v2=mise(mu.est.tie.haar.blip.3.v2,mu.blip)
mise.tie.haar.blip.3.v3=mise(mu.est.tie.haar.blip.3.v3,mu.blip)
mise.tie.haar.blip.3.v4=mise(mu.est.tie.haar.blip.3.v4,mu.blip)
mise.tie.haar.blip.3.v5=mise(mu.est.tie.haar.blip.3.v5,mu.blip)

mise.tie.haar.cor.1.v1=mise(mu.est.tie.haar.cor.1.v1,mu.cor)
mise.tie.haar.cor.1.v2=mise(mu.est.tie.haar.cor.1.v2,mu.cor)
mise.tie.haar.cor.1.v3=mise(mu.est.tie.haar.cor.1.v3,mu.cor)
mise.tie.haar.cor.1.v4=mise(mu.est.tie.haar.cor.1.v4,mu.cor)
mise.tie.haar.cor.1.v5=mise(mu.est.tie.haar.cor.1.v5,mu.cor)
mise.tie.haar.cor.3.v1=mise(mu.est.tie.haar.cor.3.v1,mu.cor)
mise.tie.haar.cor.3.v2=mise(mu.est.tie.haar.cor.3.v2,mu.cor)
mise.tie.haar.cor.3.v3=mise(mu.est.tie.haar.cor.3.v3,mu.cor)
mise.tie.haar.cor.3.v4=mise(mu.est.tie.haar.cor.3.v4,mu.cor)
mise.tie.haar.cor.3.v5=mise(mu.est.tie.haar.cor.3.v5,mu.cor)


mise.tie.s8.sp.1.v1=mise(mu.est.tie.s8.sp.1.v1,mu.sp)
mise.tie.s8.sp.1.v2=mise(mu.est.tie.s8.sp.1.v2,mu.sp)
mise.tie.s8.sp.1.v3=mise(mu.est.tie.s8.sp.1.v3,mu.sp)
mise.tie.s8.sp.1.v4=mise(mu.est.tie.s8.sp.1.v4,mu.sp)
mise.tie.s8.sp.1.v5=mise(mu.est.tie.s8.sp.1.v5,mu.sp)
mise.tie.s8.sp.3.v1=mise(mu.est.tie.s8.sp.3.v1,mu.sp)
mise.tie.s8.sp.3.v2=mise(mu.est.tie.s8.sp.3.v2,mu.sp)
mise.tie.s8.sp.3.v3=mise(mu.est.tie.s8.sp.3.v3,mu.sp)
mise.tie.s8.sp.3.v4=mise(mu.est.tie.s8.sp.3.v4,mu.sp)
mise.tie.s8.sp.3.v5=mise(mu.est.tie.s8.sp.3.v5,mu.sp)

mise.tie.s8.bump.1.v1=mise(mu.est.tie.s8.bump.1.v1,mu.bump)
mise.tie.s8.bump.1.v2=mise(mu.est.tie.s8.bump.1.v2,mu.bump)
mise.tie.s8.bump.1.v3=mise(mu.est.tie.s8.bump.1.v3,mu.bump)
mise.tie.s8.bump.1.v4=mise(mu.est.tie.s8.bump.1.v4,mu.bump)
mise.tie.s8.bump.1.v5=mise(mu.est.tie.s8.bump.1.v5,mu.bump)
mise.tie.s8.bump.3.v1=mise(mu.est.tie.s8.bump.3.v1,mu.bump)
mise.tie.s8.bump.3.v2=mise(mu.est.tie.s8.bump.3.v2,mu.bump)
mise.tie.s8.bump.3.v3=mise(mu.est.tie.s8.bump.3.v3,mu.bump)
mise.tie.s8.bump.3.v4=mise(mu.est.tie.s8.bump.3.v4,mu.bump)
mise.tie.s8.bump.3.v5=mise(mu.est.tie.s8.bump.3.v5,mu.bump)

mise.tie.s8.blk.1.v1=mise(mu.est.tie.s8.blk.1.v1,mu.blk)
mise.tie.s8.blk.1.v2=mise(mu.est.tie.s8.blk.1.v2,mu.blk)
mise.tie.s8.blk.1.v3=mise(mu.est.tie.s8.blk.1.v3,mu.blk)
mise.tie.s8.blk.1.v4=mise(mu.est.tie.s8.blk.1.v4,mu.blk)
mise.tie.s8.blk.1.v5=mise(mu.est.tie.s8.blk.1.v5,mu.blk)
mise.tie.s8.blk.3.v1=mise(mu.est.tie.s8.blk.3.v1,mu.blk)
mise.tie.s8.blk.3.v2=mise(mu.est.tie.s8.blk.3.v2,mu.blk)
mise.tie.s8.blk.3.v3=mise(mu.est.tie.s8.blk.3.v3,mu.blk)
mise.tie.s8.blk.3.v4=mise(mu.est.tie.s8.blk.3.v4,mu.blk)
mise.tie.s8.blk.3.v5=mise(mu.est.tie.s8.blk.3.v5,mu.blk)

mise.tie.s8.ang.1.v1=mise(mu.est.tie.s8.ang.1.v1,mu.ang)
mise.tie.s8.ang.1.v2=mise(mu.est.tie.s8.ang.1.v2,mu.ang)
mise.tie.s8.ang.1.v3=mise(mu.est.tie.s8.ang.1.v3,mu.ang)
mise.tie.s8.ang.1.v4=mise(mu.est.tie.s8.ang.1.v4,mu.ang)
mise.tie.s8.ang.1.v5=mise(mu.est.tie.s8.ang.1.v5,mu.ang)
mise.tie.s8.ang.3.v1=mise(mu.est.tie.s8.ang.3.v1,mu.ang)
mise.tie.s8.ang.3.v2=mise(mu.est.tie.s8.ang.3.v2,mu.ang)
mise.tie.s8.ang.3.v3=mise(mu.est.tie.s8.ang.3.v3,mu.ang)
mise.tie.s8.ang.3.v4=mise(mu.est.tie.s8.ang.3.v4,mu.ang)
mise.tie.s8.ang.3.v5=mise(mu.est.tie.s8.ang.3.v5,mu.ang)

mise.tie.s8.dop.1.v1=mise(mu.est.tie.s8.dop.1.v1,mu.dop)
mise.tie.s8.dop.1.v2=mise(mu.est.tie.s8.dop.1.v2,mu.dop)
mise.tie.s8.dop.1.v3=mise(mu.est.tie.s8.dop.1.v3,mu.dop)
mise.tie.s8.dop.1.v4=mise(mu.est.tie.s8.dop.1.v4,mu.dop)
mise.tie.s8.dop.1.v5=mise(mu.est.tie.s8.dop.1.v5,mu.dop)
mise.tie.s8.dop.3.v1=mise(mu.est.tie.s8.dop.3.v1,mu.dop)
mise.tie.s8.dop.3.v2=mise(mu.est.tie.s8.dop.3.v2,mu.dop)
mise.tie.s8.dop.3.v3=mise(mu.est.tie.s8.dop.3.v3,mu.dop)
mise.tie.s8.dop.3.v4=mise(mu.est.tie.s8.dop.3.v4,mu.dop)
mise.tie.s8.dop.3.v5=mise(mu.est.tie.s8.dop.3.v5,mu.dop)

mise.tie.s8.blip.1.v1=mise(mu.est.tie.s8.blip.1.v1,mu.blip)
mise.tie.s8.blip.1.v2=mise(mu.est.tie.s8.blip.1.v2,mu.blip)
mise.tie.s8.blip.1.v3=mise(mu.est.tie.s8.blip.1.v3,mu.blip)
mise.tie.s8.blip.1.v4=mise(mu.est.tie.s8.blip.1.v4,mu.blip)
mise.tie.s8.blip.1.v5=mise(mu.est.tie.s8.blip.1.v5,mu.blip)
mise.tie.s8.blip.3.v1=mise(mu.est.tie.s8.blip.3.v1,mu.blip)
mise.tie.s8.blip.3.v2=mise(mu.est.tie.s8.blip.3.v2,mu.blip)
mise.tie.s8.blip.3.v3=mise(mu.est.tie.s8.blip.3.v3,mu.blip)
mise.tie.s8.blip.3.v4=mise(mu.est.tie.s8.blip.3.v4,mu.blip)
mise.tie.s8.blip.3.v5=mise(mu.est.tie.s8.blip.3.v5,mu.blip)

mise.tie.s8.cor.1.v1=mise(mu.est.tie.s8.cor.1.v1,mu.cor)
mise.tie.s8.cor.1.v2=mise(mu.est.tie.s8.cor.1.v2,mu.cor)
mise.tie.s8.cor.1.v3=mise(mu.est.tie.s8.cor.1.v3,mu.cor)
mise.tie.s8.cor.1.v4=mise(mu.est.tie.s8.cor.1.v4,mu.cor)
mise.tie.s8.cor.1.v5=mise(mu.est.tie.s8.cor.1.v5,mu.cor)
mise.tie.s8.cor.3.v1=mise(mu.est.tie.s8.cor.3.v1,mu.cor)
mise.tie.s8.cor.3.v2=mise(mu.est.tie.s8.cor.3.v2,mu.cor)
mise.tie.s8.cor.3.v3=mise(mu.est.tie.s8.cor.3.v3,mu.cor)
mise.tie.s8.cor.3.v4=mise(mu.est.tie.s8.cor.3.v4,mu.cor)
mise.tie.s8.cor.3.v5=mise(mu.est.tie.s8.cor.3.v5,mu.cor)



mise.asht.haar.sp.1.v1=mise(mu.est.asht.haar.sp.1.v1,mu.sp)
mise.asht.haar.sp.1.v2=mise(mu.est.asht.haar.sp.1.v2,mu.sp)
mise.asht.haar.sp.1.v3=mise(mu.est.asht.haar.sp.1.v3,mu.sp)
mise.asht.haar.sp.1.v4=mise(mu.est.asht.haar.sp.1.v4,mu.sp)
mise.asht.haar.sp.1.v5=mise(mu.est.asht.haar.sp.1.v5,mu.sp)
mise.asht.haar.sp.3.v1=mise(mu.est.asht.haar.sp.3.v1,mu.sp)
mise.asht.haar.sp.3.v2=mise(mu.est.asht.haar.sp.3.v2,mu.sp)
mise.asht.haar.sp.3.v3=mise(mu.est.asht.haar.sp.3.v3,mu.sp)
mise.asht.haar.sp.3.v4=mise(mu.est.asht.haar.sp.3.v4,mu.sp)
mise.asht.haar.sp.3.v5=mise(mu.est.asht.haar.sp.3.v5,mu.sp)

mise.asht.haar.bump.1.v1=mise(mu.est.asht.haar.bump.1.v1,mu.bump)
mise.asht.haar.bump.1.v2=mise(mu.est.asht.haar.bump.1.v2,mu.bump)
mise.asht.haar.bump.1.v3=mise(mu.est.asht.haar.bump.1.v3,mu.bump)
mise.asht.haar.bump.1.v4=mise(mu.est.asht.haar.bump.1.v4,mu.bump)
mise.asht.haar.bump.1.v5=mise(mu.est.asht.haar.bump.1.v5,mu.bump)
mise.asht.haar.bump.3.v1=mise(mu.est.asht.haar.bump.3.v1,mu.bump)
mise.asht.haar.bump.3.v2=mise(mu.est.asht.haar.bump.3.v2,mu.bump)
mise.asht.haar.bump.3.v3=mise(mu.est.asht.haar.bump.3.v3,mu.bump)
mise.asht.haar.bump.3.v4=mise(mu.est.asht.haar.bump.3.v4,mu.bump)
mise.asht.haar.bump.3.v5=mise(mu.est.asht.haar.bump.3.v5,mu.bump)

mise.asht.haar.blk.1.v1=mise(mu.est.asht.haar.blk.1.v1,mu.blk)
mise.asht.haar.blk.1.v2=mise(mu.est.asht.haar.blk.1.v2,mu.blk)
mise.asht.haar.blk.1.v3=mise(mu.est.asht.haar.blk.1.v3,mu.blk)
mise.asht.haar.blk.1.v4=mise(mu.est.asht.haar.blk.1.v4,mu.blk)
mise.asht.haar.blk.1.v5=mise(mu.est.asht.haar.blk.1.v5,mu.blk)
mise.asht.haar.blk.3.v1=mise(mu.est.asht.haar.blk.3.v1,mu.blk)
mise.asht.haar.blk.3.v2=mise(mu.est.asht.haar.blk.3.v2,mu.blk)
mise.asht.haar.blk.3.v3=mise(mu.est.asht.haar.blk.3.v3,mu.blk)
mise.asht.haar.blk.3.v4=mise(mu.est.asht.haar.blk.3.v4,mu.blk)
mise.asht.haar.blk.3.v5=mise(mu.est.asht.haar.blk.3.v5,mu.blk)

mise.asht.haar.ang.1.v1=mise(mu.est.asht.haar.ang.1.v1,mu.ang)
mise.asht.haar.ang.1.v2=mise(mu.est.asht.haar.ang.1.v2,mu.ang)
mise.asht.haar.ang.1.v3=mise(mu.est.asht.haar.ang.1.v3,mu.ang)
mise.asht.haar.ang.1.v4=mise(mu.est.asht.haar.ang.1.v4,mu.ang)
mise.asht.haar.ang.1.v5=mise(mu.est.asht.haar.ang.1.v5,mu.ang)
mise.asht.haar.ang.3.v1=mise(mu.est.asht.haar.ang.3.v1,mu.ang)
mise.asht.haar.ang.3.v2=mise(mu.est.asht.haar.ang.3.v2,mu.ang)
mise.asht.haar.ang.3.v3=mise(mu.est.asht.haar.ang.3.v3,mu.ang)
mise.asht.haar.ang.3.v4=mise(mu.est.asht.haar.ang.3.v4,mu.ang)
mise.asht.haar.ang.3.v5=mise(mu.est.asht.haar.ang.3.v5,mu.ang)

mise.asht.haar.dop.1.v1=mise(mu.est.asht.haar.dop.1.v1,mu.dop)
mise.asht.haar.dop.1.v2=mise(mu.est.asht.haar.dop.1.v2,mu.dop)
mise.asht.haar.dop.1.v3=mise(mu.est.asht.haar.dop.1.v3,mu.dop)
mise.asht.haar.dop.1.v4=mise(mu.est.asht.haar.dop.1.v4,mu.dop)
mise.asht.haar.dop.1.v5=mise(mu.est.asht.haar.dop.1.v5,mu.dop)
mise.asht.haar.dop.3.v1=mise(mu.est.asht.haar.dop.3.v1,mu.dop)
mise.asht.haar.dop.3.v2=mise(mu.est.asht.haar.dop.3.v2,mu.dop)
mise.asht.haar.dop.3.v3=mise(mu.est.asht.haar.dop.3.v3,mu.dop)
mise.asht.haar.dop.3.v4=mise(mu.est.asht.haar.dop.3.v4,mu.dop)
mise.asht.haar.dop.3.v5=mise(mu.est.asht.haar.dop.3.v5,mu.dop)

mise.asht.haar.blip.1.v1=mise(mu.est.asht.haar.blip.1.v1,mu.blip)
mise.asht.haar.blip.1.v2=mise(mu.est.asht.haar.blip.1.v2,mu.blip)
mise.asht.haar.blip.1.v3=mise(mu.est.asht.haar.blip.1.v3,mu.blip)
mise.asht.haar.blip.1.v4=mise(mu.est.asht.haar.blip.1.v4,mu.blip)
mise.asht.haar.blip.1.v5=mise(mu.est.asht.haar.blip.1.v5,mu.blip)
mise.asht.haar.blip.3.v1=mise(mu.est.asht.haar.blip.3.v1,mu.blip)
mise.asht.haar.blip.3.v2=mise(mu.est.asht.haar.blip.3.v2,mu.blip)
mise.asht.haar.blip.3.v3=mise(mu.est.asht.haar.blip.3.v3,mu.blip)
mise.asht.haar.blip.3.v4=mise(mu.est.asht.haar.blip.3.v4,mu.blip)
mise.asht.haar.blip.3.v5=mise(mu.est.asht.haar.blip.3.v5,mu.blip)

mise.asht.haar.cor.1.v1=mise(mu.est.asht.haar.cor.1.v1,mu.cor)
mise.asht.haar.cor.1.v2=mise(mu.est.asht.haar.cor.1.v2,mu.cor)
mise.asht.haar.cor.1.v3=mise(mu.est.asht.haar.cor.1.v3,mu.cor)
mise.asht.haar.cor.1.v4=mise(mu.est.asht.haar.cor.1.v4,mu.cor)
mise.asht.haar.cor.1.v5=mise(mu.est.asht.haar.cor.1.v5,mu.cor)
mise.asht.haar.cor.3.v1=mise(mu.est.asht.haar.cor.3.v1,mu.cor)
mise.asht.haar.cor.3.v2=mise(mu.est.asht.haar.cor.3.v2,mu.cor)
mise.asht.haar.cor.3.v3=mise(mu.est.asht.haar.cor.3.v3,mu.cor)
mise.asht.haar.cor.3.v4=mise(mu.est.asht.haar.cor.3.v4,mu.cor)
mise.asht.haar.cor.3.v5=mise(mu.est.asht.haar.cor.3.v5,mu.cor)


mise.asht.s8.sp.1.v1=mise(mu.est.asht.s8.sp.1.v1,mu.sp)
mise.asht.s8.sp.1.v2=mise(mu.est.asht.s8.sp.1.v2,mu.sp)
mise.asht.s8.sp.1.v3=mise(mu.est.asht.s8.sp.1.v3,mu.sp)
mise.asht.s8.sp.1.v4=mise(mu.est.asht.s8.sp.1.v4,mu.sp)
mise.asht.s8.sp.1.v5=mise(mu.est.asht.s8.sp.1.v5,mu.sp)
mise.asht.s8.sp.3.v1=mise(mu.est.asht.s8.sp.3.v1,mu.sp)
mise.asht.s8.sp.3.v2=mise(mu.est.asht.s8.sp.3.v2,mu.sp)
mise.asht.s8.sp.3.v3=mise(mu.est.asht.s8.sp.3.v3,mu.sp)
mise.asht.s8.sp.3.v4=mise(mu.est.asht.s8.sp.3.v4,mu.sp)
mise.asht.s8.sp.3.v5=mise(mu.est.asht.s8.sp.3.v5,mu.sp)

mise.asht.s8.bump.1.v1=mise(mu.est.asht.s8.bump.1.v1,mu.bump)
mise.asht.s8.bump.1.v2=mise(mu.est.asht.s8.bump.1.v2,mu.bump)
mise.asht.s8.bump.1.v3=mise(mu.est.asht.s8.bump.1.v3,mu.bump)
mise.asht.s8.bump.1.v4=mise(mu.est.asht.s8.bump.1.v4,mu.bump)
mise.asht.s8.bump.1.v5=mise(mu.est.asht.s8.bump.1.v5,mu.bump)
mise.asht.s8.bump.3.v1=mise(mu.est.asht.s8.bump.3.v1,mu.bump)
mise.asht.s8.bump.3.v2=mise(mu.est.asht.s8.bump.3.v2,mu.bump)
mise.asht.s8.bump.3.v3=mise(mu.est.asht.s8.bump.3.v3,mu.bump)
mise.asht.s8.bump.3.v4=mise(mu.est.asht.s8.bump.3.v4,mu.bump)
mise.asht.s8.bump.3.v5=mise(mu.est.asht.s8.bump.3.v5,mu.bump)

mise.asht.s8.blk.1.v1=mise(mu.est.asht.s8.blk.1.v1,mu.blk)
mise.asht.s8.blk.1.v2=mise(mu.est.asht.s8.blk.1.v2,mu.blk)
mise.asht.s8.blk.1.v3=mise(mu.est.asht.s8.blk.1.v3,mu.blk)
mise.asht.s8.blk.1.v4=mise(mu.est.asht.s8.blk.1.v4,mu.blk)
mise.asht.s8.blk.1.v5=mise(mu.est.asht.s8.blk.1.v5,mu.blk)
mise.asht.s8.blk.3.v1=mise(mu.est.asht.s8.blk.3.v1,mu.blk)
mise.asht.s8.blk.3.v2=mise(mu.est.asht.s8.blk.3.v2,mu.blk)
mise.asht.s8.blk.3.v3=mise(mu.est.asht.s8.blk.3.v3,mu.blk)
mise.asht.s8.blk.3.v4=mise(mu.est.asht.s8.blk.3.v4,mu.blk)
mise.asht.s8.blk.3.v5=mise(mu.est.asht.s8.blk.3.v5,mu.blk)

mise.asht.s8.ang.1.v1=mise(mu.est.asht.s8.ang.1.v1,mu.ang)
mise.asht.s8.ang.1.v2=mise(mu.est.asht.s8.ang.1.v2,mu.ang)
mise.asht.s8.ang.1.v3=mise(mu.est.asht.s8.ang.1.v3,mu.ang)
mise.asht.s8.ang.1.v4=mise(mu.est.asht.s8.ang.1.v4,mu.ang)
mise.asht.s8.ang.1.v5=mise(mu.est.asht.s8.ang.1.v5,mu.ang)
mise.asht.s8.ang.3.v1=mise(mu.est.asht.s8.ang.3.v1,mu.ang)
mise.asht.s8.ang.3.v2=mise(mu.est.asht.s8.ang.3.v2,mu.ang)
mise.asht.s8.ang.3.v3=mise(mu.est.asht.s8.ang.3.v3,mu.ang)
mise.asht.s8.ang.3.v4=mise(mu.est.asht.s8.ang.3.v4,mu.ang)
mise.asht.s8.ang.3.v5=mise(mu.est.asht.s8.ang.3.v5,mu.ang)

mise.asht.s8.dop.1.v1=mise(mu.est.asht.s8.dop.1.v1,mu.dop)
mise.asht.s8.dop.1.v2=mise(mu.est.asht.s8.dop.1.v2,mu.dop)
mise.asht.s8.dop.1.v3=mise(mu.est.asht.s8.dop.1.v3,mu.dop)
mise.asht.s8.dop.1.v4=mise(mu.est.asht.s8.dop.1.v4,mu.dop)
mise.asht.s8.dop.1.v5=mise(mu.est.asht.s8.dop.1.v5,mu.dop)
mise.asht.s8.dop.3.v1=mise(mu.est.asht.s8.dop.3.v1,mu.dop)
mise.asht.s8.dop.3.v2=mise(mu.est.asht.s8.dop.3.v2,mu.dop)
mise.asht.s8.dop.3.v3=mise(mu.est.asht.s8.dop.3.v3,mu.dop)
mise.asht.s8.dop.3.v4=mise(mu.est.asht.s8.dop.3.v4,mu.dop)
mise.asht.s8.dop.3.v5=mise(mu.est.asht.s8.dop.3.v5,mu.dop)

mise.asht.s8.blip.1.v1=mise(mu.est.asht.s8.blip.1.v1,mu.blip)
mise.asht.s8.blip.1.v2=mise(mu.est.asht.s8.blip.1.v2,mu.blip)
mise.asht.s8.blip.1.v3=mise(mu.est.asht.s8.blip.1.v3,mu.blip)
mise.asht.s8.blip.1.v4=mise(mu.est.asht.s8.blip.1.v4,mu.blip)
mise.asht.s8.blip.1.v5=mise(mu.est.asht.s8.blip.1.v5,mu.blip)
mise.asht.s8.blip.3.v1=mise(mu.est.asht.s8.blip.3.v1,mu.blip)
mise.asht.s8.blip.3.v2=mise(mu.est.asht.s8.blip.3.v2,mu.blip)
mise.asht.s8.blip.3.v3=mise(mu.est.asht.s8.blip.3.v3,mu.blip)
mise.asht.s8.blip.3.v4=mise(mu.est.asht.s8.blip.3.v4,mu.blip)
mise.asht.s8.blip.3.v5=mise(mu.est.asht.s8.blip.3.v5,mu.blip)

mise.asht.s8.cor.1.v1=mise(mu.est.asht.s8.cor.1.v1,mu.cor)
mise.asht.s8.cor.1.v2=mise(mu.est.asht.s8.cor.1.v2,mu.cor)
mise.asht.s8.cor.1.v3=mise(mu.est.asht.s8.cor.1.v3,mu.cor)
mise.asht.s8.cor.1.v4=mise(mu.est.asht.s8.cor.1.v4,mu.cor)
mise.asht.s8.cor.1.v5=mise(mu.est.asht.s8.cor.1.v5,mu.cor)
mise.asht.s8.cor.3.v1=mise(mu.est.asht.s8.cor.3.v1,mu.cor)
mise.asht.s8.cor.3.v2=mise(mu.est.asht.s8.cor.3.v2,mu.cor)
mise.asht.s8.cor.3.v3=mise(mu.est.asht.s8.cor.3.v3,mu.cor)
mise.asht.s8.cor.3.v4=mise(mu.est.asht.s8.cor.3.v4,mu.cor)
mise.asht.s8.cor.3.v5=mise(mu.est.asht.s8.cor.3.v5,mu.cor)



mise.tit.haar.sp.1.v1=mise(mu.est.tit.haar.sp.1.v1,mu.sp)
mise.tit.haar.sp.1.v2=mise(mu.est.tit.haar.sp.1.v2,mu.sp)
mise.tit.haar.sp.1.v3=mise(mu.est.tit.haar.sp.1.v3,mu.sp)
mise.tit.haar.sp.1.v4=mise(mu.est.tit.haar.sp.1.v4,mu.sp)
mise.tit.haar.sp.1.v5=mise(mu.est.tit.haar.sp.1.v5,mu.sp)
mise.tit.haar.sp.3.v1=mise(mu.est.tit.haar.sp.3.v1,mu.sp)
mise.tit.haar.sp.3.v2=mise(mu.est.tit.haar.sp.3.v2,mu.sp)
mise.tit.haar.sp.3.v3=mise(mu.est.tit.haar.sp.3.v3,mu.sp)
mise.tit.haar.sp.3.v4=mise(mu.est.tit.haar.sp.3.v4,mu.sp)
mise.tit.haar.sp.3.v5=mise(mu.est.tit.haar.sp.3.v5,mu.sp)

mise.tit.haar.bump.1.v1=mise(mu.est.tit.haar.bump.1.v1,mu.bump)
mise.tit.haar.bump.1.v2=mise(mu.est.tit.haar.bump.1.v2,mu.bump)
mise.tit.haar.bump.1.v3=mise(mu.est.tit.haar.bump.1.v3,mu.bump)
mise.tit.haar.bump.1.v4=mise(mu.est.tit.haar.bump.1.v4,mu.bump)
mise.tit.haar.bump.1.v5=mise(mu.est.tit.haar.bump.1.v5,mu.bump)
mise.tit.haar.bump.3.v1=mise(mu.est.tit.haar.bump.3.v1,mu.bump)
mise.tit.haar.bump.3.v2=mise(mu.est.tit.haar.bump.3.v2,mu.bump)
mise.tit.haar.bump.3.v3=mise(mu.est.tit.haar.bump.3.v3,mu.bump)
mise.tit.haar.bump.3.v4=mise(mu.est.tit.haar.bump.3.v4,mu.bump)
mise.tit.haar.bump.3.v5=mise(mu.est.tit.haar.bump.3.v5,mu.bump)

mise.tit.haar.blk.1.v1=mise(mu.est.tit.haar.blk.1.v1,mu.blk)
mise.tit.haar.blk.1.v2=mise(mu.est.tit.haar.blk.1.v2,mu.blk)
mise.tit.haar.blk.1.v3=mise(mu.est.tit.haar.blk.1.v3,mu.blk)
mise.tit.haar.blk.1.v4=mise(mu.est.tit.haar.blk.1.v4,mu.blk)
mise.tit.haar.blk.1.v5=mise(mu.est.tit.haar.blk.1.v5,mu.blk)
mise.tit.haar.blk.3.v1=mise(mu.est.tit.haar.blk.3.v1,mu.blk)
mise.tit.haar.blk.3.v2=mise(mu.est.tit.haar.blk.3.v2,mu.blk)
mise.tit.haar.blk.3.v3=mise(mu.est.tit.haar.blk.3.v3,mu.blk)
mise.tit.haar.blk.3.v4=mise(mu.est.tit.haar.blk.3.v4,mu.blk)
mise.tit.haar.blk.3.v5=mise(mu.est.tit.haar.blk.3.v5,mu.blk)

mise.tit.haar.ang.1.v1=mise(mu.est.tit.haar.ang.1.v1,mu.ang)
mise.tit.haar.ang.1.v2=mise(mu.est.tit.haar.ang.1.v2,mu.ang)
mise.tit.haar.ang.1.v3=mise(mu.est.tit.haar.ang.1.v3,mu.ang)
mise.tit.haar.ang.1.v4=mise(mu.est.tit.haar.ang.1.v4,mu.ang)
mise.tit.haar.ang.1.v5=mise(mu.est.tit.haar.ang.1.v5,mu.ang)
mise.tit.haar.ang.3.v1=mise(mu.est.tit.haar.ang.3.v1,mu.ang)
mise.tit.haar.ang.3.v2=mise(mu.est.tit.haar.ang.3.v2,mu.ang)
mise.tit.haar.ang.3.v3=mise(mu.est.tit.haar.ang.3.v3,mu.ang)
mise.tit.haar.ang.3.v4=mise(mu.est.tit.haar.ang.3.v4,mu.ang)
mise.tit.haar.ang.3.v5=mise(mu.est.tit.haar.ang.3.v5,mu.ang)

mise.tit.haar.dop.1.v1=mise(mu.est.tit.haar.dop.1.v1,mu.dop)
mise.tit.haar.dop.1.v2=mise(mu.est.tit.haar.dop.1.v2,mu.dop)
mise.tit.haar.dop.1.v3=mise(mu.est.tit.haar.dop.1.v3,mu.dop)
mise.tit.haar.dop.1.v4=mise(mu.est.tit.haar.dop.1.v4,mu.dop)
mise.tit.haar.dop.1.v5=mise(mu.est.tit.haar.dop.1.v5,mu.dop)
mise.tit.haar.dop.3.v1=mise(mu.est.tit.haar.dop.3.v1,mu.dop)
mise.tit.haar.dop.3.v2=mise(mu.est.tit.haar.dop.3.v2,mu.dop)
mise.tit.haar.dop.3.v3=mise(mu.est.tit.haar.dop.3.v3,mu.dop)
mise.tit.haar.dop.3.v4=mise(mu.est.tit.haar.dop.3.v4,mu.dop)
mise.tit.haar.dop.3.v5=mise(mu.est.tit.haar.dop.3.v5,mu.dop)

mise.tit.haar.blip.1.v1=mise(mu.est.tit.haar.blip.1.v1,mu.blip)
mise.tit.haar.blip.1.v2=mise(mu.est.tit.haar.blip.1.v2,mu.blip)
mise.tit.haar.blip.1.v3=mise(mu.est.tit.haar.blip.1.v3,mu.blip)
mise.tit.haar.blip.1.v4=mise(mu.est.tit.haar.blip.1.v4,mu.blip)
mise.tit.haar.blip.1.v5=mise(mu.est.tit.haar.blip.1.v5,mu.blip)
mise.tit.haar.blip.3.v1=mise(mu.est.tit.haar.blip.3.v1,mu.blip)
mise.tit.haar.blip.3.v2=mise(mu.est.tit.haar.blip.3.v2,mu.blip)
mise.tit.haar.blip.3.v3=mise(mu.est.tit.haar.blip.3.v3,mu.blip)
mise.tit.haar.blip.3.v4=mise(mu.est.tit.haar.blip.3.v4,mu.blip)
mise.tit.haar.blip.3.v5=mise(mu.est.tit.haar.blip.3.v5,mu.blip)

mise.tit.haar.cor.1.v1=mise(mu.est.tit.haar.cor.1.v1,mu.cor)
mise.tit.haar.cor.1.v2=mise(mu.est.tit.haar.cor.1.v2,mu.cor)
mise.tit.haar.cor.1.v3=mise(mu.est.tit.haar.cor.1.v3,mu.cor)
mise.tit.haar.cor.1.v4=mise(mu.est.tit.haar.cor.1.v4,mu.cor)
mise.tit.haar.cor.1.v5=mise(mu.est.tit.haar.cor.1.v5,mu.cor)
mise.tit.haar.cor.3.v1=mise(mu.est.tit.haar.cor.3.v1,mu.cor)
mise.tit.haar.cor.3.v2=mise(mu.est.tit.haar.cor.3.v2,mu.cor)
mise.tit.haar.cor.3.v3=mise(mu.est.tit.haar.cor.3.v3,mu.cor)
mise.tit.haar.cor.3.v4=mise(mu.est.tit.haar.cor.3.v4,mu.cor)
mise.tit.haar.cor.3.v5=mise(mu.est.tit.haar.cor.3.v5,mu.cor)


mise.tit.s8.sp.1.v1=mise(mu.est.tit.s8.sp.1.v1,mu.sp)
mise.tit.s8.sp.1.v2=mise(mu.est.tit.s8.sp.1.v2,mu.sp)
mise.tit.s8.sp.1.v3=mise(mu.est.tit.s8.sp.1.v3,mu.sp)
mise.tit.s8.sp.1.v4=mise(mu.est.tit.s8.sp.1.v4,mu.sp)
mise.tit.s8.sp.1.v5=mise(mu.est.tit.s8.sp.1.v5,mu.sp)
mise.tit.s8.sp.3.v1=mise(mu.est.tit.s8.sp.3.v1,mu.sp)
mise.tit.s8.sp.3.v2=mise(mu.est.tit.s8.sp.3.v2,mu.sp)
mise.tit.s8.sp.3.v3=mise(mu.est.tit.s8.sp.3.v3,mu.sp)
mise.tit.s8.sp.3.v4=mise(mu.est.tit.s8.sp.3.v4,mu.sp)
mise.tit.s8.sp.3.v5=mise(mu.est.tit.s8.sp.3.v5,mu.sp)

mise.tit.s8.bump.1.v1=mise(mu.est.tit.s8.bump.1.v1,mu.bump)
mise.tit.s8.bump.1.v2=mise(mu.est.tit.s8.bump.1.v2,mu.bump)
mise.tit.s8.bump.1.v3=mise(mu.est.tit.s8.bump.1.v3,mu.bump)
mise.tit.s8.bump.1.v4=mise(mu.est.tit.s8.bump.1.v4,mu.bump)
mise.tit.s8.bump.1.v5=mise(mu.est.tit.s8.bump.1.v5,mu.bump)
mise.tit.s8.bump.3.v1=mise(mu.est.tit.s8.bump.3.v1,mu.bump)
mise.tit.s8.bump.3.v2=mise(mu.est.tit.s8.bump.3.v2,mu.bump)
mise.tit.s8.bump.3.v3=mise(mu.est.tit.s8.bump.3.v3,mu.bump)
mise.tit.s8.bump.3.v4=mise(mu.est.tit.s8.bump.3.v4,mu.bump)
mise.tit.s8.bump.3.v5=mise(mu.est.tit.s8.bump.3.v5,mu.bump)

mise.tit.s8.blk.1.v1=mise(mu.est.tit.s8.blk.1.v1,mu.blk)
mise.tit.s8.blk.1.v2=mise(mu.est.tit.s8.blk.1.v2,mu.blk)
mise.tit.s8.blk.1.v3=mise(mu.est.tit.s8.blk.1.v3,mu.blk)
mise.tit.s8.blk.1.v4=mise(mu.est.tit.s8.blk.1.v4,mu.blk)
mise.tit.s8.blk.1.v5=mise(mu.est.tit.s8.blk.1.v5,mu.blk)
mise.tit.s8.blk.3.v1=mise(mu.est.tit.s8.blk.3.v1,mu.blk)
mise.tit.s8.blk.3.v2=mise(mu.est.tit.s8.blk.3.v2,mu.blk)
mise.tit.s8.blk.3.v3=mise(mu.est.tit.s8.blk.3.v3,mu.blk)
mise.tit.s8.blk.3.v4=mise(mu.est.tit.s8.blk.3.v4,mu.blk)
mise.tit.s8.blk.3.v5=mise(mu.est.tit.s8.blk.3.v5,mu.blk)

mise.tit.s8.ang.1.v1=mise(mu.est.tit.s8.ang.1.v1,mu.ang)
mise.tit.s8.ang.1.v2=mise(mu.est.tit.s8.ang.1.v2,mu.ang)
mise.tit.s8.ang.1.v3=mise(mu.est.tit.s8.ang.1.v3,mu.ang)
mise.tit.s8.ang.1.v4=mise(mu.est.tit.s8.ang.1.v4,mu.ang)
mise.tit.s8.ang.1.v5=mise(mu.est.tit.s8.ang.1.v5,mu.ang)
mise.tit.s8.ang.3.v1=mise(mu.est.tit.s8.ang.3.v1,mu.ang)
mise.tit.s8.ang.3.v2=mise(mu.est.tit.s8.ang.3.v2,mu.ang)
mise.tit.s8.ang.3.v3=mise(mu.est.tit.s8.ang.3.v3,mu.ang)
mise.tit.s8.ang.3.v4=mise(mu.est.tit.s8.ang.3.v4,mu.ang)
mise.tit.s8.ang.3.v5=mise(mu.est.tit.s8.ang.3.v5,mu.ang)

mise.tit.s8.dop.1.v1=mise(mu.est.tit.s8.dop.1.v1,mu.dop)
mise.tit.s8.dop.1.v2=mise(mu.est.tit.s8.dop.1.v2,mu.dop)
mise.tit.s8.dop.1.v3=mise(mu.est.tit.s8.dop.1.v3,mu.dop)
mise.tit.s8.dop.1.v4=mise(mu.est.tit.s8.dop.1.v4,mu.dop)
mise.tit.s8.dop.1.v5=mise(mu.est.tit.s8.dop.1.v5,mu.dop)
mise.tit.s8.dop.3.v1=mise(mu.est.tit.s8.dop.3.v1,mu.dop)
mise.tit.s8.dop.3.v2=mise(mu.est.tit.s8.dop.3.v2,mu.dop)
mise.tit.s8.dop.3.v3=mise(mu.est.tit.s8.dop.3.v3,mu.dop)
mise.tit.s8.dop.3.v4=mise(mu.est.tit.s8.dop.3.v4,mu.dop)
mise.tit.s8.dop.3.v5=mise(mu.est.tit.s8.dop.3.v5,mu.dop)

mise.tit.s8.blip.1.v1=mise(mu.est.tit.s8.blip.1.v1,mu.blip)
mise.tit.s8.blip.1.v2=mise(mu.est.tit.s8.blip.1.v2,mu.blip)
mise.tit.s8.blip.1.v3=mise(mu.est.tit.s8.blip.1.v3,mu.blip)
mise.tit.s8.blip.1.v4=mise(mu.est.tit.s8.blip.1.v4,mu.blip)
mise.tit.s8.blip.1.v5=mise(mu.est.tit.s8.blip.1.v5,mu.blip)
mise.tit.s8.blip.3.v1=mise(mu.est.tit.s8.blip.3.v1,mu.blip)
mise.tit.s8.blip.3.v2=mise(mu.est.tit.s8.blip.3.v2,mu.blip)
mise.tit.s8.blip.3.v3=mise(mu.est.tit.s8.blip.3.v3,mu.blip)
mise.tit.s8.blip.3.v4=mise(mu.est.tit.s8.blip.3.v4,mu.blip)
mise.tit.s8.blip.3.v5=mise(mu.est.tit.s8.blip.3.v5,mu.blip)

mise.tit.s8.cor.1.v1=mise(mu.est.tit.s8.cor.1.v1,mu.cor)
mise.tit.s8.cor.1.v2=mise(mu.est.tit.s8.cor.1.v2,mu.cor)
mise.tit.s8.cor.1.v3=mise(mu.est.tit.s8.cor.1.v3,mu.cor)
mise.tit.s8.cor.1.v4=mise(mu.est.tit.s8.cor.1.v4,mu.cor)
mise.tit.s8.cor.1.v5=mise(mu.est.tit.s8.cor.1.v5,mu.cor)
mise.tit.s8.cor.3.v1=mise(mu.est.tit.s8.cor.3.v1,mu.cor)
mise.tit.s8.cor.3.v2=mise(mu.est.tit.s8.cor.3.v2,mu.cor)
mise.tit.s8.cor.3.v3=mise(mu.est.tit.s8.cor.3.v3,mu.cor)
mise.tit.s8.cor.3.v4=mise(mu.est.tit.s8.cor.3.v4,mu.cor)
mise.tit.s8.cor.3.v5=mise(mu.est.tit.s8.cor.3.v5,mu.cor)



mise.tieb.haar.sp.1.v1=mise(mu.est.tieb.haar.sp.1.v1,mu.sp)
mise.tieb.haar.sp.1.v2=mise(mu.est.tieb.haar.sp.1.v2,mu.sp)
mise.tieb.haar.sp.1.v3=mise(mu.est.tieb.haar.sp.1.v3,mu.sp)
mise.tieb.haar.sp.1.v4=mise(mu.est.tieb.haar.sp.1.v4,mu.sp)
mise.tieb.haar.sp.1.v5=mise(mu.est.tieb.haar.sp.1.v5,mu.sp)
mise.tieb.haar.sp.3.v1=mise(mu.est.tieb.haar.sp.3.v1,mu.sp)
mise.tieb.haar.sp.3.v2=mise(mu.est.tieb.haar.sp.3.v2,mu.sp)
mise.tieb.haar.sp.3.v3=mise(mu.est.tieb.haar.sp.3.v3,mu.sp)
mise.tieb.haar.sp.3.v4=mise(mu.est.tieb.haar.sp.3.v4,mu.sp)
mise.tieb.haar.sp.3.v5=mise(mu.est.tieb.haar.sp.3.v5,mu.sp)

mise.tieb.haar.bump.1.v1=mise(mu.est.tieb.haar.bump.1.v1,mu.bump)
mise.tieb.haar.bump.1.v2=mise(mu.est.tieb.haar.bump.1.v2,mu.bump)
mise.tieb.haar.bump.1.v3=mise(mu.est.tieb.haar.bump.1.v3,mu.bump)
mise.tieb.haar.bump.1.v4=mise(mu.est.tieb.haar.bump.1.v4,mu.bump)
mise.tieb.haar.bump.1.v5=mise(mu.est.tieb.haar.bump.1.v5,mu.bump)
mise.tieb.haar.bump.3.v1=mise(mu.est.tieb.haar.bump.3.v1,mu.bump)
mise.tieb.haar.bump.3.v2=mise(mu.est.tieb.haar.bump.3.v2,mu.bump)
mise.tieb.haar.bump.3.v3=mise(mu.est.tieb.haar.bump.3.v3,mu.bump)
mise.tieb.haar.bump.3.v4=mise(mu.est.tieb.haar.bump.3.v4,mu.bump)
mise.tieb.haar.bump.3.v5=mise(mu.est.tieb.haar.bump.3.v5,mu.bump)

mise.tieb.haar.blk.1.v1=mise(mu.est.tieb.haar.blk.1.v1,mu.blk)
mise.tieb.haar.blk.1.v2=mise(mu.est.tieb.haar.blk.1.v2,mu.blk)
mise.tieb.haar.blk.1.v3=mise(mu.est.tieb.haar.blk.1.v3,mu.blk)
mise.tieb.haar.blk.1.v4=mise(mu.est.tieb.haar.blk.1.v4,mu.blk)
mise.tieb.haar.blk.1.v5=mise(mu.est.tieb.haar.blk.1.v5,mu.blk)
mise.tieb.haar.blk.3.v1=mise(mu.est.tieb.haar.blk.3.v1,mu.blk)
mise.tieb.haar.blk.3.v2=mise(mu.est.tieb.haar.blk.3.v2,mu.blk)
mise.tieb.haar.blk.3.v3=mise(mu.est.tieb.haar.blk.3.v3,mu.blk)
mise.tieb.haar.blk.3.v4=mise(mu.est.tieb.haar.blk.3.v4,mu.blk)
mise.tieb.haar.blk.3.v5=mise(mu.est.tieb.haar.blk.3.v5,mu.blk)

mise.tieb.haar.ang.1.v1=mise(mu.est.tieb.haar.ang.1.v1,mu.ang)
mise.tieb.haar.ang.1.v2=mise(mu.est.tieb.haar.ang.1.v2,mu.ang)
mise.tieb.haar.ang.1.v3=mise(mu.est.tieb.haar.ang.1.v3,mu.ang)
mise.tieb.haar.ang.1.v4=mise(mu.est.tieb.haar.ang.1.v4,mu.ang)
mise.tieb.haar.ang.1.v5=mise(mu.est.tieb.haar.ang.1.v5,mu.ang)
mise.tieb.haar.ang.3.v1=mise(mu.est.tieb.haar.ang.3.v1,mu.ang)
mise.tieb.haar.ang.3.v2=mise(mu.est.tieb.haar.ang.3.v2,mu.ang)
mise.tieb.haar.ang.3.v3=mise(mu.est.tieb.haar.ang.3.v3,mu.ang)
mise.tieb.haar.ang.3.v4=mise(mu.est.tieb.haar.ang.3.v4,mu.ang)
mise.tieb.haar.ang.3.v5=mise(mu.est.tieb.haar.ang.3.v5,mu.ang)

mise.tieb.haar.dop.1.v1=mise(mu.est.tieb.haar.dop.1.v1,mu.dop)
mise.tieb.haar.dop.1.v2=mise(mu.est.tieb.haar.dop.1.v2,mu.dop)
mise.tieb.haar.dop.1.v3=mise(mu.est.tieb.haar.dop.1.v3,mu.dop)
mise.tieb.haar.dop.1.v4=mise(mu.est.tieb.haar.dop.1.v4,mu.dop)
mise.tieb.haar.dop.1.v5=mise(mu.est.tieb.haar.dop.1.v5,mu.dop)
mise.tieb.haar.dop.3.v1=mise(mu.est.tieb.haar.dop.3.v1,mu.dop)
mise.tieb.haar.dop.3.v2=mise(mu.est.tieb.haar.dop.3.v2,mu.dop)
mise.tieb.haar.dop.3.v3=mise(mu.est.tieb.haar.dop.3.v3,mu.dop)
mise.tieb.haar.dop.3.v4=mise(mu.est.tieb.haar.dop.3.v4,mu.dop)
mise.tieb.haar.dop.3.v5=mise(mu.est.tieb.haar.dop.3.v5,mu.dop)

mise.tieb.haar.blip.1.v1=mise(mu.est.tieb.haar.blip.1.v1,mu.blip)
mise.tieb.haar.blip.1.v2=mise(mu.est.tieb.haar.blip.1.v2,mu.blip)
mise.tieb.haar.blip.1.v3=mise(mu.est.tieb.haar.blip.1.v3,mu.blip)
mise.tieb.haar.blip.1.v4=mise(mu.est.tieb.haar.blip.1.v4,mu.blip)
mise.tieb.haar.blip.1.v5=mise(mu.est.tieb.haar.blip.1.v5,mu.blip)
mise.tieb.haar.blip.3.v1=mise(mu.est.tieb.haar.blip.3.v1,mu.blip)
mise.tieb.haar.blip.3.v2=mise(mu.est.tieb.haar.blip.3.v2,mu.blip)
mise.tieb.haar.blip.3.v3=mise(mu.est.tieb.haar.blip.3.v3,mu.blip)
mise.tieb.haar.blip.3.v4=mise(mu.est.tieb.haar.blip.3.v4,mu.blip)
mise.tieb.haar.blip.3.v5=mise(mu.est.tieb.haar.blip.3.v5,mu.blip)

mise.tieb.haar.cor.1.v1=mise(mu.est.tieb.haar.cor.1.v1,mu.cor)
mise.tieb.haar.cor.1.v2=mise(mu.est.tieb.haar.cor.1.v2,mu.cor)
mise.tieb.haar.cor.1.v3=mise(mu.est.tieb.haar.cor.1.v3,mu.cor)
mise.tieb.haar.cor.1.v4=mise(mu.est.tieb.haar.cor.1.v4,mu.cor)
mise.tieb.haar.cor.1.v5=mise(mu.est.tieb.haar.cor.1.v5,mu.cor)
mise.tieb.haar.cor.3.v1=mise(mu.est.tieb.haar.cor.3.v1,mu.cor)
mise.tieb.haar.cor.3.v2=mise(mu.est.tieb.haar.cor.3.v2,mu.cor)
mise.tieb.haar.cor.3.v3=mise(mu.est.tieb.haar.cor.3.v3,mu.cor)
mise.tieb.haar.cor.3.v4=mise(mu.est.tieb.haar.cor.3.v4,mu.cor)
mise.tieb.haar.cor.3.v5=mise(mu.est.tieb.haar.cor.3.v5,mu.cor)


mise.tieb.s8.sp.1.v1=mise(mu.est.tieb.s8.sp.1.v1,mu.sp)
mise.tieb.s8.sp.1.v2=mise(mu.est.tieb.s8.sp.1.v2,mu.sp)
mise.tieb.s8.sp.1.v3=mise(mu.est.tieb.s8.sp.1.v3,mu.sp)
mise.tieb.s8.sp.1.v4=mise(mu.est.tieb.s8.sp.1.v4,mu.sp)
mise.tieb.s8.sp.1.v5=mise(mu.est.tieb.s8.sp.1.v5,mu.sp)
mise.tieb.s8.sp.3.v1=mise(mu.est.tieb.s8.sp.3.v1,mu.sp)
mise.tieb.s8.sp.3.v2=mise(mu.est.tieb.s8.sp.3.v2,mu.sp)
mise.tieb.s8.sp.3.v3=mise(mu.est.tieb.s8.sp.3.v3,mu.sp)
mise.tieb.s8.sp.3.v4=mise(mu.est.tieb.s8.sp.3.v4,mu.sp)
mise.tieb.s8.sp.3.v5=mise(mu.est.tieb.s8.sp.3.v5,mu.sp)

mise.tieb.s8.bump.1.v1=mise(mu.est.tieb.s8.bump.1.v1,mu.bump)
mise.tieb.s8.bump.1.v2=mise(mu.est.tieb.s8.bump.1.v2,mu.bump)
mise.tieb.s8.bump.1.v3=mise(mu.est.tieb.s8.bump.1.v3,mu.bump)
mise.tieb.s8.bump.1.v4=mise(mu.est.tieb.s8.bump.1.v4,mu.bump)
mise.tieb.s8.bump.1.v5=mise(mu.est.tieb.s8.bump.1.v5,mu.bump)
mise.tieb.s8.bump.3.v1=mise(mu.est.tieb.s8.bump.3.v1,mu.bump)
mise.tieb.s8.bump.3.v2=mise(mu.est.tieb.s8.bump.3.v2,mu.bump)
mise.tieb.s8.bump.3.v3=mise(mu.est.tieb.s8.bump.3.v3,mu.bump)
mise.tieb.s8.bump.3.v4=mise(mu.est.tieb.s8.bump.3.v4,mu.bump)
mise.tieb.s8.bump.3.v5=mise(mu.est.tieb.s8.bump.3.v5,mu.bump)

mise.tieb.s8.blk.1.v1=mise(mu.est.tieb.s8.blk.1.v1,mu.blk)
mise.tieb.s8.blk.1.v2=mise(mu.est.tieb.s8.blk.1.v2,mu.blk)
mise.tieb.s8.blk.1.v3=mise(mu.est.tieb.s8.blk.1.v3,mu.blk)
mise.tieb.s8.blk.1.v4=mise(mu.est.tieb.s8.blk.1.v4,mu.blk)
mise.tieb.s8.blk.1.v5=mise(mu.est.tieb.s8.blk.1.v5,mu.blk)
mise.tieb.s8.blk.3.v1=mise(mu.est.tieb.s8.blk.3.v1,mu.blk)
mise.tieb.s8.blk.3.v2=mise(mu.est.tieb.s8.blk.3.v2,mu.blk)
mise.tieb.s8.blk.3.v3=mise(mu.est.tieb.s8.blk.3.v3,mu.blk)
mise.tieb.s8.blk.3.v4=mise(mu.est.tieb.s8.blk.3.v4,mu.blk)
mise.tieb.s8.blk.3.v5=mise(mu.est.tieb.s8.blk.3.v5,mu.blk)

mise.tieb.s8.ang.1.v1=mise(mu.est.tieb.s8.ang.1.v1,mu.ang)
mise.tieb.s8.ang.1.v2=mise(mu.est.tieb.s8.ang.1.v2,mu.ang)
mise.tieb.s8.ang.1.v3=mise(mu.est.tieb.s8.ang.1.v3,mu.ang)
mise.tieb.s8.ang.1.v4=mise(mu.est.tieb.s8.ang.1.v4,mu.ang)
mise.tieb.s8.ang.1.v5=mise(mu.est.tieb.s8.ang.1.v5,mu.ang)
mise.tieb.s8.ang.3.v1=mise(mu.est.tieb.s8.ang.3.v1,mu.ang)
mise.tieb.s8.ang.3.v2=mise(mu.est.tieb.s8.ang.3.v2,mu.ang)
mise.tieb.s8.ang.3.v3=mise(mu.est.tieb.s8.ang.3.v3,mu.ang)
mise.tieb.s8.ang.3.v4=mise(mu.est.tieb.s8.ang.3.v4,mu.ang)
mise.tieb.s8.ang.3.v5=mise(mu.est.tieb.s8.ang.3.v5,mu.ang)

mise.tieb.s8.dop.1.v1=mise(mu.est.tieb.s8.dop.1.v1,mu.dop)
mise.tieb.s8.dop.1.v2=mise(mu.est.tieb.s8.dop.1.v2,mu.dop)
mise.tieb.s8.dop.1.v3=mise(mu.est.tieb.s8.dop.1.v3,mu.dop)
mise.tieb.s8.dop.1.v4=mise(mu.est.tieb.s8.dop.1.v4,mu.dop)
mise.tieb.s8.dop.1.v5=mise(mu.est.tieb.s8.dop.1.v5,mu.dop)
mise.tieb.s8.dop.3.v1=mise(mu.est.tieb.s8.dop.3.v1,mu.dop)
mise.tieb.s8.dop.3.v2=mise(mu.est.tieb.s8.dop.3.v2,mu.dop)
mise.tieb.s8.dop.3.v3=mise(mu.est.tieb.s8.dop.3.v3,mu.dop)
mise.tieb.s8.dop.3.v4=mise(mu.est.tieb.s8.dop.3.v4,mu.dop)
mise.tieb.s8.dop.3.v5=mise(mu.est.tieb.s8.dop.3.v5,mu.dop)

mise.tieb.s8.blip.1.v1=mise(mu.est.tieb.s8.blip.1.v1,mu.blip)
mise.tieb.s8.blip.1.v2=mise(mu.est.tieb.s8.blip.1.v2,mu.blip)
mise.tieb.s8.blip.1.v3=mise(mu.est.tieb.s8.blip.1.v3,mu.blip)
mise.tieb.s8.blip.1.v4=mise(mu.est.tieb.s8.blip.1.v4,mu.blip)
mise.tieb.s8.blip.1.v5=mise(mu.est.tieb.s8.blip.1.v5,mu.blip)
mise.tieb.s8.blip.3.v1=mise(mu.est.tieb.s8.blip.3.v1,mu.blip)
mise.tieb.s8.blip.3.v2=mise(mu.est.tieb.s8.blip.3.v2,mu.blip)
mise.tieb.s8.blip.3.v3=mise(mu.est.tieb.s8.blip.3.v3,mu.blip)
mise.tieb.s8.blip.3.v4=mise(mu.est.tieb.s8.blip.3.v4,mu.blip)
mise.tieb.s8.blip.3.v5=mise(mu.est.tieb.s8.blip.3.v5,mu.blip)

mise.tieb.s8.cor.1.v1=mise(mu.est.tieb.s8.cor.1.v1,mu.cor)
mise.tieb.s8.cor.1.v2=mise(mu.est.tieb.s8.cor.1.v2,mu.cor)
mise.tieb.s8.cor.1.v3=mise(mu.est.tieb.s8.cor.1.v3,mu.cor)
mise.tieb.s8.cor.1.v4=mise(mu.est.tieb.s8.cor.1.v4,mu.cor)
mise.tieb.s8.cor.1.v5=mise(mu.est.tieb.s8.cor.1.v5,mu.cor)
mise.tieb.s8.cor.3.v1=mise(mu.est.tieb.s8.cor.3.v1,mu.cor)
mise.tieb.s8.cor.3.v2=mise(mu.est.tieb.s8.cor.3.v2,mu.cor)
mise.tieb.s8.cor.3.v3=mise(mu.est.tieb.s8.cor.3.v3,mu.cor)
mise.tieb.s8.cor.3.v4=mise(mu.est.tieb.s8.cor.3.v4,mu.cor)
mise.tieb.s8.cor.3.v5=mise(mu.est.tieb.s8.cor.3.v5,mu.cor)



mise.ebayes.sp.1.v1=mise(mu.est.ebayes.sp.1.v1,mu.sp)
mise.ebayes.sp.1.v2=mise(mu.est.ebayes.sp.1.v2,mu.sp)
mise.ebayes.sp.1.v3=mise(mu.est.ebayes.sp.1.v3,mu.sp)
mise.ebayes.sp.1.v4=mise(mu.est.ebayes.sp.1.v4,mu.sp)
mise.ebayes.sp.1.v5=mise(mu.est.ebayes.sp.1.v5,mu.sp)
mise.ebayes.sp.3.v1=mise(mu.est.ebayes.sp.3.v1,mu.sp)
mise.ebayes.sp.3.v2=mise(mu.est.ebayes.sp.3.v2,mu.sp)
mise.ebayes.sp.3.v3=mise(mu.est.ebayes.sp.3.v3,mu.sp)
mise.ebayes.sp.3.v4=mise(mu.est.ebayes.sp.3.v4,mu.sp)
mise.ebayes.sp.3.v5=mise(mu.est.ebayes.sp.3.v5,mu.sp)

mise.ebayes.bump.1.v1=mise(mu.est.ebayes.bump.1.v1,mu.bump)
mise.ebayes.bump.1.v2=mise(mu.est.ebayes.bump.1.v2,mu.bump)
mise.ebayes.bump.1.v3=mise(mu.est.ebayes.bump.1.v3,mu.bump)
mise.ebayes.bump.1.v4=mise(mu.est.ebayes.bump.1.v4,mu.bump)
mise.ebayes.bump.1.v5=mise(mu.est.ebayes.bump.1.v5,mu.bump)
mise.ebayes.bump.3.v1=mise(mu.est.ebayes.bump.3.v1,mu.bump)
mise.ebayes.bump.3.v2=mise(mu.est.ebayes.bump.3.v2,mu.bump)
mise.ebayes.bump.3.v3=mise(mu.est.ebayes.bump.3.v3,mu.bump)
mise.ebayes.bump.3.v4=mise(mu.est.ebayes.bump.3.v4,mu.bump)
mise.ebayes.bump.3.v5=mise(mu.est.ebayes.bump.3.v5,mu.bump)

mise.ebayes.blk.1.v1=mise(mu.est.ebayes.blk.1.v1,mu.blk)
mise.ebayes.blk.1.v2=mise(mu.est.ebayes.blk.1.v2,mu.blk)
mise.ebayes.blk.1.v3=mise(mu.est.ebayes.blk.1.v3,mu.blk)
mise.ebayes.blk.1.v4=mise(mu.est.ebayes.blk.1.v4,mu.blk)
mise.ebayes.blk.1.v5=mise(mu.est.ebayes.blk.1.v5,mu.blk)
mise.ebayes.blk.3.v1=mise(mu.est.ebayes.blk.3.v1,mu.blk)
mise.ebayes.blk.3.v2=mise(mu.est.ebayes.blk.3.v2,mu.blk)
mise.ebayes.blk.3.v3=mise(mu.est.ebayes.blk.3.v3,mu.blk)
mise.ebayes.blk.3.v4=mise(mu.est.ebayes.blk.3.v4,mu.blk)
mise.ebayes.blk.3.v5=mise(mu.est.ebayes.blk.3.v5,mu.blk)

mise.ebayes.ang.1.v1=mise(mu.est.ebayes.ang.1.v1,mu.ang)
mise.ebayes.ang.1.v2=mise(mu.est.ebayes.ang.1.v2,mu.ang)
mise.ebayes.ang.1.v3=mise(mu.est.ebayes.ang.1.v3,mu.ang)
mise.ebayes.ang.1.v4=mise(mu.est.ebayes.ang.1.v4,mu.ang)
mise.ebayes.ang.1.v5=mise(mu.est.ebayes.ang.1.v5,mu.ang)
mise.ebayes.ang.3.v1=mise(mu.est.ebayes.ang.3.v1,mu.ang)
mise.ebayes.ang.3.v2=mise(mu.est.ebayes.ang.3.v2,mu.ang)
mise.ebayes.ang.3.v3=mise(mu.est.ebayes.ang.3.v3,mu.ang)
mise.ebayes.ang.3.v4=mise(mu.est.ebayes.ang.3.v4,mu.ang)
mise.ebayes.ang.3.v5=mise(mu.est.ebayes.ang.3.v5,mu.ang)

mise.ebayes.dop.1.v1=mise(mu.est.ebayes.dop.1.v1,mu.dop)
mise.ebayes.dop.1.v2=mise(mu.est.ebayes.dop.1.v2,mu.dop)
mise.ebayes.dop.1.v3=mise(mu.est.ebayes.dop.1.v3,mu.dop)
mise.ebayes.dop.1.v4=mise(mu.est.ebayes.dop.1.v4,mu.dop)
mise.ebayes.dop.1.v5=mise(mu.est.ebayes.dop.1.v5,mu.dop)
mise.ebayes.dop.3.v1=mise(mu.est.ebayes.dop.3.v1,mu.dop)
mise.ebayes.dop.3.v2=mise(mu.est.ebayes.dop.3.v2,mu.dop)
mise.ebayes.dop.3.v3=mise(mu.est.ebayes.dop.3.v3,mu.dop)
mise.ebayes.dop.3.v4=mise(mu.est.ebayes.dop.3.v4,mu.dop)
mise.ebayes.dop.3.v5=mise(mu.est.ebayes.dop.3.v5,mu.dop)

mise.ebayes.blip.1.v1=mise(mu.est.ebayes.blip.1.v1,mu.blip)
mise.ebayes.blip.1.v2=mise(mu.est.ebayes.blip.1.v2,mu.blip)
mise.ebayes.blip.1.v3=mise(mu.est.ebayes.blip.1.v3,mu.blip)
mise.ebayes.blip.1.v4=mise(mu.est.ebayes.blip.1.v4,mu.blip)
mise.ebayes.blip.1.v5=mise(mu.est.ebayes.blip.1.v5,mu.blip)
mise.ebayes.blip.3.v1=mise(mu.est.ebayes.blip.3.v1,mu.blip)
mise.ebayes.blip.3.v2=mise(mu.est.ebayes.blip.3.v2,mu.blip)
mise.ebayes.blip.3.v3=mise(mu.est.ebayes.blip.3.v3,mu.blip)
mise.ebayes.blip.3.v4=mise(mu.est.ebayes.blip.3.v4,mu.blip)
mise.ebayes.blip.3.v5=mise(mu.est.ebayes.blip.3.v5,mu.blip)

mise.ebayes.cor.1.v1=mise(mu.est.ebayes.cor.1.v1,mu.cor)
mise.ebayes.cor.1.v2=mise(mu.est.ebayes.cor.1.v2,mu.cor)
mise.ebayes.cor.1.v3=mise(mu.est.ebayes.cor.1.v3,mu.cor)
mise.ebayes.cor.1.v4=mise(mu.est.ebayes.cor.1.v4,mu.cor)
mise.ebayes.cor.1.v5=mise(mu.est.ebayes.cor.1.v5,mu.cor)
mise.ebayes.cor.3.v1=mise(mu.est.ebayes.cor.3.v1,mu.cor)
mise.ebayes.cor.3.v2=mise(mu.est.ebayes.cor.3.v2,mu.cor)
mise.ebayes.cor.3.v3=mise(mu.est.ebayes.cor.3.v3,mu.cor)
mise.ebayes.cor.3.v4=mise(mu.est.ebayes.cor.3.v4,mu.cor)
mise.ebayes.cor.3.v5=mise(mu.est.ebayes.cor.3.v5,mu.cor)




mise.ashe.s8.j.sp.1.v1=mise(mu.est.ashe.s8.j.sp.1.v1,mu.sp)
mise.ashe.s8.j.sp.1.v2=mise(mu.est.ashe.s8.j.sp.1.v2,mu.sp)
mise.ashe.s8.j.sp.1.v3=mise(mu.est.ashe.s8.j.sp.1.v3,mu.sp)
mise.ashe.s8.j.sp.1.v4=mise(mu.est.ashe.s8.j.sp.1.v4,mu.sp)
mise.ashe.s8.j.sp.1.v5=mise(mu.est.ashe.s8.j.sp.1.v5,mu.sp)
mise.ashe.s8.j.sp.3.v1=mise(mu.est.ashe.s8.j.sp.3.v1,mu.sp)
mise.ashe.s8.j.sp.3.v2=mise(mu.est.ashe.s8.j.sp.3.v2,mu.sp)
mise.ashe.s8.j.sp.3.v3=mise(mu.est.ashe.s8.j.sp.3.v3,mu.sp)
mise.ashe.s8.j.sp.3.v4=mise(mu.est.ashe.s8.j.sp.3.v4,mu.sp)
mise.ashe.s8.j.sp.3.v5=mise(mu.est.ashe.s8.j.sp.3.v5,mu.sp)

mise.ashe.s8.j.bump.1.v1=mise(mu.est.ashe.s8.j.bump.1.v1,mu.bump)
mise.ashe.s8.j.bump.1.v2=mise(mu.est.ashe.s8.j.bump.1.v2,mu.bump)
mise.ashe.s8.j.bump.1.v3=mise(mu.est.ashe.s8.j.bump.1.v3,mu.bump)
mise.ashe.s8.j.bump.1.v4=mise(mu.est.ashe.s8.j.bump.1.v4,mu.bump)
mise.ashe.s8.j.bump.1.v5=mise(mu.est.ashe.s8.j.bump.1.v5,mu.bump)
mise.ashe.s8.j.bump.3.v1=mise(mu.est.ashe.s8.j.bump.3.v1,mu.bump)
mise.ashe.s8.j.bump.3.v2=mise(mu.est.ashe.s8.j.bump.3.v2,mu.bump)
mise.ashe.s8.j.bump.3.v3=mise(mu.est.ashe.s8.j.bump.3.v3,mu.bump)
mise.ashe.s8.j.bump.3.v4=mise(mu.est.ashe.s8.j.bump.3.v4,mu.bump)
mise.ashe.s8.j.bump.3.v5=mise(mu.est.ashe.s8.j.bump.3.v5,mu.bump)

mise.ashe.s8.j.blk.1.v1=mise(mu.est.ashe.s8.j.blk.1.v1,mu.blk)
mise.ashe.s8.j.blk.1.v2=mise(mu.est.ashe.s8.j.blk.1.v2,mu.blk)
mise.ashe.s8.j.blk.1.v3=mise(mu.est.ashe.s8.j.blk.1.v3,mu.blk)
mise.ashe.s8.j.blk.1.v4=mise(mu.est.ashe.s8.j.blk.1.v4,mu.blk)
mise.ashe.s8.j.blk.1.v5=mise(mu.est.ashe.s8.j.blk.1.v5,mu.blk)
mise.ashe.s8.j.blk.3.v1=mise(mu.est.ashe.s8.j.blk.3.v1,mu.blk)
mise.ashe.s8.j.blk.3.v2=mise(mu.est.ashe.s8.j.blk.3.v2,mu.blk)
mise.ashe.s8.j.blk.3.v3=mise(mu.est.ashe.s8.j.blk.3.v3,mu.blk)
mise.ashe.s8.j.blk.3.v4=mise(mu.est.ashe.s8.j.blk.3.v4,mu.blk)
mise.ashe.s8.j.blk.3.v5=mise(mu.est.ashe.s8.j.blk.3.v5,mu.blk)

mise.ashe.s8.j.ang.1.v1=mise(mu.est.ashe.s8.j.ang.1.v1,mu.ang)
mise.ashe.s8.j.ang.1.v2=mise(mu.est.ashe.s8.j.ang.1.v2,mu.ang)
mise.ashe.s8.j.ang.1.v3=mise(mu.est.ashe.s8.j.ang.1.v3,mu.ang)
mise.ashe.s8.j.ang.1.v4=mise(mu.est.ashe.s8.j.ang.1.v4,mu.ang)
mise.ashe.s8.j.ang.1.v5=mise(mu.est.ashe.s8.j.ang.1.v5,mu.ang)
mise.ashe.s8.j.ang.3.v1=mise(mu.est.ashe.s8.j.ang.3.v1,mu.ang)
mise.ashe.s8.j.ang.3.v2=mise(mu.est.ashe.s8.j.ang.3.v2,mu.ang)
mise.ashe.s8.j.ang.3.v3=mise(mu.est.ashe.s8.j.ang.3.v3,mu.ang)
mise.ashe.s8.j.ang.3.v4=mise(mu.est.ashe.s8.j.ang.3.v4,mu.ang)
mise.ashe.s8.j.ang.3.v5=mise(mu.est.ashe.s8.j.ang.3.v5,mu.ang)

mise.ashe.s8.j.dop.1.v1=mise(mu.est.ashe.s8.j.dop.1.v1,mu.dop)
mise.ashe.s8.j.dop.1.v2=mise(mu.est.ashe.s8.j.dop.1.v2,mu.dop)
mise.ashe.s8.j.dop.1.v3=mise(mu.est.ashe.s8.j.dop.1.v3,mu.dop)
mise.ashe.s8.j.dop.1.v4=mise(mu.est.ashe.s8.j.dop.1.v4,mu.dop)
mise.ashe.s8.j.dop.1.v5=mise(mu.est.ashe.s8.j.dop.1.v5,mu.dop)
mise.ashe.s8.j.dop.3.v1=mise(mu.est.ashe.s8.j.dop.3.v1,mu.dop)
mise.ashe.s8.j.dop.3.v2=mise(mu.est.ashe.s8.j.dop.3.v2,mu.dop)
mise.ashe.s8.j.dop.3.v3=mise(mu.est.ashe.s8.j.dop.3.v3,mu.dop)
mise.ashe.s8.j.dop.3.v4=mise(mu.est.ashe.s8.j.dop.3.v4,mu.dop)
mise.ashe.s8.j.dop.3.v5=mise(mu.est.ashe.s8.j.dop.3.v5,mu.dop)

mise.ashe.s8.j.blip.1.v1=mise(mu.est.ashe.s8.j.blip.1.v1,mu.blip)
mise.ashe.s8.j.blip.1.v2=mise(mu.est.ashe.s8.j.blip.1.v2,mu.blip)
mise.ashe.s8.j.blip.1.v3=mise(mu.est.ashe.s8.j.blip.1.v3,mu.blip)
mise.ashe.s8.j.blip.1.v4=mise(mu.est.ashe.s8.j.blip.1.v4,mu.blip)
mise.ashe.s8.j.blip.1.v5=mise(mu.est.ashe.s8.j.blip.1.v5,mu.blip)
mise.ashe.s8.j.blip.3.v1=mise(mu.est.ashe.s8.j.blip.3.v1,mu.blip)
mise.ashe.s8.j.blip.3.v2=mise(mu.est.ashe.s8.j.blip.3.v2,mu.blip)
mise.ashe.s8.j.blip.3.v3=mise(mu.est.ashe.s8.j.blip.3.v3,mu.blip)
mise.ashe.s8.j.blip.3.v4=mise(mu.est.ashe.s8.j.blip.3.v4,mu.blip)
mise.ashe.s8.j.blip.3.v5=mise(mu.est.ashe.s8.j.blip.3.v5,mu.blip)

mise.ashe.s8.j.cor.1.v1=mise(mu.est.ashe.s8.j.cor.1.v1,mu.cor)
mise.ashe.s8.j.cor.1.v2=mise(mu.est.ashe.s8.j.cor.1.v2,mu.cor)
mise.ashe.s8.j.cor.1.v3=mise(mu.est.ashe.s8.j.cor.1.v3,mu.cor)
mise.ashe.s8.j.cor.1.v4=mise(mu.est.ashe.s8.j.cor.1.v4,mu.cor)
mise.ashe.s8.j.cor.1.v5=mise(mu.est.ashe.s8.j.cor.1.v5,mu.cor)
mise.ashe.s8.j.cor.3.v1=mise(mu.est.ashe.s8.j.cor.3.v1,mu.cor)
mise.ashe.s8.j.cor.3.v2=mise(mu.est.ashe.s8.j.cor.3.v2,mu.cor)
mise.ashe.s8.j.cor.3.v3=mise(mu.est.ashe.s8.j.cor.3.v3,mu.cor)
mise.ashe.s8.j.cor.3.v4=mise(mu.est.ashe.s8.j.cor.3.v4,mu.cor)
mise.ashe.s8.j.cor.3.v5=mise(mu.est.ashe.s8.j.cor.3.v5,mu.cor)






################################################





#1-ti,homo
#2-js,homo
#3-bams,hoomo
#4-nblk,homo
#5-sure,homo
#6-postmean,homo
#7-ebayes,homo
#8-ash_haar,homo
#9-ash_s8,homo
#10-ash_haar,est
#11-ash_s8,est
#12-ash_jash,est
#13-ti_haar,est
#14-ti_s8,est
#15-ti_haar,ash_est
#16-ti_s8,ash_est
#17-ash_haar,true
#18-ash_s8,true
#19-ti_haar,true
#20-ti_s8,true



mise.sp.1.v1=c(mise.ti.sp.1.v1,
mise.js.sp.1.v1,
mise.bams.sp.1.v1,
mise.nblk.sp.1.v1,
mise.sure.sp.1.v1,
mise.postmean.sp.1.v1,
mise.ebayes.sp.1.v1,
mise.ash.haar.sp.1.v1,
mise.ash.s8.sp.1.v1,
mise.ashe.haar.sp.1.v1,
mise.ashe.s8.sp.1.v1,
mise.ashe.s8.j.sp.1.v1,
mise.tie.haar.sp.1.v1,
mise.tie.s8.sp.1.v1,
mise.tieb.haar.sp.1.v1,
mise.tieb.s8.sp.1.v1,
mise.asht.haar.sp.1.v1,
mise.asht.s8.sp.1.v1,
mise.tit.haar.sp.1.v1,
mise.tit.s8.sp.1.v1
)

names(mise.sp.1.v1)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)


mise.sp.1.v2=c(mise.ti.sp.1.v2,
mise.js.sp.1.v2,
mise.bams.sp.1.v2,
mise.nblk.sp.1.v2,
mise.sure.sp.1.v2,
mise.postmean.sp.1.v2,
mise.ebayes.sp.1.v2,
mise.ash.haar.sp.1.v2,
mise.ash.s8.sp.1.v2,
mise.ashe.haar.sp.1.v2,
mise.ashe.s8.sp.1.v2,
mise.ashe.s8.j.sp.1.v2,
mise.tie.haar.sp.1.v2,
mise.tie.s8.sp.1.v2,
mise.tieb.haar.sp.1.v2,
mise.tieb.s8.sp.1.v2,
mise.asht.haar.sp.1.v2,
mise.asht.s8.sp.1.v2,
mise.tit.haar.sp.1.v2,
mise.tit.s8.sp.1.v2
)

names(mise.sp.1.v2)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)




mise.sp.1.v3=c(mise.ti.sp.1.v3,
mise.js.sp.1.v3,
mise.bams.sp.1.v3,
mise.nblk.sp.1.v3,
mise.sure.sp.1.v3,
mise.postmean.sp.1.v3,
mise.ebayes.sp.1.v3,
mise.ash.haar.sp.1.v3,
mise.ash.s8.sp.1.v3,
mise.ashe.haar.sp.1.v3,
mise.ashe.s8.sp.1.v3,
mise.ashe.s8.j.sp.1.v3,
mise.tie.haar.sp.1.v3,
mise.tie.s8.sp.1.v3,
mise.tieb.haar.sp.1.v3,
mise.tieb.s8.sp.1.v3,
mise.asht.haar.sp.1.v3,
mise.asht.s8.sp.1.v3,
mise.tit.haar.sp.1.v3,
mise.tit.s8.sp.1.v3
)

names(mise.sp.1.v3)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)




mise.sp.1.v4=c(mise.ti.sp.1.v4,
mise.js.sp.1.v4,
mise.bams.sp.1.v4,
mise.nblk.sp.1.v4,
mise.sure.sp.1.v4,
mise.postmean.sp.1.v4,
mise.ebayes.sp.1.v4,
mise.ash.haar.sp.1.v4,
mise.ash.s8.sp.1.v4,
mise.ashe.haar.sp.1.v4,
mise.ashe.s8.sp.1.v4,
mise.ashe.s8.j.sp.1.v4,
mise.tie.haar.sp.1.v4,
mise.tie.s8.sp.1.v4,
mise.tieb.haar.sp.1.v4,
mise.tieb.s8.sp.1.v4,
mise.asht.haar.sp.1.v4,
mise.asht.s8.sp.1.v4,
mise.tit.haar.sp.1.v4,
mise.tit.s8.sp.1.v4
)

names(mise.sp.1.v4)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)





mise.sp.1.v5=c(mise.ti.sp.1.v5,
mise.js.sp.1.v5,
mise.bams.sp.1.v5,
mise.nblk.sp.1.v5,
mise.sure.sp.1.v5,
mise.postmean.sp.1.v5,
mise.ebayes.sp.1.v5,
mise.ash.haar.sp.1.v5,
mise.ash.s8.sp.1.v5,
mise.ashe.haar.sp.1.v5,
mise.ashe.s8.sp.1.v5,
mise.ashe.s8.j.sp.1.v5,
mise.tie.haar.sp.1.v5,
mise.tie.s8.sp.1.v5,
mise.tieb.haar.sp.1.v5,
mise.tieb.s8.sp.1.v5,
mise.asht.haar.sp.1.v5,
mise.asht.s8.sp.1.v5,
mise.tit.haar.sp.1.v5,
mise.tit.s8.sp.1.v5
)

names(mise.sp.1.v5)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)



mise.sp.3.v1=c(mise.ti.sp.3.v1,
mise.js.sp.3.v1,
mise.bams.sp.3.v1,
mise.nblk.sp.3.v1,
mise.sure.sp.3.v1,
mise.postmean.sp.3.v1,
mise.ebayes.sp.3.v1,
mise.ash.haar.sp.3.v1,
mise.ash.s8.sp.3.v1,
mise.ashe.haar.sp.3.v1,
mise.ashe.s8.sp.3.v1,
mise.ashe.s8.j.sp.3.v1,
mise.tie.haar.sp.3.v1,
mise.tie.s8.sp.3.v1,
mise.tieb.haar.sp.3.v1,
mise.tieb.s8.sp.3.v1,
mise.asht.haar.sp.3.v1,
mise.asht.s8.sp.3.v1,
mise.tit.haar.sp.3.v1,
mise.tit.s8.sp.3.v1
)

names(mise.sp.3.v1)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)


mise.sp.3.v2=c(mise.ti.sp.3.v2,
mise.js.sp.3.v2,
mise.bams.sp.3.v2,
mise.nblk.sp.3.v2,
mise.sure.sp.3.v2,
mise.postmean.sp.3.v2,
mise.ebayes.sp.3.v2,
mise.ash.haar.sp.3.v2,
mise.ash.s8.sp.3.v2,
mise.ashe.haar.sp.3.v2,
mise.ashe.s8.sp.3.v2,
mise.ashe.s8.j.sp.3.v2,
mise.tie.haar.sp.3.v2,
mise.tie.s8.sp.3.v2,
mise.tieb.haar.sp.3.v2,
mise.tieb.s8.sp.3.v2,
mise.asht.haar.sp.3.v2,
mise.asht.s8.sp.3.v2,
mise.tit.haar.sp.3.v2,
mise.tit.s8.sp.3.v2
)

names(mise.sp.3.v2)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)




mise.sp.3.v3=c(mise.ti.sp.3.v3,
mise.js.sp.3.v3,
mise.bams.sp.3.v3,
mise.nblk.sp.3.v3,
mise.sure.sp.3.v3,
mise.postmean.sp.3.v3,
mise.ebayes.sp.3.v3,
mise.ash.haar.sp.3.v3,
mise.ash.s8.sp.3.v3,
mise.ashe.haar.sp.3.v3,
mise.ashe.s8.sp.3.v3,
mise.ashe.s8.j.sp.3.v3,
mise.tie.haar.sp.3.v3,
mise.tie.s8.sp.3.v3,
mise.tieb.haar.sp.3.v3,
mise.tieb.s8.sp.3.v3,
mise.asht.haar.sp.3.v3,
mise.asht.s8.sp.3.v3,
mise.tit.haar.sp.3.v3,
mise.tit.s8.sp.3.v3
)

names(mise.sp.3.v3)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)




mise.sp.3.v4=c(mise.ti.sp.3.v4,
mise.js.sp.3.v4,
mise.bams.sp.3.v4,
mise.nblk.sp.3.v4,
mise.sure.sp.3.v4,
mise.postmean.sp.3.v4,
mise.ebayes.sp.3.v4,
mise.ash.haar.sp.3.v4,
mise.ash.s8.sp.3.v4,
mise.ashe.haar.sp.3.v4,
mise.ashe.s8.sp.3.v4,
mise.ashe.s8.j.sp.3.v4,
mise.tie.haar.sp.3.v4,
mise.tie.s8.sp.3.v4,
mise.tieb.haar.sp.3.v4,
mise.tieb.s8.sp.3.v4,
mise.asht.haar.sp.3.v4,
mise.asht.s8.sp.3.v4,
mise.tit.haar.sp.3.v4,
mise.tit.s8.sp.3.v4
)

names(mise.sp.3.v4)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)





mise.sp.3.v5=c(mise.ti.sp.3.v5,
mise.js.sp.3.v5,
mise.bams.sp.3.v5,
mise.nblk.sp.3.v5,
mise.sure.sp.3.v5,
mise.postmean.sp.3.v5,
mise.ebayes.sp.3.v5,
mise.ash.haar.sp.3.v5,
mise.ash.s8.sp.3.v5,
mise.ashe.haar.sp.3.v5,
mise.ashe.s8.sp.3.v5,
mise.ashe.s8.j.sp.3.v5,
mise.tie.haar.sp.3.v5,
mise.tie.s8.sp.3.v5,
mise.tieb.haar.sp.3.v5,
mise.tieb.s8.sp.3.v5,
mise.asht.haar.sp.3.v5,
mise.asht.s8.sp.3.v5,
mise.tit.haar.sp.3.v5,
mise.tit.s8.sp.3.v5
)

names(mise.sp.3.v5)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)





mise.bump.1.v1=c(mise.ti.bump.1.v1,
mise.js.bump.1.v1,
mise.bams.bump.1.v1,
mise.nblk.bump.1.v1,
mise.sure.bump.1.v1,
mise.postmean.bump.1.v1,
mise.ebayes.bump.1.v1,
mise.ash.haar.bump.1.v1,
mise.ash.s8.bump.1.v1,
mise.ashe.haar.bump.1.v1,
mise.ashe.s8.bump.1.v1,
mise.ashe.s8.j.bump.1.v1,
mise.tie.haar.bump.1.v1,
mise.tie.s8.bump.1.v1,
mise.tieb.haar.bump.1.v1,
mise.tieb.s8.bump.1.v1,
mise.asht.haar.bump.1.v1,
mise.asht.s8.bump.1.v1,
mise.tit.haar.bump.1.v1,
mise.tit.s8.bump.1.v1
)

names(mise.bump.1.v1)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)


mise.bump.1.v2=c(mise.ti.bump.1.v2,
mise.js.bump.1.v2,
mise.bams.bump.1.v2,
mise.nblk.bump.1.v2,
mise.sure.bump.1.v2,
mise.postmean.bump.1.v2,
mise.ebayes.bump.1.v2,
mise.ash.haar.bump.1.v2,
mise.ash.s8.bump.1.v2,
mise.ashe.haar.bump.1.v2,
mise.ashe.s8.bump.1.v2,
mise.ashe.s8.j.bump.1.v2,
mise.tie.haar.bump.1.v2,
mise.tie.s8.bump.1.v2,
mise.tieb.haar.bump.1.v2,
mise.tieb.s8.bump.1.v2,
mise.asht.haar.bump.1.v2,
mise.asht.s8.bump.1.v2,
mise.tit.haar.bump.1.v2,
mise.tit.s8.bump.1.v2
)

names(mise.bump.1.v2)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)




mise.bump.1.v3=c(mise.ti.bump.1.v3,
mise.js.bump.1.v3,
mise.bams.bump.1.v3,
mise.nblk.bump.1.v3,
mise.sure.bump.1.v3,
mise.postmean.bump.1.v3,
mise.ebayes.bump.1.v3,
mise.ash.haar.bump.1.v3,
mise.ash.s8.bump.1.v3,
mise.ashe.haar.bump.1.v3,
mise.ashe.s8.bump.1.v3,
mise.ashe.s8.j.bump.1.v3,
mise.tie.haar.bump.1.v3,
mise.tie.s8.bump.1.v3,
mise.tieb.haar.bump.1.v3,
mise.tieb.s8.bump.1.v3,
mise.asht.haar.bump.1.v3,
mise.asht.s8.bump.1.v3,
mise.tit.haar.bump.1.v3,
mise.tit.s8.bump.1.v3
)

names(mise.bump.1.v3)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)




mise.bump.1.v4=c(mise.ti.bump.1.v4,
mise.js.bump.1.v4,
mise.bams.bump.1.v4,
mise.nblk.bump.1.v4,
mise.sure.bump.1.v4,
mise.postmean.bump.1.v4,
mise.ebayes.bump.1.v4,
mise.ash.haar.bump.1.v4,
mise.ash.s8.bump.1.v4,
mise.ashe.haar.bump.1.v4,
mise.ashe.s8.bump.1.v4,
mise.ashe.s8.j.bump.1.v4,
mise.tie.haar.bump.1.v4,
mise.tie.s8.bump.1.v4,
mise.tieb.haar.bump.1.v4,
mise.tieb.s8.bump.1.v4,
mise.asht.haar.bump.1.v4,
mise.asht.s8.bump.1.v4,
mise.tit.haar.bump.1.v4,
mise.tit.s8.bump.1.v4
)

names(mise.bump.1.v4)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)





mise.bump.1.v5=c(mise.ti.bump.1.v5,
mise.js.bump.1.v5,
mise.bams.bump.1.v5,
mise.nblk.bump.1.v5,
mise.sure.bump.1.v5,
mise.postmean.bump.1.v5,
mise.ebayes.bump.1.v5,
mise.ash.haar.bump.1.v5,
mise.ash.s8.bump.1.v5,
mise.ashe.haar.bump.1.v5,
mise.ashe.s8.bump.1.v5,
mise.ashe.s8.j.bump.1.v5,
mise.tie.haar.bump.1.v5,
mise.tie.s8.bump.1.v5,
mise.tieb.haar.bump.1.v5,
mise.tieb.s8.bump.1.v5,
mise.asht.haar.bump.1.v5,
mise.asht.s8.bump.1.v5,
mise.tit.haar.bump.1.v5,
mise.tit.s8.bump.1.v5
)

names(mise.bump.1.v5)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)



mise.bump.3.v1=c(mise.ti.bump.3.v1,
mise.js.bump.3.v1,
mise.bams.bump.3.v1,
mise.nblk.bump.3.v1,
mise.sure.bump.3.v1,
mise.postmean.bump.3.v1,
mise.ebayes.bump.3.v1,
mise.ash.haar.bump.3.v1,
mise.ash.s8.bump.3.v1,
mise.ashe.haar.bump.3.v1,
mise.ashe.s8.bump.3.v1,
mise.ashe.s8.j.bump.3.v1,
mise.tie.haar.bump.3.v1,
mise.tie.s8.bump.3.v1,
mise.tieb.haar.bump.3.v1,
mise.tieb.s8.bump.3.v1,
mise.asht.haar.bump.3.v1,
mise.asht.s8.bump.3.v1,
mise.tit.haar.bump.3.v1,
mise.tit.s8.bump.3.v1
)

names(mise.bump.3.v1)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)


mise.bump.3.v2=c(mise.ti.bump.3.v2,
mise.js.bump.3.v2,
mise.bams.bump.3.v2,
mise.nblk.bump.3.v2,
mise.sure.bump.3.v2,
mise.postmean.bump.3.v2,
mise.ebayes.bump.3.v2,
mise.ash.haar.bump.3.v2,
mise.ash.s8.bump.3.v2,
mise.ashe.haar.bump.3.v2,
mise.ashe.s8.bump.3.v2,
mise.ashe.s8.j.bump.3.v2,
mise.tie.haar.bump.3.v2,
mise.tie.s8.bump.3.v2,
mise.tieb.haar.bump.3.v2,
mise.tieb.s8.bump.3.v2,
mise.asht.haar.bump.3.v2,
mise.asht.s8.bump.3.v2,
mise.tit.haar.bump.3.v2,
mise.tit.s8.bump.3.v2
)

names(mise.bump.3.v2)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)




mise.bump.3.v3=c(mise.ti.bump.3.v3,
mise.js.bump.3.v3,
mise.bams.bump.3.v3,
mise.nblk.bump.3.v3,
mise.sure.bump.3.v3,
mise.postmean.bump.3.v3,
mise.ebayes.bump.3.v3,
mise.ash.haar.bump.3.v3,
mise.ash.s8.bump.3.v3,
mise.ashe.haar.bump.3.v3,
mise.ashe.s8.bump.3.v3,
mise.ashe.s8.j.bump.3.v3,
mise.tie.haar.bump.3.v3,
mise.tie.s8.bump.3.v3,
mise.tieb.haar.bump.3.v3,
mise.tieb.s8.bump.3.v3,
mise.asht.haar.bump.3.v3,
mise.asht.s8.bump.3.v3,
mise.tit.haar.bump.3.v3,
mise.tit.s8.bump.3.v3
)

names(mise.bump.3.v3)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)




mise.bump.3.v4=c(mise.ti.bump.3.v4,
mise.js.bump.3.v4,
mise.bams.bump.3.v4,
mise.nblk.bump.3.v4,
mise.sure.bump.3.v4,
mise.postmean.bump.3.v4,
mise.ebayes.bump.3.v4,
mise.ash.haar.bump.3.v4,
mise.ash.s8.bump.3.v4,
mise.ashe.haar.bump.3.v4,
mise.ashe.s8.bump.3.v4,
mise.ashe.s8.j.bump.3.v4,
mise.tie.haar.bump.3.v4,
mise.tie.s8.bump.3.v4,
mise.tieb.haar.bump.3.v4,
mise.tieb.s8.bump.3.v4,
mise.asht.haar.bump.3.v4,
mise.asht.s8.bump.3.v4,
mise.tit.haar.bump.3.v4,
mise.tit.s8.bump.3.v4
)

names(mise.bump.3.v4)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)





mise.bump.3.v5=c(mise.ti.bump.3.v5,
mise.js.bump.3.v5,
mise.bams.bump.3.v5,
mise.nblk.bump.3.v5,
mise.sure.bump.3.v5,
mise.postmean.bump.3.v5,
mise.ebayes.bump.3.v5,
mise.ash.haar.bump.3.v5,
mise.ash.s8.bump.3.v5,
mise.ashe.haar.bump.3.v5,
mise.ashe.s8.bump.3.v5,
mise.ashe.s8.j.bump.3.v5,
mise.tie.haar.bump.3.v5,
mise.tie.s8.bump.3.v5,
mise.tieb.haar.bump.3.v5,
mise.tieb.s8.bump.3.v5,
mise.asht.haar.bump.3.v5,
mise.asht.s8.bump.3.v5,
mise.tit.haar.bump.3.v5,
mise.tit.s8.bump.3.v5
)

names(mise.bump.3.v5)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)







mise.blk.1.v1=c(mise.ti.blk.1.v1,
mise.js.blk.1.v1,
mise.bams.blk.1.v1,
mise.nblk.blk.1.v1,
mise.sure.blk.1.v1,
mise.postmean.blk.1.v1,
mise.ebayes.blk.1.v1,
mise.ash.haar.blk.1.v1,
mise.ash.s8.blk.1.v1,
mise.ashe.haar.blk.1.v1,
mise.ashe.s8.blk.1.v1,
mise.ashe.s8.j.blk.1.v1,
mise.tie.haar.blk.1.v1,
mise.tie.s8.blk.1.v1,
mise.tieb.haar.blk.1.v1,
mise.tieb.s8.blk.1.v1,
mise.asht.haar.blk.1.v1,
mise.asht.s8.blk.1.v1,
mise.tit.haar.blk.1.v1,
mise.tit.s8.blk.1.v1
)

names(mise.blk.1.v1)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)


mise.blk.1.v2=c(mise.ti.blk.1.v2,
mise.js.blk.1.v2,
mise.bams.blk.1.v2,
mise.nblk.blk.1.v2,
mise.sure.blk.1.v2,
mise.postmean.blk.1.v2,
mise.ebayes.blk.1.v2,
mise.ash.haar.blk.1.v2,
mise.ash.s8.blk.1.v2,
mise.ashe.haar.blk.1.v2,
mise.ashe.s8.blk.1.v2,
mise.ashe.s8.j.blk.1.v2,
mise.tie.haar.blk.1.v2,
mise.tie.s8.blk.1.v2,
mise.tieb.haar.blk.1.v2,
mise.tieb.s8.blk.1.v2,
mise.asht.haar.blk.1.v2,
mise.asht.s8.blk.1.v2,
mise.tit.haar.blk.1.v2,
mise.tit.s8.blk.1.v2
)

names(mise.blk.1.v2)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)




mise.blk.1.v3=c(mise.ti.blk.1.v3,
mise.js.blk.1.v3,
mise.bams.blk.1.v3,
mise.nblk.blk.1.v3,
mise.sure.blk.1.v3,
mise.postmean.blk.1.v3,
mise.ebayes.blk.1.v3,
mise.ash.haar.blk.1.v3,
mise.ash.s8.blk.1.v3,
mise.ashe.haar.blk.1.v3,
mise.ashe.s8.blk.1.v3,
mise.ashe.s8.j.blk.1.v3,
mise.tie.haar.blk.1.v3,
mise.tie.s8.blk.1.v3,
mise.tieb.haar.blk.1.v3,
mise.tieb.s8.blk.1.v3,
mise.asht.haar.blk.1.v3,
mise.asht.s8.blk.1.v3,
mise.tit.haar.blk.1.v3,
mise.tit.s8.blk.1.v3
)

names(mise.blk.1.v3)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)




mise.blk.1.v4=c(mise.ti.blk.1.v4,
mise.js.blk.1.v4,
mise.bams.blk.1.v4,
mise.nblk.blk.1.v4,
mise.sure.blk.1.v4,
mise.postmean.blk.1.v4,
mise.ebayes.blk.1.v4,
mise.ash.haar.blk.1.v4,
mise.ash.s8.blk.1.v4,
mise.ashe.haar.blk.1.v4,
mise.ashe.s8.blk.1.v4,
mise.ashe.s8.j.blk.1.v4,
mise.tie.haar.blk.1.v4,
mise.tie.s8.blk.1.v4,
mise.tieb.haar.blk.1.v4,
mise.tieb.s8.blk.1.v4,
mise.asht.haar.blk.1.v4,
mise.asht.s8.blk.1.v4,
mise.tit.haar.blk.1.v4,
mise.tit.s8.blk.1.v4
)

names(mise.blk.1.v4)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)





mise.blk.1.v5=c(mise.ti.blk.1.v5,
mise.js.blk.1.v5,
mise.bams.blk.1.v5,
mise.nblk.blk.1.v5,
mise.sure.blk.1.v5,
mise.postmean.blk.1.v5,
mise.ebayes.blk.1.v5,
mise.ash.haar.blk.1.v5,
mise.ash.s8.blk.1.v5,
mise.ashe.haar.blk.1.v5,
mise.ashe.s8.blk.1.v5,
mise.ashe.s8.j.blk.1.v5,
mise.tie.haar.blk.1.v5,
mise.tie.s8.blk.1.v5,
mise.tieb.haar.blk.1.v5,
mise.tieb.s8.blk.1.v5,
mise.asht.haar.blk.1.v5,
mise.asht.s8.blk.1.v5,
mise.tit.haar.blk.1.v5,
mise.tit.s8.blk.1.v5
)

names(mise.blk.1.v5)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)



mise.blk.3.v1=c(mise.ti.blk.3.v1,
mise.js.blk.3.v1,
mise.bams.blk.3.v1,
mise.nblk.blk.3.v1,
mise.sure.blk.3.v1,
mise.postmean.blk.3.v1,
mise.ebayes.blk.3.v1,
mise.ash.haar.blk.3.v1,
mise.ash.s8.blk.3.v1,
mise.ashe.haar.blk.3.v1,
mise.ashe.s8.blk.3.v1,
mise.ashe.s8.j.blk.3.v1,
mise.tie.haar.blk.3.v1,
mise.tie.s8.blk.3.v1,
mise.tieb.haar.blk.3.v1,
mise.tieb.s8.blk.3.v1,
mise.asht.haar.blk.3.v1,
mise.asht.s8.blk.3.v1,
mise.tit.haar.blk.3.v1,
mise.tit.s8.blk.3.v1
)

names(mise.blk.3.v1)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)


mise.blk.3.v2=c(mise.ti.blk.3.v2,
mise.js.blk.3.v2,
mise.bams.blk.3.v2,
mise.nblk.blk.3.v2,
mise.sure.blk.3.v2,
mise.postmean.blk.3.v2,
mise.ebayes.blk.3.v2,
mise.ash.haar.blk.3.v2,
mise.ash.s8.blk.3.v2,
mise.ashe.haar.blk.3.v2,
mise.ashe.s8.blk.3.v2,
mise.ashe.s8.j.blk.3.v2,
mise.tie.haar.blk.3.v2,
mise.tie.s8.blk.3.v2,
mise.tieb.haar.blk.3.v2,
mise.tieb.s8.blk.3.v2,
mise.asht.haar.blk.3.v2,
mise.asht.s8.blk.3.v2,
mise.tit.haar.blk.3.v2,
mise.tit.s8.blk.3.v2
)

names(mise.blk.3.v2)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)




mise.blk.3.v3=c(mise.ti.blk.3.v3,
mise.js.blk.3.v3,
mise.bams.blk.3.v3,
mise.nblk.blk.3.v3,
mise.sure.blk.3.v3,
mise.postmean.blk.3.v3,
mise.ebayes.blk.3.v3,
mise.ash.haar.blk.3.v3,
mise.ash.s8.blk.3.v3,
mise.ashe.haar.blk.3.v3,
mise.ashe.s8.blk.3.v3,
mise.ashe.s8.j.blk.3.v3,
mise.tie.haar.blk.3.v3,
mise.tie.s8.blk.3.v3,
mise.tieb.haar.blk.3.v3,
mise.tieb.s8.blk.3.v3,
mise.asht.haar.blk.3.v3,
mise.asht.s8.blk.3.v3,
mise.tit.haar.blk.3.v3,
mise.tit.s8.blk.3.v3
)

names(mise.blk.3.v3)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)




mise.blk.3.v4=c(mise.ti.blk.3.v4,
mise.js.blk.3.v4,
mise.bams.blk.3.v4,
mise.nblk.blk.3.v4,
mise.sure.blk.3.v4,
mise.postmean.blk.3.v4,
mise.ebayes.blk.3.v4,
mise.ash.haar.blk.3.v4,
mise.ash.s8.blk.3.v4,
mise.ashe.haar.blk.3.v4,
mise.ashe.s8.blk.3.v4,
mise.ashe.s8.j.blk.3.v4,
mise.tie.haar.blk.3.v4,
mise.tie.s8.blk.3.v4,
mise.tieb.haar.blk.3.v4,
mise.tieb.s8.blk.3.v4,
mise.asht.haar.blk.3.v4,
mise.asht.s8.blk.3.v4,
mise.tit.haar.blk.3.v4,
mise.tit.s8.blk.3.v4
)

names(mise.blk.3.v4)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)





mise.blk.3.v5=c(mise.ti.blk.3.v5,
mise.js.blk.3.v5,
mise.bams.blk.3.v5,
mise.nblk.blk.3.v5,
mise.sure.blk.3.v5,
mise.postmean.blk.3.v5,
mise.ebayes.blk.3.v5,
mise.ash.haar.blk.3.v5,
mise.ash.s8.blk.3.v5,
mise.ashe.haar.blk.3.v5,
mise.ashe.s8.blk.3.v5,
mise.ashe.s8.j.blk.3.v5,
mise.tie.haar.blk.3.v5,
mise.tie.s8.blk.3.v5,
mise.tieb.haar.blk.3.v5,
mise.tieb.s8.blk.3.v5,
mise.asht.haar.blk.3.v5,
mise.asht.s8.blk.3.v5,
mise.tit.haar.blk.3.v5,
mise.tit.s8.blk.3.v5
)

names(mise.blk.3.v5)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)







mise.ang.1.v1=c(mise.ti.ang.1.v1,
mise.js.ang.1.v1,
mise.bams.ang.1.v1,
mise.nblk.ang.1.v1,
mise.sure.ang.1.v1,
mise.postmean.ang.1.v1,
mise.ebayes.ang.1.v1,
mise.ash.haar.ang.1.v1,
mise.ash.s8.ang.1.v1,
mise.ashe.haar.ang.1.v1,
mise.ashe.s8.ang.1.v1,
mise.ashe.s8.j.ang.1.v1,
mise.tie.haar.ang.1.v1,
mise.tie.s8.ang.1.v1,
mise.tieb.haar.ang.1.v1,
mise.tieb.s8.ang.1.v1,
mise.asht.haar.ang.1.v1,
mise.asht.s8.ang.1.v1,
mise.tit.haar.ang.1.v1,
mise.tit.s8.ang.1.v1
)

names(mise.ang.1.v1)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)


mise.ang.1.v2=c(mise.ti.ang.1.v2,
mise.js.ang.1.v2,
mise.bams.ang.1.v2,
mise.nblk.ang.1.v2,
mise.sure.ang.1.v2,
mise.postmean.ang.1.v2,
mise.ebayes.ang.1.v2,
mise.ash.haar.ang.1.v2,
mise.ash.s8.ang.1.v2,
mise.ashe.haar.ang.1.v2,
mise.ashe.s8.ang.1.v2,
mise.ashe.s8.j.ang.1.v2,
mise.tie.haar.ang.1.v2,
mise.tie.s8.ang.1.v2,
mise.tieb.haar.ang.1.v2,
mise.tieb.s8.ang.1.v2,
mise.asht.haar.ang.1.v2,
mise.asht.s8.ang.1.v2,
mise.tit.haar.ang.1.v2,
mise.tit.s8.ang.1.v2
)

names(mise.ang.1.v2)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)




mise.ang.1.v3=c(mise.ti.ang.1.v3,
mise.js.ang.1.v3,
mise.bams.ang.1.v3,
mise.nblk.ang.1.v3,
mise.sure.ang.1.v3,
mise.postmean.ang.1.v3,
mise.ebayes.ang.1.v3,
mise.ash.haar.ang.1.v3,
mise.ash.s8.ang.1.v3,
mise.ashe.haar.ang.1.v3,
mise.ashe.s8.ang.1.v3,
mise.ashe.s8.j.ang.1.v3,
mise.tie.haar.ang.1.v3,
mise.tie.s8.ang.1.v3,
mise.tieb.haar.ang.1.v3,
mise.tieb.s8.ang.1.v3,
mise.asht.haar.ang.1.v3,
mise.asht.s8.ang.1.v3,
mise.tit.haar.ang.1.v3,
mise.tit.s8.ang.1.v3
)

names(mise.ang.1.v3)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)




mise.ang.1.v4=c(mise.ti.ang.1.v4,
mise.js.ang.1.v4,
mise.bams.ang.1.v4,
mise.nblk.ang.1.v4,
mise.sure.ang.1.v4,
mise.postmean.ang.1.v4,
mise.ebayes.ang.1.v4,
mise.ash.haar.ang.1.v4,
mise.ash.s8.ang.1.v4,
mise.ashe.haar.ang.1.v4,
mise.ashe.s8.ang.1.v4,
mise.ashe.s8.j.ang.1.v4,
mise.tie.haar.ang.1.v4,
mise.tie.s8.ang.1.v4,
mise.tieb.haar.ang.1.v4,
mise.tieb.s8.ang.1.v4,
mise.asht.haar.ang.1.v4,
mise.asht.s8.ang.1.v4,
mise.tit.haar.ang.1.v4,
mise.tit.s8.ang.1.v4
)

names(mise.ang.1.v4)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)





mise.ang.1.v5=c(mise.ti.ang.1.v5,
mise.js.ang.1.v5,
mise.bams.ang.1.v5,
mise.nblk.ang.1.v5,
mise.sure.ang.1.v5,
mise.postmean.ang.1.v5,
mise.ebayes.ang.1.v5,
mise.ash.haar.ang.1.v5,
mise.ash.s8.ang.1.v5,
mise.ashe.haar.ang.1.v5,
mise.ashe.s8.ang.1.v5,
mise.ashe.s8.j.ang.1.v5,
mise.tie.haar.ang.1.v5,
mise.tie.s8.ang.1.v5,
mise.tieb.haar.ang.1.v5,
mise.tieb.s8.ang.1.v5,
mise.asht.haar.ang.1.v5,
mise.asht.s8.ang.1.v5,
mise.tit.haar.ang.1.v5,
mise.tit.s8.ang.1.v5
)

names(mise.ang.1.v5)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)



mise.ang.3.v1=c(mise.ti.ang.3.v1,
mise.js.ang.3.v1,
mise.bams.ang.3.v1,
mise.nblk.ang.3.v1,
mise.sure.ang.3.v1,
mise.postmean.ang.3.v1,
mise.ebayes.ang.3.v1,
mise.ash.haar.ang.3.v1,
mise.ash.s8.ang.3.v1,
mise.ashe.haar.ang.3.v1,
mise.ashe.s8.ang.3.v1,
mise.ashe.s8.j.ang.3.v1,
mise.tie.haar.ang.3.v1,
mise.tie.s8.ang.3.v1,
mise.tieb.haar.ang.3.v1,
mise.tieb.s8.ang.3.v1,
mise.asht.haar.ang.3.v1,
mise.asht.s8.ang.3.v1,
mise.tit.haar.ang.3.v1,
mise.tit.s8.ang.3.v1
)

names(mise.ang.3.v1)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)


mise.ang.3.v2=c(mise.ti.ang.3.v2,
mise.js.ang.3.v2,
mise.bams.ang.3.v2,
mise.nblk.ang.3.v2,
mise.sure.ang.3.v2,
mise.postmean.ang.3.v2,
mise.ebayes.ang.3.v2,
mise.ash.haar.ang.3.v2,
mise.ash.s8.ang.3.v2,
mise.ashe.haar.ang.3.v2,
mise.ashe.s8.ang.3.v2,
mise.ashe.s8.j.ang.3.v2,
mise.tie.haar.ang.3.v2,
mise.tie.s8.ang.3.v2,
mise.tieb.haar.ang.3.v2,
mise.tieb.s8.ang.3.v2,
mise.asht.haar.ang.3.v2,
mise.asht.s8.ang.3.v2,
mise.tit.haar.ang.3.v2,
mise.tit.s8.ang.3.v2
)

names(mise.ang.3.v2)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)




mise.ang.3.v3=c(mise.ti.ang.3.v3,
mise.js.ang.3.v3,
mise.bams.ang.3.v3,
mise.nblk.ang.3.v3,
mise.sure.ang.3.v3,
mise.postmean.ang.3.v3,
mise.ebayes.ang.3.v3,
mise.ash.haar.ang.3.v3,
mise.ash.s8.ang.3.v3,
mise.ashe.haar.ang.3.v3,
mise.ashe.s8.ang.3.v3,
mise.ashe.s8.j.ang.3.v3,
mise.tie.haar.ang.3.v3,
mise.tie.s8.ang.3.v3,
mise.tieb.haar.ang.3.v3,
mise.tieb.s8.ang.3.v3,
mise.asht.haar.ang.3.v3,
mise.asht.s8.ang.3.v3,
mise.tit.haar.ang.3.v3,
mise.tit.s8.ang.3.v3
)

names(mise.ang.3.v3)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)




mise.ang.3.v4=c(mise.ti.ang.3.v4,
mise.js.ang.3.v4,
mise.bams.ang.3.v4,
mise.nblk.ang.3.v4,
mise.sure.ang.3.v4,
mise.postmean.ang.3.v4,
mise.ebayes.ang.3.v4,
mise.ash.haar.ang.3.v4,
mise.ash.s8.ang.3.v4,
mise.ashe.haar.ang.3.v4,
mise.ashe.s8.ang.3.v4,
mise.ashe.s8.j.ang.3.v4,
mise.tie.haar.ang.3.v4,
mise.tie.s8.ang.3.v4,
mise.tieb.haar.ang.3.v4,
mise.tieb.s8.ang.3.v4,
mise.asht.haar.ang.3.v4,
mise.asht.s8.ang.3.v4,
mise.tit.haar.ang.3.v4,
mise.tit.s8.ang.3.v4
)

names(mise.ang.3.v4)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)





mise.ang.3.v5=c(mise.ti.ang.3.v5,
mise.js.ang.3.v5,
mise.bams.ang.3.v5,
mise.nblk.ang.3.v5,
mise.sure.ang.3.v5,
mise.postmean.ang.3.v5,
mise.ebayes.ang.3.v5,
mise.ash.haar.ang.3.v5,
mise.ash.s8.ang.3.v5,
mise.ashe.haar.ang.3.v5,
mise.ashe.s8.ang.3.v5,
mise.ashe.s8.j.ang.3.v5,
mise.tie.haar.ang.3.v5,
mise.tie.s8.ang.3.v5,
mise.tieb.haar.ang.3.v5,
mise.tieb.s8.ang.3.v5,
mise.asht.haar.ang.3.v5,
mise.asht.s8.ang.3.v5,
mise.tit.haar.ang.3.v5,
mise.tit.s8.ang.3.v5
)

names(mise.ang.3.v5)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)







mise.dop.1.v1=c(mise.ti.dop.1.v1,
mise.js.dop.1.v1,
mise.bams.dop.1.v1,
mise.nblk.dop.1.v1,
mise.sure.dop.1.v1,
mise.postmean.dop.1.v1,
mise.ebayes.dop.1.v1,
mise.ash.haar.dop.1.v1,
mise.ash.s8.dop.1.v1,
mise.ashe.haar.dop.1.v1,
mise.ashe.s8.dop.1.v1,
mise.ashe.s8.j.dop.1.v1,
mise.tie.haar.dop.1.v1,
mise.tie.s8.dop.1.v1,
mise.tieb.haar.dop.1.v1,
mise.tieb.s8.dop.1.v1,
mise.asht.haar.dop.1.v1,
mise.asht.s8.dop.1.v1,
mise.tit.haar.dop.1.v1,
mise.tit.s8.dop.1.v1
)

names(mise.dop.1.v1)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)


mise.dop.1.v2=c(mise.ti.dop.1.v2,
mise.js.dop.1.v2,
mise.bams.dop.1.v2,
mise.nblk.dop.1.v2,
mise.sure.dop.1.v2,
mise.postmean.dop.1.v2,
mise.ebayes.dop.1.v2,
mise.ash.haar.dop.1.v2,
mise.ash.s8.dop.1.v2,
mise.ashe.haar.dop.1.v2,
mise.ashe.s8.dop.1.v2,
mise.ashe.s8.j.dop.1.v2,
mise.tie.haar.dop.1.v2,
mise.tie.s8.dop.1.v2,
mise.tieb.haar.dop.1.v2,
mise.tieb.s8.dop.1.v2,
mise.asht.haar.dop.1.v2,
mise.asht.s8.dop.1.v2,
mise.tit.haar.dop.1.v2,
mise.tit.s8.dop.1.v2
)

names(mise.dop.1.v2)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)




mise.dop.1.v3=c(mise.ti.dop.1.v3,
mise.js.dop.1.v3,
mise.bams.dop.1.v3,
mise.nblk.dop.1.v3,
mise.sure.dop.1.v3,
mise.postmean.dop.1.v3,
mise.ebayes.dop.1.v3,
mise.ash.haar.dop.1.v3,
mise.ash.s8.dop.1.v3,
mise.ashe.haar.dop.1.v3,
mise.ashe.s8.dop.1.v3,
mise.ashe.s8.j.dop.1.v3,
mise.tie.haar.dop.1.v3,
mise.tie.s8.dop.1.v3,
mise.tieb.haar.dop.1.v3,
mise.tieb.s8.dop.1.v3,
mise.asht.haar.dop.1.v3,
mise.asht.s8.dop.1.v3,
mise.tit.haar.dop.1.v3,
mise.tit.s8.dop.1.v3
)

names(mise.dop.1.v3)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)




mise.dop.1.v4=c(mise.ti.dop.1.v4,
mise.js.dop.1.v4,
mise.bams.dop.1.v4,
mise.nblk.dop.1.v4,
mise.sure.dop.1.v4,
mise.postmean.dop.1.v4,
mise.ebayes.dop.1.v4,
mise.ash.haar.dop.1.v4,
mise.ash.s8.dop.1.v4,
mise.ashe.haar.dop.1.v4,
mise.ashe.s8.dop.1.v4,
mise.ashe.s8.j.dop.1.v4,
mise.tie.haar.dop.1.v4,
mise.tie.s8.dop.1.v4,
mise.tieb.haar.dop.1.v4,
mise.tieb.s8.dop.1.v4,
mise.asht.haar.dop.1.v4,
mise.asht.s8.dop.1.v4,
mise.tit.haar.dop.1.v4,
mise.tit.s8.dop.1.v4
)

names(mise.dop.1.v4)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)





mise.dop.1.v5=c(mise.ti.dop.1.v5,
mise.js.dop.1.v5,
mise.bams.dop.1.v5,
mise.nblk.dop.1.v5,
mise.sure.dop.1.v5,
mise.postmean.dop.1.v5,
mise.ebayes.dop.1.v5,
mise.ash.haar.dop.1.v5,
mise.ash.s8.dop.1.v5,
mise.ashe.haar.dop.1.v5,
mise.ashe.s8.dop.1.v5,
mise.ashe.s8.j.dop.1.v5,
mise.tie.haar.dop.1.v5,
mise.tie.s8.dop.1.v5,
mise.tieb.haar.dop.1.v5,
mise.tieb.s8.dop.1.v5,
mise.asht.haar.dop.1.v5,
mise.asht.s8.dop.1.v5,
mise.tit.haar.dop.1.v5,
mise.tit.s8.dop.1.v5
)

names(mise.dop.1.v5)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)



mise.dop.3.v1=c(mise.ti.dop.3.v1,
mise.js.dop.3.v1,
mise.bams.dop.3.v1,
mise.nblk.dop.3.v1,
mise.sure.dop.3.v1,
mise.postmean.dop.3.v1,
mise.ebayes.dop.3.v1,
mise.ash.haar.dop.3.v1,
mise.ash.s8.dop.3.v1,
mise.ashe.haar.dop.3.v1,
mise.ashe.s8.dop.3.v1,
mise.ashe.s8.j.dop.3.v1,
mise.tie.haar.dop.3.v1,
mise.tie.s8.dop.3.v1,
mise.tieb.haar.dop.3.v1,
mise.tieb.s8.dop.3.v1,
mise.asht.haar.dop.3.v1,
mise.asht.s8.dop.3.v1,
mise.tit.haar.dop.3.v1,
mise.tit.s8.dop.3.v1
)

names(mise.dop.3.v1)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)


mise.dop.3.v2=c(mise.ti.dop.3.v2,
mise.js.dop.3.v2,
mise.bams.dop.3.v2,
mise.nblk.dop.3.v2,
mise.sure.dop.3.v2,
mise.postmean.dop.3.v2,
mise.ebayes.dop.3.v2,
mise.ash.haar.dop.3.v2,
mise.ash.s8.dop.3.v2,
mise.ashe.haar.dop.3.v2,
mise.ashe.s8.dop.3.v2,
mise.ashe.s8.j.dop.3.v2,
mise.tie.haar.dop.3.v2,
mise.tie.s8.dop.3.v2,
mise.tieb.haar.dop.3.v2,
mise.tieb.s8.dop.3.v2,
mise.asht.haar.dop.3.v2,
mise.asht.s8.dop.3.v2,
mise.tit.haar.dop.3.v2,
mise.tit.s8.dop.3.v2
)

names(mise.dop.3.v2)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)




mise.dop.3.v3=c(mise.ti.dop.3.v3,
mise.js.dop.3.v3,
mise.bams.dop.3.v3,
mise.nblk.dop.3.v3,
mise.sure.dop.3.v3,
mise.postmean.dop.3.v3,
mise.ebayes.dop.3.v3,
mise.ash.haar.dop.3.v3,
mise.ash.s8.dop.3.v3,
mise.ashe.haar.dop.3.v3,
mise.ashe.s8.dop.3.v3,
mise.ashe.s8.j.dop.3.v3,
mise.tie.haar.dop.3.v3,
mise.tie.s8.dop.3.v3,
mise.tieb.haar.dop.3.v3,
mise.tieb.s8.dop.3.v3,
mise.asht.haar.dop.3.v3,
mise.asht.s8.dop.3.v3,
mise.tit.haar.dop.3.v3,
mise.tit.s8.dop.3.v3
)

names(mise.dop.3.v3)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)




mise.dop.3.v4=c(mise.ti.dop.3.v4,
mise.js.dop.3.v4,
mise.bams.dop.3.v4,
mise.nblk.dop.3.v4,
mise.sure.dop.3.v4,
mise.postmean.dop.3.v4,
mise.ebayes.dop.3.v4,
mise.ash.haar.dop.3.v4,
mise.ash.s8.dop.3.v4,
mise.ashe.haar.dop.3.v4,
mise.ashe.s8.dop.3.v4,
mise.ashe.s8.j.dop.3.v4,
mise.tie.haar.dop.3.v4,
mise.tie.s8.dop.3.v4,
mise.tieb.haar.dop.3.v4,
mise.tieb.s8.dop.3.v4,
mise.asht.haar.dop.3.v4,
mise.asht.s8.dop.3.v4,
mise.tit.haar.dop.3.v4,
mise.tit.s8.dop.3.v4
)

names(mise.dop.3.v4)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)





mise.dop.3.v5=c(mise.ti.dop.3.v5,
mise.js.dop.3.v5,
mise.bams.dop.3.v5,
mise.nblk.dop.3.v5,
mise.sure.dop.3.v5,
mise.postmean.dop.3.v5,
mise.ebayes.dop.3.v5,
mise.ash.haar.dop.3.v5,
mise.ash.s8.dop.3.v5,
mise.ashe.haar.dop.3.v5,
mise.ashe.s8.dop.3.v5,
mise.ashe.s8.j.dop.3.v5,
mise.tie.haar.dop.3.v5,
mise.tie.s8.dop.3.v5,
mise.tieb.haar.dop.3.v5,
mise.tieb.s8.dop.3.v5,
mise.asht.haar.dop.3.v5,
mise.asht.s8.dop.3.v5,
mise.tit.haar.dop.3.v5,
mise.tit.s8.dop.3.v5
)

names(mise.dop.3.v5)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)








mise.blip.1.v1=c(mise.ti.blip.1.v1,
mise.js.blip.1.v1,
mise.bams.blip.1.v1,
mise.nblk.blip.1.v1,
mise.sure.blip.1.v1,
mise.postmean.blip.1.v1,
mise.ebayes.blip.1.v1,
mise.ash.haar.blip.1.v1,
mise.ash.s8.blip.1.v1,
mise.ashe.haar.blip.1.v1,
mise.ashe.s8.blip.1.v1,
mise.ashe.s8.j.blip.1.v1,
mise.tie.haar.blip.1.v1,
mise.tie.s8.blip.1.v1,
mise.tieb.haar.blip.1.v1,
mise.tieb.s8.blip.1.v1,
mise.asht.haar.blip.1.v1,
mise.asht.s8.blip.1.v1,
mise.tit.haar.blip.1.v1,
mise.tit.s8.blip.1.v1
)

names(mise.blip.1.v1)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)


mise.blip.1.v2=c(mise.ti.blip.1.v2,
mise.js.blip.1.v2,
mise.bams.blip.1.v2,
mise.nblk.blip.1.v2,
mise.sure.blip.1.v2,
mise.postmean.blip.1.v2,
mise.ebayes.blip.1.v2,
mise.ash.haar.blip.1.v2,
mise.ash.s8.blip.1.v2,
mise.ashe.haar.blip.1.v2,
mise.ashe.s8.blip.1.v2,
mise.ashe.s8.j.blip.1.v2,
mise.tie.haar.blip.1.v2,
mise.tie.s8.blip.1.v2,
mise.tieb.haar.blip.1.v2,
mise.tieb.s8.blip.1.v2,
mise.asht.haar.blip.1.v2,
mise.asht.s8.blip.1.v2,
mise.tit.haar.blip.1.v2,
mise.tit.s8.blip.1.v2
)

names(mise.blip.1.v2)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)




mise.blip.1.v3=c(mise.ti.blip.1.v3,
mise.js.blip.1.v3,
mise.bams.blip.1.v3,
mise.nblk.blip.1.v3,
mise.sure.blip.1.v3,
mise.postmean.blip.1.v3,
mise.ebayes.blip.1.v3,
mise.ash.haar.blip.1.v3,
mise.ash.s8.blip.1.v3,
mise.ashe.haar.blip.1.v3,
mise.ashe.s8.blip.1.v3,
mise.ashe.s8.j.blip.1.v3,
mise.tie.haar.blip.1.v3,
mise.tie.s8.blip.1.v3,
mise.tieb.haar.blip.1.v3,
mise.tieb.s8.blip.1.v3,
mise.asht.haar.blip.1.v3,
mise.asht.s8.blip.1.v3,
mise.tit.haar.blip.1.v3,
mise.tit.s8.blip.1.v3
)

names(mise.blip.1.v3)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)




mise.blip.1.v4=c(mise.ti.blip.1.v4,
mise.js.blip.1.v4,
mise.bams.blip.1.v4,
mise.nblk.blip.1.v4,
mise.sure.blip.1.v4,
mise.postmean.blip.1.v4,
mise.ebayes.blip.1.v4,
mise.ash.haar.blip.1.v4,
mise.ash.s8.blip.1.v4,
mise.ashe.haar.blip.1.v4,
mise.ashe.s8.blip.1.v4,
mise.ashe.s8.j.blip.1.v4,
mise.tie.haar.blip.1.v4,
mise.tie.s8.blip.1.v4,
mise.tieb.haar.blip.1.v4,
mise.tieb.s8.blip.1.v4,
mise.asht.haar.blip.1.v4,
mise.asht.s8.blip.1.v4,
mise.tit.haar.blip.1.v4,
mise.tit.s8.blip.1.v4
)

names(mise.blip.1.v4)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)





mise.blip.1.v5=c(mise.ti.blip.1.v5,
mise.js.blip.1.v5,
mise.bams.blip.1.v5,
mise.nblk.blip.1.v5,
mise.sure.blip.1.v5,
mise.postmean.blip.1.v5,
mise.ebayes.blip.1.v5,
mise.ash.haar.blip.1.v5,
mise.ash.s8.blip.1.v5,
mise.ashe.haar.blip.1.v5,
mise.ashe.s8.blip.1.v5,
mise.ashe.s8.j.blip.1.v5,
mise.tie.haar.blip.1.v5,
mise.tie.s8.blip.1.v5,
mise.tieb.haar.blip.1.v5,
mise.tieb.s8.blip.1.v5,
mise.asht.haar.blip.1.v5,
mise.asht.s8.blip.1.v5,
mise.tit.haar.blip.1.v5,
mise.tit.s8.blip.1.v5
)

names(mise.blip.1.v5)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)



mise.blip.3.v1=c(mise.ti.blip.3.v1,
mise.js.blip.3.v1,
mise.bams.blip.3.v1,
mise.nblk.blip.3.v1,
mise.sure.blip.3.v1,
mise.postmean.blip.3.v1,
mise.ebayes.blip.3.v1,
mise.ash.haar.blip.3.v1,
mise.ash.s8.blip.3.v1,
mise.ashe.haar.blip.3.v1,
mise.ashe.s8.blip.3.v1,
mise.ashe.s8.j.blip.3.v1,
mise.tie.haar.blip.3.v1,
mise.tie.s8.blip.3.v1,
mise.tieb.haar.blip.3.v1,
mise.tieb.s8.blip.3.v1,
mise.asht.haar.blip.3.v1,
mise.asht.s8.blip.3.v1,
mise.tit.haar.blip.3.v1,
mise.tit.s8.blip.3.v1
)

names(mise.blip.3.v1)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)


mise.blip.3.v2=c(mise.ti.blip.3.v2,
mise.js.blip.3.v2,
mise.bams.blip.3.v2,
mise.nblk.blip.3.v2,
mise.sure.blip.3.v2,
mise.postmean.blip.3.v2,
mise.ebayes.blip.3.v2,
mise.ash.haar.blip.3.v2,
mise.ash.s8.blip.3.v2,
mise.ashe.haar.blip.3.v2,
mise.ashe.s8.blip.3.v2,
mise.ashe.s8.j.blip.3.v2,
mise.tie.haar.blip.3.v2,
mise.tie.s8.blip.3.v2,
mise.tieb.haar.blip.3.v2,
mise.tieb.s8.blip.3.v2,
mise.asht.haar.blip.3.v2,
mise.asht.s8.blip.3.v2,
mise.tit.haar.blip.3.v2,
mise.tit.s8.blip.3.v2
)

names(mise.blip.3.v2)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)




mise.blip.3.v3=c(mise.ti.blip.3.v3,
mise.js.blip.3.v3,
mise.bams.blip.3.v3,
mise.nblk.blip.3.v3,
mise.sure.blip.3.v3,
mise.postmean.blip.3.v3,
mise.ebayes.blip.3.v3,
mise.ash.haar.blip.3.v3,
mise.ash.s8.blip.3.v3,
mise.ashe.haar.blip.3.v3,
mise.ashe.s8.blip.3.v3,
mise.ashe.s8.j.blip.3.v3,
mise.tie.haar.blip.3.v3,
mise.tie.s8.blip.3.v3,
mise.tieb.haar.blip.3.v3,
mise.tieb.s8.blip.3.v3,
mise.asht.haar.blip.3.v3,
mise.asht.s8.blip.3.v3,
mise.tit.haar.blip.3.v3,
mise.tit.s8.blip.3.v3
)

names(mise.blip.3.v3)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)




mise.blip.3.v4=c(mise.ti.blip.3.v4,
mise.js.blip.3.v4,
mise.bams.blip.3.v4,
mise.nblk.blip.3.v4,
mise.sure.blip.3.v4,
mise.postmean.blip.3.v4,
mise.ebayes.blip.3.v4,
mise.ash.haar.blip.3.v4,
mise.ash.s8.blip.3.v4,
mise.ashe.haar.blip.3.v4,
mise.ashe.s8.blip.3.v4,
mise.ashe.s8.j.blip.3.v4,
mise.tie.haar.blip.3.v4,
mise.tie.s8.blip.3.v4,
mise.tieb.haar.blip.3.v4,
mise.tieb.s8.blip.3.v4,
mise.asht.haar.blip.3.v4,
mise.asht.s8.blip.3.v4,
mise.tit.haar.blip.3.v4,
mise.tit.s8.blip.3.v4
)

names(mise.blip.3.v4)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)





mise.blip.3.v5=c(mise.ti.blip.3.v5,
mise.js.blip.3.v5,
mise.bams.blip.3.v5,
mise.nblk.blip.3.v5,
mise.sure.blip.3.v5,
mise.postmean.blip.3.v5,
mise.ebayes.blip.3.v5,
mise.ash.haar.blip.3.v5,
mise.ash.s8.blip.3.v5,
mise.ashe.haar.blip.3.v5,
mise.ashe.s8.blip.3.v5,
mise.ashe.s8.j.blip.3.v5,
mise.tie.haar.blip.3.v5,
mise.tie.s8.blip.3.v5,
mise.tieb.haar.blip.3.v5,
mise.tieb.s8.blip.3.v5,
mise.asht.haar.blip.3.v5,
mise.asht.s8.blip.3.v5,
mise.tit.haar.blip.3.v5,
mise.tit.s8.blip.3.v5
)

names(mise.blip.3.v5)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)








mise.cor.1.v1=c(mise.ti.cor.1.v1,
mise.js.cor.1.v1,
mise.bams.cor.1.v1,
mise.nblk.cor.1.v1,
mise.sure.cor.1.v1,
mise.postmean.cor.1.v1,
mise.ebayes.cor.1.v1,
mise.ash.haar.cor.1.v1,
mise.ash.s8.cor.1.v1,
mise.ashe.haar.cor.1.v1,
mise.ashe.s8.cor.1.v1,
mise.ashe.s8.j.cor.1.v1,
mise.tie.haar.cor.1.v1,
mise.tie.s8.cor.1.v1,
mise.tieb.haar.cor.1.v1,
mise.tieb.s8.cor.1.v1,
mise.asht.haar.cor.1.v1,
mise.asht.s8.cor.1.v1,
mise.tit.haar.cor.1.v1,
mise.tit.s8.cor.1.v1
)

names(mise.cor.1.v1)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)


mise.cor.1.v2=c(mise.ti.cor.1.v2,
mise.js.cor.1.v2,
mise.bams.cor.1.v2,
mise.nblk.cor.1.v2,
mise.sure.cor.1.v2,
mise.postmean.cor.1.v2,
mise.ebayes.cor.1.v2,
mise.ash.haar.cor.1.v2,
mise.ash.s8.cor.1.v2,
mise.ashe.haar.cor.1.v2,
mise.ashe.s8.cor.1.v2,
mise.ashe.s8.j.cor.1.v2,
mise.tie.haar.cor.1.v2,
mise.tie.s8.cor.1.v2,
mise.tieb.haar.cor.1.v2,
mise.tieb.s8.cor.1.v2,
mise.asht.haar.cor.1.v2,
mise.asht.s8.cor.1.v2,
mise.tit.haar.cor.1.v2,
mise.tit.s8.cor.1.v2
)

names(mise.cor.1.v2)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)




mise.cor.1.v3=c(mise.ti.cor.1.v3,
mise.js.cor.1.v3,
mise.bams.cor.1.v3,
mise.nblk.cor.1.v3,
mise.sure.cor.1.v3,
mise.postmean.cor.1.v3,
mise.ebayes.cor.1.v3,
mise.ash.haar.cor.1.v3,
mise.ash.s8.cor.1.v3,
mise.ashe.haar.cor.1.v3,
mise.ashe.s8.cor.1.v3,
mise.ashe.s8.j.cor.1.v3,
mise.tie.haar.cor.1.v3,
mise.tie.s8.cor.1.v3,
mise.tieb.haar.cor.1.v3,
mise.tieb.s8.cor.1.v3,
mise.asht.haar.cor.1.v3,
mise.asht.s8.cor.1.v3,
mise.tit.haar.cor.1.v3,
mise.tit.s8.cor.1.v3
)

names(mise.cor.1.v3)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)




mise.cor.1.v4=c(mise.ti.cor.1.v4,
mise.js.cor.1.v4,
mise.bams.cor.1.v4,
mise.nblk.cor.1.v4,
mise.sure.cor.1.v4,
mise.postmean.cor.1.v4,
mise.ebayes.cor.1.v4,
mise.ash.haar.cor.1.v4,
mise.ash.s8.cor.1.v4,
mise.ashe.haar.cor.1.v4,
mise.ashe.s8.cor.1.v4,
mise.ashe.s8.j.cor.1.v4,
mise.tie.haar.cor.1.v4,
mise.tie.s8.cor.1.v4,
mise.tieb.haar.cor.1.v4,
mise.tieb.s8.cor.1.v4,
mise.asht.haar.cor.1.v4,
mise.asht.s8.cor.1.v4,
mise.tit.haar.cor.1.v4,
mise.tit.s8.cor.1.v4
)

names(mise.cor.1.v4)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)





mise.cor.1.v5=c(mise.ti.cor.1.v5,
mise.js.cor.1.v5,
mise.bams.cor.1.v5,
mise.nblk.cor.1.v5,
mise.sure.cor.1.v5,
mise.postmean.cor.1.v5,
mise.ebayes.cor.1.v5,
mise.ash.haar.cor.1.v5,
mise.ash.s8.cor.1.v5,
mise.ashe.haar.cor.1.v5,
mise.ashe.s8.cor.1.v5,
mise.ashe.s8.j.cor.1.v5,
mise.tie.haar.cor.1.v5,
mise.tie.s8.cor.1.v5,
mise.tieb.haar.cor.1.v5,
mise.tieb.s8.cor.1.v5,
mise.asht.haar.cor.1.v5,
mise.asht.s8.cor.1.v5,
mise.tit.haar.cor.1.v5,
mise.tit.s8.cor.1.v5
)

names(mise.cor.1.v5)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)



mise.cor.3.v1=c(mise.ti.cor.3.v1,
mise.js.cor.3.v1,
mise.bams.cor.3.v1,
mise.nblk.cor.3.v1,
mise.sure.cor.3.v1,
mise.postmean.cor.3.v1,
mise.ebayes.cor.3.v1,
mise.ash.haar.cor.3.v1,
mise.ash.s8.cor.3.v1,
mise.ashe.haar.cor.3.v1,
mise.ashe.s8.cor.3.v1,
mise.ashe.s8.j.cor.3.v1,
mise.tie.haar.cor.3.v1,
mise.tie.s8.cor.3.v1,
mise.tieb.haar.cor.3.v1,
mise.tieb.s8.cor.3.v1,
mise.asht.haar.cor.3.v1,
mise.asht.s8.cor.3.v1,
mise.tit.haar.cor.3.v1,
mise.tit.s8.cor.3.v1
)

names(mise.cor.3.v1)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)


mise.cor.3.v2=c(mise.ti.cor.3.v2,
mise.js.cor.3.v2,
mise.bams.cor.3.v2,
mise.nblk.cor.3.v2,
mise.sure.cor.3.v2,
mise.postmean.cor.3.v2,
mise.ebayes.cor.3.v2,
mise.ash.haar.cor.3.v2,
mise.ash.s8.cor.3.v2,
mise.ashe.haar.cor.3.v2,
mise.ashe.s8.cor.3.v2,
mise.ashe.s8.j.cor.3.v2,
mise.tie.haar.cor.3.v2,
mise.tie.s8.cor.3.v2,
mise.tieb.haar.cor.3.v2,
mise.tieb.s8.cor.3.v2,
mise.asht.haar.cor.3.v2,
mise.asht.s8.cor.3.v2,
mise.tit.haar.cor.3.v2,
mise.tit.s8.cor.3.v2
)

names(mise.cor.3.v2)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)




mise.cor.3.v3=c(mise.ti.cor.3.v3,
mise.js.cor.3.v3,
mise.bams.cor.3.v3,
mise.nblk.cor.3.v3,
mise.sure.cor.3.v3,
mise.postmean.cor.3.v3,
mise.ebayes.cor.3.v3,
mise.ash.haar.cor.3.v3,
mise.ash.s8.cor.3.v3,
mise.ashe.haar.cor.3.v3,
mise.ashe.s8.cor.3.v3,
mise.ashe.s8.j.cor.3.v3,
mise.tie.haar.cor.3.v3,
mise.tie.s8.cor.3.v3,
mise.tieb.haar.cor.3.v3,
mise.tieb.s8.cor.3.v3,
mise.asht.haar.cor.3.v3,
mise.asht.s8.cor.3.v3,
mise.tit.haar.cor.3.v3,
mise.tit.s8.cor.3.v3
)

names(mise.cor.3.v3)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)




mise.cor.3.v4=c(mise.ti.cor.3.v4,
mise.js.cor.3.v4,
mise.bams.cor.3.v4,
mise.nblk.cor.3.v4,
mise.sure.cor.3.v4,
mise.postmean.cor.3.v4,
mise.ebayes.cor.3.v4,
mise.ash.haar.cor.3.v4,
mise.ash.s8.cor.3.v4,
mise.ashe.haar.cor.3.v4,
mise.ashe.s8.cor.3.v4,
mise.ashe.s8.j.cor.3.v4,
mise.tie.haar.cor.3.v4,
mise.tie.s8.cor.3.v4,
mise.tieb.haar.cor.3.v4,
mise.tieb.s8.cor.3.v4,
mise.asht.haar.cor.3.v4,
mise.asht.s8.cor.3.v4,
mise.tit.haar.cor.3.v4,
mise.tit.s8.cor.3.v4
)

names(mise.cor.3.v4)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)





mise.cor.3.v5=c(mise.ti.cor.3.v5,
mise.js.cor.3.v5,
mise.bams.cor.3.v5,
mise.nblk.cor.3.v5,
mise.sure.cor.3.v5,
mise.postmean.cor.3.v5,
mise.ebayes.cor.3.v5,
mise.ash.haar.cor.3.v5,
mise.ash.s8.cor.3.v5,
mise.ashe.haar.cor.3.v5,
mise.ashe.s8.cor.3.v5,
mise.ashe.s8.j.cor.3.v5,
mise.tie.haar.cor.3.v5,
mise.tie.s8.cor.3.v5,
mise.tieb.haar.cor.3.v5,
mise.tieb.s8.cor.3.v5,
mise.asht.haar.cor.3.v5,
mise.asht.s8.cor.3.v5,
mise.tit.haar.cor.3.v5,
mise.tit.s8.cor.3.v5
)

names(mise.cor.3.v5)=c("TI_homo",
"blockJS_homo",
"BAMS_homo",
"NeighBlk_homo",
"Sure_homo",
"SingPostMean_homo",
"Ebayes_homo",
"ash_haar_homo",
"ash_s8_homo",
"ash_haar_est",
"ash_s8_est",
"ash_jash_est",
"TI_haar_est",
"TI_s8_est",
"TI_haar_ash_est",
"TI_s8_ash_est",
"ash_haar_true",
"ash_s8_true",
"TI_haar_true",
"TI_s8_true"
)





save.image("D:/Grad School/Spring 2013/multiscale_ash/simulation_1d_g/n_1024/res_sim_1024.RData")


round(sort(mise.sp.1.v1),digits=1)
round(sort(mise.sp.1.v2),digits=1)
round(sort(mise.sp.1.v3),digits=1)
round(sort(mise.sp.1.v4),digits=1)
round(sort(mise.sp.1.v5),digits=1)
round(sort(mise.sp.3.v1),digits=1)
round(sort(mise.sp.3.v2),digits=1)
round(sort(mise.sp.3.v3),digits=1)
round(sort(mise.sp.3.v4),digits=1)
round(sort(mise.sp.3.v5),digits=1)


round(sort(mise.bump.1.v1),digits=1)
round(sort(mise.bump.1.v2),digits=1)
round(sort(mise.bump.1.v3),digits=1)
round(sort(mise.bump.1.v4),digits=1)
round(sort(mise.bump.1.v5),digits=1)
round(sort(mise.bump.3.v1),digits=1)
round(sort(mise.bump.3.v2),digits=1)
round(sort(mise.bump.3.v3),digits=1)
round(sort(mise.bump.3.v4),digits=1)
round(sort(mise.bump.3.v5),digits=1)


round(sort(mise.blk.1.v1),digits=1)
round(sort(mise.blk.1.v2),digits=1)
round(sort(mise.blk.1.v3),digits=1)
round(sort(mise.blk.1.v4),digits=1)
round(sort(mise.blk.1.v5),digits=1)
round(sort(mise.blk.3.v1),digits=1)
round(sort(mise.blk.3.v2),digits=1)
round(sort(mise.blk.3.v3),digits=1)
round(sort(mise.blk.3.v4),digits=1)
round(sort(mise.blk.3.v5),digits=1)


round(sort(mise.ang.1.v1),digits=1)
round(sort(mise.ang.1.v2),digits=1)
round(sort(mise.ang.1.v3),digits=1)
round(sort(mise.ang.1.v4),digits=1)
round(sort(mise.ang.1.v5),digits=1)
round(sort(mise.ang.3.v1),digits=1)
round(sort(mise.ang.3.v2),digits=1)
round(sort(mise.ang.3.v3),digits=1)
round(sort(mise.ang.3.v4),digits=1)
round(sort(mise.ang.3.v5),digits=1)


round(sort(mise.dop.1.v1),digits=1)
round(sort(mise.dop.1.v2),digits=1)
round(sort(mise.dop.1.v3),digits=1)
round(sort(mise.dop.1.v4),digits=1)
round(sort(mise.dop.1.v5),digits=1)
round(sort(mise.dop.3.v1),digits=1)
round(sort(mise.dop.3.v2),digits=1)
round(sort(mise.dop.3.v3),digits=1)
round(sort(mise.dop.3.v4),digits=1)
round(sort(mise.dop.3.v5),digits=1)


round(sort(mise.blip.1.v1),digits=1)
round(sort(mise.blip.1.v2),digits=1)
round(sort(mise.blip.1.v3),digits=1)
round(sort(mise.blip.1.v4),digits=1)
round(sort(mise.blip.1.v5),digits=1)
round(sort(mise.blip.3.v1),digits=1)
round(sort(mise.blip.3.v2),digits=1)
round(sort(mise.blip.3.v3),digits=1)
round(sort(mise.blip.3.v4),digits=1)
round(sort(mise.blip.3.v5),digits=1)


round(sort(mise.cor.1.v1),digits=1)
round(sort(mise.cor.1.v2),digits=1)
round(sort(mise.cor.1.v3),digits=1)
round(sort(mise.cor.1.v4),digits=1)
round(sort(mise.cor.1.v5),digits=1)
round(sort(mise.cor.3.v1),digits=1)
round(sort(mise.cor.3.v2),digits=1)
round(sort(mise.cor.3.v3),digits=1)
round(sort(mise.cor.3.v4),digits=1)
round(sort(mise.cor.3.v5),digits=1)


############################################################
############################################################
############################################################

summary(
c(which(order(mise.sp.1.v1)==20),
which(order(mise.sp.1.v2)==20),
which(order(mise.sp.1.v3)==20),
which(order(mise.sp.1.v4)==20),
which(order(mise.sp.1.v5)==20),
which(order(mise.sp.3.v1)==20),
which(order(mise.sp.3.v2)==20),
which(order(mise.sp.3.v3)==20),
which(order(mise.sp.3.v4)==20),
which(order(mise.sp.3.v5)==20),


which(order(mise.bump.1.v1)==20),
which(order(mise.bump.1.v2)==20),
which(order(mise.bump.1.v3)==20),
which(order(mise.bump.1.v4)==20),
which(order(mise.bump.1.v5)==20),
which(order(mise.bump.3.v1)==20),
which(order(mise.bump.3.v2)==20),
which(order(mise.bump.3.v3)==20),
which(order(mise.bump.3.v4)==20),
which(order(mise.bump.3.v5)==20),


which(order(mise.blk.1.v1)==20),
which(order(mise.blk.1.v2)==20),
which(order(mise.blk.1.v3)==20),
which(order(mise.blk.1.v4)==20),
which(order(mise.blk.1.v5)==20),
which(order(mise.blk.3.v1)==20),
which(order(mise.blk.3.v2)==20),
which(order(mise.blk.3.v3)==20),
which(order(mise.blk.3.v4)==20),
which(order(mise.blk.3.v5)==20),


which(order(mise.ang.1.v1)==20),
which(order(mise.ang.1.v2)==20),
which(order(mise.ang.1.v3)==20),
which(order(mise.ang.1.v4)==20),
which(order(mise.ang.1.v5)==20),
which(order(mise.ang.3.v1)==20),
which(order(mise.ang.3.v2)==20),
which(order(mise.ang.3.v3)==20),
which(order(mise.ang.3.v4)==20),
which(order(mise.ang.3.v5)==20),


which(order(mise.dop.1.v1)==20),
which(order(mise.dop.1.v2)==20),
which(order(mise.dop.1.v3)==20),
which(order(mise.dop.1.v4)==20),
which(order(mise.dop.1.v5)==20),
which(order(mise.dop.3.v1)==20),
which(order(mise.dop.3.v2)==20),
which(order(mise.dop.3.v3)==20),
which(order(mise.dop.3.v4)==20),
which(order(mise.dop.3.v5)==20),


which(order(mise.blip.1.v1)==20),
which(order(mise.blip.1.v2)==20),
which(order(mise.blip.1.v3)==20),
which(order(mise.blip.1.v4)==20),
which(order(mise.blip.1.v5)==20),
which(order(mise.blip.3.v1)==20),
which(order(mise.blip.3.v2)==20),
which(order(mise.blip.3.v3)==20),
which(order(mise.blip.3.v4)==20),
which(order(mise.blip.3.v5)==20),


which(order(mise.cor.1.v1)==20),
which(order(mise.cor.1.v2)==20),
which(order(mise.cor.1.v3)==20),
which(order(mise.cor.1.v4)==20),
which(order(mise.cor.1.v5)==20),
which(order(mise.cor.3.v1)==20),
which(order(mise.cor.3.v2)==20),
which(order(mise.cor.3.v3)==20),
which(order(mise.cor.3.v4)==20),
which(order(mise.cor.3.v5)==20)
),breaks=100,xlab="",main="")

