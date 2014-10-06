library(haarfisz)

##generate signal
mse=function(x,y) mean((x-y)^2)
l2norm=function(x) sum(x^2)
mise=function(x,y) 10000*mean(apply(x-matrix(rep(y,times=100),nrow=100,byrow=T),1,l2norm)/l2norm(y))

n=1024
t=1:n/n

spike.f=function(x) (0.75*exp(-500*(x-0.23)^2)+1.5*exp(-2000*(x-0.33)^2)+3*exp(-8000*(x-0.47)^2)+2.25*exp(-16000*(x-0.69)^2)+0.5*exp(-32000*(x-0.83)^2))
mu.s=spike.f(t)

dop.f=function(x) sqrt(x*(1-x))*sin((2*pi*1.05)/(x+0.05))
mu.ang=dop.f(t)

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

I_1 = exp(-(abs(t-0.2)/0.01)^1.2)*(t<=0.2) + exp(-(abs(t-0.2)/0.03)^1.2)*(t>0.2);
I_2 = exp(-(abs(t-0.3)/0.01)^1.2)*(t<=0.3) + exp(-(abs(t-0.3)/0.03)^1.2)*(t>0.3);
I_3 = exp(-(abs(t-0.4)/0.01)^1.2)*(t<=0.4) + exp(-(abs(t-0.4)/0.03)^1.2)*(t>0.4);
mu.bur = 2.99/4.51804*(4*I_1+3*I_2+4.5*I_3);


pos=c(.1, .13, .15, .23, .25, .40, .44, .65, .76, .78, .81)
hgt = 2.88/5*c(4, (-5), 3, (-4), 5, (-4.2), 2.1, 4.3, (-3.1), 2.1, (-4.2))
mu.cb = rep(0,n)
for(j in 1:length(pos)){
  mu.cb = mu.cb + (1 + sign(t-pos[j]))*(hgt[j]/2)
}
mu.cb[mu.cb<0]=0

pos = c(.1, .13, .15, .23, .25, .40, .44, .65, .76, .78, .81)
hgt = 2.97/5*c(4, 5, 3, 4, 5, 4.2, 2.1, 4.3, 3.1, 5.1, 4.2)
wth = c(.005, .005, .006, .01, .01, .03, .01, .01, .005, .008, .005)
mu.b = rep(0,n)
for(j in 1:length(pos)){
  mu.b = mu.b + hgt[j]/(( 1 + (abs(t - pos[j])/wth[j]))^4)
}





mu.s.1=1*(0.01+mu.s)
mu.s.8=8/3*(3/64+mu.s)
mu.s.128=128/3*(3/128^2+mu.s)
mu.s.nz=20*(1+mu.s)

#mu.ang.1=0.01-2.99*0.5/max(mu.ang)*min(mu.ang)+2.99*0.5/max(mu.ang)*mu.ang
#mu.ang.8=1/8-7.785*0.5/max(mu.ang)*min(mu.ang)+7.785*0.5/max(mu.ang)*mu.ang
#mu.ang.128=1/128-128*0.5/max(mu.ang)*min(mu.ang)+128*0.5/max(mu.ang)*mu.ang
#mu.ang.nz=20*(1-3*min(mu.ang)+3*mu.ang)

mu.ang.1=1*(0.01+mu.ang)
mu.ang.8=8/3*(3/64+mu.ang)
mu.ang.128=1/128+128/3*mu.ang
mu.ang.nz=20*(1+mu.ang)

mu.hs.1=1*(0.01+mu.hs)
mu.hs.8=8/3*(3/64+mu.hs)
mu.hs.128=1/128+128/3*mu.hs
mu.hs.nz=20*(1+mu.hs)

mu.bur.1=1*(0.01+mu.bur)
mu.bur.8=8/3*(3/64+mu.bur)
mu.bur.128=1/128+128/3*mu.bur
mu.bur.nz=20*(1+mu.bur)

mu.cb.1=1*(0.01+mu.cb)
mu.cb.8=8/3*(3/64+mu.cb)
mu.cb.128=128/3*(3/128^2+mu.cb)
mu.cb.nz=20*(1+mu.cb)

mu.b.1=1*(0.01+mu.b)
mu.b.8=8/3*(3/64+mu.b)
mu.b.128=128/3*(3/128^2+mu.b)
mu.b.nz=20*(1+mu.b)


##generate data
set.seed(1002)
sim.m.s.1=matrix(rpois(100*n,mu.s.1),nrow=100,ncol=n,byrow=T)
set.seed(1002)
sim.m.s.8=matrix(rpois(100*n,mu.s.8),nrow=100,ncol=n,byrow=T)
set.seed(1002)
sim.m.s.128=matrix(rpois(100*n,mu.s.128),nrow=100,ncol=n,byrow=T)
set.seed(1002)
sim.m.s.nz=matrix(rpois(100*n,mu.s.nz),nrow=100,ncol=n,byrow=T)

set.seed(1002)
sim.m.ang.1=matrix(rpois(100*n,mu.ang.1),nrow=100,ncol=n,byrow=T)
set.seed(1002)
sim.m.ang.8=matrix(rpois(100*n,mu.ang.8),nrow=100,ncol=n,byrow=T)
set.seed(1002)
sim.m.ang.128=matrix(rpois(100*n,mu.ang.128),nrow=100,ncol=n,byrow=T)
set.seed(1002)
sim.m.ang.nz=matrix(rpois(100*n,mu.ang.nz),nrow=100,ncol=n,byrow=T)

set.seed(1002)
sim.m.hs.1=matrix(rpois(100*n,mu.hs.1),nrow=100,ncol=n,byrow=T)
set.seed(1002)
sim.m.hs.8=matrix(rpois(100*n,mu.hs.8),nrow=100,ncol=n,byrow=T)
set.seed(1002)
sim.m.hs.128=matrix(rpois(100*n,mu.hs.128),nrow=100,ncol=n,byrow=T)
set.seed(1002)
sim.m.hs.nz=matrix(rpois(100*n,mu.hs.nz),nrow=100,ncol=n,byrow=T)

set.seed(1002)
sim.m.bur.1=matrix(rpois(100*n,mu.bur.1),nrow=100,ncol=n,byrow=T)
set.seed(1002)
sim.m.bur.8=matrix(rpois(100*n,mu.bur.8),nrow=100,ncol=n,byrow=T)
set.seed(1002)
sim.m.bur.128=matrix(rpois(100*n,mu.bur.128),nrow=100,ncol=n,byrow=T)
set.seed(1002)
sim.m.bur.nz=matrix(rpois(100*n,mu.bur.nz),nrow=100,ncol=n,byrow=T)

set.seed(1002)
sim.m.cb.1=matrix(rpois(100*n,mu.cb.1),nrow=100,ncol=n,byrow=T)
set.seed(1002)
sim.m.cb.8=matrix(rpois(100*n,mu.cb.8),nrow=100,ncol=n,byrow=T)
set.seed(1002)
sim.m.cb.128=matrix(rpois(100*n,mu.cb.128),nrow=100,ncol=n,byrow=T)
set.seed(1002)
sim.m.cb.nz=matrix(rpois(100*n,mu.cb.nz),nrow=100,ncol=n,byrow=T)

set.seed(1002)
sim.m.b.1=matrix(rpois(100*n,mu.b.1),nrow=100,ncol=n,byrow=T)
set.seed(1002)
sim.m.b.8=matrix(rpois(100*n,mu.b.8),nrow=100,ncol=n,byrow=T)
set.seed(1002)
sim.m.b.128=matrix(rpois(100*n,mu.b.128),nrow=100,ncol=n,byrow=T)
set.seed(1002)
sim.m.b.nz=matrix(rpois(100*n,mu.b.nz),nrow=100,ncol=n,byrow=T)



#library(multiseq)
source("~/ashwave/Rcode/binnorm.R")
#source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/PoissonBinomial.funcs.R"))  
#source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/deltamethod.R"))  


source("~/ashwave/Rcode/binshrink.R")



est.ash1.s.1=apply(sim.m.s.1,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=2)
est.ash2.s.1=apply(sim.m.s.1,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,mixsd=c(0,1.815))
est.ash3.s.1=apply(sim.m.s.1,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,nullcheck=FALSE,VB=TRUE)
est.ash4.s.1=apply(sim.m.s.1,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=0)
est.ash5.s.1=apply(sim.m.s.1,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=2)
est.ash6.s.1=apply(sim.m.s.1,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,mixsd=c(0,1.815))
est.ash7.s.1=apply(sim.m.s.1,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,nullcheck=FALSE,VB=TRUE)
est.ash8.s.1=apply(sim.m.s.1,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=0)

est.ash1.ang.1=apply(sim.m.ang.1,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=2)
est.ash2.ang.1=apply(sim.m.ang.1,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,mixsd=c(0,1.815))
est.ash3.ang.1=apply(sim.m.ang.1,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,nullcheck=FALSE,VB=TRUE)
est.ash4.ang.1=apply(sim.m.ang.1,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=0)
est.ash5.ang.1=apply(sim.m.ang.1,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=2)
est.ash6.ang.1=apply(sim.m.ang.1,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,mixsd=c(0,1.815))
est.ash7.ang.1=apply(sim.m.ang.1,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,nullcheck=FALSE,VB=TRUE)
est.ash8.ang.1=apply(sim.m.ang.1,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=0)


est.ash1.hs.1=apply(sim.m.hs.1,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=2)
est.ash2.hs.1=apply(sim.m.hs.1,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,mixsd=c(0,1.815))
est.ash3.hs.1=apply(sim.m.hs.1,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,nullcheck=FALSE,VB=TRUE)
est.ash4.hs.1=apply(sim.m.hs.1,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=0)
est.ash5.hs.1=apply(sim.m.hs.1,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=2)
est.ash6.hs.1=apply(sim.m.hs.1,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,mixsd=c(0,1.815))
est.ash7.hs.1=apply(sim.m.hs.1,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,nullcheck=FALSE,VB=TRUE)
est.ash8.hs.1=apply(sim.m.hs.1,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=0)


est.ash1.bur.1=apply(sim.m.bur.1,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=2)
est.ash2.bur.1=apply(sim.m.bur.1,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,mixsd=c(0,1.815))
est.ash3.bur.1=apply(sim.m.bur.1,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,nullcheck=FALSE,VB=TRUE)
est.ash4.bur.1=apply(sim.m.bur.1,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=0)
est.ash5.bur.1=apply(sim.m.bur.1,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=2)
est.ash6.bur.1=apply(sim.m.bur.1,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,mixsd=c(0,1.815))
est.ash7.bur.1=apply(sim.m.bur.1,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,nullcheck=FALSE,VB=TRUE)
est.ash8.bur.1=apply(sim.m.bur.1,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=0)


est.ash1.cb.1=apply(sim.m.cb.1,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=2)
est.ash2.cb.1=apply(sim.m.cb.1,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,mixsd=c(0,1.815))
est.ash3.cb.1=apply(sim.m.cb.1,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,nullcheck=FALSE,VB=TRUE)
est.ash4.cb.1=apply(sim.m.cb.1,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=0)
est.ash5.cb.1=apply(sim.m.cb.1,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=2)
est.ash6.cb.1=apply(sim.m.cb.1,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,mixsd=c(0,1.815))
est.ash7.cb.1=apply(sim.m.cb.1,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,nullcheck=FALSE,VB=TRUE)
est.ash8.cb.1=apply(sim.m.cb.1,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=0)


est.ash1.b.1=apply(sim.m.b.1,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=2)
est.ash2.b.1=apply(sim.m.b.1,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,mixsd=c(0,1.815))
est.ash3.b.1=apply(sim.m.b.1,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,nullcheck=FALSE,VB=TRUE)
est.ash4.b.1=apply(sim.m.b.1,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=0)
est.ash5.b.1=apply(sim.m.b.1,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=2)
est.ash6.b.1=apply(sim.m.b.1,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,mixsd=c(0,1.815))
est.ash7.b.1=apply(sim.m.b.1,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,nullcheck=FALSE,VB=TRUE)
est.ash8.b.1=apply(sim.m.b.1,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=0)




est.ash1.s.8=apply(sim.m.s.8,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=2)
est.ash2.s.8=apply(sim.m.s.8,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,mixsd=c(0,1.815))
est.ash3.s.8=apply(sim.m.s.8,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,nullcheck=FALSE,VB=TRUE)
est.ash4.s.8=apply(sim.m.s.8,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=0)
est.ash5.s.8=apply(sim.m.s.8,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=2)
est.ash6.s.8=apply(sim.m.s.8,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,mixsd=c(0,1.815))
est.ash7.s.8=apply(sim.m.s.8,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,nullcheck=FALSE,VB=TRUE)
est.ash8.s.8=apply(sim.m.s.8,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=0)

est.ash1.ang.8=apply(sim.m.ang.8,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=2)
est.ash2.ang.8=apply(sim.m.ang.8,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,mixsd=c(0,1.815))
est.ash3.ang.8=apply(sim.m.ang.8,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,nullcheck=FALSE,VB=TRUE)
est.ash4.ang.8=apply(sim.m.ang.8,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=0)
est.ash5.ang.8=apply(sim.m.ang.8,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=2)
est.ash6.ang.8=apply(sim.m.ang.8,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,mixsd=c(0,1.815))
est.ash7.ang.8=apply(sim.m.ang.8,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,nullcheck=FALSE,VB=TRUE)
est.ash8.ang.8=apply(sim.m.ang.8,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=0)


est.ash1.hs.8=apply(sim.m.hs.8,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=2)
est.ash2.hs.8=apply(sim.m.hs.8,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,mixsd=c(0,1.815))
est.ash3.hs.8=apply(sim.m.hs.8,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,nullcheck=FALSE,VB=TRUE)
est.ash4.hs.8=apply(sim.m.hs.8,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=0)
est.ash5.hs.8=apply(sim.m.hs.8,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=2)
est.ash6.hs.8=apply(sim.m.hs.8,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,mixsd=c(0,1.815))
est.ash7.hs.8=apply(sim.m.hs.8,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,nullcheck=FALSE,VB=TRUE)
est.ash8.hs.8=apply(sim.m.hs.8,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=0)


est.ash1.bur.8=apply(sim.m.bur.8,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=2)
est.ash2.bur.8=apply(sim.m.bur.8,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,mixsd=c(0,1.815))
est.ash3.bur.8=apply(sim.m.bur.8,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,nullcheck=FALSE,VB=TRUE)
est.ash4.bur.8=apply(sim.m.bur.8,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=0)
est.ash5.bur.8=apply(sim.m.bur.8,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=2)
est.ash6.bur.8=apply(sim.m.bur.8,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,mixsd=c(0,1.815))
est.ash7.bur.8=apply(sim.m.bur.8,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,nullcheck=FALSE,VB=TRUE)
est.ash8.bur.8=apply(sim.m.bur.8,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=0)


est.ash1.cb.8=apply(sim.m.cb.8,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=2)
est.ash2.cb.8=apply(sim.m.cb.8,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,mixsd=c(0,1.815))
est.ash3.cb.8=apply(sim.m.cb.8,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,nullcheck=FALSE,VB=TRUE)
est.ash4.cb.8=apply(sim.m.cb.8,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=0)
est.ash5.cb.8=apply(sim.m.cb.8,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=2)
est.ash6.cb.8=apply(sim.m.cb.8,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,mixsd=c(0,1.815))
est.ash7.cb.8=apply(sim.m.cb.8,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,nullcheck=FALSE,VB=TRUE)
est.ash8.cb.8=apply(sim.m.cb.8,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=0)


est.ash1.b.8=apply(sim.m.b.8,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=2)
est.ash2.b.8=apply(sim.m.b.8,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,mixsd=c(0,1.815))
est.ash3.b.8=apply(sim.m.b.8,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,nullcheck=FALSE,VB=TRUE)
est.ash4.b.8=apply(sim.m.b.8,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=0)
est.ash5.b.8=apply(sim.m.b.8,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=2)
est.ash6.b.8=apply(sim.m.b.8,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,mixsd=c(0,1.815))
est.ash7.b.8=apply(sim.m.b.8,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,nullcheck=FALSE,VB=TRUE)
est.ash8.b.8=apply(sim.m.b.8,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=0)




est.ash1.s.128=apply(sim.m.s.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=2)
est.ash2.s.128=apply(sim.m.s.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,mixsd=c(0,1.815))
est.ash3.s.128=apply(sim.m.s.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,nullcheck=FALSE,VB=TRUE)
est.ash4.s.128=apply(sim.m.s.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=0)
est.ash5.s.128=apply(sim.m.s.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=2)
est.ash6.s.128=apply(sim.m.s.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,mixsd=c(0,1.815))
est.ash7.s.128=apply(sim.m.s.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,nullcheck=FALSE,VB=TRUE)
est.ash8.s.128=apply(sim.m.s.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=0)

est.ash1.ang.128=apply(sim.m.ang.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=2)
est.ash2.ang.128=apply(sim.m.ang.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,mixsd=c(0,1.815))
est.ash3.ang.128=apply(sim.m.ang.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,nullcheck=FALSE,VB=TRUE)
est.ash4.ang.128=apply(sim.m.ang.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=0)
est.ash5.ang.128=apply(sim.m.ang.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=2)
est.ash6.ang.128=apply(sim.m.ang.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,mixsd=c(0,1.815))
est.ash7.ang.128=apply(sim.m.ang.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,nullcheck=FALSE,VB=TRUE)
est.ash8.ang.128=apply(sim.m.ang.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=0)


est.ash1.hs.128=apply(sim.m.hs.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=2)
est.ash2.hs.128=apply(sim.m.hs.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,mixsd=c(0,1.815))
est.ash3.hs.128=apply(sim.m.hs.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,nullcheck=FALSE,VB=TRUE)
est.ash4.hs.128=apply(sim.m.hs.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=0)
est.ash5.hs.128=apply(sim.m.hs.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=2)
est.ash6.hs.128=apply(sim.m.hs.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,mixsd=c(0,1.815))
est.ash7.hs.128=apply(sim.m.hs.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,nullcheck=FALSE,VB=TRUE)
est.ash8.hs.128=apply(sim.m.hs.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=0)


est.ash1.bur.128=apply(sim.m.bur.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=2)
est.ash2.bur.128=apply(sim.m.bur.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,mixsd=c(0,1.815))
est.ash3.bur.128=apply(sim.m.bur.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,nullcheck=FALSE,VB=TRUE)
est.ash4.bur.128=apply(sim.m.bur.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=0)
est.ash5.bur.128=apply(sim.m.bur.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=2)
est.ash6.bur.128=apply(sim.m.bur.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,mixsd=c(0,1.815))
est.ash7.bur.128=apply(sim.m.bur.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,nullcheck=FALSE,VB=TRUE)
est.ash8.bur.128=apply(sim.m.bur.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=0)


est.ash1.cb.128=apply(sim.m.cb.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=2)
est.ash2.cb.128=apply(sim.m.cb.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,mixsd=c(0,1.815))
est.ash3.cb.128=apply(sim.m.cb.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,nullcheck=FALSE,VB=TRUE)
est.ash4.cb.128=apply(sim.m.cb.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=0)
est.ash5.cb.128=apply(sim.m.cb.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=2)
est.ash6.cb.128=apply(sim.m.cb.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,mixsd=c(0,1.815))
est.ash7.cb.128=apply(sim.m.cb.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,nullcheck=FALSE,VB=TRUE)
est.ash8.cb.128=apply(sim.m.cb.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=0)


est.ash1.b.128=apply(sim.m.b.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=2)
est.ash2.b.128=apply(sim.m.b.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,mixsd=c(0,1.815))
est.ash3.b.128=apply(sim.m.b.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,nullcheck=FALSE,VB=TRUE)
est.ash4.b.128=apply(sim.m.b.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=0)
est.ash5.b.128=apply(sim.m.b.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=2)
est.ash6.b.128=apply(sim.m.b.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,mixsd=c(0,1.815))
est.ash7.b.128=apply(sim.m.b.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,nullcheck=FALSE,VB=TRUE)
est.ash8.b.128=apply(sim.m.b.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=0)



est.ash1.s.128=apply(sim.m.s.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=2)
est.ash2.s.128=apply(sim.m.s.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,mixsd=c(0,1.815))
est.ash3.s.128=apply(sim.m.s.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,nullcheck=FALSE,VB=TRUE)
est.ash4.s.128=apply(sim.m.s.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=0)
est.ash5.s.128=apply(sim.m.s.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=2)
est.ash6.s.128=apply(sim.m.s.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,mixsd=c(0,1.815))
est.ash7.s.128=apply(sim.m.s.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,nullcheck=FALSE,VB=TRUE)
est.ash8.s.128=apply(sim.m.s.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=0)

est.ash1.ang.128=apply(sim.m.ang.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=2)
est.ash2.ang.128=apply(sim.m.ang.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,mixsd=c(0,1.815))
est.ash3.ang.128=apply(sim.m.ang.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,nullcheck=FALSE,VB=TRUE)
est.ash4.ang.128=apply(sim.m.ang.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=0)
est.ash5.ang.128=apply(sim.m.ang.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=2)
est.ash6.ang.128=apply(sim.m.ang.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,mixsd=c(0,1.815))
est.ash7.ang.128=apply(sim.m.ang.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,nullcheck=FALSE,VB=TRUE)
est.ash8.ang.128=apply(sim.m.ang.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=0)


est.ash1.hs.128=apply(sim.m.hs.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=2)
est.ash2.hs.128=apply(sim.m.hs.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,mixsd=c(0,1.815))
est.ash3.hs.128=apply(sim.m.hs.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,nullcheck=FALSE,VB=TRUE)
est.ash4.hs.128=apply(sim.m.hs.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=0)
est.ash5.hs.128=apply(sim.m.hs.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=2)
est.ash6.hs.128=apply(sim.m.hs.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,mixsd=c(0,1.815))
est.ash7.hs.128=apply(sim.m.hs.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,nullcheck=FALSE,VB=TRUE)
est.ash8.hs.128=apply(sim.m.hs.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=0)


est.ash1.bur.128=apply(sim.m.bur.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=2)
est.ash2.bur.128=apply(sim.m.bur.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,mixsd=c(0,1.815))
est.ash3.bur.128=apply(sim.m.bur.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,nullcheck=FALSE,VB=TRUE)
est.ash4.bur.128=apply(sim.m.bur.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=0)
est.ash5.bur.128=apply(sim.m.bur.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=2)
est.ash6.bur.128=apply(sim.m.bur.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,mixsd=c(0,1.815))
est.ash7.bur.128=apply(sim.m.bur.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,nullcheck=FALSE,VB=TRUE)
est.ash8.bur.128=apply(sim.m.bur.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=0)


est.ash1.cb.128=apply(sim.m.cb.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=2)
est.ash2.cb.128=apply(sim.m.cb.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,mixsd=c(0,1.815))
est.ash3.cb.128=apply(sim.m.cb.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,nullcheck=FALSE,VB=TRUE)
est.ash4.cb.128=apply(sim.m.cb.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=0)
est.ash5.cb.128=apply(sim.m.cb.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=2)
est.ash6.cb.128=apply(sim.m.cb.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,mixsd=c(0,1.815))
est.ash7.cb.128=apply(sim.m.cb.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,nullcheck=FALSE,VB=TRUE)
est.ash8.cb.128=apply(sim.m.cb.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=0)


est.ash1.b.128=apply(sim.m.b.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=2)
est.ash2.b.128=apply(sim.m.b.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,mixsd=c(0,1.815))
est.ash3.b.128=apply(sim.m.b.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,nullcheck=FALSE,VB=TRUE)
est.ash4.b.128=apply(sim.m.b.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=0)
est.ash5.b.128=apply(sim.m.b.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=2)
est.ash6.b.128=apply(sim.m.b.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,mixsd=c(0,1.815))
est.ash7.b.128=apply(sim.m.b.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,nullcheck=FALSE,VB=TRUE)
est.ash8.b.128=apply(sim.m.b.128,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=0)

est.ash1.s.nz=apply(sim.m.s.nz,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=2)
est.ash2.s.nz=apply(sim.m.s.nz,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,mixsd=c(0,1.815))
est.ash3.s.nz=apply(sim.m.s.nz,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,nullcheck=FALSE,VB=TRUE)
est.ash4.s.nz=apply(sim.m.s.nz,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=0)
est.ash5.s.nz=apply(sim.m.s.nz,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=2)
est.ash6.s.nz=apply(sim.m.s.nz,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,mixsd=c(0,1.815))
est.ash7.s.nz=apply(sim.m.s.nz,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,nullcheck=FALSE,VB=TRUE)
est.ash8.s.nz=apply(sim.m.s.nz,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=0)

est.ash1.ang.nz=apply(sim.m.ang.nz,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=2)
est.ash2.ang.nz=apply(sim.m.ang.nz,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,mixsd=c(0,1.815))
est.ash3.ang.nz=apply(sim.m.ang.nz,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,nullcheck=FALSE,VB=TRUE)
est.ash4.ang.nz=apply(sim.m.ang.nz,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=0)
est.ash5.ang.nz=apply(sim.m.ang.nz,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=2)
est.ash6.ang.nz=apply(sim.m.ang.nz,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,mixsd=c(0,1.815))
est.ash7.ang.nz=apply(sim.m.ang.nz,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,nullcheck=FALSE,VB=TRUE)
est.ash8.ang.nz=apply(sim.m.ang.nz,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=0)


est.ash1.hs.nz=apply(sim.m.hs.nz,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=2)
est.ash2.hs.nz=apply(sim.m.hs.nz,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,mixsd=c(0,1.815))
est.ash3.hs.nz=apply(sim.m.hs.nz,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,nullcheck=FALSE,VB=TRUE)
est.ash4.hs.nz=apply(sim.m.hs.nz,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=0)
est.ash5.hs.nz=apply(sim.m.hs.nz,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=2)
est.ash6.hs.nz=apply(sim.m.hs.nz,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,mixsd=c(0,1.815))
est.ash7.hs.nz=apply(sim.m.hs.nz,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,nullcheck=FALSE,VB=TRUE)
est.ash8.hs.nz=apply(sim.m.hs.nz,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=0)


est.ash1.bur.nz=apply(sim.m.bur.nz,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=2)
est.ash2.bur.nz=apply(sim.m.bur.nz,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,mixsd=c(0,1.815))
est.ash3.bur.nz=apply(sim.m.bur.nz,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,nullcheck=FALSE,VB=TRUE)
est.ash4.bur.nz=apply(sim.m.bur.nz,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=0)
est.ash5.bur.nz=apply(sim.m.bur.nz,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=2)
est.ash6.bur.nz=apply(sim.m.bur.nz,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,mixsd=c(0,1.815))
est.ash7.bur.nz=apply(sim.m.bur.nz,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,nullcheck=FALSE,VB=TRUE)
est.ash8.bur.nz=apply(sim.m.bur.nz,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=0)


est.ash1.cb.nz=apply(sim.m.cb.nz,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=2)
est.ash2.cb.nz=apply(sim.m.cb.nz,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,mixsd=c(0,1.815))
est.ash3.cb.nz=apply(sim.m.cb.nz,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,nullcheck=FALSE,VB=TRUE)
est.ash4.cb.nz=apply(sim.m.cb.nz,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=0)
est.ash5.cb.nz=apply(sim.m.cb.nz,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=2)
est.ash6.cb.nz=apply(sim.m.cb.nz,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,mixsd=c(0,1.815))
est.ash7.cb.nz=apply(sim.m.cb.nz,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,nullcheck=FALSE,VB=TRUE)
est.ash8.cb.nz=apply(sim.m.cb.nz,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=0)


est.ash1.b.nz=apply(sim.m.b.nz,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=2)
est.ash2.b.nz=apply(sim.m.b.nz,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,mixsd=c(0,1.815))
est.ash3.b.nz=apply(sim.m.b.nz,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,nullcheck=FALSE,VB=TRUE)
est.ash4.b.nz=apply(sim.m.b.nz,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=FALSE,gridmult=0)
est.ash5.b.nz=apply(sim.m.b.nz,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=2)
est.ash6.b.nz=apply(sim.m.b.nz,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,mixsd=c(0,1.815))
est.ash7.b.nz=apply(sim.m.b.nz,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,nullcheck=FALSE,VB=TRUE)
est.ash8.b.nz=apply(sim.m.b.nz,1,binnorm.smooth,log=FALSE,lev=0,return.est=TRUE,pseudocounts=0.5,all=TRUE,gridmult=0)





q=c(1/5*2^(seq(-1,5,2)),12.8,51.2)
est.BMSMshrink1.s.1=apply(sim.m.s.1,1,binshrink,1,1,return.est=TRUE)
est.BMSMshrink2.s.1=apply(sim.m.s.1,1,binshrink,q,2,return.est=TRUE)

est.BMSMshrink1.ang.1=apply(sim.m.ang.1,1,binshrink,1,1,return.est=TRUE)
est.BMSMshrink2.ang.1=apply(sim.m.ang.1,1,binshrink,q,2,return.est=TRUE)

est.BMSMshrink1.hs.1=apply(sim.m.hs.1,1,binshrink,1,1,return.est=TRUE)
est.BMSMshrink2.hs.1=apply(sim.m.hs.1,1,binshrink,q,2,return.est=TRUE)

est.BMSMshrink1.bur.1=apply(sim.m.bur.1,1,binshrink,1,1,return.est=TRUE)
est.BMSMshrink2.bur.1=apply(sim.m.bur.1,1,binshrink,q,2,return.est=TRUE)

est.BMSMshrink1.cb.1=apply(sim.m.cb.1,1,binshrink,1,1,return.est=TRUE)
est.BMSMshrink2.cb.1=apply(sim.m.cb.1,1,binshrink,q,2,return.est=TRUE)

est.BMSMshrink1.b.1=apply(sim.m.b.1,1,binshrink,1,1,return.est=TRUE)
est.BMSMshrink2.b.1=apply(sim.m.b.1,1,binshrink,q,2,return.est=TRUE)


est.BMSMshrink1.s.8=apply(sim.m.s.8,1,binshrink,1,1,return.est=TRUE)
est.BMSMshrink2.s.8=apply(sim.m.s.8,1,binshrink,q,2,return.est=TRUE)

est.BMSMshrink1.ang.8=apply(sim.m.ang.8,1,binshrink,1,1,return.est=TRUE)
est.BMSMshrink2.ang.8=apply(sim.m.ang.8,1,binshrink,q,2,return.est=TRUE)

est.BMSMshrink1.hs.8=apply(sim.m.hs.8,1,binshrink,1,1,return.est=TRUE)
est.BMSMshrink2.hs.8=apply(sim.m.hs.8,1,binshrink,q,2,return.est=TRUE)

est.BMSMshrink1.bur.8=apply(sim.m.bur.8,1,binshrink,1,1,return.est=TRUE)
est.BMSMshrink2.bur.8=apply(sim.m.bur.8,1,binshrink,q,2,return.est=TRUE)

est.BMSMshrink1.cb.8=apply(sim.m.cb.8,1,binshrink,1,1,return.est=TRUE)
est.BMSMshrink2.cb.8=apply(sim.m.cb.8,1,binshrink,q,2,return.est=TRUE)

est.BMSMshrink1.b.8=apply(sim.m.b.8,1,binshrink,1,1,return.est=TRUE)
est.BMSMshrink2.b.8=apply(sim.m.b.8,1,binshrink,q,2,return.est=TRUE)


est.BMSMshrink1.s.128=apply(sim.m.s.128,1,binshrink,1,1,return.est=TRUE)
est.BMSMshrink2.s.128=apply(sim.m.s.128,1,binshrink,q,2,return.est=TRUE)

est.BMSMshrink1.ang.128=apply(sim.m.ang.128,1,binshrink,1,1,return.est=TRUE)
est.BMSMshrink2.ang.128=apply(sim.m.ang.128,1,binshrink,q,2,return.est=TRUE)

est.BMSMshrink1.hs.128=apply(sim.m.hs.128,1,binshrink,1,1,return.est=TRUE)
est.BMSMshrink2.hs.128=apply(sim.m.hs.128,1,binshrink,q,2,return.est=TRUE)

est.BMSMshrink1.bur.128=apply(sim.m.bur.128,1,binshrink,1,1,return.est=TRUE)
est.BMSMshrink2.bur.128=apply(sim.m.bur.128,1,binshrink,q,2,return.est=TRUE)

est.BMSMshrink1.cb.128=apply(sim.m.cb.128,1,binshrink,1,1,return.est=TRUE)
est.BMSMshrink2.cb.128=apply(sim.m.cb.128,1,binshrink,q,2,return.est=TRUE)

est.BMSMshrink1.b.128=apply(sim.m.b.128,1,binshrink,1,1,return.est=TRUE)
est.BMSMshrink2.b.128=apply(sim.m.b.128,1,binshrink,q,2,return.est=TRUE)



est.BMSMshrink1.s.nz=apply(sim.m.s.nz,1,binshrink,1,1,return.est=TRUE)
est.BMSMshrink2.s.nz=apply(sim.m.s.nz,1,binshrink,q,2,return.est=TRUE)

est.BMSMshrink1.ang.nz=apply(sim.m.ang.nz,1,binshrink,1,1,return.est=TRUE)
est.BMSMshrink2.ang.nz=apply(sim.m.ang.nz,1,binshrink,q,2,return.est=TRUE)

est.BMSMshrink1.hs.nz=apply(sim.m.hs.nz,1,binshrink,1,1,return.est=TRUE)
est.BMSMshrink2.hs.nz=apply(sim.m.hs.nz,1,binshrink,q,2,return.est=TRUE)

est.BMSMshrink1.bur.nz=apply(sim.m.bur.nz,1,binshrink,1,1,return.est=TRUE)
est.BMSMshrink2.bur.nz=apply(sim.m.bur.nz,1,binshrink,q,2,return.est=TRUE)

est.BMSMshrink1.cb.nz=apply(sim.m.cb.nz,1,binshrink,1,1,return.est=TRUE)
est.BMSMshrink2.cb.nz=apply(sim.m.cb.nz,1,binshrink,q,2,return.est=TRUE)

est.BMSMshrink1.b.nz=apply(sim.m.b.nz,1,binshrink,1,1,return.est=TRUE)
est.BMSMshrink2.b.nz=apply(sim.m.b.nz,1,binshrink,q,2,return.est=TRUE)



hf.la10.ti2=function (x, filter.number = 10, family = "DaubLeAsymm", min.level = 2, noise.level = 1) 
{
    TT <- length(x)
    thresh <- noise.level * sqrt(2 * log(TT))
    x.w <- wd(x, filter.number, family, type = "station")
    x.w.t <- threshold(x.w, levels = (min.level):(x.w$nlevels - 1), policy = "manual", value = thresh, type = "hard")
    x.w.t.r <- AvBasis(convert(x.w.t))
    return(x.w.t.r)
}


hf.la10.ti3=function (x, filter.number = 10, family = "DaubLeAsymm", min.level = 3, noise.level = 1) 
{
    TT <- length(x)
    thresh <- noise.level * sqrt(2 * log(TT))
    x.w <- wd(x, filter.number, family, type = "station")
    x.w.t <- threshold(x.w, levels = (min.level):(x.w$nlevels - 1), policy = "manual", value = thresh, type = "hard")
    x.w.t.r <- AvBasis(convert(x.w.t))
    return(x.w.t.r)
}

hf.la10.ti4=function (x, filter.number = 10, family = "DaubLeAsymm", min.level = 4, noise.level = 1) 
{
    TT <- length(x)
    thresh <- noise.level * sqrt(2 * log(TT))
    x.w <- wd(x, filter.number, family, type = "station")
    x.w.t <- threshold(x.w, levels = (min.level):(x.w$nlevels - 1), policy = "manual", value = thresh, type = "hard")
    x.w.t.r <- AvBasis(convert(x.w.t))
    return(x.w.t.r)
}

hf.la10.ti5=function (x, filter.number = 10, family = "DaubLeAsymm", min.level = 5, noise.level = 1) 
{
    TT <- length(x)
    thresh <- noise.level * sqrt(2 * log(TT))
    x.w <- wd(x, filter.number, family, type = "station")
    x.w.t <- threshold(x.w, levels = (min.level):(x.w$nlevels - 1), policy = "manual", value = thresh, type = "hard")
    x.w.t.r <- AvBasis(convert(x.w.t))
    return(x.w.t.r)
}

hf.la10.ti6=function (x, filter.number = 10, family = "DaubLeAsymm", min.level = 6, noise.level = 1) 
{
    TT <- length(x)
    thresh <- noise.level * sqrt(2 * log(TT))
    x.w <- wd(x, filter.number, family, type = "station")
    x.w.t <- threshold(x.w, levels = (min.level):(x.w$nlevels - 1), policy = "manual", value = thresh, type = "hard")
    x.w.t.r <- AvBasis(convert(x.w.t))
    return(x.w.t.r)
}

hf.la10.ti7=function (x, filter.number = 10, family = "DaubLeAsymm", min.level = 7, noise.level = 1) 
{
    TT <- length(x)
    thresh <- noise.level * sqrt(2 * log(TT))
    x.w <- wd(x, filter.number, family, type = "station")
    x.w.t <- threshold(x.w, levels = (min.level):(x.w$nlevels - 1), policy = "manual", value = thresh, type = "hard")
    x.w.t.r <- AvBasis(convert(x.w.t))
    return(x.w.t.r)
}


est.hf.ti.r.2.s.1=apply(sim.m.s.1,1,denoise.poisson,meth.1=hf.la10.ti2,cs.1=50,hybrid=F)
est.hf.ti.r.3.s.1=apply(sim.m.s.1,1,denoise.poisson,meth.1=hf.la10.ti3,cs.1=50,hybrid=F)
est.hf.ti.r.4.s.1=apply(sim.m.s.1,1,denoise.poisson,meth.1=hf.la10.ti4,cs.1=50,hybrid=F)
est.hf.ti.r.5.s.1=apply(sim.m.s.1,1,denoise.poisson,meth.1=hf.la10.ti5,cs.1=50,hybrid=F)
est.hf.ti.r.6.s.1=apply(sim.m.s.1,1,denoise.poisson,meth.1=hf.la10.ti6,cs.1=50,hybrid=F)
est.hf.ti.r.7.s.1=apply(sim.m.s.1,1,denoise.poisson,meth.1=hf.la10.ti7,cs.1=50,hybrid=F)

est.hf.ti.r.2.ang.1=apply(sim.m.ang.1,1,denoise.poisson,meth.1=hf.la10.ti2,cs.1=50,hybrid=F)
est.hf.ti.r.3.ang.1=apply(sim.m.ang.1,1,denoise.poisson,meth.1=hf.la10.ti3,cs.1=50,hybrid=F)
est.hf.ti.r.4.ang.1=apply(sim.m.ang.1,1,denoise.poisson,meth.1=hf.la10.ti4,cs.1=50,hybrid=F)
est.hf.ti.r.5.ang.1=apply(sim.m.ang.1,1,denoise.poisson,meth.1=hf.la10.ti5,cs.1=50,hybrid=F)
est.hf.ti.r.6.ang.1=apply(sim.m.ang.1,1,denoise.poisson,meth.1=hf.la10.ti6,cs.1=50,hybrid=F)
est.hf.ti.r.7.ang.1=apply(sim.m.ang.1,1,denoise.poisson,meth.1=hf.la10.ti7,cs.1=50,hybrid=F)

est.hf.ti.r.2.hs.1=apply(sim.m.hs.1,1,denoise.poisson,meth.1=hf.la10.ti2,cs.1=50,hybrid=F)
est.hf.ti.r.3.hs.1=apply(sim.m.hs.1,1,denoise.poisson,meth.1=hf.la10.ti3,cs.1=50,hybrid=F)
est.hf.ti.r.4.hs.1=apply(sim.m.hs.1,1,denoise.poisson,meth.1=hf.la10.ti4,cs.1=50,hybrid=F)
est.hf.ti.r.5.hs.1=apply(sim.m.hs.1,1,denoise.poisson,meth.1=hf.la10.ti5,cs.1=50,hybrid=F)
est.hf.ti.r.6.hs.1=apply(sim.m.hs.1,1,denoise.poisson,meth.1=hf.la10.ti6,cs.1=50,hybrid=F)
est.hf.ti.r.7.hs.1=apply(sim.m.hs.1,1,denoise.poisson,meth.1=hf.la10.ti7,cs.1=50,hybrid=F)

est.hf.ti.r.2.bur.1=apply(sim.m.bur.1,1,denoise.poisson,meth.1=hf.la10.ti2,cs.1=50,hybrid=F)
est.hf.ti.r.3.bur.1=apply(sim.m.bur.1,1,denoise.poisson,meth.1=hf.la10.ti3,cs.1=50,hybrid=F)
est.hf.ti.r.4.bur.1=apply(sim.m.bur.1,1,denoise.poisson,meth.1=hf.la10.ti4,cs.1=50,hybrid=F)
est.hf.ti.r.5.bur.1=apply(sim.m.bur.1,1,denoise.poisson,meth.1=hf.la10.ti5,cs.1=50,hybrid=F)
est.hf.ti.r.6.bur.1=apply(sim.m.bur.1,1,denoise.poisson,meth.1=hf.la10.ti6,cs.1=50,hybrid=F)
est.hf.ti.r.7.bur.1=apply(sim.m.bur.1,1,denoise.poisson,meth.1=hf.la10.ti7,cs.1=50,hybrid=F)

est.hf.ti.r.2.cb.1=apply(sim.m.cb.1,1,denoise.poisson,meth.1=hf.la10.ti2,cs.1=50,hybrid=F)
est.hf.ti.r.3.cb.1=apply(sim.m.cb.1,1,denoise.poisson,meth.1=hf.la10.ti3,cs.1=50,hybrid=F)
est.hf.ti.r.4.cb.1=apply(sim.m.cb.1,1,denoise.poisson,meth.1=hf.la10.ti4,cs.1=50,hybrid=F)
est.hf.ti.r.5.cb.1=apply(sim.m.cb.1,1,denoise.poisson,meth.1=hf.la10.ti5,cs.1=50,hybrid=F)
est.hf.ti.r.6.cb.1=apply(sim.m.cb.1,1,denoise.poisson,meth.1=hf.la10.ti6,cs.1=50,hybrid=F)
est.hf.ti.r.7.cb.1=apply(sim.m.cb.1,1,denoise.poisson,meth.1=hf.la10.ti7,cs.1=50,hybrid=F)

est.hf.ti.r.2.b.1=apply(sim.m.b.1,1,denoise.poisson,meth.1=hf.la10.ti2,cs.1=50,hybrid=F)
est.hf.ti.r.3.b.1=apply(sim.m.b.1,1,denoise.poisson,meth.1=hf.la10.ti3,cs.1=50,hybrid=F)
est.hf.ti.r.4.b.1=apply(sim.m.b.1,1,denoise.poisson,meth.1=hf.la10.ti4,cs.1=50,hybrid=F)
est.hf.ti.r.5.b.1=apply(sim.m.b.1,1,denoise.poisson,meth.1=hf.la10.ti5,cs.1=50,hybrid=F)
est.hf.ti.r.6.b.1=apply(sim.m.b.1,1,denoise.poisson,meth.1=hf.la10.ti6,cs.1=50,hybrid=F)
est.hf.ti.r.7.b.1=apply(sim.m.b.1,1,denoise.poisson,meth.1=hf.la10.ti7,cs.1=50,hybrid=F)


est.hf.ti.r.2.s.8=apply(sim.m.s.8,1,denoise.poisson,meth.1=hf.la10.ti2,cs.1=50,hybrid=F)
est.hf.ti.r.3.s.8=apply(sim.m.s.8,1,denoise.poisson,meth.1=hf.la10.ti3,cs.1=50,hybrid=F)
est.hf.ti.r.4.s.8=apply(sim.m.s.8,1,denoise.poisson,meth.1=hf.la10.ti4,cs.1=50,hybrid=F)
est.hf.ti.r.5.s.8=apply(sim.m.s.8,1,denoise.poisson,meth.1=hf.la10.ti5,cs.1=50,hybrid=F)
est.hf.ti.r.6.s.8=apply(sim.m.s.8,1,denoise.poisson,meth.1=hf.la10.ti6,cs.1=50,hybrid=F)
est.hf.ti.r.7.s.8=apply(sim.m.s.8,1,denoise.poisson,meth.1=hf.la10.ti7,cs.1=50,hybrid=F)

est.hf.ti.r.2.ang.8=apply(sim.m.ang.8,1,denoise.poisson,meth.1=hf.la10.ti2,cs.1=50,hybrid=F)
est.hf.ti.r.3.ang.8=apply(sim.m.ang.8,1,denoise.poisson,meth.1=hf.la10.ti3,cs.1=50,hybrid=F)
est.hf.ti.r.4.ang.8=apply(sim.m.ang.8,1,denoise.poisson,meth.1=hf.la10.ti4,cs.1=50,hybrid=F)
est.hf.ti.r.5.ang.8=apply(sim.m.ang.8,1,denoise.poisson,meth.1=hf.la10.ti5,cs.1=50,hybrid=F)
est.hf.ti.r.6.ang.8=apply(sim.m.ang.8,1,denoise.poisson,meth.1=hf.la10.ti6,cs.1=50,hybrid=F)
est.hf.ti.r.7.ang.8=apply(sim.m.ang.8,1,denoise.poisson,meth.1=hf.la10.ti7,cs.1=50,hybrid=F)

est.hf.ti.r.2.hs.8=apply(sim.m.hs.8,1,denoise.poisson,meth.1=hf.la10.ti2,cs.1=50,hybrid=F)
est.hf.ti.r.3.hs.8=apply(sim.m.hs.8,1,denoise.poisson,meth.1=hf.la10.ti3,cs.1=50,hybrid=F)
est.hf.ti.r.4.hs.8=apply(sim.m.hs.8,1,denoise.poisson,meth.1=hf.la10.ti4,cs.1=50,hybrid=F)
est.hf.ti.r.5.hs.8=apply(sim.m.hs.8,1,denoise.poisson,meth.1=hf.la10.ti5,cs.1=50,hybrid=F)
est.hf.ti.r.6.hs.8=apply(sim.m.hs.8,1,denoise.poisson,meth.1=hf.la10.ti6,cs.1=50,hybrid=F)
est.hf.ti.r.7.hs.8=apply(sim.m.hs.8,1,denoise.poisson,meth.1=hf.la10.ti7,cs.1=50,hybrid=F)

est.hf.ti.r.2.bur.8=apply(sim.m.bur.8,1,denoise.poisson,meth.1=hf.la10.ti2,cs.1=50,hybrid=F)
est.hf.ti.r.3.bur.8=apply(sim.m.bur.8,1,denoise.poisson,meth.1=hf.la10.ti3,cs.1=50,hybrid=F)
est.hf.ti.r.4.bur.8=apply(sim.m.bur.8,1,denoise.poisson,meth.1=hf.la10.ti4,cs.1=50,hybrid=F)
est.hf.ti.r.5.bur.8=apply(sim.m.bur.8,1,denoise.poisson,meth.1=hf.la10.ti5,cs.1=50,hybrid=F)
est.hf.ti.r.6.bur.8=apply(sim.m.bur.8,1,denoise.poisson,meth.1=hf.la10.ti6,cs.1=50,hybrid=F)
est.hf.ti.r.7.bur.8=apply(sim.m.bur.8,1,denoise.poisson,meth.1=hf.la10.ti7,cs.1=50,hybrid=F)

est.hf.ti.r.2.cb.8=apply(sim.m.cb.8,1,denoise.poisson,meth.1=hf.la10.ti2,cs.1=50,hybrid=F)
est.hf.ti.r.3.cb.8=apply(sim.m.cb.8,1,denoise.poisson,meth.1=hf.la10.ti3,cs.1=50,hybrid=F)
est.hf.ti.r.4.cb.8=apply(sim.m.cb.8,1,denoise.poisson,meth.1=hf.la10.ti4,cs.1=50,hybrid=F)
est.hf.ti.r.5.cb.8=apply(sim.m.cb.8,1,denoise.poisson,meth.1=hf.la10.ti5,cs.1=50,hybrid=F)
est.hf.ti.r.6.cb.8=apply(sim.m.cb.8,1,denoise.poisson,meth.1=hf.la10.ti6,cs.1=50,hybrid=F)
est.hf.ti.r.7.cb.8=apply(sim.m.cb.8,1,denoise.poisson,meth.1=hf.la10.ti7,cs.1=50,hybrid=F)

est.hf.ti.r.2.b.8=apply(sim.m.b.8,1,denoise.poisson,meth.1=hf.la10.ti2,cs.1=50,hybrid=F)
est.hf.ti.r.3.b.8=apply(sim.m.b.8,1,denoise.poisson,meth.1=hf.la10.ti3,cs.1=50,hybrid=F)
est.hf.ti.r.4.b.8=apply(sim.m.b.8,1,denoise.poisson,meth.1=hf.la10.ti4,cs.1=50,hybrid=F)
est.hf.ti.r.5.b.8=apply(sim.m.b.8,1,denoise.poisson,meth.1=hf.la10.ti5,cs.1=50,hybrid=F)
est.hf.ti.r.6.b.8=apply(sim.m.b.8,1,denoise.poisson,meth.1=hf.la10.ti6,cs.1=50,hybrid=F)
est.hf.ti.r.7.b.8=apply(sim.m.b.8,1,denoise.poisson,meth.1=hf.la10.ti7,cs.1=50,hybrid=F)


est.hf.ti.r.2.s.128=apply(sim.m.s.128,1,denoise.poisson,meth.1=hf.la10.ti2,cs.1=50,hybrid=F)
est.hf.ti.r.3.s.128=apply(sim.m.s.128,1,denoise.poisson,meth.1=hf.la10.ti3,cs.1=50,hybrid=F)
est.hf.ti.r.4.s.128=apply(sim.m.s.128,1,denoise.poisson,meth.1=hf.la10.ti4,cs.1=50,hybrid=F)
est.hf.ti.r.5.s.128=apply(sim.m.s.128,1,denoise.poisson,meth.1=hf.la10.ti5,cs.1=50,hybrid=F)
est.hf.ti.r.6.s.128=apply(sim.m.s.128,1,denoise.poisson,meth.1=hf.la10.ti6,cs.1=50,hybrid=F)
est.hf.ti.r.7.s.128=apply(sim.m.s.128,1,denoise.poisson,meth.1=hf.la10.ti7,cs.1=50,hybrid=F)

est.hf.ti.r.2.ang.128=apply(sim.m.ang.128,1,denoise.poisson,meth.1=hf.la10.ti2,cs.1=50,hybrid=F)
est.hf.ti.r.3.ang.128=apply(sim.m.ang.128,1,denoise.poisson,meth.1=hf.la10.ti3,cs.1=50,hybrid=F)
est.hf.ti.r.4.ang.128=apply(sim.m.ang.128,1,denoise.poisson,meth.1=hf.la10.ti4,cs.1=50,hybrid=F)
est.hf.ti.r.5.ang.128=apply(sim.m.ang.128,1,denoise.poisson,meth.1=hf.la10.ti5,cs.1=50,hybrid=F)
est.hf.ti.r.6.ang.128=apply(sim.m.ang.128,1,denoise.poisson,meth.1=hf.la10.ti6,cs.1=50,hybrid=F)
est.hf.ti.r.7.ang.128=apply(sim.m.ang.128,1,denoise.poisson,meth.1=hf.la10.ti7,cs.1=50,hybrid=F)

est.hf.ti.r.2.hs.128=apply(sim.m.hs.128,1,denoise.poisson,meth.1=hf.la10.ti2,cs.1=50,hybrid=F)
est.hf.ti.r.3.hs.128=apply(sim.m.hs.128,1,denoise.poisson,meth.1=hf.la10.ti3,cs.1=50,hybrid=F)
est.hf.ti.r.4.hs.128=apply(sim.m.hs.128,1,denoise.poisson,meth.1=hf.la10.ti4,cs.1=50,hybrid=F)
est.hf.ti.r.5.hs.128=apply(sim.m.hs.128,1,denoise.poisson,meth.1=hf.la10.ti5,cs.1=50,hybrid=F)
est.hf.ti.r.6.hs.128=apply(sim.m.hs.128,1,denoise.poisson,meth.1=hf.la10.ti6,cs.1=50,hybrid=F)
est.hf.ti.r.7.hs.128=apply(sim.m.hs.128,1,denoise.poisson,meth.1=hf.la10.ti7,cs.1=50,hybrid=F)

est.hf.ti.r.2.bur.128=apply(sim.m.bur.128,1,denoise.poisson,meth.1=hf.la10.ti2,cs.1=50,hybrid=F)
est.hf.ti.r.3.bur.128=apply(sim.m.bur.128,1,denoise.poisson,meth.1=hf.la10.ti3,cs.1=50,hybrid=F)
est.hf.ti.r.4.bur.128=apply(sim.m.bur.128,1,denoise.poisson,meth.1=hf.la10.ti4,cs.1=50,hybrid=F)
est.hf.ti.r.5.bur.128=apply(sim.m.bur.128,1,denoise.poisson,meth.1=hf.la10.ti5,cs.1=50,hybrid=F)
est.hf.ti.r.6.bur.128=apply(sim.m.bur.128,1,denoise.poisson,meth.1=hf.la10.ti6,cs.1=50,hybrid=F)
est.hf.ti.r.7.bur.128=apply(sim.m.bur.128,1,denoise.poisson,meth.1=hf.la10.ti7,cs.1=50,hybrid=F)

est.hf.ti.r.2.cb.128=apply(sim.m.cb.128,1,denoise.poisson,meth.1=hf.la10.ti2,cs.1=50,hybrid=F)
est.hf.ti.r.3.cb.128=apply(sim.m.cb.128,1,denoise.poisson,meth.1=hf.la10.ti3,cs.1=50,hybrid=F)
est.hf.ti.r.4.cb.128=apply(sim.m.cb.128,1,denoise.poisson,meth.1=hf.la10.ti4,cs.1=50,hybrid=F)
est.hf.ti.r.5.cb.128=apply(sim.m.cb.128,1,denoise.poisson,meth.1=hf.la10.ti5,cs.1=50,hybrid=F)
est.hf.ti.r.6.cb.128=apply(sim.m.cb.128,1,denoise.poisson,meth.1=hf.la10.ti6,cs.1=50,hybrid=F)
est.hf.ti.r.7.cb.128=apply(sim.m.cb.128,1,denoise.poisson,meth.1=hf.la10.ti7,cs.1=50,hybrid=F)

est.hf.ti.r.2.b.128=apply(sim.m.b.128,1,denoise.poisson,meth.1=hf.la10.ti2,cs.1=50,hybrid=F)
est.hf.ti.r.3.b.128=apply(sim.m.b.128,1,denoise.poisson,meth.1=hf.la10.ti3,cs.1=50,hybrid=F)
est.hf.ti.r.4.b.128=apply(sim.m.b.128,1,denoise.poisson,meth.1=hf.la10.ti4,cs.1=50,hybrid=F)
est.hf.ti.r.5.b.128=apply(sim.m.b.128,1,denoise.poisson,meth.1=hf.la10.ti5,cs.1=50,hybrid=F)
est.hf.ti.r.6.b.128=apply(sim.m.b.128,1,denoise.poisson,meth.1=hf.la10.ti6,cs.1=50,hybrid=F)
est.hf.ti.r.7.b.128=apply(sim.m.b.128,1,denoise.poisson,meth.1=hf.la10.ti7,cs.1=50,hybrid=F)


est.hf.ti.r.2.s.nz=apply(sim.m.s.nz,1,denoise.poisson,meth.1=hf.la10.ti2,cs.1=50,hybrid=F)
est.hf.ti.r.3.s.nz=apply(sim.m.s.nz,1,denoise.poisson,meth.1=hf.la10.ti3,cs.1=50,hybrid=F)
est.hf.ti.r.4.s.nz=apply(sim.m.s.nz,1,denoise.poisson,meth.1=hf.la10.ti4,cs.1=50,hybrid=F)
est.hf.ti.r.5.s.nz=apply(sim.m.s.nz,1,denoise.poisson,meth.1=hf.la10.ti5,cs.1=50,hybrid=F)
est.hf.ti.r.6.s.nz=apply(sim.m.s.nz,1,denoise.poisson,meth.1=hf.la10.ti6,cs.1=50,hybrid=F)
est.hf.ti.r.7.s.nz=apply(sim.m.s.nz,1,denoise.poisson,meth.1=hf.la10.ti7,cs.1=50,hybrid=F)

est.hf.ti.r.2.ang.nz=apply(sim.m.ang.nz,1,denoise.poisson,meth.1=hf.la10.ti2,cs.1=50,hybrid=F)
est.hf.ti.r.3.ang.nz=apply(sim.m.ang.nz,1,denoise.poisson,meth.1=hf.la10.ti3,cs.1=50,hybrid=F)
est.hf.ti.r.4.ang.nz=apply(sim.m.ang.nz,1,denoise.poisson,meth.1=hf.la10.ti4,cs.1=50,hybrid=F)
est.hf.ti.r.5.ang.nz=apply(sim.m.ang.nz,1,denoise.poisson,meth.1=hf.la10.ti5,cs.1=50,hybrid=F)
est.hf.ti.r.6.ang.nz=apply(sim.m.ang.nz,1,denoise.poisson,meth.1=hf.la10.ti6,cs.1=50,hybrid=F)
est.hf.ti.r.7.ang.nz=apply(sim.m.ang.nz,1,denoise.poisson,meth.1=hf.la10.ti7,cs.1=50,hybrid=F)

est.hf.ti.r.2.hs.nz=apply(sim.m.hs.nz,1,denoise.poisson,meth.1=hf.la10.ti2,cs.1=50,hybrid=F)
est.hf.ti.r.3.hs.nz=apply(sim.m.hs.nz,1,denoise.poisson,meth.1=hf.la10.ti3,cs.1=50,hybrid=F)
est.hf.ti.r.4.hs.nz=apply(sim.m.hs.nz,1,denoise.poisson,meth.1=hf.la10.ti4,cs.1=50,hybrid=F)
est.hf.ti.r.5.hs.nz=apply(sim.m.hs.nz,1,denoise.poisson,meth.1=hf.la10.ti5,cs.1=50,hybrid=F)
est.hf.ti.r.6.hs.nz=apply(sim.m.hs.nz,1,denoise.poisson,meth.1=hf.la10.ti6,cs.1=50,hybrid=F)
est.hf.ti.r.7.hs.nz=apply(sim.m.hs.nz,1,denoise.poisson,meth.1=hf.la10.ti7,cs.1=50,hybrid=F)

est.hf.ti.r.2.bur.nz=apply(sim.m.bur.nz,1,denoise.poisson,meth.1=hf.la10.ti2,cs.1=50,hybrid=F)
est.hf.ti.r.3.bur.nz=apply(sim.m.bur.nz,1,denoise.poisson,meth.1=hf.la10.ti3,cs.1=50,hybrid=F)
est.hf.ti.r.4.bur.nz=apply(sim.m.bur.nz,1,denoise.poisson,meth.1=hf.la10.ti4,cs.1=50,hybrid=F)
est.hf.ti.r.5.bur.nz=apply(sim.m.bur.nz,1,denoise.poisson,meth.1=hf.la10.ti5,cs.1=50,hybrid=F)
est.hf.ti.r.6.bur.nz=apply(sim.m.bur.nz,1,denoise.poisson,meth.1=hf.la10.ti6,cs.1=50,hybrid=F)
est.hf.ti.r.7.bur.nz=apply(sim.m.bur.nz,1,denoise.poisson,meth.1=hf.la10.ti7,cs.1=50,hybrid=F)

est.hf.ti.r.2.cb.nz=apply(sim.m.cb.nz,1,denoise.poisson,meth.1=hf.la10.ti2,cs.1=50,hybrid=F)
est.hf.ti.r.3.cb.nz=apply(sim.m.cb.nz,1,denoise.poisson,meth.1=hf.la10.ti3,cs.1=50,hybrid=F)
est.hf.ti.r.4.cb.nz=apply(sim.m.cb.nz,1,denoise.poisson,meth.1=hf.la10.ti4,cs.1=50,hybrid=F)
est.hf.ti.r.5.cb.nz=apply(sim.m.cb.nz,1,denoise.poisson,meth.1=hf.la10.ti5,cs.1=50,hybrid=F)
est.hf.ti.r.6.cb.nz=apply(sim.m.cb.nz,1,denoise.poisson,meth.1=hf.la10.ti6,cs.1=50,hybrid=F)
est.hf.ti.r.7.cb.nz=apply(sim.m.cb.nz,1,denoise.poisson,meth.1=hf.la10.ti7,cs.1=50,hybrid=F)

est.hf.ti.r.2.b.nz=apply(sim.m.b.nz,1,denoise.poisson,meth.1=hf.la10.ti2,cs.1=50,hybrid=F)
est.hf.ti.r.3.b.nz=apply(sim.m.b.nz,1,denoise.poisson,meth.1=hf.la10.ti3,cs.1=50,hybrid=F)
est.hf.ti.r.4.b.nz=apply(sim.m.b.nz,1,denoise.poisson,meth.1=hf.la10.ti4,cs.1=50,hybrid=F)
est.hf.ti.r.5.b.nz=apply(sim.m.b.nz,1,denoise.poisson,meth.1=hf.la10.ti5,cs.1=50,hybrid=F)
est.hf.ti.r.6.b.nz=apply(sim.m.b.nz,1,denoise.poisson,meth.1=hf.la10.ti6,cs.1=50,hybrid=F)
est.hf.ti.r.7.b.nz=apply(sim.m.b.nz,1,denoise.poisson,meth.1=hf.la10.ti7,cs.1=50,hybrid=F)


save.image("res_sim_pois.RData")

