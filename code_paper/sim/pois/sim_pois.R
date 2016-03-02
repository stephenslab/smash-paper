library(smash)
library(haarfisz)

#define functions
mse=function(x,y) mean((x-y)^2)
l2norm=function(x) sum(x^2)
mise=function(x,y) 10000*apply(x-matrix(rep(y,times=100),nrow=100,byrow=T),1,l2norm)/l2norm(y)


#generate test functions
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



#generate different intensities
mu.s.1=1*(0.01+mu.s)
mu.s.8=8/3*(3/64+mu.s)
mu.s.128=128/3*(3/128^2+mu.s)

mu.ang.1=1*(0.01+mu.ang)
mu.ang.8=8/3*(3/64+mu.ang)
mu.ang.128=1/128+128/3*mu.ang

mu.hs.1=1*(0.01+mu.hs)
mu.hs.8=8/3*(3/64+mu.hs)
mu.hs.128=1/128+128/3*mu.hs

mu.bur.1=1*(0.01+mu.bur)
mu.bur.8=8/3*(3/64+mu.bur)
mu.bur.128=1/128+128/3*mu.bur

mu.cb.1=1*(0.01+mu.cb)
mu.cb.8=8/3*(3/64+mu.cb)
mu.cb.128=128/3*(3/128^2+mu.cb)

mu.b.1=1*(0.01+mu.b)
mu.b.8=8/3*(3/64+mu.b)
mu.b.128=128/3*(3/128^2+mu.b)


##generate data
set.seed(1002)
sim.m.s.1=matrix(rpois(100*n,mu.s.1),nrow=100,ncol=n,byrow=T)
set.seed(1002)
sim.m.s.8=matrix(rpois(100*n,mu.s.8),nrow=100,ncol=n,byrow=T)
set.seed(1002)
sim.m.s.128=matrix(rpois(100*n,mu.s.128),nrow=100,ncol=n,byrow=T)

set.seed(1002)
sim.m.ang.1=matrix(rpois(100*n,mu.ang.1),nrow=100,ncol=n,byrow=T)
set.seed(1002)
sim.m.ang.8=matrix(rpois(100*n,mu.ang.8),nrow=100,ncol=n,byrow=T)
set.seed(1002)
sim.m.ang.128=matrix(rpois(100*n,mu.ang.128),nrow=100,ncol=n,byrow=T)

set.seed(1002)
sim.m.hs.1=matrix(rpois(100*n,mu.hs.1),nrow=100,ncol=n,byrow=T)
set.seed(1002)
sim.m.hs.8=matrix(rpois(100*n,mu.hs.8),nrow=100,ncol=n,byrow=T)
set.seed(1002)
sim.m.hs.128=matrix(rpois(100*n,mu.hs.128),nrow=100,ncol=n,byrow=T)

set.seed(1002)
sim.m.bur.1=matrix(rpois(100*n,mu.bur.1),nrow=100,ncol=n,byrow=T)
set.seed(1002)
sim.m.bur.8=matrix(rpois(100*n,mu.bur.8),nrow=100,ncol=n,byrow=T)
set.seed(1002)
sim.m.bur.128=matrix(rpois(100*n,mu.bur.128),nrow=100,ncol=n,byrow=T)

set.seed(1002)
sim.m.cb.1=matrix(rpois(100*n,mu.cb.1),nrow=100,ncol=n,byrow=T)
set.seed(1002)
sim.m.cb.8=matrix(rpois(100*n,mu.cb.8),nrow=100,ncol=n,byrow=T)
set.seed(1002)
sim.m.cb.128=matrix(rpois(100*n,mu.cb.128),nrow=100,ncol=n,byrow=T)

set.seed(1002)
sim.m.b.1=matrix(rpois(100*n,mu.b.1),nrow=100,ncol=n,byrow=T)
set.seed(1002)
sim.m.b.8=matrix(rpois(100*n,mu.b.8),nrow=100,ncol=n,byrow=T)
set.seed(1002)
sim.m.b.128=matrix(rpois(100*n,mu.b.128),nrow=100,ncol=n,byrow=T)



#write data to be used in matlab
write(t(sim.m.s.1),'sim_s_1.txt',ncolumns=n)
write(t(sim.m.s.8),'sim_s_8.txt',ncolumns=n)
write(t(sim.m.s.128),'sim_s_128.txt',ncolumns=n)
write(t(sim.m.ang.1),'sim_ang_1.txt',ncolumns=n)
write(t(sim.m.ang.8),'sim_ang_8.txt',ncolumns=n)
write(t(sim.m.ang.128),'sim_ang_128.txt',ncolumns=n)
write(t(sim.m.hs.1),'sim_hs_1.txt',ncolumns=n)
write(t(sim.m.hs.8),'sim_hs_8.txt',ncolumns=n)
write(t(sim.m.hs.128),'sim_hs_128.txt',ncolumns=n)
write(t(sim.m.bur.1),'sim_bur_1.txt',ncolumns=n)
write(t(sim.m.bur.8),'sim_bur_8.txt',ncolumns=n)
write(t(sim.m.bur.128),'sim_bur_128.txt',ncolumns=n)
write(t(sim.m.cb.1),'sim_cb_1.txt',ncolumns=n)
write(t(sim.m.cb.8),'sim_cb_8.txt',ncolumns=n)
write(t(sim.m.cb.128),'sim_cb_128.txt',ncolumns=n)
write(t(sim.m.b.1),'sim_b_1.txt',ncolumns=n)
write(t(sim.m.b.8),'sim_b_8.txt',ncolumns=n)
write(t(sim.m.b.128),'sim_b_128.txt',ncolumns=n)




#run smash on data
est.ash.s.1=apply(sim.m.s.1,1,ashsmooth.pois)

est.ash.ang.1=apply(sim.m.ang.1,1,ashsmooth.pois)

est.ash.hs.1=apply(sim.m.hs.1,1,ashsmooth.pois)

est.ash.bur.1=apply(sim.m.bur.1,1,ashsmooth.pois)

est.ash.cb.1=apply(sim.m.cb.1,1,ashsmooth.pois)

est.ash.b.1=apply(sim.m.b.1,1,ashsmooth.pois)



est.ash.s.8=apply(sim.m.s.8,1,ashsmooth.pois)

est.ash.ang.8=apply(sim.m.ang.8,1,ashsmooth.pois)

est.ash.hs.8=apply(sim.m.hs.8,1,ashsmooth.pois)

est.ash.bur.8=apply(sim.m.bur.8,1,ashsmooth.pois)

est.ash.cb.8=apply(sim.m.cb.8,1,ashsmooth.pois)

est.ash.b.8=apply(sim.m.b.8,1,ashsmooth.pois)





est.ash.s.128=apply(sim.m.s.128,1,ashsmooth.pois)

est.ash.ang.128=apply(sim.m.ang.128,1,ashsmooth.pois)

est.ash.hs.128=apply(sim.m.hs.128,1,ashsmooth.pois)

est.ash.bur.128=apply(sim.m.bur.128,1,ashsmooth.pois)

est.ash.cb.128=apply(sim.m.cb.128,1,ashsmooth.pois)

est.ash.b.128=apply(sim.m.b.128,1,ashsmooth.pois)





#define gaussian denoising steps for haarfisz
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

#run haarfisz
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




#load in results from matlab runs (anscombe - cv, universal, haar thresholding, BMSM, BMIE, l1 penalty )
est.an.ti.cvh.2.s.1=data.matrix(read.csv('est_an_ti_cvh_2_s_1.csv',header=F))
est.an.ti.cvh.3.s.1=data.matrix(read.csv('est_an_ti_cvh_3_s_1.csv',header=F))
est.an.ti.cvh.4.s.1=data.matrix(read.csv('est_an_ti_cvh_4_s_1.csv',header=F))
est.an.ti.cvh.5.s.1=data.matrix(read.csv('est_an_ti_cvh_5_s_1.csv',header=F))
est.an.ti.cvh.6.s.1=data.matrix(read.csv('est_an_ti_cvh_6_s_1.csv',header=F))
est.an.ti.cvh.7.s.1=data.matrix(read.csv('est_an_ti_cvh_7_s_1.csv',header=F))
est.an.ti.cvh.2.s.8=data.matrix(read.csv('est_an_ti_cvh_2_s_8.csv',header=F))
est.an.ti.cvh.3.s.8=data.matrix(read.csv('est_an_ti_cvh_3_s_8.csv',header=F))
est.an.ti.cvh.4.s.8=data.matrix(read.csv('est_an_ti_cvh_4_s_8.csv',header=F))
est.an.ti.cvh.5.s.8=data.matrix(read.csv('est_an_ti_cvh_5_s_8.csv',header=F))
est.an.ti.cvh.6.s.8=data.matrix(read.csv('est_an_ti_cvh_6_s_8.csv',header=F))
est.an.ti.cvh.7.s.8=data.matrix(read.csv('est_an_ti_cvh_7_s_8.csv',header=F))
est.an.ti.cvh.2.s.128=data.matrix(read.csv('est_an_ti_cvh_2_s_128.csv',header=F))
est.an.ti.cvh.3.s.128=data.matrix(read.csv('est_an_ti_cvh_3_s_128.csv',header=F))
est.an.ti.cvh.4.s.128=data.matrix(read.csv('est_an_ti_cvh_4_s_128.csv',header=F))
est.an.ti.cvh.5.s.128=data.matrix(read.csv('est_an_ti_cvh_5_s_128.csv',header=F))
est.an.ti.cvh.6.s.128=data.matrix(read.csv('est_an_ti_cvh_6_s_128.csv',header=F))
est.an.ti.cvh.7.s.128=data.matrix(read.csv('est_an_ti_cvh_7_s_128.csv',header=F))



est.an.ti.cvh.2.ang.1=data.matrix(read.csv('est_an_ti_cvh_2_ang_1.csv',header=F))
est.an.ti.cvh.3.ang.1=data.matrix(read.csv('est_an_ti_cvh_3_ang_1.csv',header=F))
est.an.ti.cvh.4.ang.1=data.matrix(read.csv('est_an_ti_cvh_4_ang_1.csv',header=F))
est.an.ti.cvh.5.ang.1=data.matrix(read.csv('est_an_ti_cvh_5_ang_1.csv',header=F))
est.an.ti.cvh.6.ang.1=data.matrix(read.csv('est_an_ti_cvh_6_ang_1.csv',header=F))
est.an.ti.cvh.7.ang.1=data.matrix(read.csv('est_an_ti_cvh_7_ang_1.csv',header=F))
est.an.ti.cvh.2.ang.8=data.matrix(read.csv('est_an_ti_cvh_2_ang_8.csv',header=F))
est.an.ti.cvh.3.ang.8=data.matrix(read.csv('est_an_ti_cvh_3_ang_8.csv',header=F))
est.an.ti.cvh.4.ang.8=data.matrix(read.csv('est_an_ti_cvh_4_ang_8.csv',header=F))
est.an.ti.cvh.5.ang.8=data.matrix(read.csv('est_an_ti_cvh_5_ang_8.csv',header=F))
est.an.ti.cvh.6.ang.8=data.matrix(read.csv('est_an_ti_cvh_6_ang_8.csv',header=F))
est.an.ti.cvh.7.ang.8=data.matrix(read.csv('est_an_ti_cvh_7_ang_8.csv',header=F))
est.an.ti.cvh.2.ang.128=data.matrix(read.csv('est_an_ti_cvh_2_ang_128.csv',header=F))
est.an.ti.cvh.3.ang.128=data.matrix(read.csv('est_an_ti_cvh_3_ang_128.csv',header=F))
est.an.ti.cvh.4.ang.128=data.matrix(read.csv('est_an_ti_cvh_4_ang_128.csv',header=F))
est.an.ti.cvh.5.ang.128=data.matrix(read.csv('est_an_ti_cvh_5_ang_128.csv',header=F))
est.an.ti.cvh.6.ang.128=data.matrix(read.csv('est_an_ti_cvh_6_ang_128.csv',header=F))
est.an.ti.cvh.7.ang.128=data.matrix(read.csv('est_an_ti_cvh_7_ang_128.csv',header=F))



est.an.ti.cvh.2.hs.1=data.matrix(read.csv('est_an_ti_cvh_2_hs_1.csv',header=F))
est.an.ti.cvh.3.hs.1=data.matrix(read.csv('est_an_ti_cvh_3_hs_1.csv',header=F))
est.an.ti.cvh.4.hs.1=data.matrix(read.csv('est_an_ti_cvh_4_hs_1.csv',header=F))
est.an.ti.cvh.5.hs.1=data.matrix(read.csv('est_an_ti_cvh_5_hs_1.csv',header=F))
est.an.ti.cvh.6.hs.1=data.matrix(read.csv('est_an_ti_cvh_6_hs_1.csv',header=F))
est.an.ti.cvh.7.hs.1=data.matrix(read.csv('est_an_ti_cvh_7_hs_1.csv',header=F))
est.an.ti.cvh.2.hs.8=data.matrix(read.csv('est_an_ti_cvh_2_hs_8.csv',header=F))
est.an.ti.cvh.3.hs.8=data.matrix(read.csv('est_an_ti_cvh_3_hs_8.csv',header=F))
est.an.ti.cvh.4.hs.8=data.matrix(read.csv('est_an_ti_cvh_4_hs_8.csv',header=F))
est.an.ti.cvh.5.hs.8=data.matrix(read.csv('est_an_ti_cvh_5_hs_8.csv',header=F))
est.an.ti.cvh.6.hs.8=data.matrix(read.csv('est_an_ti_cvh_6_hs_8.csv',header=F))
est.an.ti.cvh.7.hs.8=data.matrix(read.csv('est_an_ti_cvh_7_hs_8.csv',header=F))
est.an.ti.cvh.2.hs.128=data.matrix(read.csv('est_an_ti_cvh_2_hs_128.csv',header=F))
est.an.ti.cvh.3.hs.128=data.matrix(read.csv('est_an_ti_cvh_3_hs_128.csv',header=F))
est.an.ti.cvh.4.hs.128=data.matrix(read.csv('est_an_ti_cvh_4_hs_128.csv',header=F))
est.an.ti.cvh.5.hs.128=data.matrix(read.csv('est_an_ti_cvh_5_hs_128.csv',header=F))
est.an.ti.cvh.6.hs.128=data.matrix(read.csv('est_an_ti_cvh_6_hs_128.csv',header=F))
est.an.ti.cvh.7.hs.128=data.matrix(read.csv('est_an_ti_cvh_7_hs_128.csv',header=F))



est.an.ti.cvh.2.bur.1=data.matrix(read.csv('est_an_ti_cvh_2_bur_1.csv',header=F))
est.an.ti.cvh.3.bur.1=data.matrix(read.csv('est_an_ti_cvh_3_bur_1.csv',header=F))
est.an.ti.cvh.4.bur.1=data.matrix(read.csv('est_an_ti_cvh_4_bur_1.csv',header=F))
est.an.ti.cvh.5.bur.1=data.matrix(read.csv('est_an_ti_cvh_5_bur_1.csv',header=F))
est.an.ti.cvh.6.bur.1=data.matrix(read.csv('est_an_ti_cvh_6_bur_1.csv',header=F))
est.an.ti.cvh.7.bur.1=data.matrix(read.csv('est_an_ti_cvh_7_bur_1.csv',header=F))
est.an.ti.cvh.2.bur.8=data.matrix(read.csv('est_an_ti_cvh_2_bur_8.csv',header=F))
est.an.ti.cvh.3.bur.8=data.matrix(read.csv('est_an_ti_cvh_3_bur_8.csv',header=F))
est.an.ti.cvh.4.bur.8=data.matrix(read.csv('est_an_ti_cvh_4_bur_8.csv',header=F))
est.an.ti.cvh.5.bur.8=data.matrix(read.csv('est_an_ti_cvh_5_bur_8.csv',header=F))
est.an.ti.cvh.6.bur.8=data.matrix(read.csv('est_an_ti_cvh_6_bur_8.csv',header=F))
est.an.ti.cvh.7.bur.8=data.matrix(read.csv('est_an_ti_cvh_7_bur_8.csv',header=F))
est.an.ti.cvh.2.bur.128=data.matrix(read.csv('est_an_ti_cvh_2_bur_128.csv',header=F))
est.an.ti.cvh.3.bur.128=data.matrix(read.csv('est_an_ti_cvh_3_bur_128.csv',header=F))
est.an.ti.cvh.4.bur.128=data.matrix(read.csv('est_an_ti_cvh_4_bur_128.csv',header=F))
est.an.ti.cvh.5.bur.128=data.matrix(read.csv('est_an_ti_cvh_5_bur_128.csv',header=F))
est.an.ti.cvh.6.bur.128=data.matrix(read.csv('est_an_ti_cvh_6_bur_128.csv',header=F))
est.an.ti.cvh.7.bur.128=data.matrix(read.csv('est_an_ti_cvh_7_bur_128.csv',header=F))


est.an.ti.cvh.2.cb.1=data.matrix(read.csv('est_an_ti_cvh_2_cb_1.csv',header=F))
est.an.ti.cvh.3.cb.1=data.matrix(read.csv('est_an_ti_cvh_3_cb_1.csv',header=F))
est.an.ti.cvh.4.cb.1=data.matrix(read.csv('est_an_ti_cvh_4_cb_1.csv',header=F))
est.an.ti.cvh.5.cb.1=data.matrix(read.csv('est_an_ti_cvh_5_cb_1.csv',header=F))
est.an.ti.cvh.6.cb.1=data.matrix(read.csv('est_an_ti_cvh_6_cb_1.csv',header=F))
est.an.ti.cvh.7.cb.1=data.matrix(read.csv('est_an_ti_cvh_7_cb_1.csv',header=F))
est.an.ti.cvh.2.cb.8=data.matrix(read.csv('est_an_ti_cvh_2_cb_8.csv',header=F))
est.an.ti.cvh.3.cb.8=data.matrix(read.csv('est_an_ti_cvh_3_cb_8.csv',header=F))
est.an.ti.cvh.4.cb.8=data.matrix(read.csv('est_an_ti_cvh_4_cb_8.csv',header=F))
est.an.ti.cvh.5.cb.8=data.matrix(read.csv('est_an_ti_cvh_5_cb_8.csv',header=F))
est.an.ti.cvh.6.cb.8=data.matrix(read.csv('est_an_ti_cvh_6_cb_8.csv',header=F))
est.an.ti.cvh.7.cb.8=data.matrix(read.csv('est_an_ti_cvh_7_cb_8.csv',header=F))
est.an.ti.cvh.2.cb.128=data.matrix(read.csv('est_an_ti_cvh_2_cb_128.csv',header=F))
est.an.ti.cvh.3.cb.128=data.matrix(read.csv('est_an_ti_cvh_3_cb_128.csv',header=F))
est.an.ti.cvh.4.cb.128=data.matrix(read.csv('est_an_ti_cvh_4_cb_128.csv',header=F))
est.an.ti.cvh.5.cb.128=data.matrix(read.csv('est_an_ti_cvh_5_cb_128.csv',header=F))
est.an.ti.cvh.6.cb.128=data.matrix(read.csv('est_an_ti_cvh_6_cb_128.csv',header=F))
est.an.ti.cvh.7.cb.128=data.matrix(read.csv('est_an_ti_cvh_7_cb_128.csv',header=F))


est.an.ti.cvh.2.b.1=data.matrix(read.csv('est_an_ti_cvh_2_b_1.csv',header=F))
est.an.ti.cvh.3.b.1=data.matrix(read.csv('est_an_ti_cvh_3_b_1.csv',header=F))
est.an.ti.cvh.4.b.1=data.matrix(read.csv('est_an_ti_cvh_4_b_1.csv',header=F))
est.an.ti.cvh.5.b.1=data.matrix(read.csv('est_an_ti_cvh_5_b_1.csv',header=F))
est.an.ti.cvh.6.b.1=data.matrix(read.csv('est_an_ti_cvh_6_b_1.csv',header=F))
est.an.ti.cvh.7.b.1=data.matrix(read.csv('est_an_ti_cvh_7_b_1.csv',header=F))
est.an.ti.cvh.2.b.8=data.matrix(read.csv('est_an_ti_cvh_2_b_8.csv',header=F))
est.an.ti.cvh.3.b.8=data.matrix(read.csv('est_an_ti_cvh_3_b_8.csv',header=F))
est.an.ti.cvh.4.b.8=data.matrix(read.csv('est_an_ti_cvh_4_b_8.csv',header=F))
est.an.ti.cvh.5.b.8=data.matrix(read.csv('est_an_ti_cvh_5_b_8.csv',header=F))
est.an.ti.cvh.6.b.8=data.matrix(read.csv('est_an_ti_cvh_6_b_8.csv',header=F))
est.an.ti.cvh.7.b.8=data.matrix(read.csv('est_an_ti_cvh_7_b_8.csv',header=F))
est.an.ti.cvh.2.b.128=data.matrix(read.csv('est_an_ti_cvh_2_b_128.csv',header=F))
est.an.ti.cvh.3.b.128=data.matrix(read.csv('est_an_ti_cvh_3_b_128.csv',header=F))
est.an.ti.cvh.4.b.128=data.matrix(read.csv('est_an_ti_cvh_4_b_128.csv',header=F))
est.an.ti.cvh.5.b.128=data.matrix(read.csv('est_an_ti_cvh_5_b_128.csv',header=F))
est.an.ti.cvh.6.b.128=data.matrix(read.csv('est_an_ti_cvh_6_b_128.csv',header=F))
est.an.ti.cvh.7.b.128=data.matrix(read.csv('est_an_ti_cvh_7_b_128.csv',header=F))



est.an.ti.uh.2.s.1=data.matrix(read.csv('est_an_ti_uh_2_s_1.csv',header=F))
est.an.ti.uh.3.s.1=data.matrix(read.csv('est_an_ti_uh_3_s_1.csv',header=F))
est.an.ti.uh.4.s.1=data.matrix(read.csv('est_an_ti_uh_4_s_1.csv',header=F))
est.an.ti.uh.5.s.1=data.matrix(read.csv('est_an_ti_uh_5_s_1.csv',header=F))
est.an.ti.uh.6.s.1=data.matrix(read.csv('est_an_ti_uh_6_s_1.csv',header=F))
est.an.ti.uh.7.s.1=data.matrix(read.csv('est_an_ti_uh_7_s_1.csv',header=F))
est.an.ti.uh.2.s.8=data.matrix(read.csv('est_an_ti_uh_2_s_8.csv',header=F))
est.an.ti.uh.3.s.8=data.matrix(read.csv('est_an_ti_uh_3_s_8.csv',header=F))
est.an.ti.uh.4.s.8=data.matrix(read.csv('est_an_ti_uh_4_s_8.csv',header=F))
est.an.ti.uh.5.s.8=data.matrix(read.csv('est_an_ti_uh_5_s_8.csv',header=F))
est.an.ti.uh.6.s.8=data.matrix(read.csv('est_an_ti_uh_6_s_8.csv',header=F))
est.an.ti.uh.7.s.8=data.matrix(read.csv('est_an_ti_uh_7_s_8.csv',header=F))
est.an.ti.uh.2.s.128=data.matrix(read.csv('est_an_ti_uh_2_s_128.csv',header=F))
est.an.ti.uh.3.s.128=data.matrix(read.csv('est_an_ti_uh_3_s_128.csv',header=F))
est.an.ti.uh.4.s.128=data.matrix(read.csv('est_an_ti_uh_4_s_128.csv',header=F))
est.an.ti.uh.5.s.128=data.matrix(read.csv('est_an_ti_uh_5_s_128.csv',header=F))
est.an.ti.uh.6.s.128=data.matrix(read.csv('est_an_ti_uh_6_s_128.csv',header=F))
est.an.ti.uh.7.s.128=data.matrix(read.csv('est_an_ti_uh_7_s_128.csv',header=F))



est.an.ti.uh.2.ang.1=data.matrix(read.csv('est_an_ti_uh_2_ang_1.csv',header=F))
est.an.ti.uh.3.ang.1=data.matrix(read.csv('est_an_ti_uh_3_ang_1.csv',header=F))
est.an.ti.uh.4.ang.1=data.matrix(read.csv('est_an_ti_uh_4_ang_1.csv',header=F))
est.an.ti.uh.5.ang.1=data.matrix(read.csv('est_an_ti_uh_5_ang_1.csv',header=F))
est.an.ti.uh.6.ang.1=data.matrix(read.csv('est_an_ti_uh_6_ang_1.csv',header=F))
est.an.ti.uh.7.ang.1=data.matrix(read.csv('est_an_ti_uh_7_ang_1.csv',header=F))
est.an.ti.uh.2.ang.8=data.matrix(read.csv('est_an_ti_uh_2_ang_8.csv',header=F))
est.an.ti.uh.3.ang.8=data.matrix(read.csv('est_an_ti_uh_3_ang_8.csv',header=F))
est.an.ti.uh.4.ang.8=data.matrix(read.csv('est_an_ti_uh_4_ang_8.csv',header=F))
est.an.ti.uh.5.ang.8=data.matrix(read.csv('est_an_ti_uh_5_ang_8.csv',header=F))
est.an.ti.uh.6.ang.8=data.matrix(read.csv('est_an_ti_uh_6_ang_8.csv',header=F))
est.an.ti.uh.7.ang.8=data.matrix(read.csv('est_an_ti_uh_7_ang_8.csv',header=F))
est.an.ti.uh.2.ang.128=data.matrix(read.csv('est_an_ti_uh_2_ang_128.csv',header=F))
est.an.ti.uh.3.ang.128=data.matrix(read.csv('est_an_ti_uh_3_ang_128.csv',header=F))
est.an.ti.uh.4.ang.128=data.matrix(read.csv('est_an_ti_uh_4_ang_128.csv',header=F))
est.an.ti.uh.5.ang.128=data.matrix(read.csv('est_an_ti_uh_5_ang_128.csv',header=F))
est.an.ti.uh.6.ang.128=data.matrix(read.csv('est_an_ti_uh_6_ang_128.csv',header=F))
est.an.ti.uh.7.ang.128=data.matrix(read.csv('est_an_ti_uh_7_ang_128.csv',header=F))



est.an.ti.uh.2.hs.1=data.matrix(read.csv('est_an_ti_uh_2_hs_1.csv',header=F))
est.an.ti.uh.3.hs.1=data.matrix(read.csv('est_an_ti_uh_3_hs_1.csv',header=F))
est.an.ti.uh.4.hs.1=data.matrix(read.csv('est_an_ti_uh_4_hs_1.csv',header=F))
est.an.ti.uh.5.hs.1=data.matrix(read.csv('est_an_ti_uh_5_hs_1.csv',header=F))
est.an.ti.uh.6.hs.1=data.matrix(read.csv('est_an_ti_uh_6_hs_1.csv',header=F))
est.an.ti.uh.7.hs.1=data.matrix(read.csv('est_an_ti_uh_7_hs_1.csv',header=F))
est.an.ti.uh.2.hs.8=data.matrix(read.csv('est_an_ti_uh_2_hs_8.csv',header=F))
est.an.ti.uh.3.hs.8=data.matrix(read.csv('est_an_ti_uh_3_hs_8.csv',header=F))
est.an.ti.uh.4.hs.8=data.matrix(read.csv('est_an_ti_uh_4_hs_8.csv',header=F))
est.an.ti.uh.5.hs.8=data.matrix(read.csv('est_an_ti_uh_5_hs_8.csv',header=F))
est.an.ti.uh.6.hs.8=data.matrix(read.csv('est_an_ti_uh_6_hs_8.csv',header=F))
est.an.ti.uh.7.hs.8=data.matrix(read.csv('est_an_ti_uh_7_hs_8.csv',header=F))
est.an.ti.uh.2.hs.128=data.matrix(read.csv('est_an_ti_uh_2_hs_128.csv',header=F))
est.an.ti.uh.3.hs.128=data.matrix(read.csv('est_an_ti_uh_3_hs_128.csv',header=F))
est.an.ti.uh.4.hs.128=data.matrix(read.csv('est_an_ti_uh_4_hs_128.csv',header=F))
est.an.ti.uh.5.hs.128=data.matrix(read.csv('est_an_ti_uh_5_hs_128.csv',header=F))
est.an.ti.uh.6.hs.128=data.matrix(read.csv('est_an_ti_uh_6_hs_128.csv',header=F))
est.an.ti.uh.7.hs.128=data.matrix(read.csv('est_an_ti_uh_7_hs_128.csv',header=F))



est.an.ti.uh.2.bur.1=data.matrix(read.csv('est_an_ti_uh_2_bur_1.csv',header=F))
est.an.ti.uh.3.bur.1=data.matrix(read.csv('est_an_ti_uh_3_bur_1.csv',header=F))
est.an.ti.uh.4.bur.1=data.matrix(read.csv('est_an_ti_uh_4_bur_1.csv',header=F))
est.an.ti.uh.5.bur.1=data.matrix(read.csv('est_an_ti_uh_5_bur_1.csv',header=F))
est.an.ti.uh.6.bur.1=data.matrix(read.csv('est_an_ti_uh_6_bur_1.csv',header=F))
est.an.ti.uh.7.bur.1=data.matrix(read.csv('est_an_ti_uh_7_bur_1.csv',header=F))
est.an.ti.uh.2.bur.8=data.matrix(read.csv('est_an_ti_uh_2_bur_8.csv',header=F))
est.an.ti.uh.3.bur.8=data.matrix(read.csv('est_an_ti_uh_3_bur_8.csv',header=F))
est.an.ti.uh.4.bur.8=data.matrix(read.csv('est_an_ti_uh_4_bur_8.csv',header=F))
est.an.ti.uh.5.bur.8=data.matrix(read.csv('est_an_ti_uh_5_bur_8.csv',header=F))
est.an.ti.uh.6.bur.8=data.matrix(read.csv('est_an_ti_uh_6_bur_8.csv',header=F))
est.an.ti.uh.7.bur.8=data.matrix(read.csv('est_an_ti_uh_7_bur_8.csv',header=F))
est.an.ti.uh.2.bur.128=data.matrix(read.csv('est_an_ti_uh_2_bur_128.csv',header=F))
est.an.ti.uh.3.bur.128=data.matrix(read.csv('est_an_ti_uh_3_bur_128.csv',header=F))
est.an.ti.uh.4.bur.128=data.matrix(read.csv('est_an_ti_uh_4_bur_128.csv',header=F))
est.an.ti.uh.5.bur.128=data.matrix(read.csv('est_an_ti_uh_5_bur_128.csv',header=F))
est.an.ti.uh.6.bur.128=data.matrix(read.csv('est_an_ti_uh_6_bur_128.csv',header=F))
est.an.ti.uh.7.bur.128=data.matrix(read.csv('est_an_ti_uh_7_bur_128.csv',header=F))


est.an.ti.uh.2.cb.1=data.matrix(read.csv('est_an_ti_uh_2_cb_1.csv',header=F))
est.an.ti.uh.3.cb.1=data.matrix(read.csv('est_an_ti_uh_3_cb_1.csv',header=F))
est.an.ti.uh.4.cb.1=data.matrix(read.csv('est_an_ti_uh_4_cb_1.csv',header=F))
est.an.ti.uh.5.cb.1=data.matrix(read.csv('est_an_ti_uh_5_cb_1.csv',header=F))
est.an.ti.uh.6.cb.1=data.matrix(read.csv('est_an_ti_uh_6_cb_1.csv',header=F))
est.an.ti.uh.7.cb.1=data.matrix(read.csv('est_an_ti_uh_7_cb_1.csv',header=F))
est.an.ti.uh.2.cb.8=data.matrix(read.csv('est_an_ti_uh_2_cb_8.csv',header=F))
est.an.ti.uh.3.cb.8=data.matrix(read.csv('est_an_ti_uh_3_cb_8.csv',header=F))
est.an.ti.uh.4.cb.8=data.matrix(read.csv('est_an_ti_uh_4_cb_8.csv',header=F))
est.an.ti.uh.5.cb.8=data.matrix(read.csv('est_an_ti_uh_5_cb_8.csv',header=F))
est.an.ti.uh.6.cb.8=data.matrix(read.csv('est_an_ti_uh_6_cb_8.csv',header=F))
est.an.ti.uh.7.cb.8=data.matrix(read.csv('est_an_ti_uh_7_cb_8.csv',header=F))
est.an.ti.uh.2.cb.128=data.matrix(read.csv('est_an_ti_uh_2_cb_128.csv',header=F))
est.an.ti.uh.3.cb.128=data.matrix(read.csv('est_an_ti_uh_3_cb_128.csv',header=F))
est.an.ti.uh.4.cb.128=data.matrix(read.csv('est_an_ti_uh_4_cb_128.csv',header=F))
est.an.ti.uh.5.cb.128=data.matrix(read.csv('est_an_ti_uh_5_cb_128.csv',header=F))
est.an.ti.uh.6.cb.128=data.matrix(read.csv('est_an_ti_uh_6_cb_128.csv',header=F))
est.an.ti.uh.7.cb.128=data.matrix(read.csv('est_an_ti_uh_7_cb_128.csv',header=F))


est.an.ti.uh.2.b.1=data.matrix(read.csv('est_an_ti_uh_2_b_1.csv',header=F))
est.an.ti.uh.3.b.1=data.matrix(read.csv('est_an_ti_uh_3_b_1.csv',header=F))
est.an.ti.uh.4.b.1=data.matrix(read.csv('est_an_ti_uh_4_b_1.csv',header=F))
est.an.ti.uh.5.b.1=data.matrix(read.csv('est_an_ti_uh_5_b_1.csv',header=F))
est.an.ti.uh.6.b.1=data.matrix(read.csv('est_an_ti_uh_6_b_1.csv',header=F))
est.an.ti.uh.7.b.1=data.matrix(read.csv('est_an_ti_uh_7_b_1.csv',header=F))
est.an.ti.uh.2.b.8=data.matrix(read.csv('est_an_ti_uh_2_b_8.csv',header=F))
est.an.ti.uh.3.b.8=data.matrix(read.csv('est_an_ti_uh_3_b_8.csv',header=F))
est.an.ti.uh.4.b.8=data.matrix(read.csv('est_an_ti_uh_4_b_8.csv',header=F))
est.an.ti.uh.5.b.8=data.matrix(read.csv('est_an_ti_uh_5_b_8.csv',header=F))
est.an.ti.uh.6.b.8=data.matrix(read.csv('est_an_ti_uh_6_b_8.csv',header=F))
est.an.ti.uh.7.b.8=data.matrix(read.csv('est_an_ti_uh_7_b_8.csv',header=F))
est.an.ti.uh.2.b.128=data.matrix(read.csv('est_an_ti_uh_2_b_128.csv',header=F))
est.an.ti.uh.3.b.128=data.matrix(read.csv('est_an_ti_uh_3_b_128.csv',header=F))
est.an.ti.uh.4.b.128=data.matrix(read.csv('est_an_ti_uh_4_b_128.csv',header=F))
est.an.ti.uh.5.b.128=data.matrix(read.csv('est_an_ti_uh_5_b_128.csv',header=F))
est.an.ti.uh.6.b.128=data.matrix(read.csv('est_an_ti_uh_6_b_128.csv',header=F))
est.an.ti.uh.7.b.128=data.matrix(read.csv('est_an_ti_uh_7_b_128.csv',header=F))




est.an.ti.us.2.s.1=data.matrix(read.csv('est_an_ti_us_2_s_1.csv',header=F))
est.an.ti.us.3.s.1=data.matrix(read.csv('est_an_ti_us_3_s_1.csv',header=F))
est.an.ti.us.4.s.1=data.matrix(read.csv('est_an_ti_us_4_s_1.csv',header=F))
est.an.ti.us.5.s.1=data.matrix(read.csv('est_an_ti_us_5_s_1.csv',header=F))
est.an.ti.us.6.s.1=data.matrix(read.csv('est_an_ti_us_6_s_1.csv',header=F))
est.an.ti.us.7.s.1=data.matrix(read.csv('est_an_ti_us_7_s_1.csv',header=F))
est.an.ti.us.2.s.8=data.matrix(read.csv('est_an_ti_us_2_s_8.csv',header=F))
est.an.ti.us.3.s.8=data.matrix(read.csv('est_an_ti_us_3_s_8.csv',header=F))
est.an.ti.us.4.s.8=data.matrix(read.csv('est_an_ti_us_4_s_8.csv',header=F))
est.an.ti.us.5.s.8=data.matrix(read.csv('est_an_ti_us_5_s_8.csv',header=F))
est.an.ti.us.6.s.8=data.matrix(read.csv('est_an_ti_us_6_s_8.csv',header=F))
est.an.ti.us.7.s.8=data.matrix(read.csv('est_an_ti_us_7_s_8.csv',header=F))
est.an.ti.us.2.s.128=data.matrix(read.csv('est_an_ti_us_2_s_128.csv',header=F))
est.an.ti.us.3.s.128=data.matrix(read.csv('est_an_ti_us_3_s_128.csv',header=F))
est.an.ti.us.4.s.128=data.matrix(read.csv('est_an_ti_us_4_s_128.csv',header=F))
est.an.ti.us.5.s.128=data.matrix(read.csv('est_an_ti_us_5_s_128.csv',header=F))
est.an.ti.us.6.s.128=data.matrix(read.csv('est_an_ti_us_6_s_128.csv',header=F))
est.an.ti.us.7.s.128=data.matrix(read.csv('est_an_ti_us_7_s_128.csv',header=F))



est.an.ti.us.2.ang.1=data.matrix(read.csv('est_an_ti_us_2_ang_1.csv',header=F))
est.an.ti.us.3.ang.1=data.matrix(read.csv('est_an_ti_us_3_ang_1.csv',header=F))
est.an.ti.us.4.ang.1=data.matrix(read.csv('est_an_ti_us_4_ang_1.csv',header=F))
est.an.ti.us.5.ang.1=data.matrix(read.csv('est_an_ti_us_5_ang_1.csv',header=F))
est.an.ti.us.6.ang.1=data.matrix(read.csv('est_an_ti_us_6_ang_1.csv',header=F))
est.an.ti.us.7.ang.1=data.matrix(read.csv('est_an_ti_us_7_ang_1.csv',header=F))
est.an.ti.us.2.ang.8=data.matrix(read.csv('est_an_ti_us_2_ang_8.csv',header=F))
est.an.ti.us.3.ang.8=data.matrix(read.csv('est_an_ti_us_3_ang_8.csv',header=F))
est.an.ti.us.4.ang.8=data.matrix(read.csv('est_an_ti_us_4_ang_8.csv',header=F))
est.an.ti.us.5.ang.8=data.matrix(read.csv('est_an_ti_us_5_ang_8.csv',header=F))
est.an.ti.us.6.ang.8=data.matrix(read.csv('est_an_ti_us_6_ang_8.csv',header=F))
est.an.ti.us.7.ang.8=data.matrix(read.csv('est_an_ti_us_7_ang_8.csv',header=F))
est.an.ti.us.2.ang.128=data.matrix(read.csv('est_an_ti_us_2_ang_128.csv',header=F))
est.an.ti.us.3.ang.128=data.matrix(read.csv('est_an_ti_us_3_ang_128.csv',header=F))
est.an.ti.us.4.ang.128=data.matrix(read.csv('est_an_ti_us_4_ang_128.csv',header=F))
est.an.ti.us.5.ang.128=data.matrix(read.csv('est_an_ti_us_5_ang_128.csv',header=F))
est.an.ti.us.6.ang.128=data.matrix(read.csv('est_an_ti_us_6_ang_128.csv',header=F))
est.an.ti.us.7.ang.128=data.matrix(read.csv('est_an_ti_us_7_ang_128.csv',header=F))



est.an.ti.us.2.hs.1=data.matrix(read.csv('est_an_ti_us_2_hs_1.csv',header=F))
est.an.ti.us.3.hs.1=data.matrix(read.csv('est_an_ti_us_3_hs_1.csv',header=F))
est.an.ti.us.4.hs.1=data.matrix(read.csv('est_an_ti_us_4_hs_1.csv',header=F))
est.an.ti.us.5.hs.1=data.matrix(read.csv('est_an_ti_us_5_hs_1.csv',header=F))
est.an.ti.us.6.hs.1=data.matrix(read.csv('est_an_ti_us_6_hs_1.csv',header=F))
est.an.ti.us.7.hs.1=data.matrix(read.csv('est_an_ti_us_7_hs_1.csv',header=F))
est.an.ti.us.2.hs.8=data.matrix(read.csv('est_an_ti_us_2_hs_8.csv',header=F))
est.an.ti.us.3.hs.8=data.matrix(read.csv('est_an_ti_us_3_hs_8.csv',header=F))
est.an.ti.us.4.hs.8=data.matrix(read.csv('est_an_ti_us_4_hs_8.csv',header=F))
est.an.ti.us.5.hs.8=data.matrix(read.csv('est_an_ti_us_5_hs_8.csv',header=F))
est.an.ti.us.6.hs.8=data.matrix(read.csv('est_an_ti_us_6_hs_8.csv',header=F))
est.an.ti.us.7.hs.8=data.matrix(read.csv('est_an_ti_us_7_hs_8.csv',header=F))
est.an.ti.us.2.hs.128=data.matrix(read.csv('est_an_ti_us_2_hs_128.csv',header=F))
est.an.ti.us.3.hs.128=data.matrix(read.csv('est_an_ti_us_3_hs_128.csv',header=F))
est.an.ti.us.4.hs.128=data.matrix(read.csv('est_an_ti_us_4_hs_128.csv',header=F))
est.an.ti.us.5.hs.128=data.matrix(read.csv('est_an_ti_us_5_hs_128.csv',header=F))
est.an.ti.us.6.hs.128=data.matrix(read.csv('est_an_ti_us_6_hs_128.csv',header=F))
est.an.ti.us.7.hs.128=data.matrix(read.csv('est_an_ti_us_7_hs_128.csv',header=F))



est.an.ti.us.2.bur.1=data.matrix(read.csv('est_an_ti_us_2_bur_1.csv',header=F))
est.an.ti.us.3.bur.1=data.matrix(read.csv('est_an_ti_us_3_bur_1.csv',header=F))
est.an.ti.us.4.bur.1=data.matrix(read.csv('est_an_ti_us_4_bur_1.csv',header=F))
est.an.ti.us.5.bur.1=data.matrix(read.csv('est_an_ti_us_5_bur_1.csv',header=F))
est.an.ti.us.6.bur.1=data.matrix(read.csv('est_an_ti_us_6_bur_1.csv',header=F))
est.an.ti.us.7.bur.1=data.matrix(read.csv('est_an_ti_us_7_bur_1.csv',header=F))
est.an.ti.us.2.bur.8=data.matrix(read.csv('est_an_ti_us_2_bur_8.csv',header=F))
est.an.ti.us.3.bur.8=data.matrix(read.csv('est_an_ti_us_3_bur_8.csv',header=F))
est.an.ti.us.4.bur.8=data.matrix(read.csv('est_an_ti_us_4_bur_8.csv',header=F))
est.an.ti.us.5.bur.8=data.matrix(read.csv('est_an_ti_us_5_bur_8.csv',header=F))
est.an.ti.us.6.bur.8=data.matrix(read.csv('est_an_ti_us_6_bur_8.csv',header=F))
est.an.ti.us.7.bur.8=data.matrix(read.csv('est_an_ti_us_7_bur_8.csv',header=F))
est.an.ti.us.2.bur.128=data.matrix(read.csv('est_an_ti_us_2_bur_128.csv',header=F))
est.an.ti.us.3.bur.128=data.matrix(read.csv('est_an_ti_us_3_bur_128.csv',header=F))
est.an.ti.us.4.bur.128=data.matrix(read.csv('est_an_ti_us_4_bur_128.csv',header=F))
est.an.ti.us.5.bur.128=data.matrix(read.csv('est_an_ti_us_5_bur_128.csv',header=F))
est.an.ti.us.6.bur.128=data.matrix(read.csv('est_an_ti_us_6_bur_128.csv',header=F))
est.an.ti.us.7.bur.128=data.matrix(read.csv('est_an_ti_us_7_bur_128.csv',header=F))


est.an.ti.us.2.cb.1=data.matrix(read.csv('est_an_ti_us_2_cb_1.csv',header=F))
est.an.ti.us.3.cb.1=data.matrix(read.csv('est_an_ti_us_3_cb_1.csv',header=F))
est.an.ti.us.4.cb.1=data.matrix(read.csv('est_an_ti_us_4_cb_1.csv',header=F))
est.an.ti.us.5.cb.1=data.matrix(read.csv('est_an_ti_us_5_cb_1.csv',header=F))
est.an.ti.us.6.cb.1=data.matrix(read.csv('est_an_ti_us_6_cb_1.csv',header=F))
est.an.ti.us.7.cb.1=data.matrix(read.csv('est_an_ti_us_7_cb_1.csv',header=F))
est.an.ti.us.2.cb.8=data.matrix(read.csv('est_an_ti_us_2_cb_8.csv',header=F))
est.an.ti.us.3.cb.8=data.matrix(read.csv('est_an_ti_us_3_cb_8.csv',header=F))
est.an.ti.us.4.cb.8=data.matrix(read.csv('est_an_ti_us_4_cb_8.csv',header=F))
est.an.ti.us.5.cb.8=data.matrix(read.csv('est_an_ti_us_5_cb_8.csv',header=F))
est.an.ti.us.6.cb.8=data.matrix(read.csv('est_an_ti_us_6_cb_8.csv',header=F))
est.an.ti.us.7.cb.8=data.matrix(read.csv('est_an_ti_us_7_cb_8.csv',header=F))
est.an.ti.us.2.cb.128=data.matrix(read.csv('est_an_ti_us_2_cb_128.csv',header=F))
est.an.ti.us.3.cb.128=data.matrix(read.csv('est_an_ti_us_3_cb_128.csv',header=F))
est.an.ti.us.4.cb.128=data.matrix(read.csv('est_an_ti_us_4_cb_128.csv',header=F))
est.an.ti.us.5.cb.128=data.matrix(read.csv('est_an_ti_us_5_cb_128.csv',header=F))
est.an.ti.us.6.cb.128=data.matrix(read.csv('est_an_ti_us_6_cb_128.csv',header=F))
est.an.ti.us.7.cb.128=data.matrix(read.csv('est_an_ti_us_7_cb_128.csv',header=F))


est.an.ti.us.2.b.1=data.matrix(read.csv('est_an_ti_us_2_b_1.csv',header=F))
est.an.ti.us.3.b.1=data.matrix(read.csv('est_an_ti_us_3_b_1.csv',header=F))
est.an.ti.us.4.b.1=data.matrix(read.csv('est_an_ti_us_4_b_1.csv',header=F))
est.an.ti.us.5.b.1=data.matrix(read.csv('est_an_ti_us_5_b_1.csv',header=F))
est.an.ti.us.6.b.1=data.matrix(read.csv('est_an_ti_us_6_b_1.csv',header=F))
est.an.ti.us.7.b.1=data.matrix(read.csv('est_an_ti_us_7_b_1.csv',header=F))
est.an.ti.us.2.b.8=data.matrix(read.csv('est_an_ti_us_2_b_8.csv',header=F))
est.an.ti.us.3.b.8=data.matrix(read.csv('est_an_ti_us_3_b_8.csv',header=F))
est.an.ti.us.4.b.8=data.matrix(read.csv('est_an_ti_us_4_b_8.csv',header=F))
est.an.ti.us.5.b.8=data.matrix(read.csv('est_an_ti_us_5_b_8.csv',header=F))
est.an.ti.us.6.b.8=data.matrix(read.csv('est_an_ti_us_6_b_8.csv',header=F))
est.an.ti.us.7.b.8=data.matrix(read.csv('est_an_ti_us_7_b_8.csv',header=F))
est.an.ti.us.2.b.128=data.matrix(read.csv('est_an_ti_us_2_b_128.csv',header=F))
est.an.ti.us.3.b.128=data.matrix(read.csv('est_an_ti_us_3_b_128.csv',header=F))
est.an.ti.us.4.b.128=data.matrix(read.csv('est_an_ti_us_4_b_128.csv',header=F))
est.an.ti.us.5.b.128=data.matrix(read.csv('est_an_ti_us_5_b_128.csv',header=F))
est.an.ti.us.6.b.128=data.matrix(read.csv('est_an_ti_us_6_b_128.csv',header=F))
est.an.ti.us.7.b.128=data.matrix(read.csv('est_an_ti_us_7_b_128.csv',header=F))





est.haar.h.2.s.1=data.matrix(read.csv('est_haar_h_2_s_1.csv',header=F))
est.haar.h.3.s.1=data.matrix(read.csv('est_haar_h_3_s_1.csv',header=F))
est.haar.h.4.s.1=data.matrix(read.csv('est_haar_h_4_s_1.csv',header=F))
est.haar.h.5.s.1=data.matrix(read.csv('est_haar_h_5_s_1.csv',header=F))
est.haar.h.6.s.1=data.matrix(read.csv('est_haar_h_6_s_1.csv',header=F))
est.haar.h.7.s.1=data.matrix(read.csv('est_haar_h_7_s_1.csv',header=F))
est.haar.h.2.s.8=data.matrix(read.csv('est_haar_h_2_s_8.csv',header=F))
est.haar.h.3.s.8=data.matrix(read.csv('est_haar_h_3_s_8.csv',header=F))
est.haar.h.4.s.8=data.matrix(read.csv('est_haar_h_4_s_8.csv',header=F))
est.haar.h.5.s.8=data.matrix(read.csv('est_haar_h_5_s_8.csv',header=F))
est.haar.h.6.s.8=data.matrix(read.csv('est_haar_h_6_s_8.csv',header=F))
est.haar.h.7.s.8=data.matrix(read.csv('est_haar_h_7_s_8.csv',header=F))
est.haar.h.2.s.128=data.matrix(read.csv('est_haar_h_2_s_128.csv',header=F))
est.haar.h.3.s.128=data.matrix(read.csv('est_haar_h_3_s_128.csv',header=F))
est.haar.h.4.s.128=data.matrix(read.csv('est_haar_h_4_s_128.csv',header=F))
est.haar.h.5.s.128=data.matrix(read.csv('est_haar_h_5_s_128.csv',header=F))
est.haar.h.6.s.128=data.matrix(read.csv('est_haar_h_6_s_128.csv',header=F))
est.haar.h.7.s.128=data.matrix(read.csv('est_haar_h_7_s_128.csv',header=F))


est.haar.h.2.ang.1=data.matrix(read.csv('est_haar_h_2_ang_1.csv',header=F))
est.haar.h.3.ang.1=data.matrix(read.csv('est_haar_h_3_ang_1.csv',header=F))
est.haar.h.4.ang.1=data.matrix(read.csv('est_haar_h_4_ang_1.csv',header=F))
est.haar.h.5.ang.1=data.matrix(read.csv('est_haar_h_5_ang_1.csv',header=F))
est.haar.h.6.ang.1=data.matrix(read.csv('est_haar_h_6_ang_1.csv',header=F))
est.haar.h.7.ang.1=data.matrix(read.csv('est_haar_h_7_ang_1.csv',header=F))
est.haar.h.2.ang.8=data.matrix(read.csv('est_haar_h_2_ang_8.csv',header=F))
est.haar.h.3.ang.8=data.matrix(read.csv('est_haar_h_3_ang_8.csv',header=F))
est.haar.h.4.ang.8=data.matrix(read.csv('est_haar_h_4_ang_8.csv',header=F))
est.haar.h.5.ang.8=data.matrix(read.csv('est_haar_h_5_ang_8.csv',header=F))
est.haar.h.6.ang.8=data.matrix(read.csv('est_haar_h_6_ang_8.csv',header=F))
est.haar.h.7.ang.8=data.matrix(read.csv('est_haar_h_7_ang_8.csv',header=F))
est.haar.h.2.ang.128=data.matrix(read.csv('est_haar_h_2_ang_128.csv',header=F))
est.haar.h.3.ang.128=data.matrix(read.csv('est_haar_h_3_ang_128.csv',header=F))
est.haar.h.4.ang.128=data.matrix(read.csv('est_haar_h_4_ang_128.csv',header=F))
est.haar.h.5.ang.128=data.matrix(read.csv('est_haar_h_5_ang_128.csv',header=F))
est.haar.h.6.ang.128=data.matrix(read.csv('est_haar_h_6_ang_128.csv',header=F))
est.haar.h.7.ang.128=data.matrix(read.csv('est_haar_h_7_ang_128.csv',header=F))


est.haar.h.2.hs.1=data.matrix(read.csv('est_haar_h_2_hs_1.csv',header=F))
est.haar.h.3.hs.1=data.matrix(read.csv('est_haar_h_3_hs_1.csv',header=F))
est.haar.h.4.hs.1=data.matrix(read.csv('est_haar_h_4_hs_1.csv',header=F))
est.haar.h.5.hs.1=data.matrix(read.csv('est_haar_h_5_hs_1.csv',header=F))
est.haar.h.6.hs.1=data.matrix(read.csv('est_haar_h_6_hs_1.csv',header=F))
est.haar.h.7.hs.1=data.matrix(read.csv('est_haar_h_7_hs_1.csv',header=F))
est.haar.h.2.hs.8=data.matrix(read.csv('est_haar_h_2_hs_8.csv',header=F))
est.haar.h.3.hs.8=data.matrix(read.csv('est_haar_h_3_hs_8.csv',header=F))
est.haar.h.4.hs.8=data.matrix(read.csv('est_haar_h_4_hs_8.csv',header=F))
est.haar.h.5.hs.8=data.matrix(read.csv('est_haar_h_5_hs_8.csv',header=F))
est.haar.h.6.hs.8=data.matrix(read.csv('est_haar_h_6_hs_8.csv',header=F))
est.haar.h.7.hs.8=data.matrix(read.csv('est_haar_h_7_hs_8.csv',header=F))
est.haar.h.2.hs.128=data.matrix(read.csv('est_haar_h_2_hs_128.csv',header=F))
est.haar.h.3.hs.128=data.matrix(read.csv('est_haar_h_3_hs_128.csv',header=F))
est.haar.h.4.hs.128=data.matrix(read.csv('est_haar_h_4_hs_128.csv',header=F))
est.haar.h.5.hs.128=data.matrix(read.csv('est_haar_h_5_hs_128.csv',header=F))
est.haar.h.6.hs.128=data.matrix(read.csv('est_haar_h_6_hs_128.csv',header=F))
est.haar.h.7.hs.128=data.matrix(read.csv('est_haar_h_7_hs_128.csv',header=F))


est.haar.h.2.bur.1=data.matrix(read.csv('est_haar_h_2_bur_1.csv',header=F))
est.haar.h.3.bur.1=data.matrix(read.csv('est_haar_h_3_bur_1.csv',header=F))
est.haar.h.4.bur.1=data.matrix(read.csv('est_haar_h_4_bur_1.csv',header=F))
est.haar.h.5.bur.1=data.matrix(read.csv('est_haar_h_5_bur_1.csv',header=F))
est.haar.h.6.bur.1=data.matrix(read.csv('est_haar_h_6_bur_1.csv',header=F))
est.haar.h.7.bur.1=data.matrix(read.csv('est_haar_h_7_bur_1.csv',header=F))
est.haar.h.2.bur.8=data.matrix(read.csv('est_haar_h_2_bur_8.csv',header=F))
est.haar.h.3.bur.8=data.matrix(read.csv('est_haar_h_3_bur_8.csv',header=F))
est.haar.h.4.bur.8=data.matrix(read.csv('est_haar_h_4_bur_8.csv',header=F))
est.haar.h.5.bur.8=data.matrix(read.csv('est_haar_h_5_bur_8.csv',header=F))
est.haar.h.6.bur.8=data.matrix(read.csv('est_haar_h_6_bur_8.csv',header=F))
est.haar.h.7.bur.8=data.matrix(read.csv('est_haar_h_7_bur_8.csv',header=F))
est.haar.h.2.bur.128=data.matrix(read.csv('est_haar_h_2_bur_128.csv',header=F))
est.haar.h.3.bur.128=data.matrix(read.csv('est_haar_h_3_bur_128.csv',header=F))
est.haar.h.4.bur.128=data.matrix(read.csv('est_haar_h_4_bur_128.csv',header=F))
est.haar.h.5.bur.128=data.matrix(read.csv('est_haar_h_5_bur_128.csv',header=F))
est.haar.h.6.bur.128=data.matrix(read.csv('est_haar_h_6_bur_128.csv',header=F))
est.haar.h.7.bur.128=data.matrix(read.csv('est_haar_h_7_bur_128.csv',header=F))


est.haar.h.2.cb.1=data.matrix(read.csv('est_haar_h_2_cb_1.csv',header=F))
est.haar.h.3.cb.1=data.matrix(read.csv('est_haar_h_3_cb_1.csv',header=F))
est.haar.h.4.cb.1=data.matrix(read.csv('est_haar_h_4_cb_1.csv',header=F))
est.haar.h.5.cb.1=data.matrix(read.csv('est_haar_h_5_cb_1.csv',header=F))
est.haar.h.6.cb.1=data.matrix(read.csv('est_haar_h_6_cb_1.csv',header=F))
est.haar.h.7.cb.1=data.matrix(read.csv('est_haar_h_7_cb_1.csv',header=F))
est.haar.h.2.cb.8=data.matrix(read.csv('est_haar_h_2_cb_8.csv',header=F))
est.haar.h.3.cb.8=data.matrix(read.csv('est_haar_h_3_cb_8.csv',header=F))
est.haar.h.4.cb.8=data.matrix(read.csv('est_haar_h_4_cb_8.csv',header=F))
est.haar.h.5.cb.8=data.matrix(read.csv('est_haar_h_5_cb_8.csv',header=F))
est.haar.h.6.cb.8=data.matrix(read.csv('est_haar_h_6_cb_8.csv',header=F))
est.haar.h.7.cb.8=data.matrix(read.csv('est_haar_h_7_cb_8.csv',header=F))
est.haar.h.2.cb.128=data.matrix(read.csv('est_haar_h_2_cb_128.csv',header=F))
est.haar.h.3.cb.128=data.matrix(read.csv('est_haar_h_3_cb_128.csv',header=F))
est.haar.h.4.cb.128=data.matrix(read.csv('est_haar_h_4_cb_128.csv',header=F))
est.haar.h.5.cb.128=data.matrix(read.csv('est_haar_h_5_cb_128.csv',header=F))
est.haar.h.6.cb.128=data.matrix(read.csv('est_haar_h_6_cb_128.csv',header=F))
est.haar.h.7.cb.128=data.matrix(read.csv('est_haar_h_7_cb_128.csv',header=F))


est.haar.h.2.b.1=data.matrix(read.csv('est_haar_h_2_b_1.csv',header=F))
est.haar.h.3.b.1=data.matrix(read.csv('est_haar_h_3_b_1.csv',header=F))
est.haar.h.4.b.1=data.matrix(read.csv('est_haar_h_4_b_1.csv',header=F))
est.haar.h.5.b.1=data.matrix(read.csv('est_haar_h_5_b_1.csv',header=F))
est.haar.h.6.b.1=data.matrix(read.csv('est_haar_h_6_b_1.csv',header=F))
est.haar.h.7.b.1=data.matrix(read.csv('est_haar_h_7_b_1.csv',header=F))
est.haar.h.2.b.8=data.matrix(read.csv('est_haar_h_2_b_8.csv',header=F))
est.haar.h.3.b.8=data.matrix(read.csv('est_haar_h_3_b_8.csv',header=F))
est.haar.h.4.b.8=data.matrix(read.csv('est_haar_h_4_b_8.csv',header=F))
est.haar.h.5.b.8=data.matrix(read.csv('est_haar_h_5_b_8.csv',header=F))
est.haar.h.6.b.8=data.matrix(read.csv('est_haar_h_6_b_8.csv',header=F))
est.haar.h.7.b.8=data.matrix(read.csv('est_haar_h_7_b_8.csv',header=F))
est.haar.h.2.b.128=data.matrix(read.csv('est_haar_h_2_b_128.csv',header=F))
est.haar.h.3.b.128=data.matrix(read.csv('est_haar_h_3_b_128.csv',header=F))
est.haar.h.4.b.128=data.matrix(read.csv('est_haar_h_4_b_128.csv',header=F))
est.haar.h.5.b.128=data.matrix(read.csv('est_haar_h_5_b_128.csv',header=F))
est.haar.h.6.b.128=data.matrix(read.csv('est_haar_h_6_b_128.csv',header=F))
est.haar.h.7.b.128=data.matrix(read.csv('est_haar_h_7_b_128.csv',header=F))




est.pen.s.1=data.matrix(read.csv('est_pen_s_1.csv',header=F))
est.pen.s.8=data.matrix(read.csv('est_pen_s_8.csv',header=F))
est.pen.s.128=data.matrix(read.csv('est_pen_s_128.csv',header=F))

est.pen.ang.1=data.matrix(read.csv('est_pen_ang_1.csv',header=F))
est.pen.ang.8=data.matrix(read.csv('est_pen_ang_8.csv',header=F))
est.pen.ang.128=data.matrix(read.csv('est_pen_ang_128.csv',header=F))

est.pen.hs.1=data.matrix(read.csv('est_pen_hs_1.csv',header=F))
est.pen.hs.8=data.matrix(read.csv('est_pen_hs_8.csv',header=F))
est.pen.hs.128=data.matrix(read.csv('est_pen_hs_128.csv',header=F))

est.pen.bur.1=data.matrix(read.csv('est_pen_bur_1.csv',header=F))
est.pen.bur.8=data.matrix(read.csv('est_pen_bur_8.csv',header=F))
est.pen.bur.128=data.matrix(read.csv('est_pen_bur_128.csv',header=F))

est.pen.cb.1=data.matrix(read.csv('est_pen_cb_1.csv',header=F))
est.pen.cb.8=data.matrix(read.csv('est_pen_cb_8.csv',header=F))
est.pen.cb.128=data.matrix(read.csv('est_pen_cb_128.csv',header=F))

est.pen.b.1=data.matrix(read.csv('est_pen_b_1.csv',header=F))
est.pen.b.8=data.matrix(read.csv('est_pen_b_8.csv',header=F))
est.pen.b.128=data.matrix(read.csv('est_pen_b_128.csv',header=F))



est.BMSM.s.1=data.matrix(read.csv('est_BMSM_s_1.csv',header=F))
est.BMSM.s.8=data.matrix(read.csv('est_BMSM_s_8.csv',header=F))
est.BMSM.s.128=data.matrix(read.csv('est_BMSM_s_128.csv',header=F))

est.BMSM.ang.1=data.matrix(read.csv('est_BMSM_ang_1.csv',header=F))
est.BMSM.ang.8=data.matrix(read.csv('est_BMSM_ang_8.csv',header=F))
est.BMSM.ang.128=data.matrix(read.csv('est_BMSM_ang_128.csv',header=F))

est.BMSM.hs.1=data.matrix(read.csv('est_BMSM_hs_1.csv',header=F))
est.BMSM.hs.8=data.matrix(read.csv('est_BMSM_hs_8.csv',header=F))
est.BMSM.hs.128=data.matrix(read.csv('est_BMSM_hs_128.csv',header=F))

est.BMSM.bur.1=data.matrix(read.csv('est_BMSM_bur_1.csv',header=F))
est.BMSM.bur.8=data.matrix(read.csv('est_BMSM_bur_8.csv',header=F))
est.BMSM.bur.128=data.matrix(read.csv('est_BMSM_bur_128.csv',header=F))

est.BMSM.cb.1=data.matrix(read.csv('est_BMSM_cb_1.csv',header=F))
est.BMSM.cb.8=data.matrix(read.csv('est_BMSM_cb_8.csv',header=F))
est.BMSM.cb.128=data.matrix(read.csv('est_BMSM_cb_128.csv',header=F))

est.BMSM.b.1=data.matrix(read.csv('est_BMSM_b_1.csv',header=F))
est.BMSM.b.8=data.matrix(read.csv('est_BMSM_b_8.csv',header=F))
est.BMSM.b.128=data.matrix(read.csv('est_BMSM_b_128.csv',header=F))



est.BMIE.0.s.1=data.matrix(read.csv('est_BMIE_0_s_1.csv',header=F))
est.BMIE.2.s.1=data.matrix(read.csv('est_BMIE_2_s_1.csv',header=F))
est.BMIE.3.s.1=data.matrix(read.csv('est_BMIE_3_s_1.csv',header=F))
est.BMIE.4.s.1=data.matrix(read.csv('est_BMIE_4_s_1.csv',header=F))
est.BMIE.5.s.1=data.matrix(read.csv('est_BMIE_5_s_1.csv',header=F))
est.BMIE.6.s.1=data.matrix(read.csv('est_BMIE_6_s_1.csv',header=F))
est.BMIE.7.s.1=data.matrix(read.csv('est_BMIE_7_s_1.csv',header=F))

est.BMIE.0.s.8=data.matrix(read.csv('est_BMIE_0_s_8.csv',header=F))
est.BMIE.2.s.8=data.matrix(read.csv('est_BMIE_2_s_8.csv',header=F))
est.BMIE.3.s.8=data.matrix(read.csv('est_BMIE_3_s_8.csv',header=F))
est.BMIE.4.s.8=data.matrix(read.csv('est_BMIE_4_s_8.csv',header=F))
est.BMIE.5.s.8=data.matrix(read.csv('est_BMIE_5_s_8.csv',header=F))
est.BMIE.6.s.8=data.matrix(read.csv('est_BMIE_6_s_8.csv',header=F))
est.BMIE.7.s.8=data.matrix(read.csv('est_BMIE_7_s_8.csv',header=F))

est.BMIE.0.s.128=data.matrix(read.csv('est_BMIE_0_s_128.csv',header=F))
est.BMIE.2.s.128=data.matrix(read.csv('est_BMIE_2_s_128.csv',header=F))
est.BMIE.3.s.128=data.matrix(read.csv('est_BMIE_3_s_128.csv',header=F))
est.BMIE.4.s.128=data.matrix(read.csv('est_BMIE_4_s_128.csv',header=F))
est.BMIE.5.s.128=data.matrix(read.csv('est_BMIE_5_s_128.csv',header=F))
est.BMIE.6.s.128=data.matrix(read.csv('est_BMIE_6_s_128.csv',header=F))
est.BMIE.7.s.128=data.matrix(read.csv('est_BMIE_7_s_128.csv',header=F))



est.BMIE.0.ang.1=data.matrix(read.csv('est_BMIE_0_ang_1.csv',header=F))
est.BMIE.2.ang.1=data.matrix(read.csv('est_BMIE_2_ang_1.csv',header=F))
est.BMIE.3.ang.1=data.matrix(read.csv('est_BMIE_3_ang_1.csv',header=F))
est.BMIE.4.ang.1=data.matrix(read.csv('est_BMIE_4_ang_1.csv',header=F))
est.BMIE.5.ang.1=data.matrix(read.csv('est_BMIE_5_ang_1.csv',header=F))
est.BMIE.6.ang.1=data.matrix(read.csv('est_BMIE_6_ang_1.csv',header=F))
est.BMIE.7.ang.1=data.matrix(read.csv('est_BMIE_7_ang_1.csv',header=F))

est.BMIE.0.ang.8=data.matrix(read.csv('est_BMIE_0_ang_8.csv',header=F))
est.BMIE.2.ang.8=data.matrix(read.csv('est_BMIE_2_ang_8.csv',header=F))
est.BMIE.3.ang.8=data.matrix(read.csv('est_BMIE_3_ang_8.csv',header=F))
est.BMIE.4.ang.8=data.matrix(read.csv('est_BMIE_4_ang_8.csv',header=F))
est.BMIE.5.ang.8=data.matrix(read.csv('est_BMIE_5_ang_8.csv',header=F))
est.BMIE.6.ang.8=data.matrix(read.csv('est_BMIE_6_ang_8.csv',header=F))
est.BMIE.7.ang.8=data.matrix(read.csv('est_BMIE_7_ang_8.csv',header=F))

est.BMIE.0.ang.128=data.matrix(read.csv('est_BMIE_0_ang_128.csv',header=F))
est.BMIE.2.ang.128=data.matrix(read.csv('est_BMIE_2_ang_128.csv',header=F))
est.BMIE.3.ang.128=data.matrix(read.csv('est_BMIE_3_ang_128.csv',header=F))
est.BMIE.4.ang.128=data.matrix(read.csv('est_BMIE_4_ang_128.csv',header=F))
est.BMIE.5.ang.128=data.matrix(read.csv('est_BMIE_5_ang_128.csv',header=F))
est.BMIE.6.ang.128=data.matrix(read.csv('est_BMIE_6_ang_128.csv',header=F))
est.BMIE.7.ang.128=data.matrix(read.csv('est_BMIE_7_ang_128.csv',header=F))


est.BMIE.0.hs.1=data.matrix(read.csv('est_BMIE_0_hs_1.csv',header=F))
est.BMIE.2.hs.1=data.matrix(read.csv('est_BMIE_2_hs_1.csv',header=F))
est.BMIE.3.hs.1=data.matrix(read.csv('est_BMIE_3_hs_1.csv',header=F))
est.BMIE.4.hs.1=data.matrix(read.csv('est_BMIE_4_hs_1.csv',header=F))
est.BMIE.5.hs.1=data.matrix(read.csv('est_BMIE_5_hs_1.csv',header=F))
est.BMIE.6.hs.1=data.matrix(read.csv('est_BMIE_6_hs_1.csv',header=F))
est.BMIE.7.hs.1=data.matrix(read.csv('est_BMIE_7_hs_1.csv',header=F))

est.BMIE.0.hs.8=data.matrix(read.csv('est_BMIE_0_hs_8.csv',header=F))
est.BMIE.2.hs.8=data.matrix(read.csv('est_BMIE_2_hs_8.csv',header=F))
est.BMIE.3.hs.8=data.matrix(read.csv('est_BMIE_3_hs_8.csv',header=F))
est.BMIE.4.hs.8=data.matrix(read.csv('est_BMIE_4_hs_8.csv',header=F))
est.BMIE.5.hs.8=data.matrix(read.csv('est_BMIE_5_hs_8.csv',header=F))
est.BMIE.6.hs.8=data.matrix(read.csv('est_BMIE_6_hs_8.csv',header=F))
est.BMIE.7.hs.8=data.matrix(read.csv('est_BMIE_7_hs_8.csv',header=F))

est.BMIE.0.hs.128=data.matrix(read.csv('est_BMIE_0_hs_128.csv',header=F))
est.BMIE.2.hs.128=data.matrix(read.csv('est_BMIE_2_hs_128.csv',header=F))
est.BMIE.3.hs.128=data.matrix(read.csv('est_BMIE_3_hs_128.csv',header=F))
est.BMIE.4.hs.128=data.matrix(read.csv('est_BMIE_4_hs_128.csv',header=F))
est.BMIE.5.hs.128=data.matrix(read.csv('est_BMIE_5_hs_128.csv',header=F))
est.BMIE.6.hs.128=data.matrix(read.csv('est_BMIE_6_hs_128.csv',header=F))
est.BMIE.7.hs.128=data.matrix(read.csv('est_BMIE_7_hs_128.csv',header=F))



est.BMIE.0.bur.1=data.matrix(read.csv('est_BMIE_0_bur_1.csv',header=F))
est.BMIE.2.bur.1=data.matrix(read.csv('est_BMIE_2_bur_1.csv',header=F))
est.BMIE.3.bur.1=data.matrix(read.csv('est_BMIE_3_bur_1.csv',header=F))
est.BMIE.4.bur.1=data.matrix(read.csv('est_BMIE_4_bur_1.csv',header=F))
est.BMIE.5.bur.1=data.matrix(read.csv('est_BMIE_5_bur_1.csv',header=F))
est.BMIE.6.bur.1=data.matrix(read.csv('est_BMIE_6_bur_1.csv',header=F))
est.BMIE.7.bur.1=data.matrix(read.csv('est_BMIE_7_bur_1.csv',header=F))

est.BMIE.0.bur.8=data.matrix(read.csv('est_BMIE_0_bur_8.csv',header=F))
est.BMIE.2.bur.8=data.matrix(read.csv('est_BMIE_2_bur_8.csv',header=F))
est.BMIE.3.bur.8=data.matrix(read.csv('est_BMIE_3_bur_8.csv',header=F))
est.BMIE.4.bur.8=data.matrix(read.csv('est_BMIE_4_bur_8.csv',header=F))
est.BMIE.5.bur.8=data.matrix(read.csv('est_BMIE_5_bur_8.csv',header=F))
est.BMIE.6.bur.8=data.matrix(read.csv('est_BMIE_6_bur_8.csv',header=F))
est.BMIE.7.bur.8=data.matrix(read.csv('est_BMIE_7_bur_8.csv',header=F))

est.BMIE.0.bur.128=data.matrix(read.csv('est_BMIE_0_bur_128.csv',header=F))
est.BMIE.2.bur.128=data.matrix(read.csv('est_BMIE_2_bur_128.csv',header=F))
est.BMIE.3.bur.128=data.matrix(read.csv('est_BMIE_3_bur_128.csv',header=F))
est.BMIE.4.bur.128=data.matrix(read.csv('est_BMIE_4_bur_128.csv',header=F))
est.BMIE.5.bur.128=data.matrix(read.csv('est_BMIE_5_bur_128.csv',header=F))
est.BMIE.6.bur.128=data.matrix(read.csv('est_BMIE_6_bur_128.csv',header=F))
est.BMIE.7.bur.128=data.matrix(read.csv('est_BMIE_7_bur_128.csv',header=F))


est.BMIE.0.cb.1=data.matrix(read.csv('est_BMIE_0_cb_1.csv',header=F))
est.BMIE.2.cb.1=data.matrix(read.csv('est_BMIE_2_cb_1.csv',header=F))
est.BMIE.3.cb.1=data.matrix(read.csv('est_BMIE_3_cb_1.csv',header=F))
est.BMIE.4.cb.1=data.matrix(read.csv('est_BMIE_4_cb_1.csv',header=F))
est.BMIE.5.cb.1=data.matrix(read.csv('est_BMIE_5_cb_1.csv',header=F))
est.BMIE.6.cb.1=data.matrix(read.csv('est_BMIE_6_cb_1.csv',header=F))
est.BMIE.7.cb.1=data.matrix(read.csv('est_BMIE_7_cb_1.csv',header=F))

est.BMIE.0.cb.8=data.matrix(read.csv('est_BMIE_0_cb_8.csv',header=F))
est.BMIE.2.cb.8=data.matrix(read.csv('est_BMIE_2_cb_8.csv',header=F))
est.BMIE.3.cb.8=data.matrix(read.csv('est_BMIE_3_cb_8.csv',header=F))
est.BMIE.4.cb.8=data.matrix(read.csv('est_BMIE_4_cb_8.csv',header=F))
est.BMIE.5.cb.8=data.matrix(read.csv('est_BMIE_5_cb_8.csv',header=F))
est.BMIE.6.cb.8=data.matrix(read.csv('est_BMIE_6_cb_8.csv',header=F))
est.BMIE.7.cb.8=data.matrix(read.csv('est_BMIE_7_cb_8.csv',header=F))

est.BMIE.0.cb.128=data.matrix(read.csv('est_BMIE_0_cb_128.csv',header=F))
est.BMIE.2.cb.128=data.matrix(read.csv('est_BMIE_2_cb_128.csv',header=F))
est.BMIE.3.cb.128=data.matrix(read.csv('est_BMIE_3_cb_128.csv',header=F))
est.BMIE.4.cb.128=data.matrix(read.csv('est_BMIE_4_cb_128.csv',header=F))
est.BMIE.5.cb.128=data.matrix(read.csv('est_BMIE_5_cb_128.csv',header=F))
est.BMIE.6.cb.128=data.matrix(read.csv('est_BMIE_6_cb_128.csv',header=F))
est.BMIE.7.cb.128=data.matrix(read.csv('est_BMIE_7_cb_128.csv',header=F))


est.BMIE.0.b.1=data.matrix(read.csv('est_BMIE_0_b_1.csv',header=F))
est.BMIE.2.b.1=data.matrix(read.csv('est_BMIE_2_b_1.csv',header=F))
est.BMIE.3.b.1=data.matrix(read.csv('est_BMIE_3_b_1.csv',header=F))
est.BMIE.4.b.1=data.matrix(read.csv('est_BMIE_4_b_1.csv',header=F))
est.BMIE.5.b.1=data.matrix(read.csv('est_BMIE_5_b_1.csv',header=F))
est.BMIE.6.b.1=data.matrix(read.csv('est_BMIE_6_b_1.csv',header=F))
est.BMIE.7.b.1=data.matrix(read.csv('est_BMIE_7_b_1.csv',header=F))

est.BMIE.0.b.8=data.matrix(read.csv('est_BMIE_0_b_8.csv',header=F))
est.BMIE.2.b.8=data.matrix(read.csv('est_BMIE_2_b_8.csv',header=F))
est.BMIE.3.b.8=data.matrix(read.csv('est_BMIE_3_b_8.csv',header=F))
est.BMIE.4.b.8=data.matrix(read.csv('est_BMIE_4_b_8.csv',header=F))
est.BMIE.5.b.8=data.matrix(read.csv('est_BMIE_5_b_8.csv',header=F))
est.BMIE.6.b.8=data.matrix(read.csv('est_BMIE_6_b_8.csv',header=F))
est.BMIE.7.b.8=data.matrix(read.csv('est_BMIE_7_b_8.csv',header=F))

est.BMIE.0.b.128=data.matrix(read.csv('est_BMIE_0_b_128.csv',header=F))
est.BMIE.2.b.128=data.matrix(read.csv('est_BMIE_2_b_128.csv',header=F))
est.BMIE.3.b.128=data.matrix(read.csv('est_BMIE_3_b_128.csv',header=F))
est.BMIE.4.b.128=data.matrix(read.csv('est_BMIE_4_b_128.csv',header=F))
est.BMIE.5.b.128=data.matrix(read.csv('est_BMIE_5_b_128.csv',header=F))
est.BMIE.6.b.128=data.matrix(read.csv('est_BMIE_6_b_128.csv',header=F))
est.BMIE.7.b.128=data.matrix(read.csv('est_BMIE_7_b_128.csv',header=F))




#compute mises

mise.ash.s.1=mise(t(est.ash.s.1),mu.s.1)


mise.ash.s.8=mise(t(est.ash.s.8),mu.s.8)


mise.ash.s.128=mise(t(est.ash.s.128),mu.s.128)




mise.ash.ang.1=mise(t(est.ash.ang.1),mu.ang.1)


mise.ash.ang.8=mise(t(est.ash.ang.8),mu.ang.8)


mise.ash.ang.128=mise(t(est.ash.ang.128),mu.ang.128)




mise.ash.hs.1=mise(t(est.ash.hs.1),mu.hs.1)


mise.ash.hs.8=mise(t(est.ash.hs.8),mu.hs.8)


mise.ash.hs.128=mise(t(est.ash.hs.128),mu.hs.128)




mise.ash.bur.1=mise(t(est.ash.bur.1),mu.bur.1)


mise.ash.bur.8=mise(t(est.ash.bur.8),mu.bur.8)


mise.ash.bur.128=mise(t(est.ash.bur.128),mu.bur.128)





mise.ash.cb.1=mise(t(est.ash.cb.1),mu.cb.1)


mise.ash.cb.8=mise(t(est.ash.cb.8),mu.cb.8)


mise.ash.cb.128=mise(t(est.ash.cb.128),mu.cb.128)





mise.ash.b.1=mise(t(est.ash.b.1),mu.b.1)


mise.ash.b.8=mise(t(est.ash.b.8),mu.b.8)


mise.ash.b.128=mise(t(est.ash.b.128),mu.b.128)


mise.hf.ti.r.2.s.1=mise(t(est.hf.ti.r.2.s.1),mu.s.1)
mise.hf.ti.r.3.s.1=mise(t(est.hf.ti.r.3.s.1),mu.s.1)
mise.hf.ti.r.4.s.1=mise(t(est.hf.ti.r.4.s.1),mu.s.1)
mise.hf.ti.r.5.s.1=mise(t(est.hf.ti.r.5.s.1),mu.s.1)
mise.hf.ti.r.6.s.1=mise(t(est.hf.ti.r.6.s.1),mu.s.1)
mise.hf.ti.r.7.s.1=mise(t(est.hf.ti.r.7.s.1),mu.s.1)

mise.hf.ti.r.2.s.8=mise(t(est.hf.ti.r.2.s.8),mu.s.8)
mise.hf.ti.r.3.s.8=mise(t(est.hf.ti.r.3.s.8),mu.s.8)
mise.hf.ti.r.4.s.8=mise(t(est.hf.ti.r.4.s.8),mu.s.8)
mise.hf.ti.r.5.s.8=mise(t(est.hf.ti.r.5.s.8),mu.s.8)
mise.hf.ti.r.6.s.8=mise(t(est.hf.ti.r.6.s.8),mu.s.8)
mise.hf.ti.r.7.s.8=mise(t(est.hf.ti.r.7.s.8),mu.s.8)

mise.hf.ti.r.2.s.128=mise(t(est.hf.ti.r.2.s.128),mu.s.128)
mise.hf.ti.r.3.s.128=mise(t(est.hf.ti.r.3.s.128),mu.s.128)
mise.hf.ti.r.4.s.128=mise(t(est.hf.ti.r.4.s.128),mu.s.128)
mise.hf.ti.r.5.s.128=mise(t(est.hf.ti.r.5.s.128),mu.s.128)
mise.hf.ti.r.6.s.128=mise(t(est.hf.ti.r.6.s.128),mu.s.128)
mise.hf.ti.r.7.s.128=mise(t(est.hf.ti.r.7.s.128),mu.s.128)





mise.hf.ti.r.2.ang.1=mise(t(est.hf.ti.r.2.ang.1),mu.ang.1)
mise.hf.ti.r.3.ang.1=mise(t(est.hf.ti.r.3.ang.1),mu.ang.1)
mise.hf.ti.r.4.ang.1=mise(t(est.hf.ti.r.4.ang.1),mu.ang.1)
mise.hf.ti.r.5.ang.1=mise(t(est.hf.ti.r.5.ang.1),mu.ang.1)
mise.hf.ti.r.6.ang.1=mise(t(est.hf.ti.r.6.ang.1),mu.ang.1)
mise.hf.ti.r.7.ang.1=mise(t(est.hf.ti.r.7.ang.1),mu.ang.1)

mise.hf.ti.r.2.ang.8=mise(t(est.hf.ti.r.2.ang.8),mu.ang.8)
mise.hf.ti.r.3.ang.8=mise(t(est.hf.ti.r.3.ang.8),mu.ang.8)
mise.hf.ti.r.4.ang.8=mise(t(est.hf.ti.r.4.ang.8),mu.ang.8)
mise.hf.ti.r.5.ang.8=mise(t(est.hf.ti.r.5.ang.8),mu.ang.8)
mise.hf.ti.r.6.ang.8=mise(t(est.hf.ti.r.6.ang.8),mu.ang.8)
mise.hf.ti.r.7.ang.8=mise(t(est.hf.ti.r.7.ang.8),mu.ang.8)

mise.hf.ti.r.2.ang.128=mise(t(est.hf.ti.r.2.ang.128),mu.ang.128)
mise.hf.ti.r.3.ang.128=mise(t(est.hf.ti.r.3.ang.128),mu.ang.128)
mise.hf.ti.r.4.ang.128=mise(t(est.hf.ti.r.4.ang.128),mu.ang.128)
mise.hf.ti.r.5.ang.128=mise(t(est.hf.ti.r.5.ang.128),mu.ang.128)
mise.hf.ti.r.6.ang.128=mise(t(est.hf.ti.r.6.ang.128),mu.ang.128)
mise.hf.ti.r.7.ang.128=mise(t(est.hf.ti.r.7.ang.128),mu.ang.128)




mise.hf.ti.r.2.hs.1=mise(t(est.hf.ti.r.2.hs.1),mu.hs.1)
mise.hf.ti.r.3.hs.1=mise(t(est.hf.ti.r.3.hs.1),mu.hs.1)
mise.hf.ti.r.4.hs.1=mise(t(est.hf.ti.r.4.hs.1),mu.hs.1)
mise.hf.ti.r.5.hs.1=mise(t(est.hf.ti.r.5.hs.1),mu.hs.1)
mise.hf.ti.r.6.hs.1=mise(t(est.hf.ti.r.6.hs.1),mu.hs.1)
mise.hf.ti.r.7.hs.1=mise(t(est.hf.ti.r.7.hs.1),mu.hs.1)

mise.hf.ti.r.2.hs.8=mise(t(est.hf.ti.r.2.hs.8),mu.hs.8)
mise.hf.ti.r.3.hs.8=mise(t(est.hf.ti.r.3.hs.8),mu.hs.8)
mise.hf.ti.r.4.hs.8=mise(t(est.hf.ti.r.4.hs.8),mu.hs.8)
mise.hf.ti.r.5.hs.8=mise(t(est.hf.ti.r.5.hs.8),mu.hs.8)
mise.hf.ti.r.6.hs.8=mise(t(est.hf.ti.r.6.hs.8),mu.hs.8)
mise.hf.ti.r.7.hs.8=mise(t(est.hf.ti.r.7.hs.8),mu.hs.8)

mise.hf.ti.r.2.hs.128=mise(t(est.hf.ti.r.2.hs.128),mu.hs.128)
mise.hf.ti.r.3.hs.128=mise(t(est.hf.ti.r.3.hs.128),mu.hs.128)
mise.hf.ti.r.4.hs.128=mise(t(est.hf.ti.r.4.hs.128),mu.hs.128)
mise.hf.ti.r.5.hs.128=mise(t(est.hf.ti.r.5.hs.128),mu.hs.128)
mise.hf.ti.r.6.hs.128=mise(t(est.hf.ti.r.6.hs.128),mu.hs.128)
mise.hf.ti.r.7.hs.128=mise(t(est.hf.ti.r.7.hs.128),mu.hs.128)




mise.hf.ti.r.2.bur.1=mise(t(est.hf.ti.r.2.bur.1),mu.bur.1)
mise.hf.ti.r.3.bur.1=mise(t(est.hf.ti.r.3.bur.1),mu.bur.1)
mise.hf.ti.r.4.bur.1=mise(t(est.hf.ti.r.4.bur.1),mu.bur.1)
mise.hf.ti.r.5.bur.1=mise(t(est.hf.ti.r.5.bur.1),mu.bur.1)
mise.hf.ti.r.6.bur.1=mise(t(est.hf.ti.r.6.bur.1),mu.bur.1)
mise.hf.ti.r.7.bur.1=mise(t(est.hf.ti.r.7.bur.1),mu.bur.1)

mise.hf.ti.r.2.bur.8=mise(t(est.hf.ti.r.2.bur.8),mu.bur.8)
mise.hf.ti.r.3.bur.8=mise(t(est.hf.ti.r.3.bur.8),mu.bur.8)
mise.hf.ti.r.4.bur.8=mise(t(est.hf.ti.r.4.bur.8),mu.bur.8)
mise.hf.ti.r.5.bur.8=mise(t(est.hf.ti.r.5.bur.8),mu.bur.8)
mise.hf.ti.r.6.bur.8=mise(t(est.hf.ti.r.6.bur.8),mu.bur.8)
mise.hf.ti.r.7.bur.8=mise(t(est.hf.ti.r.7.bur.8),mu.bur.8)

mise.hf.ti.r.2.bur.128=mise(t(est.hf.ti.r.2.bur.128),mu.bur.128)
mise.hf.ti.r.3.bur.128=mise(t(est.hf.ti.r.3.bur.128),mu.bur.128)
mise.hf.ti.r.4.bur.128=mise(t(est.hf.ti.r.4.bur.128),mu.bur.128)
mise.hf.ti.r.5.bur.128=mise(t(est.hf.ti.r.5.bur.128),mu.bur.128)
mise.hf.ti.r.6.bur.128=mise(t(est.hf.ti.r.6.bur.128),mu.bur.128)
mise.hf.ti.r.7.bur.128=mise(t(est.hf.ti.r.7.bur.128),mu.bur.128)




mise.hf.ti.r.2.cb.1=mise(t(est.hf.ti.r.2.cb.1),mu.cb.1)
mise.hf.ti.r.3.cb.1=mise(t(est.hf.ti.r.3.cb.1),mu.cb.1)
mise.hf.ti.r.4.cb.1=mise(t(est.hf.ti.r.4.cb.1),mu.cb.1)
mise.hf.ti.r.5.cb.1=mise(t(est.hf.ti.r.5.cb.1),mu.cb.1)
mise.hf.ti.r.6.cb.1=mise(t(est.hf.ti.r.6.cb.1),mu.cb.1)
mise.hf.ti.r.7.cb.1=mise(t(est.hf.ti.r.7.cb.1),mu.cb.1)

mise.hf.ti.r.2.cb.8=mise(t(est.hf.ti.r.2.cb.8),mu.cb.8)
mise.hf.ti.r.3.cb.8=mise(t(est.hf.ti.r.3.cb.8),mu.cb.8)
mise.hf.ti.r.4.cb.8=mise(t(est.hf.ti.r.4.cb.8),mu.cb.8)
mise.hf.ti.r.5.cb.8=mise(t(est.hf.ti.r.5.cb.8),mu.cb.8)
mise.hf.ti.r.6.cb.8=mise(t(est.hf.ti.r.6.cb.8),mu.cb.8)
mise.hf.ti.r.7.cb.8=mise(t(est.hf.ti.r.7.cb.8),mu.cb.8)

mise.hf.ti.r.2.cb.128=mise(t(est.hf.ti.r.2.cb.128),mu.cb.128)
mise.hf.ti.r.3.cb.128=mise(t(est.hf.ti.r.3.cb.128),mu.cb.128)
mise.hf.ti.r.4.cb.128=mise(t(est.hf.ti.r.4.cb.128),mu.cb.128)
mise.hf.ti.r.5.cb.128=mise(t(est.hf.ti.r.5.cb.128),mu.cb.128)
mise.hf.ti.r.6.cb.128=mise(t(est.hf.ti.r.6.cb.128),mu.cb.128)
mise.hf.ti.r.7.cb.128=mise(t(est.hf.ti.r.7.cb.128),mu.cb.128)




mise.hf.ti.r.2.b.1=mise(t(est.hf.ti.r.2.b.1),mu.b.1)
mise.hf.ti.r.3.b.1=mise(t(est.hf.ti.r.3.b.1),mu.b.1)
mise.hf.ti.r.4.b.1=mise(t(est.hf.ti.r.4.b.1),mu.b.1)
mise.hf.ti.r.5.b.1=mise(t(est.hf.ti.r.5.b.1),mu.b.1)
mise.hf.ti.r.6.b.1=mise(t(est.hf.ti.r.6.b.1),mu.b.1)
mise.hf.ti.r.7.b.1=mise(t(est.hf.ti.r.7.b.1),mu.b.1)

mise.hf.ti.r.2.b.8=mise(t(est.hf.ti.r.2.b.8),mu.b.8)
mise.hf.ti.r.3.b.8=mise(t(est.hf.ti.r.3.b.8),mu.b.8)
mise.hf.ti.r.4.b.8=mise(t(est.hf.ti.r.4.b.8),mu.b.8)
mise.hf.ti.r.5.b.8=mise(t(est.hf.ti.r.5.b.8),mu.b.8)
mise.hf.ti.r.6.b.8=mise(t(est.hf.ti.r.6.b.8),mu.b.8)
mise.hf.ti.r.7.b.8=mise(t(est.hf.ti.r.7.b.8),mu.b.8)

mise.hf.ti.r.2.b.128=mise(t(est.hf.ti.r.2.b.128),mu.b.128)
mise.hf.ti.r.3.b.128=mise(t(est.hf.ti.r.3.b.128),mu.b.128)
mise.hf.ti.r.4.b.128=mise(t(est.hf.ti.r.4.b.128),mu.b.128)
mise.hf.ti.r.5.b.128=mise(t(est.hf.ti.r.5.b.128),mu.b.128)
mise.hf.ti.r.6.b.128=mise(t(est.hf.ti.r.6.b.128),mu.b.128)
mise.hf.ti.r.7.b.128=mise(t(est.hf.ti.r.7.b.128),mu.b.128)







mise.an.ti.cvh.2.s.1=mise((est.an.ti.cvh.2.s.1),mu.s.1)
mise.an.ti.cvh.3.s.1=mise((est.an.ti.cvh.3.s.1),mu.s.1)
mise.an.ti.cvh.4.s.1=mise((est.an.ti.cvh.4.s.1),mu.s.1)
mise.an.ti.cvh.5.s.1=mise((est.an.ti.cvh.5.s.1),mu.s.1)
mise.an.ti.cvh.6.s.1=mise((est.an.ti.cvh.6.s.1),mu.s.1)
mise.an.ti.cvh.7.s.1=mise((est.an.ti.cvh.7.s.1),mu.s.1)

mise.an.ti.cvh.2.s.8=mise((est.an.ti.cvh.2.s.8),mu.s.8)
mise.an.ti.cvh.3.s.8=mise((est.an.ti.cvh.3.s.8),mu.s.8)
mise.an.ti.cvh.4.s.8=mise((est.an.ti.cvh.4.s.8),mu.s.8)
mise.an.ti.cvh.5.s.8=mise((est.an.ti.cvh.5.s.8),mu.s.8)
mise.an.ti.cvh.6.s.8=mise((est.an.ti.cvh.6.s.8),mu.s.8)
mise.an.ti.cvh.7.s.8=mise((est.an.ti.cvh.7.s.8),mu.s.8)

mise.an.ti.cvh.2.s.128=mise((est.an.ti.cvh.2.s.128),mu.s.128)
mise.an.ti.cvh.3.s.128=mise((est.an.ti.cvh.3.s.128),mu.s.128)
mise.an.ti.cvh.4.s.128=mise((est.an.ti.cvh.4.s.128),mu.s.128)
mise.an.ti.cvh.5.s.128=mise((est.an.ti.cvh.5.s.128),mu.s.128)
mise.an.ti.cvh.6.s.128=mise((est.an.ti.cvh.6.s.128),mu.s.128)
mise.an.ti.cvh.7.s.128=mise((est.an.ti.cvh.7.s.128),mu.s.128)





mise.an.ti.cvh.2.ang.1=mise((est.an.ti.cvh.2.ang.1),mu.ang.1)
mise.an.ti.cvh.3.ang.1=mise((est.an.ti.cvh.3.ang.1),mu.ang.1)
mise.an.ti.cvh.4.ang.1=mise((est.an.ti.cvh.4.ang.1),mu.ang.1)
mise.an.ti.cvh.5.ang.1=mise((est.an.ti.cvh.5.ang.1),mu.ang.1)
mise.an.ti.cvh.6.ang.1=mise((est.an.ti.cvh.6.ang.1),mu.ang.1)
mise.an.ti.cvh.7.ang.1=mise((est.an.ti.cvh.7.ang.1),mu.ang.1)

mise.an.ti.cvh.2.ang.8=mise((est.an.ti.cvh.2.ang.8),mu.ang.8)
mise.an.ti.cvh.3.ang.8=mise((est.an.ti.cvh.3.ang.8),mu.ang.8)
mise.an.ti.cvh.4.ang.8=mise((est.an.ti.cvh.4.ang.8),mu.ang.8)
mise.an.ti.cvh.5.ang.8=mise((est.an.ti.cvh.5.ang.8),mu.ang.8)
mise.an.ti.cvh.6.ang.8=mise((est.an.ti.cvh.6.ang.8),mu.ang.8)
mise.an.ti.cvh.7.ang.8=mise((est.an.ti.cvh.7.ang.8),mu.ang.8)

mise.an.ti.cvh.2.ang.128=mise((est.an.ti.cvh.2.ang.128),mu.ang.128)
mise.an.ti.cvh.3.ang.128=mise((est.an.ti.cvh.3.ang.128),mu.ang.128)
mise.an.ti.cvh.4.ang.128=mise((est.an.ti.cvh.4.ang.128),mu.ang.128)
mise.an.ti.cvh.5.ang.128=mise((est.an.ti.cvh.5.ang.128),mu.ang.128)
mise.an.ti.cvh.6.ang.128=mise((est.an.ti.cvh.6.ang.128),mu.ang.128)
mise.an.ti.cvh.7.ang.128=mise((est.an.ti.cvh.7.ang.128),mu.ang.128)



mise.an.ti.cvh.2.hs.1=mise((est.an.ti.cvh.2.hs.1),mu.hs.1)
mise.an.ti.cvh.3.hs.1=mise((est.an.ti.cvh.3.hs.1),mu.hs.1)
mise.an.ti.cvh.4.hs.1=mise((est.an.ti.cvh.4.hs.1),mu.hs.1)
mise.an.ti.cvh.5.hs.1=mise((est.an.ti.cvh.5.hs.1),mu.hs.1)
mise.an.ti.cvh.6.hs.1=mise((est.an.ti.cvh.6.hs.1),mu.hs.1)
mise.an.ti.cvh.7.hs.1=mise((est.an.ti.cvh.7.hs.1),mu.hs.1)

mise.an.ti.cvh.2.hs.8=mise((est.an.ti.cvh.2.hs.8),mu.hs.8)
mise.an.ti.cvh.3.hs.8=mise((est.an.ti.cvh.3.hs.8),mu.hs.8)
mise.an.ti.cvh.4.hs.8=mise((est.an.ti.cvh.4.hs.8),mu.hs.8)
mise.an.ti.cvh.5.hs.8=mise((est.an.ti.cvh.5.hs.8),mu.hs.8)
mise.an.ti.cvh.6.hs.8=mise((est.an.ti.cvh.6.hs.8),mu.hs.8)
mise.an.ti.cvh.7.hs.8=mise((est.an.ti.cvh.7.hs.8),mu.hs.8)

mise.an.ti.cvh.2.hs.128=mise((est.an.ti.cvh.2.hs.128),mu.hs.128)
mise.an.ti.cvh.3.hs.128=mise((est.an.ti.cvh.3.hs.128),mu.hs.128)
mise.an.ti.cvh.4.hs.128=mise((est.an.ti.cvh.4.hs.128),mu.hs.128)
mise.an.ti.cvh.5.hs.128=mise((est.an.ti.cvh.5.hs.128),mu.hs.128)
mise.an.ti.cvh.6.hs.128=mise((est.an.ti.cvh.6.hs.128),mu.hs.128)
mise.an.ti.cvh.7.hs.128=mise((est.an.ti.cvh.7.hs.128),mu.hs.128)





mise.an.ti.cvh.2.bur.1=mise((est.an.ti.cvh.2.bur.1),mu.bur.1)
mise.an.ti.cvh.3.bur.1=mise((est.an.ti.cvh.3.bur.1),mu.bur.1)
mise.an.ti.cvh.4.bur.1=mise((est.an.ti.cvh.4.bur.1),mu.bur.1)
mise.an.ti.cvh.5.bur.1=mise((est.an.ti.cvh.5.bur.1),mu.bur.1)
mise.an.ti.cvh.6.bur.1=mise((est.an.ti.cvh.6.bur.1),mu.bur.1)
mise.an.ti.cvh.7.bur.1=mise((est.an.ti.cvh.7.bur.1),mu.bur.1)

mise.an.ti.cvh.2.bur.8=mise((est.an.ti.cvh.2.bur.8),mu.bur.8)
mise.an.ti.cvh.3.bur.8=mise((est.an.ti.cvh.3.bur.8),mu.bur.8)
mise.an.ti.cvh.4.bur.8=mise((est.an.ti.cvh.4.bur.8),mu.bur.8)
mise.an.ti.cvh.5.bur.8=mise((est.an.ti.cvh.5.bur.8),mu.bur.8)
mise.an.ti.cvh.6.bur.8=mise((est.an.ti.cvh.6.bur.8),mu.bur.8)
mise.an.ti.cvh.7.bur.8=mise((est.an.ti.cvh.7.bur.8),mu.bur.8)

mise.an.ti.cvh.2.bur.128=mise((est.an.ti.cvh.2.bur.128),mu.bur.128)
mise.an.ti.cvh.3.bur.128=mise((est.an.ti.cvh.3.bur.128),mu.bur.128)
mise.an.ti.cvh.4.bur.128=mise((est.an.ti.cvh.4.bur.128),mu.bur.128)
mise.an.ti.cvh.5.bur.128=mise((est.an.ti.cvh.5.bur.128),mu.bur.128)
mise.an.ti.cvh.6.bur.128=mise((est.an.ti.cvh.6.bur.128),mu.bur.128)
mise.an.ti.cvh.7.bur.128=mise((est.an.ti.cvh.7.bur.128),mu.bur.128)




mise.an.ti.cvh.2.cb.1=mise((est.an.ti.cvh.2.cb.1),mu.cb.1)
mise.an.ti.cvh.3.cb.1=mise((est.an.ti.cvh.3.cb.1),mu.cb.1)
mise.an.ti.cvh.4.cb.1=mise((est.an.ti.cvh.4.cb.1),mu.cb.1)
mise.an.ti.cvh.5.cb.1=mise((est.an.ti.cvh.5.cb.1),mu.cb.1)
mise.an.ti.cvh.6.cb.1=mise((est.an.ti.cvh.6.cb.1),mu.cb.1)
mise.an.ti.cvh.7.cb.1=mise((est.an.ti.cvh.7.cb.1),mu.cb.1)

mise.an.ti.cvh.2.cb.8=mise((est.an.ti.cvh.2.cb.8),mu.cb.8)
mise.an.ti.cvh.3.cb.8=mise((est.an.ti.cvh.3.cb.8),mu.cb.8)
mise.an.ti.cvh.4.cb.8=mise((est.an.ti.cvh.4.cb.8),mu.cb.8)
mise.an.ti.cvh.5.cb.8=mise((est.an.ti.cvh.5.cb.8),mu.cb.8)
mise.an.ti.cvh.6.cb.8=mise((est.an.ti.cvh.6.cb.8),mu.cb.8)
mise.an.ti.cvh.7.cb.8=mise((est.an.ti.cvh.7.cb.8),mu.cb.8)

mise.an.ti.cvh.2.cb.128=mise((est.an.ti.cvh.2.cb.128),mu.cb.128)
mise.an.ti.cvh.3.cb.128=mise((est.an.ti.cvh.3.cb.128),mu.cb.128)
mise.an.ti.cvh.4.cb.128=mise((est.an.ti.cvh.4.cb.128),mu.cb.128)
mise.an.ti.cvh.5.cb.128=mise((est.an.ti.cvh.5.cb.128),mu.cb.128)
mise.an.ti.cvh.6.cb.128=mise((est.an.ti.cvh.6.cb.128),mu.cb.128)
mise.an.ti.cvh.7.cb.128=mise((est.an.ti.cvh.7.cb.128),mu.cb.128)



mise.an.ti.cvh.2.b.1=mise((est.an.ti.cvh.2.b.1),mu.b.1)
mise.an.ti.cvh.3.b.1=mise((est.an.ti.cvh.3.b.1),mu.b.1)
mise.an.ti.cvh.4.b.1=mise((est.an.ti.cvh.4.b.1),mu.b.1)
mise.an.ti.cvh.5.b.1=mise((est.an.ti.cvh.5.b.1),mu.b.1)
mise.an.ti.cvh.6.b.1=mise((est.an.ti.cvh.6.b.1),mu.b.1)
mise.an.ti.cvh.7.b.1=mise((est.an.ti.cvh.7.b.1),mu.b.1)

mise.an.ti.cvh.2.b.8=mise((est.an.ti.cvh.2.b.8),mu.b.8)
mise.an.ti.cvh.3.b.8=mise((est.an.ti.cvh.3.b.8),mu.b.8)
mise.an.ti.cvh.4.b.8=mise((est.an.ti.cvh.4.b.8),mu.b.8)
mise.an.ti.cvh.5.b.8=mise((est.an.ti.cvh.5.b.8),mu.b.8)
mise.an.ti.cvh.6.b.8=mise((est.an.ti.cvh.6.b.8),mu.b.8)
mise.an.ti.cvh.7.b.8=mise((est.an.ti.cvh.7.b.8),mu.b.8)

mise.an.ti.cvh.2.b.128=mise((est.an.ti.cvh.2.b.128),mu.b.128)
mise.an.ti.cvh.3.b.128=mise((est.an.ti.cvh.3.b.128),mu.b.128)
mise.an.ti.cvh.4.b.128=mise((est.an.ti.cvh.4.b.128),mu.b.128)
mise.an.ti.cvh.5.b.128=mise((est.an.ti.cvh.5.b.128),mu.b.128)
mise.an.ti.cvh.6.b.128=mise((est.an.ti.cvh.6.b.128),mu.b.128)
mise.an.ti.cvh.7.b.128=mise((est.an.ti.cvh.7.b.128),mu.b.128)






mise.an.ti.uh.2.s.1=mise((est.an.ti.uh.2.s.1),mu.s.1)
mise.an.ti.uh.3.s.1=mise((est.an.ti.uh.3.s.1),mu.s.1)
mise.an.ti.uh.4.s.1=mise((est.an.ti.uh.4.s.1),mu.s.1)
mise.an.ti.uh.5.s.1=mise((est.an.ti.uh.5.s.1),mu.s.1)
mise.an.ti.uh.6.s.1=mise((est.an.ti.uh.6.s.1),mu.s.1)
mise.an.ti.uh.7.s.1=mise((est.an.ti.uh.7.s.1),mu.s.1)

mise.an.ti.uh.2.s.8=mise((est.an.ti.uh.2.s.8),mu.s.8)
mise.an.ti.uh.3.s.8=mise((est.an.ti.uh.3.s.8),mu.s.8)
mise.an.ti.uh.4.s.8=mise((est.an.ti.uh.4.s.8),mu.s.8)
mise.an.ti.uh.5.s.8=mise((est.an.ti.uh.5.s.8),mu.s.8)
mise.an.ti.uh.6.s.8=mise((est.an.ti.uh.6.s.8),mu.s.8)
mise.an.ti.uh.7.s.8=mise((est.an.ti.uh.7.s.8),mu.s.8)

mise.an.ti.uh.2.s.128=mise((est.an.ti.uh.2.s.128),mu.s.128)
mise.an.ti.uh.3.s.128=mise((est.an.ti.uh.3.s.128),mu.s.128)
mise.an.ti.uh.4.s.128=mise((est.an.ti.uh.4.s.128),mu.s.128)
mise.an.ti.uh.5.s.128=mise((est.an.ti.uh.5.s.128),mu.s.128)
mise.an.ti.uh.6.s.128=mise((est.an.ti.uh.6.s.128),mu.s.128)
mise.an.ti.uh.7.s.128=mise((est.an.ti.uh.7.s.128),mu.s.128)




mise.an.ti.uh.2.ang.1=mise((est.an.ti.uh.2.ang.1),mu.ang.1)
mise.an.ti.uh.3.ang.1=mise((est.an.ti.uh.3.ang.1),mu.ang.1)
mise.an.ti.uh.4.ang.1=mise((est.an.ti.uh.4.ang.1),mu.ang.1)
mise.an.ti.uh.5.ang.1=mise((est.an.ti.uh.5.ang.1),mu.ang.1)
mise.an.ti.uh.6.ang.1=mise((est.an.ti.uh.6.ang.1),mu.ang.1)
mise.an.ti.uh.7.ang.1=mise((est.an.ti.uh.7.ang.1),mu.ang.1)

mise.an.ti.uh.2.ang.8=mise((est.an.ti.uh.2.ang.8),mu.ang.8)
mise.an.ti.uh.3.ang.8=mise((est.an.ti.uh.3.ang.8),mu.ang.8)
mise.an.ti.uh.4.ang.8=mise((est.an.ti.uh.4.ang.8),mu.ang.8)
mise.an.ti.uh.5.ang.8=mise((est.an.ti.uh.5.ang.8),mu.ang.8)
mise.an.ti.uh.6.ang.8=mise((est.an.ti.uh.6.ang.8),mu.ang.8)
mise.an.ti.uh.7.ang.8=mise((est.an.ti.uh.7.ang.8),mu.ang.8)

mise.an.ti.uh.2.ang.128=mise((est.an.ti.uh.2.ang.128),mu.ang.128)
mise.an.ti.uh.3.ang.128=mise((est.an.ti.uh.3.ang.128),mu.ang.128)
mise.an.ti.uh.4.ang.128=mise((est.an.ti.uh.4.ang.128),mu.ang.128)
mise.an.ti.uh.5.ang.128=mise((est.an.ti.uh.5.ang.128),mu.ang.128)
mise.an.ti.uh.6.ang.128=mise((est.an.ti.uh.6.ang.128),mu.ang.128)
mise.an.ti.uh.7.ang.128=mise((est.an.ti.uh.7.ang.128),mu.ang.128)




mise.an.ti.uh.2.hs.1=mise((est.an.ti.uh.2.hs.1),mu.hs.1)
mise.an.ti.uh.3.hs.1=mise((est.an.ti.uh.3.hs.1),mu.hs.1)
mise.an.ti.uh.4.hs.1=mise((est.an.ti.uh.4.hs.1),mu.hs.1)
mise.an.ti.uh.5.hs.1=mise((est.an.ti.uh.5.hs.1),mu.hs.1)
mise.an.ti.uh.6.hs.1=mise((est.an.ti.uh.6.hs.1),mu.hs.1)
mise.an.ti.uh.7.hs.1=mise((est.an.ti.uh.7.hs.1),mu.hs.1)

mise.an.ti.uh.2.hs.8=mise((est.an.ti.uh.2.hs.8),mu.hs.8)
mise.an.ti.uh.3.hs.8=mise((est.an.ti.uh.3.hs.8),mu.hs.8)
mise.an.ti.uh.4.hs.8=mise((est.an.ti.uh.4.hs.8),mu.hs.8)
mise.an.ti.uh.5.hs.8=mise((est.an.ti.uh.5.hs.8),mu.hs.8)
mise.an.ti.uh.6.hs.8=mise((est.an.ti.uh.6.hs.8),mu.hs.8)
mise.an.ti.uh.7.hs.8=mise((est.an.ti.uh.7.hs.8),mu.hs.8)

mise.an.ti.uh.2.hs.128=mise((est.an.ti.uh.2.hs.128),mu.hs.128)
mise.an.ti.uh.3.hs.128=mise((est.an.ti.uh.3.hs.128),mu.hs.128)
mise.an.ti.uh.4.hs.128=mise((est.an.ti.uh.4.hs.128),mu.hs.128)
mise.an.ti.uh.5.hs.128=mise((est.an.ti.uh.5.hs.128),mu.hs.128)
mise.an.ti.uh.6.hs.128=mise((est.an.ti.uh.6.hs.128),mu.hs.128)
mise.an.ti.uh.7.hs.128=mise((est.an.ti.uh.7.hs.128),mu.hs.128)



mise.an.ti.uh.2.bur.1=mise((est.an.ti.uh.2.bur.1),mu.bur.1)
mise.an.ti.uh.3.bur.1=mise((est.an.ti.uh.3.bur.1),mu.bur.1)
mise.an.ti.uh.4.bur.1=mise((est.an.ti.uh.4.bur.1),mu.bur.1)
mise.an.ti.uh.5.bur.1=mise((est.an.ti.uh.5.bur.1),mu.bur.1)
mise.an.ti.uh.6.bur.1=mise((est.an.ti.uh.6.bur.1),mu.bur.1)
mise.an.ti.uh.7.bur.1=mise((est.an.ti.uh.7.bur.1),mu.bur.1)

mise.an.ti.uh.2.bur.8=mise((est.an.ti.uh.2.bur.8),mu.bur.8)
mise.an.ti.uh.3.bur.8=mise((est.an.ti.uh.3.bur.8),mu.bur.8)
mise.an.ti.uh.4.bur.8=mise((est.an.ti.uh.4.bur.8),mu.bur.8)
mise.an.ti.uh.5.bur.8=mise((est.an.ti.uh.5.bur.8),mu.bur.8)
mise.an.ti.uh.6.bur.8=mise((est.an.ti.uh.6.bur.8),mu.bur.8)
mise.an.ti.uh.7.bur.8=mise((est.an.ti.uh.7.bur.8),mu.bur.8)

mise.an.ti.uh.2.bur.128=mise((est.an.ti.uh.2.bur.128),mu.bur.128)
mise.an.ti.uh.3.bur.128=mise((est.an.ti.uh.3.bur.128),mu.bur.128)
mise.an.ti.uh.4.bur.128=mise((est.an.ti.uh.4.bur.128),mu.bur.128)
mise.an.ti.uh.5.bur.128=mise((est.an.ti.uh.5.bur.128),mu.bur.128)
mise.an.ti.uh.6.bur.128=mise((est.an.ti.uh.6.bur.128),mu.bur.128)
mise.an.ti.uh.7.bur.128=mise((est.an.ti.uh.7.bur.128),mu.bur.128)




mise.an.ti.uh.2.cb.1=mise((est.an.ti.uh.2.cb.1),mu.cb.1)
mise.an.ti.uh.3.cb.1=mise((est.an.ti.uh.3.cb.1),mu.cb.1)
mise.an.ti.uh.4.cb.1=mise((est.an.ti.uh.4.cb.1),mu.cb.1)
mise.an.ti.uh.5.cb.1=mise((est.an.ti.uh.5.cb.1),mu.cb.1)
mise.an.ti.uh.6.cb.1=mise((est.an.ti.uh.6.cb.1),mu.cb.1)
mise.an.ti.uh.7.cb.1=mise((est.an.ti.uh.7.cb.1),mu.cb.1)

mise.an.ti.uh.2.cb.8=mise((est.an.ti.uh.2.cb.8),mu.cb.8)
mise.an.ti.uh.3.cb.8=mise((est.an.ti.uh.3.cb.8),mu.cb.8)
mise.an.ti.uh.4.cb.8=mise((est.an.ti.uh.4.cb.8),mu.cb.8)
mise.an.ti.uh.5.cb.8=mise((est.an.ti.uh.5.cb.8),mu.cb.8)
mise.an.ti.uh.6.cb.8=mise((est.an.ti.uh.6.cb.8),mu.cb.8)
mise.an.ti.uh.7.cb.8=mise((est.an.ti.uh.7.cb.8),mu.cb.8)

mise.an.ti.uh.2.cb.128=mise((est.an.ti.uh.2.cb.128),mu.cb.128)
mise.an.ti.uh.3.cb.128=mise((est.an.ti.uh.3.cb.128),mu.cb.128)
mise.an.ti.uh.4.cb.128=mise((est.an.ti.uh.4.cb.128),mu.cb.128)
mise.an.ti.uh.5.cb.128=mise((est.an.ti.uh.5.cb.128),mu.cb.128)
mise.an.ti.uh.6.cb.128=mise((est.an.ti.uh.6.cb.128),mu.cb.128)
mise.an.ti.uh.7.cb.128=mise((est.an.ti.uh.7.cb.128),mu.cb.128)



mise.an.ti.uh.2.b.1=mise((est.an.ti.uh.2.b.1),mu.b.1)
mise.an.ti.uh.3.b.1=mise((est.an.ti.uh.3.b.1),mu.b.1)
mise.an.ti.uh.4.b.1=mise((est.an.ti.uh.4.b.1),mu.b.1)
mise.an.ti.uh.5.b.1=mise((est.an.ti.uh.5.b.1),mu.b.1)
mise.an.ti.uh.6.b.1=mise((est.an.ti.uh.6.b.1),mu.b.1)
mise.an.ti.uh.7.b.1=mise((est.an.ti.uh.7.b.1),mu.b.1)

mise.an.ti.uh.2.b.8=mise((est.an.ti.uh.2.b.8),mu.b.8)
mise.an.ti.uh.3.b.8=mise((est.an.ti.uh.3.b.8),mu.b.8)
mise.an.ti.uh.4.b.8=mise((est.an.ti.uh.4.b.8),mu.b.8)
mise.an.ti.uh.5.b.8=mise((est.an.ti.uh.5.b.8),mu.b.8)
mise.an.ti.uh.6.b.8=mise((est.an.ti.uh.6.b.8),mu.b.8)
mise.an.ti.uh.7.b.8=mise((est.an.ti.uh.7.b.8),mu.b.8)

mise.an.ti.uh.2.b.128=mise((est.an.ti.uh.2.b.128),mu.b.128)
mise.an.ti.uh.3.b.128=mise((est.an.ti.uh.3.b.128),mu.b.128)
mise.an.ti.uh.4.b.128=mise((est.an.ti.uh.4.b.128),mu.b.128)
mise.an.ti.uh.5.b.128=mise((est.an.ti.uh.5.b.128),mu.b.128)
mise.an.ti.uh.6.b.128=mise((est.an.ti.uh.6.b.128),mu.b.128)
mise.an.ti.uh.7.b.128=mise((est.an.ti.uh.7.b.128),mu.b.128)






mise.an.ti.us.2.s.1=mise((est.an.ti.us.2.s.1),mu.s.1)
mise.an.ti.us.3.s.1=mise((est.an.ti.us.3.s.1),mu.s.1)
mise.an.ti.us.4.s.1=mise((est.an.ti.us.4.s.1),mu.s.1)
mise.an.ti.us.5.s.1=mise((est.an.ti.us.5.s.1),mu.s.1)
mise.an.ti.us.6.s.1=mise((est.an.ti.us.6.s.1),mu.s.1)
mise.an.ti.us.7.s.1=mise((est.an.ti.us.7.s.1),mu.s.1)

mise.an.ti.us.2.s.8=mise((est.an.ti.us.2.s.8),mu.s.8)
mise.an.ti.us.3.s.8=mise((est.an.ti.us.3.s.8),mu.s.8)
mise.an.ti.us.4.s.8=mise((est.an.ti.us.4.s.8),mu.s.8)
mise.an.ti.us.5.s.8=mise((est.an.ti.us.5.s.8),mu.s.8)
mise.an.ti.us.6.s.8=mise((est.an.ti.us.6.s.8),mu.s.8)
mise.an.ti.us.7.s.8=mise((est.an.ti.us.7.s.8),mu.s.8)

mise.an.ti.us.2.s.128=mise((est.an.ti.us.2.s.128),mu.s.128)
mise.an.ti.us.3.s.128=mise((est.an.ti.us.3.s.128),mu.s.128)
mise.an.ti.us.4.s.128=mise((est.an.ti.us.4.s.128),mu.s.128)
mise.an.ti.us.5.s.128=mise((est.an.ti.us.5.s.128),mu.s.128)
mise.an.ti.us.6.s.128=mise((est.an.ti.us.6.s.128),mu.s.128)
mise.an.ti.us.7.s.128=mise((est.an.ti.us.7.s.128),mu.s.128)





mise.an.ti.us.2.ang.1=mise((est.an.ti.us.2.ang.1),mu.ang.1)
mise.an.ti.us.3.ang.1=mise((est.an.ti.us.3.ang.1),mu.ang.1)
mise.an.ti.us.4.ang.1=mise((est.an.ti.us.4.ang.1),mu.ang.1)
mise.an.ti.us.5.ang.1=mise((est.an.ti.us.5.ang.1),mu.ang.1)
mise.an.ti.us.6.ang.1=mise((est.an.ti.us.6.ang.1),mu.ang.1)
mise.an.ti.us.7.ang.1=mise((est.an.ti.us.7.ang.1),mu.ang.1)

mise.an.ti.us.2.ang.8=mise((est.an.ti.us.2.ang.8),mu.ang.8)
mise.an.ti.us.3.ang.8=mise((est.an.ti.us.3.ang.8),mu.ang.8)
mise.an.ti.us.4.ang.8=mise((est.an.ti.us.4.ang.8),mu.ang.8)
mise.an.ti.us.5.ang.8=mise((est.an.ti.us.5.ang.8),mu.ang.8)
mise.an.ti.us.6.ang.8=mise((est.an.ti.us.6.ang.8),mu.ang.8)
mise.an.ti.us.7.ang.8=mise((est.an.ti.us.7.ang.8),mu.ang.8)

mise.an.ti.us.2.ang.128=mise((est.an.ti.us.2.ang.128),mu.ang.128)
mise.an.ti.us.3.ang.128=mise((est.an.ti.us.3.ang.128),mu.ang.128)
mise.an.ti.us.4.ang.128=mise((est.an.ti.us.4.ang.128),mu.ang.128)
mise.an.ti.us.5.ang.128=mise((est.an.ti.us.5.ang.128),mu.ang.128)
mise.an.ti.us.6.ang.128=mise((est.an.ti.us.6.ang.128),mu.ang.128)
mise.an.ti.us.7.ang.128=mise((est.an.ti.us.7.ang.128),mu.ang.128)




mise.an.ti.us.2.hs.1=mise((est.an.ti.us.2.hs.1),mu.hs.1)
mise.an.ti.us.3.hs.1=mise((est.an.ti.us.3.hs.1),mu.hs.1)
mise.an.ti.us.4.hs.1=mise((est.an.ti.us.4.hs.1),mu.hs.1)
mise.an.ti.us.5.hs.1=mise((est.an.ti.us.5.hs.1),mu.hs.1)
mise.an.ti.us.6.hs.1=mise((est.an.ti.us.6.hs.1),mu.hs.1)
mise.an.ti.us.7.hs.1=mise((est.an.ti.us.7.hs.1),mu.hs.1)

mise.an.ti.us.2.hs.8=mise((est.an.ti.us.2.hs.8),mu.hs.8)
mise.an.ti.us.3.hs.8=mise((est.an.ti.us.3.hs.8),mu.hs.8)
mise.an.ti.us.4.hs.8=mise((est.an.ti.us.4.hs.8),mu.hs.8)
mise.an.ti.us.5.hs.8=mise((est.an.ti.us.5.hs.8),mu.hs.8)
mise.an.ti.us.6.hs.8=mise((est.an.ti.us.6.hs.8),mu.hs.8)
mise.an.ti.us.7.hs.8=mise((est.an.ti.us.7.hs.8),mu.hs.8)

mise.an.ti.us.2.hs.128=mise((est.an.ti.us.2.hs.128),mu.hs.128)
mise.an.ti.us.3.hs.128=mise((est.an.ti.us.3.hs.128),mu.hs.128)
mise.an.ti.us.4.hs.128=mise((est.an.ti.us.4.hs.128),mu.hs.128)
mise.an.ti.us.5.hs.128=mise((est.an.ti.us.5.hs.128),mu.hs.128)
mise.an.ti.us.6.hs.128=mise((est.an.ti.us.6.hs.128),mu.hs.128)
mise.an.ti.us.7.hs.128=mise((est.an.ti.us.7.hs.128),mu.hs.128)




mise.an.ti.us.2.bur.1=mise((est.an.ti.us.2.bur.1),mu.bur.1)
mise.an.ti.us.3.bur.1=mise((est.an.ti.us.3.bur.1),mu.bur.1)
mise.an.ti.us.4.bur.1=mise((est.an.ti.us.4.bur.1),mu.bur.1)
mise.an.ti.us.5.bur.1=mise((est.an.ti.us.5.bur.1),mu.bur.1)
mise.an.ti.us.6.bur.1=mise((est.an.ti.us.6.bur.1),mu.bur.1)
mise.an.ti.us.7.bur.1=mise((est.an.ti.us.7.bur.1),mu.bur.1)

mise.an.ti.us.2.bur.8=mise((est.an.ti.us.2.bur.8),mu.bur.8)
mise.an.ti.us.3.bur.8=mise((est.an.ti.us.3.bur.8),mu.bur.8)
mise.an.ti.us.4.bur.8=mise((est.an.ti.us.4.bur.8),mu.bur.8)
mise.an.ti.us.5.bur.8=mise((est.an.ti.us.5.bur.8),mu.bur.8)
mise.an.ti.us.6.bur.8=mise((est.an.ti.us.6.bur.8),mu.bur.8)
mise.an.ti.us.7.bur.8=mise((est.an.ti.us.7.bur.8),mu.bur.8)

mise.an.ti.us.2.bur.128=mise((est.an.ti.us.2.bur.128),mu.bur.128)
mise.an.ti.us.3.bur.128=mise((est.an.ti.us.3.bur.128),mu.bur.128)
mise.an.ti.us.4.bur.128=mise((est.an.ti.us.4.bur.128),mu.bur.128)
mise.an.ti.us.5.bur.128=mise((est.an.ti.us.5.bur.128),mu.bur.128)
mise.an.ti.us.6.bur.128=mise((est.an.ti.us.6.bur.128),mu.bur.128)
mise.an.ti.us.7.bur.128=mise((est.an.ti.us.7.bur.128),mu.bur.128)



mise.an.ti.us.2.cb.1=mise((est.an.ti.us.2.cb.1),mu.cb.1)
mise.an.ti.us.3.cb.1=mise((est.an.ti.us.3.cb.1),mu.cb.1)
mise.an.ti.us.4.cb.1=mise((est.an.ti.us.4.cb.1),mu.cb.1)
mise.an.ti.us.5.cb.1=mise((est.an.ti.us.5.cb.1),mu.cb.1)
mise.an.ti.us.6.cb.1=mise((est.an.ti.us.6.cb.1),mu.cb.1)
mise.an.ti.us.7.cb.1=mise((est.an.ti.us.7.cb.1),mu.cb.1)

mise.an.ti.us.2.cb.8=mise((est.an.ti.us.2.cb.8),mu.cb.8)
mise.an.ti.us.3.cb.8=mise((est.an.ti.us.3.cb.8),mu.cb.8)
mise.an.ti.us.4.cb.8=mise((est.an.ti.us.4.cb.8),mu.cb.8)
mise.an.ti.us.5.cb.8=mise((est.an.ti.us.5.cb.8),mu.cb.8)
mise.an.ti.us.6.cb.8=mise((est.an.ti.us.6.cb.8),mu.cb.8)
mise.an.ti.us.7.cb.8=mise((est.an.ti.us.7.cb.8),mu.cb.8)

mise.an.ti.us.2.cb.128=mise((est.an.ti.us.2.cb.128),mu.cb.128)
mise.an.ti.us.3.cb.128=mise((est.an.ti.us.3.cb.128),mu.cb.128)
mise.an.ti.us.4.cb.128=mise((est.an.ti.us.4.cb.128),mu.cb.128)
mise.an.ti.us.5.cb.128=mise((est.an.ti.us.5.cb.128),mu.cb.128)
mise.an.ti.us.6.cb.128=mise((est.an.ti.us.6.cb.128),mu.cb.128)
mise.an.ti.us.7.cb.128=mise((est.an.ti.us.7.cb.128),mu.cb.128)



mise.an.ti.us.2.b.1=mise((est.an.ti.us.2.b.1),mu.b.1)
mise.an.ti.us.3.b.1=mise((est.an.ti.us.3.b.1),mu.b.1)
mise.an.ti.us.4.b.1=mise((est.an.ti.us.4.b.1),mu.b.1)
mise.an.ti.us.5.b.1=mise((est.an.ti.us.5.b.1),mu.b.1)
mise.an.ti.us.6.b.1=mise((est.an.ti.us.6.b.1),mu.b.1)
mise.an.ti.us.7.b.1=mise((est.an.ti.us.7.b.1),mu.b.1)

mise.an.ti.us.2.b.8=mise((est.an.ti.us.2.b.8),mu.b.8)
mise.an.ti.us.3.b.8=mise((est.an.ti.us.3.b.8),mu.b.8)
mise.an.ti.us.4.b.8=mise((est.an.ti.us.4.b.8),mu.b.8)
mise.an.ti.us.5.b.8=mise((est.an.ti.us.5.b.8),mu.b.8)
mise.an.ti.us.6.b.8=mise((est.an.ti.us.6.b.8),mu.b.8)
mise.an.ti.us.7.b.8=mise((est.an.ti.us.7.b.8),mu.b.8)

mise.an.ti.us.2.b.128=mise((est.an.ti.us.2.b.128),mu.b.128)
mise.an.ti.us.3.b.128=mise((est.an.ti.us.3.b.128),mu.b.128)
mise.an.ti.us.4.b.128=mise((est.an.ti.us.4.b.128),mu.b.128)
mise.an.ti.us.5.b.128=mise((est.an.ti.us.5.b.128),mu.b.128)
mise.an.ti.us.6.b.128=mise((est.an.ti.us.6.b.128),mu.b.128)
mise.an.ti.us.7.b.128=mise((est.an.ti.us.7.b.128),mu.b.128)






mise.haar.h.2.s.1=mise((est.haar.h.2.s.1),mu.s.1)
mise.haar.h.3.s.1=mise((est.haar.h.3.s.1),mu.s.1)
mise.haar.h.4.s.1=mise((est.haar.h.4.s.1),mu.s.1)
mise.haar.h.5.s.1=mise((est.haar.h.5.s.1),mu.s.1)
mise.haar.h.6.s.1=mise((est.haar.h.6.s.1),mu.s.1)
mise.haar.h.7.s.1=mise((est.haar.h.7.s.1),mu.s.1)

mise.haar.h.2.s.8=mise((est.haar.h.2.s.8),mu.s.8)
mise.haar.h.3.s.8=mise((est.haar.h.3.s.8),mu.s.8)
mise.haar.h.4.s.8=mise((est.haar.h.4.s.8),mu.s.8)
mise.haar.h.5.s.8=mise((est.haar.h.5.s.8),mu.s.8)
mise.haar.h.6.s.8=mise((est.haar.h.6.s.8),mu.s.8)
mise.haar.h.7.s.8=mise((est.haar.h.7.s.8),mu.s.8)

mise.haar.h.2.s.128=mise((est.haar.h.2.s.128),mu.s.128)
mise.haar.h.3.s.128=mise((est.haar.h.3.s.128),mu.s.128)
mise.haar.h.4.s.128=mise((est.haar.h.4.s.128),mu.s.128)
mise.haar.h.5.s.128=mise((est.haar.h.5.s.128),mu.s.128)
mise.haar.h.6.s.128=mise((est.haar.h.6.s.128),mu.s.128)
mise.haar.h.7.s.128=mise((est.haar.h.7.s.128),mu.s.128)




mise.haar.h.2.ang.1=mise((est.haar.h.2.ang.1),mu.ang.1)
mise.haar.h.3.ang.1=mise((est.haar.h.3.ang.1),mu.ang.1)
mise.haar.h.4.ang.1=mise((est.haar.h.4.ang.1),mu.ang.1)
mise.haar.h.5.ang.1=mise((est.haar.h.5.ang.1),mu.ang.1)
mise.haar.h.6.ang.1=mise((est.haar.h.6.ang.1),mu.ang.1)
mise.haar.h.7.ang.1=mise((est.haar.h.7.ang.1),mu.ang.1)

mise.haar.h.2.ang.8=mise((est.haar.h.2.ang.8),mu.ang.8)
mise.haar.h.3.ang.8=mise((est.haar.h.3.ang.8),mu.ang.8)
mise.haar.h.4.ang.8=mise((est.haar.h.4.ang.8),mu.ang.8)
mise.haar.h.5.ang.8=mise((est.haar.h.5.ang.8),mu.ang.8)
mise.haar.h.6.ang.8=mise((est.haar.h.6.ang.8),mu.ang.8)
mise.haar.h.7.ang.8=mise((est.haar.h.7.ang.8),mu.ang.8)

mise.haar.h.2.ang.128=mise((est.haar.h.2.ang.128),mu.ang.128)
mise.haar.h.3.ang.128=mise((est.haar.h.3.ang.128),mu.ang.128)
mise.haar.h.4.ang.128=mise((est.haar.h.4.ang.128),mu.ang.128)
mise.haar.h.5.ang.128=mise((est.haar.h.5.ang.128),mu.ang.128)
mise.haar.h.6.ang.128=mise((est.haar.h.6.ang.128),mu.ang.128)
mise.haar.h.7.ang.128=mise((est.haar.h.7.ang.128),mu.ang.128)






mise.haar.h.2.hs.1=mise((est.haar.h.2.hs.1),mu.hs.1)
mise.haar.h.3.hs.1=mise((est.haar.h.3.hs.1),mu.hs.1)
mise.haar.h.4.hs.1=mise((est.haar.h.4.hs.1),mu.hs.1)
mise.haar.h.5.hs.1=mise((est.haar.h.5.hs.1),mu.hs.1)
mise.haar.h.6.hs.1=mise((est.haar.h.6.hs.1),mu.hs.1)
mise.haar.h.7.hs.1=mise((est.haar.h.7.hs.1),mu.hs.1)

mise.haar.h.2.hs.8=mise((est.haar.h.2.hs.8),mu.hs.8)
mise.haar.h.3.hs.8=mise((est.haar.h.3.hs.8),mu.hs.8)
mise.haar.h.4.hs.8=mise((est.haar.h.4.hs.8),mu.hs.8)
mise.haar.h.5.hs.8=mise((est.haar.h.5.hs.8),mu.hs.8)
mise.haar.h.6.hs.8=mise((est.haar.h.6.hs.8),mu.hs.8)
mise.haar.h.7.hs.8=mise((est.haar.h.7.hs.8),mu.hs.8)

mise.haar.h.2.hs.128=mise((est.haar.h.2.hs.128),mu.hs.128)
mise.haar.h.3.hs.128=mise((est.haar.h.3.hs.128),mu.hs.128)
mise.haar.h.4.hs.128=mise((est.haar.h.4.hs.128),mu.hs.128)
mise.haar.h.5.hs.128=mise((est.haar.h.5.hs.128),mu.hs.128)
mise.haar.h.6.hs.128=mise((est.haar.h.6.hs.128),mu.hs.128)
mise.haar.h.7.hs.128=mise((est.haar.h.7.hs.128),mu.hs.128)






mise.haar.h.2.bur.1=mise((est.haar.h.2.bur.1),mu.bur.1)
mise.haar.h.3.bur.1=mise((est.haar.h.3.bur.1),mu.bur.1)
mise.haar.h.4.bur.1=mise((est.haar.h.4.bur.1),mu.bur.1)
mise.haar.h.5.bur.1=mise((est.haar.h.5.bur.1),mu.bur.1)
mise.haar.h.6.bur.1=mise((est.haar.h.6.bur.1),mu.bur.1)
mise.haar.h.7.bur.1=mise((est.haar.h.7.bur.1),mu.bur.1)

mise.haar.h.2.bur.8=mise((est.haar.h.2.bur.8),mu.bur.8)
mise.haar.h.3.bur.8=mise((est.haar.h.3.bur.8),mu.bur.8)
mise.haar.h.4.bur.8=mise((est.haar.h.4.bur.8),mu.bur.8)
mise.haar.h.5.bur.8=mise((est.haar.h.5.bur.8),mu.bur.8)
mise.haar.h.6.bur.8=mise((est.haar.h.6.bur.8),mu.bur.8)
mise.haar.h.7.bur.8=mise((est.haar.h.7.bur.8),mu.bur.8)

mise.haar.h.2.bur.128=mise((est.haar.h.2.bur.128),mu.bur.128)
mise.haar.h.3.bur.128=mise((est.haar.h.3.bur.128),mu.bur.128)
mise.haar.h.4.bur.128=mise((est.haar.h.4.bur.128),mu.bur.128)
mise.haar.h.5.bur.128=mise((est.haar.h.5.bur.128),mu.bur.128)
mise.haar.h.6.bur.128=mise((est.haar.h.6.bur.128),mu.bur.128)
mise.haar.h.7.bur.128=mise((est.haar.h.7.bur.128),mu.bur.128)




mise.haar.h.2.cb.1=mise((est.haar.h.2.cb.1),mu.cb.1)
mise.haar.h.3.cb.1=mise((est.haar.h.3.cb.1),mu.cb.1)
mise.haar.h.4.cb.1=mise((est.haar.h.4.cb.1),mu.cb.1)
mise.haar.h.5.cb.1=mise((est.haar.h.5.cb.1),mu.cb.1)
mise.haar.h.6.cb.1=mise((est.haar.h.6.cb.1),mu.cb.1)
mise.haar.h.7.cb.1=mise((est.haar.h.7.cb.1),mu.cb.1)

mise.haar.h.2.cb.8=mise((est.haar.h.2.cb.8),mu.cb.8)
mise.haar.h.3.cb.8=mise((est.haar.h.3.cb.8),mu.cb.8)
mise.haar.h.4.cb.8=mise((est.haar.h.4.cb.8),mu.cb.8)
mise.haar.h.5.cb.8=mise((est.haar.h.5.cb.8),mu.cb.8)
mise.haar.h.6.cb.8=mise((est.haar.h.6.cb.8),mu.cb.8)
mise.haar.h.7.cb.8=mise((est.haar.h.7.cb.8),mu.cb.8)

mise.haar.h.2.cb.128=mise((est.haar.h.2.cb.128),mu.cb.128)
mise.haar.h.3.cb.128=mise((est.haar.h.3.cb.128),mu.cb.128)
mise.haar.h.4.cb.128=mise((est.haar.h.4.cb.128),mu.cb.128)
mise.haar.h.5.cb.128=mise((est.haar.h.5.cb.128),mu.cb.128)
mise.haar.h.6.cb.128=mise((est.haar.h.6.cb.128),mu.cb.128)
mise.haar.h.7.cb.128=mise((est.haar.h.7.cb.128),mu.cb.128)




mise.haar.h.2.b.1=mise((est.haar.h.2.b.1),mu.b.1)
mise.haar.h.3.b.1=mise((est.haar.h.3.b.1),mu.b.1)
mise.haar.h.4.b.1=mise((est.haar.h.4.b.1),mu.b.1)
mise.haar.h.5.b.1=mise((est.haar.h.5.b.1),mu.b.1)
mise.haar.h.6.b.1=mise((est.haar.h.6.b.1),mu.b.1)
mise.haar.h.7.b.1=mise((est.haar.h.7.b.1),mu.b.1)

mise.haar.h.2.b.8=mise((est.haar.h.2.b.8),mu.b.8)
mise.haar.h.3.b.8=mise((est.haar.h.3.b.8),mu.b.8)
mise.haar.h.4.b.8=mise((est.haar.h.4.b.8),mu.b.8)
mise.haar.h.5.b.8=mise((est.haar.h.5.b.8),mu.b.8)
mise.haar.h.6.b.8=mise((est.haar.h.6.b.8),mu.b.8)
mise.haar.h.7.b.8=mise((est.haar.h.7.b.8),mu.b.8)

mise.haar.h.2.b.128=mise((est.haar.h.2.b.128),mu.b.128)
mise.haar.h.3.b.128=mise((est.haar.h.3.b.128),mu.b.128)
mise.haar.h.4.b.128=mise((est.haar.h.4.b.128),mu.b.128)
mise.haar.h.5.b.128=mise((est.haar.h.5.b.128),mu.b.128)
mise.haar.h.6.b.128=mise((est.haar.h.6.b.128),mu.b.128)
mise.haar.h.7.b.128=mise((est.haar.h.7.b.128),mu.b.128)




mise.pen.s.1=mise((est.pen.s.1),mu.s.1)

mise.pen.ang.1=mise((est.pen.ang.1),mu.ang.1)

mise.pen.hs.1=mise((est.pen.hs.1),mu.hs.1)

mise.pen.bur.1=mise((est.pen.bur.1),mu.bur.1)

mise.pen.cb.1=mise((est.pen.cb.1),mu.cb.1)

mise.pen.b.1=mise((est.pen.b.1),mu.b.1)

mise.pen.s.8=mise((est.pen.s.8),mu.s.8)

mise.pen.ang.8=mise((est.pen.ang.8),mu.ang.8)

mise.pen.hs.8=mise((est.pen.hs.8),mu.hs.8)

mise.pen.bur.8=mise((est.pen.bur.8),mu.bur.8)

mise.pen.cb.8=mise((est.pen.cb.8),mu.cb.8)

mise.pen.b.8=mise((est.pen.b.8),mu.b.8)

mise.pen.s.128=mise((est.pen.s.128),mu.s.128)

mise.pen.ang.128=mise((est.pen.ang.128),mu.ang.128)

mise.pen.hs.128=mise((est.pen.hs.128),mu.hs.128)

mise.pen.bur.128=mise((est.pen.bur.128),mu.bur.128)

mise.pen.cb.128=mise((est.pen.cb.128),mu.cb.128)

mise.pen.b.128=mise((est.pen.b.128),mu.b.128)




mise.BMSM.s.1=mise((est.BMSM.s.1),mu.s.1)

mise.BMSM.ang.1=mise((est.BMSM.ang.1),mu.ang.1)

mise.BMSM.hs.1=mise((est.BMSM.hs.1),mu.hs.1)

mise.BMSM.bur.1=mise((est.BMSM.bur.1),mu.bur.1)

mise.BMSM.cb.1=mise((est.BMSM.cb.1),mu.cb.1)

mise.BMSM.b.1=mise((est.BMSM.b.1),mu.b.1)

mise.BMSM.s.8=mise((est.BMSM.s.8),mu.s.8)

mise.BMSM.ang.8=mise((est.BMSM.ang.8),mu.ang.8)

mise.BMSM.hs.8=mise((est.BMSM.hs.8),mu.hs.8)

mise.BMSM.bur.8=mise((est.BMSM.bur.8),mu.bur.8)

mise.BMSM.cb.8=mise((est.BMSM.cb.8),mu.cb.8)

mise.BMSM.b.8=mise((est.BMSM.b.8),mu.b.8)

mise.BMSM.s.128=mise((est.BMSM.s.128),mu.s.128)

mise.BMSM.ang.128=mise((est.BMSM.ang.128),mu.ang.128)

mise.BMSM.hs.128=mise((est.BMSM.hs.128),mu.hs.128)

mise.BMSM.bur.128=mise((est.BMSM.bur.128),mu.bur.128)

mise.BMSM.cb.128=mise((est.BMSM.cb.128),mu.cb.128)

mise.BMSM.b.128=mise((est.BMSM.b.128),mu.b.128)





mise.BMIE.0.s.1=mise((est.BMIE.0.s.1),mu.s.1)
mise.BMIE.2.s.1=mise((est.BMIE.2.s.1),mu.s.1)
mise.BMIE.3.s.1=mise((est.BMIE.3.s.1),mu.s.1)
mise.BMIE.4.s.1=mise((est.BMIE.4.s.1),mu.s.1)
mise.BMIE.5.s.1=mise((est.BMIE.5.s.1),mu.s.1)
mise.BMIE.6.s.1=mise((est.BMIE.6.s.1),mu.s.1)
mise.BMIE.7.s.1=mise((est.BMIE.7.s.1),mu.s.1)

mise.BMIE.0.s.8=mise((est.BMIE.0.s.8),mu.s.8)
mise.BMIE.2.s.8=mise((est.BMIE.2.s.8),mu.s.8)
mise.BMIE.3.s.8=mise((est.BMIE.3.s.8),mu.s.8)
mise.BMIE.4.s.8=mise((est.BMIE.4.s.8),mu.s.8)
mise.BMIE.5.s.8=mise((est.BMIE.5.s.8),mu.s.8)
mise.BMIE.6.s.8=mise((est.BMIE.6.s.8),mu.s.8)
mise.BMIE.7.s.8=mise((est.BMIE.7.s.8),mu.s.8)

mise.BMIE.0.s.128=mise((est.BMIE.0.s.128),mu.s.128)
mise.BMIE.2.s.128=mise((est.BMIE.2.s.128),mu.s.128)
mise.BMIE.3.s.128=mise((est.BMIE.3.s.128),mu.s.128)
mise.BMIE.4.s.128=mise((est.BMIE.4.s.128),mu.s.128)
mise.BMIE.5.s.128=mise((est.BMIE.5.s.128),mu.s.128)
mise.BMIE.6.s.128=mise((est.BMIE.6.s.128),mu.s.128)
mise.BMIE.7.s.128=mise((est.BMIE.7.s.128),mu.s.128)



mise.BMIE.0.ang.1=mise((est.BMIE.0.ang.1),mu.ang.1)
mise.BMIE.2.ang.1=mise((est.BMIE.2.ang.1),mu.ang.1)
mise.BMIE.3.ang.1=mise((est.BMIE.3.ang.1),mu.ang.1)
mise.BMIE.4.ang.1=mise((est.BMIE.4.ang.1),mu.ang.1)
mise.BMIE.5.ang.1=mise((est.BMIE.5.ang.1),mu.ang.1)
mise.BMIE.6.ang.1=mise((est.BMIE.6.ang.1),mu.ang.1)
mise.BMIE.7.ang.1=mise((est.BMIE.7.ang.1),mu.ang.1)

mise.BMIE.0.ang.8=mise((est.BMIE.0.ang.8),mu.ang.8)
mise.BMIE.2.ang.8=mise((est.BMIE.2.ang.8),mu.ang.8)
mise.BMIE.3.ang.8=mise((est.BMIE.3.ang.8),mu.ang.8)
mise.BMIE.4.ang.8=mise((est.BMIE.4.ang.8),mu.ang.8)
mise.BMIE.5.ang.8=mise((est.BMIE.5.ang.8),mu.ang.8)
mise.BMIE.6.ang.8=mise((est.BMIE.6.ang.8),mu.ang.8)
mise.BMIE.7.ang.8=mise((est.BMIE.7.ang.8),mu.ang.8)

mise.BMIE.0.ang.128=mise((est.BMIE.0.ang.128),mu.ang.128)
mise.BMIE.2.ang.128=mise((est.BMIE.2.ang.128),mu.ang.128)
mise.BMIE.3.ang.128=mise((est.BMIE.3.ang.128),mu.ang.128)
mise.BMIE.4.ang.128=mise((est.BMIE.4.ang.128),mu.ang.128)
mise.BMIE.5.ang.128=mise((est.BMIE.5.ang.128),mu.ang.128)
mise.BMIE.6.ang.128=mise((est.BMIE.6.ang.128),mu.ang.128)
mise.BMIE.7.ang.128=mise((est.BMIE.7.ang.128),mu.ang.128)




mise.BMIE.0.hs.1=mise((est.BMIE.0.hs.1),mu.hs.1)
mise.BMIE.2.hs.1=mise((est.BMIE.2.hs.1),mu.hs.1)
mise.BMIE.3.hs.1=mise((est.BMIE.3.hs.1),mu.hs.1)
mise.BMIE.4.hs.1=mise((est.BMIE.4.hs.1),mu.hs.1)
mise.BMIE.5.hs.1=mise((est.BMIE.5.hs.1),mu.hs.1)
mise.BMIE.6.hs.1=mise((est.BMIE.6.hs.1),mu.hs.1)
mise.BMIE.7.hs.1=mise((est.BMIE.7.hs.1),mu.hs.1)

mise.BMIE.0.hs.8=mise((est.BMIE.0.hs.8),mu.hs.8)
mise.BMIE.2.hs.8=mise((est.BMIE.2.hs.8),mu.hs.8)
mise.BMIE.3.hs.8=mise((est.BMIE.3.hs.8),mu.hs.8)
mise.BMIE.4.hs.8=mise((est.BMIE.4.hs.8),mu.hs.8)
mise.BMIE.5.hs.8=mise((est.BMIE.5.hs.8),mu.hs.8)
mise.BMIE.6.hs.8=mise((est.BMIE.6.hs.8),mu.hs.8)
mise.BMIE.7.hs.8=mise((est.BMIE.7.hs.8),mu.hs.8)

mise.BMIE.0.hs.128=mise((est.BMIE.0.hs.128),mu.hs.128)
mise.BMIE.2.hs.128=mise((est.BMIE.2.hs.128),mu.hs.128)
mise.BMIE.3.hs.128=mise((est.BMIE.3.hs.128),mu.hs.128)
mise.BMIE.4.hs.128=mise((est.BMIE.4.hs.128),mu.hs.128)
mise.BMIE.5.hs.128=mise((est.BMIE.5.hs.128),mu.hs.128)
mise.BMIE.6.hs.128=mise((est.BMIE.6.hs.128),mu.hs.128)
mise.BMIE.7.hs.128=mise((est.BMIE.7.hs.128),mu.hs.128)




mise.BMIE.0.bur.1=mise((est.BMIE.0.bur.1),mu.bur.1)
mise.BMIE.2.bur.1=mise((est.BMIE.2.bur.1),mu.bur.1)
mise.BMIE.3.bur.1=mise((est.BMIE.3.bur.1),mu.bur.1)
mise.BMIE.4.bur.1=mise((est.BMIE.4.bur.1),mu.bur.1)
mise.BMIE.5.bur.1=mise((est.BMIE.5.bur.1),mu.bur.1)
mise.BMIE.6.bur.1=mise((est.BMIE.6.bur.1),mu.bur.1)
mise.BMIE.7.bur.1=mise((est.BMIE.7.bur.1),mu.bur.1)

mise.BMIE.0.bur.8=mise((est.BMIE.0.bur.8),mu.bur.8)
mise.BMIE.2.bur.8=mise((est.BMIE.2.bur.8),mu.bur.8)
mise.BMIE.3.bur.8=mise((est.BMIE.3.bur.8),mu.bur.8)
mise.BMIE.4.bur.8=mise((est.BMIE.4.bur.8),mu.bur.8)
mise.BMIE.5.bur.8=mise((est.BMIE.5.bur.8),mu.bur.8)
mise.BMIE.6.bur.8=mise((est.BMIE.6.bur.8),mu.bur.8)
mise.BMIE.7.bur.8=mise((est.BMIE.7.bur.8),mu.bur.8)

mise.BMIE.0.bur.128=mise((est.BMIE.0.bur.128),mu.bur.128)
mise.BMIE.2.bur.128=mise((est.BMIE.2.bur.128),mu.bur.128)
mise.BMIE.3.bur.128=mise((est.BMIE.3.bur.128),mu.bur.128)
mise.BMIE.4.bur.128=mise((est.BMIE.4.bur.128),mu.bur.128)
mise.BMIE.5.bur.128=mise((est.BMIE.5.bur.128),mu.bur.128)
mise.BMIE.6.bur.128=mise((est.BMIE.6.bur.128),mu.bur.128)
mise.BMIE.7.bur.128=mise((est.BMIE.7.bur.128),mu.bur.128)


mise.BMIE.0.cb.1=mise((est.BMIE.0.cb.1),mu.cb.1)
mise.BMIE.2.cb.1=mise((est.BMIE.2.cb.1),mu.cb.1)
mise.BMIE.3.cb.1=mise((est.BMIE.3.cb.1),mu.cb.1)
mise.BMIE.4.cb.1=mise((est.BMIE.4.cb.1),mu.cb.1)
mise.BMIE.5.cb.1=mise((est.BMIE.5.cb.1),mu.cb.1)
mise.BMIE.6.cb.1=mise((est.BMIE.6.cb.1),mu.cb.1)
mise.BMIE.7.cb.1=mise((est.BMIE.7.cb.1),mu.cb.1)

mise.BMIE.0.cb.8=mise((est.BMIE.0.cb.8),mu.cb.8)
mise.BMIE.2.cb.8=mise((est.BMIE.2.cb.8),mu.cb.8)
mise.BMIE.3.cb.8=mise((est.BMIE.3.cb.8),mu.cb.8)
mise.BMIE.4.cb.8=mise((est.BMIE.4.cb.8),mu.cb.8)
mise.BMIE.5.cb.8=mise((est.BMIE.5.cb.8),mu.cb.8)
mise.BMIE.6.cb.8=mise((est.BMIE.6.cb.8),mu.cb.8)
mise.BMIE.7.cb.8=mise((est.BMIE.7.cb.8),mu.cb.8)

mise.BMIE.0.cb.128=mise((est.BMIE.0.cb.128),mu.cb.128)
mise.BMIE.2.cb.128=mise((est.BMIE.2.cb.128),mu.cb.128)
mise.BMIE.3.cb.128=mise((est.BMIE.3.cb.128),mu.cb.128)
mise.BMIE.4.cb.128=mise((est.BMIE.4.cb.128),mu.cb.128)
mise.BMIE.5.cb.128=mise((est.BMIE.5.cb.128),mu.cb.128)
mise.BMIE.6.cb.128=mise((est.BMIE.6.cb.128),mu.cb.128)
mise.BMIE.7.cb.128=mise((est.BMIE.7.cb.128),mu.cb.128)




mise.BMIE.0.b.1=mise((est.BMIE.0.b.1),mu.b.1)
mise.BMIE.2.b.1=mise((est.BMIE.2.b.1),mu.b.1)
mise.BMIE.3.b.1=mise((est.BMIE.3.b.1),mu.b.1)
mise.BMIE.4.b.1=mise((est.BMIE.4.b.1),mu.b.1)
mise.BMIE.5.b.1=mise((est.BMIE.5.b.1),mu.b.1)
mise.BMIE.6.b.1=mise((est.BMIE.6.b.1),mu.b.1)
mise.BMIE.7.b.1=mise((est.BMIE.7.b.1),mu.b.1)

mise.BMIE.0.b.8=mise((est.BMIE.0.b.8),mu.b.8)
mise.BMIE.2.b.8=mise((est.BMIE.2.b.8),mu.b.8)
mise.BMIE.3.b.8=mise((est.BMIE.3.b.8),mu.b.8)
mise.BMIE.4.b.8=mise((est.BMIE.4.b.8),mu.b.8)
mise.BMIE.5.b.8=mise((est.BMIE.5.b.8),mu.b.8)
mise.BMIE.6.b.8=mise((est.BMIE.6.b.8),mu.b.8)
mise.BMIE.7.b.8=mise((est.BMIE.7.b.8),mu.b.8)

mise.BMIE.0.b.128=mise((est.BMIE.0.b.128),mu.b.128)
mise.BMIE.2.b.128=mise((est.BMIE.2.b.128),mu.b.128)
mise.BMIE.3.b.128=mise((est.BMIE.3.b.128),mu.b.128)
mise.BMIE.4.b.128=mise((est.BMIE.4.b.128),mu.b.128)
mise.BMIE.5.b.128=mise((est.BMIE.5.b.128),mu.b.128)
mise.BMIE.6.b.128=mise((est.BMIE.6.b.128),mu.b.128)
mise.BMIE.7.b.128=mise((est.BMIE.7.b.128),mu.b.128)



#1 - ash
#2 - BMSM
#3 - BMMIM
#4 - haarfisz_R
#5 - anscombe_cv_h
#6 - anscombe_u_h
#7 - haar_h
#8 - l1_pen



#collect mises for each test function and each intensity level
mise.s.1=c(mean(mise.ash.s.1),
           mean(mise.BMSM.s.1),
           mean(c(mise.BMIE.4.s.1,mise.BMIE.5.s.1,mise.BMIE.6.s.1,mise.BMIE.7.s.1)),
           mean(c(mise.hf.ti.r.4.s.1,mise.hf.ti.r.5.s.1,mise.hf.ti.r.6.s.1,mise.hf.ti.r.7.s.1)),
           mean(c(mise.an.ti.cvh.4.s.1,mise.an.ti.cvh.5.s.1,mise.an.ti.cvh.6.s.1,mise.an.ti.cvh.7.s.1)),
           mean(c(mise.an.ti.uh.4.s.1,mise.an.ti.uh.5.s.1,mise.an.ti.uh.6.s.1,mise.an.ti.uh.7.s.1)),
           mean(c(mise.haar.h.4.s.1,mise.haar.h.5.s.1,mise.haar.h.6.s.1,mise.haar.h.7.s.1)),
           mean(mise.pen.s.1))

names(mise.s.1)=c("ash",
                  "BMSM",
                  "BMMIM",
                  "haarfisz_R",
                  "anscombe_cv_h",
                  "anscombe_u_h",
                  "haar_h",
                  "l1_pen")


mise.s.8=c(mean(mise.ash.s.8),
           mean(mise.BMSM.s.8),
           mean(c(mise.BMIE.4.s.8,mise.BMIE.5.s.8,mise.BMIE.6.s.8,mise.BMIE.7.s.8)),
           mean(c(mise.hf.ti.r.4.s.8,mise.hf.ti.r.5.s.8,mise.hf.ti.r.6.s.8,mise.hf.ti.r.7.s.8)),
           mean(c(mise.an.ti.cvh.4.s.8,mise.an.ti.cvh.5.s.8,mise.an.ti.cvh.6.s.8,mise.an.ti.cvh.7.s.8)),
           mean(c(mise.an.ti.uh.4.s.8,mise.an.ti.uh.5.s.8,mise.an.ti.uh.6.s.8,mise.an.ti.uh.7.s.8)),
           mean(c(mise.haar.h.4.s.8,mise.haar.h.5.s.8,mise.haar.h.6.s.8,mise.haar.h.7.s.8)),
           mean(mise.pen.s.8))

names(mise.s.8)=c("ash",
                  "BMSM",
                  "BMMIM",
                  "haarfisz_R",
                  "anscombe_cv_h",
                  "anscombe_u_h",
                  "haar_h",
                  "l1_pen")


mise.s.128=c(mean(mise.ash.s.128),
             mean(mise.BMSM.s.128),
             mean(c(mise.BMIE.4.s.128,mise.BMIE.5.s.128,mise.BMIE.6.s.128,mise.BMIE.7.s.128)),
             mean(c(mise.hf.ti.r.4.s.128,mise.hf.ti.r.5.s.128,mise.hf.ti.r.6.s.128,mise.hf.ti.r.7.s.128)),
             mean(c(mise.an.ti.cvh.4.s.128,mise.an.ti.cvh.5.s.128,mise.an.ti.cvh.6.s.128,mise.an.ti.cvh.7.s.128)),
             mean(c(mise.an.ti.uh.4.s.128,mise.an.ti.uh.5.s.128,mise.an.ti.uh.6.s.128,mise.an.ti.uh.7.s.128)),
             mean(c(mise.haar.h.4.s.128,mise.haar.h.5.s.128,mise.haar.h.6.s.128,mise.haar.h.7.s.128)),
             mean(mise.pen.s.128))

names(mise.s.128)=c("ash",
                    "BMSM",
                    "BMMIM",
                    "haarfisz_R",
                    "anscombe_cv_h",
                    "anscombe_u_h",
                    "haar_h",
                    "l1_pen")




mise.ang.1=c(mean(mise.ash.ang.1),
             mean(mise.BMSM.ang.1),
             mean(c(mise.BMIE.4.ang.1,mise.BMIE.5.ang.1,mise.BMIE.6.ang.1,mise.BMIE.7.ang.1)),
             mean(c(mise.hf.ti.r.4.ang.1,mise.hf.ti.r.5.ang.1,mise.hf.ti.r.6.ang.1,mise.hf.ti.r.7.ang.1)),
             mean(c(mise.an.ti.cvh.4.ang.1,mise.an.ti.cvh.5.ang.1,mise.an.ti.cvh.6.ang.1,mise.an.ti.cvh.7.ang.1)),
             mean(c(mise.an.ti.uh.4.ang.1,mise.an.ti.uh.5.ang.1,mise.an.ti.uh.6.ang.1,mise.an.ti.uh.7.ang.1)),
             mean(c(mise.haar.h.4.ang.1,mise.haar.h.5.ang.1,mise.haar.h.6.ang.1,mise.haar.h.7.ang.1)),
             mean(mise.pen.ang.1))

names(mise.ang.1)=c("ash",
                    "BMSM",
                    "BMMIM",
                    "haarfisz_R",
                    "anscombe_cv_h",
                    "anscombe_u_h",
                    "haar_h",
                    "l1_pen")


mise.ang.8=c(mean(mise.ash.ang.8),
             mean(mise.BMSM.ang.8),
             mean(c(mise.BMIE.4.ang.8,mise.BMIE.5.ang.8,mise.BMIE.6.ang.8,mise.BMIE.7.ang.8)),
             mean(c(mise.hf.ti.r.4.ang.8,mise.hf.ti.r.5.ang.8,mise.hf.ti.r.6.ang.8,mise.hf.ti.r.7.ang.8)),
             mean(c(mise.an.ti.cvh.4.ang.8,mise.an.ti.cvh.5.ang.8,mise.an.ti.cvh.6.ang.8,mise.an.ti.cvh.7.ang.8)),
             mean(c(mise.an.ti.uh.4.ang.8,mise.an.ti.uh.5.ang.8,mise.an.ti.uh.6.ang.8,mise.an.ti.uh.7.ang.8)),
             mean(c(mise.haar.h.4.ang.8,mise.haar.h.5.ang.8,mise.haar.h.6.ang.8,mise.haar.h.7.ang.8)),
             mean(mise.pen.ang.8))

names(mise.ang.8)=c("ash",
                    "BMSM",
                    "BMMIM",
                    "haarfisz_R",
                    "anscombe_cv_h",
                    "anscombe_u_h",
                    "haar_h",
                    "l1_pen")


mise.ang.128=c(mean(mise.ash.ang.128),
               mean(mise.BMSM.ang.128),
               mean(c(mise.BMIE.4.ang.128,mise.BMIE.5.ang.128,mise.BMIE.6.ang.128,mise.BMIE.7.ang.128)),
               mean(c(mise.hf.ti.r.4.ang.128,mise.hf.ti.r.5.ang.128,mise.hf.ti.r.6.ang.128,mise.hf.ti.r.7.ang.128)),
               mean(c(mise.an.ti.cvh.4.ang.128,mise.an.ti.cvh.5.ang.128,mise.an.ti.cvh.6.ang.128,mise.an.ti.cvh.7.ang.128)),
               mean(c(mise.an.ti.uh.4.ang.128,mise.an.ti.uh.5.ang.128,mise.an.ti.uh.6.ang.128,mise.an.ti.uh.7.ang.128)),
               mean(c(mise.haar.h.4.ang.128,mise.haar.h.5.ang.128,mise.haar.h.6.ang.128,mise.haar.h.7.ang.128)),
               mean(mise.pen.ang.128))

names(mise.ang.128)=c("ash",
                      "BMSM",
                      "BMMIM",
                      "haarfisz_R",
                      "anscombe_cv_h",
                      "anscombe_u_h",
                      "haar_h",
                      "l1_pen")




mise.hs.1=c(mean(mise.ash.hs.1),
            mean(mise.BMSM.hs.1),
            mean(c(mise.BMIE.4.hs.1,mise.BMIE.5.hs.1,mise.BMIE.6.hs.1,mise.BMIE.7.hs.1)),
            mean(c(mise.hf.ti.r.4.hs.1,mise.hf.ti.r.5.hs.1,mise.hf.ti.r.6.hs.1,mise.hf.ti.r.7.hs.1)),
            mean(c(mise.an.ti.cvh.4.hs.1,mise.an.ti.cvh.5.hs.1,mise.an.ti.cvh.6.hs.1,mise.an.ti.cvh.7.hs.1)),
            mean(c(mise.an.ti.uh.4.hs.1,mise.an.ti.uh.5.hs.1,mise.an.ti.uh.6.hs.1,mise.an.ti.uh.7.hs.1)),
            mean(c(mise.haar.h.4.hs.1,mise.haar.h.5.hs.1,mise.haar.h.6.hs.1,mise.haar.h.7.hs.1)),
            mean(mise.pen.hs.1))

names(mise.hs.1)=c("ash",
                   "BMSM",
                   "BMMIM",
                   "haarfisz_R",
                   "anscombe_cv_h",
                   "anscombe_u_h",
                   "haar_h",
                   "l1_pen")


mise.hs.8=c(mean(mise.ash.hs.8),
            mean(mise.BMSM.hs.8),
            mean(c(mise.BMIE.4.hs.8,mise.BMIE.5.hs.8,mise.BMIE.6.hs.8,mise.BMIE.7.hs.8)),
            mean(c(mise.hf.ti.r.4.hs.8,mise.hf.ti.r.5.hs.8,mise.hf.ti.r.6.hs.8,mise.hf.ti.r.7.hs.8)),
            mean(c(mise.an.ti.cvh.4.hs.8,mise.an.ti.cvh.5.hs.8,mise.an.ti.cvh.6.hs.8,mise.an.ti.cvh.7.hs.8)),
            mean(c(mise.an.ti.uh.4.hs.8,mise.an.ti.uh.5.hs.8,mise.an.ti.uh.6.hs.8,mise.an.ti.uh.7.hs.8)),
            mean(c(mise.haar.h.4.hs.8,mise.haar.h.5.hs.8,mise.haar.h.6.hs.8,mise.haar.h.7.hs.8)),
            mean(mise.pen.hs.8))

names(mise.hs.8)=c("ash",
                   "BMSM",
                   "BMMIM",
                   "haarfisz_R",
                   "anscombe_cv_h",
                   "anscombe_u_h",
                   "haar_h",
                   "l1_pen")


mise.hs.128=c(mean(mise.ash.hs.128),
              mean(mise.BMSM.hs.128),
              mean(c(mise.BMIE.4.hs.128,mise.BMIE.5.hs.128,mise.BMIE.6.hs.128,mise.BMIE.7.hs.128)),
              mean(c(mise.hf.ti.r.4.hs.128,mise.hf.ti.r.5.hs.128,mise.hf.ti.r.6.hs.128,mise.hf.ti.r.7.hs.128)),
              mean(c(mise.an.ti.cvh.4.hs.128,mise.an.ti.cvh.5.hs.128,mise.an.ti.cvh.6.hs.128,mise.an.ti.cvh.7.hs.128)),
              mean(c(mise.an.ti.uh.4.hs.128,mise.an.ti.uh.5.hs.128,mise.an.ti.uh.6.hs.128,mise.an.ti.uh.7.hs.128)),
              mean(c(mise.haar.h.4.hs.128,mise.haar.h.5.hs.128,mise.haar.h.6.hs.128,mise.haar.h.7.hs.128)),
              mean(mise.pen.hs.128))

names(mise.hs.128)=c("ash",
                     "BMSM",
                     "BMMIM",
                     "haarfisz_R",
                     "anscombe_cv_h",
                     "anscombe_u_h",
                     "haar_h",
                     "l1_pen")





mise.bur.1=c(mean(mise.ash.bur.1),
             mean(mise.BMSM.bur.1),
             mean(c(mise.BMIE.4.bur.1,mise.BMIE.5.bur.1,mise.BMIE.6.bur.1,mise.BMIE.7.bur.1)),
             mean(c(mise.hf.ti.r.4.bur.1,mise.hf.ti.r.5.bur.1,mise.hf.ti.r.6.bur.1,mise.hf.ti.r.7.bur.1)),
             mean(c(mise.an.ti.cvh.4.bur.1,mise.an.ti.cvh.5.bur.1,mise.an.ti.cvh.6.bur.1,mise.an.ti.cvh.7.bur.1)),
             mean(c(mise.an.ti.uh.4.bur.1,mise.an.ti.uh.5.bur.1,mise.an.ti.uh.6.bur.1,mise.an.ti.uh.7.bur.1)),
             mean(c(mise.haar.h.4.bur.1,mise.haar.h.5.bur.1,mise.haar.h.6.bur.1,mise.haar.h.7.bur.1)),
             mean(mise.pen.bur.1))

names(mise.bur.1)=c("ash",
                    "BMSM",
                    "BMMIM",
                    "haarfisz_R",
                    "anscombe_cv_h",
                    "anscombe_u_h",
                    "haar_h",
                    "l1_pen")


mise.bur.8=c(mean(mise.ash.bur.8),
             mean(mise.BMSM.bur.8),
             mean(c(mise.BMIE.4.bur.8,mise.BMIE.5.bur.8,mise.BMIE.6.bur.8,mise.BMIE.7.bur.8)),
             mean(c(mise.hf.ti.r.4.bur.8,mise.hf.ti.r.5.bur.8,mise.hf.ti.r.6.bur.8,mise.hf.ti.r.7.bur.8)),
             mean(c(mise.an.ti.cvh.4.bur.8,mise.an.ti.cvh.5.bur.8,mise.an.ti.cvh.6.bur.8,mise.an.ti.cvh.7.bur.8)),
             mean(c(mise.an.ti.uh.4.bur.8,mise.an.ti.uh.5.bur.8,mise.an.ti.uh.6.bur.8,mise.an.ti.uh.7.bur.8)),
             mean(c(mise.haar.h.4.bur.8,mise.haar.h.5.bur.8,mise.haar.h.6.bur.8,mise.haar.h.7.bur.8)),
             mean(mise.pen.bur.8))

names(mise.bur.8)=c("ash",
                    "BMSM",
                    "BMMIM",
                    "haarfisz_R",
                    "anscombe_cv_h",
                    "anscombe_u_h",
                    "haar_h",
                    "l1_pen")


mise.bur.128=c(mean(mise.ash.bur.128),
               mean(mise.BMSM.bur.128),
               mean(c(mise.BMIE.4.bur.128,mise.BMIE.5.bur.128,mise.BMIE.6.bur.128,mise.BMIE.7.bur.128)),
               mean(c(mise.hf.ti.r.4.bur.128,mise.hf.ti.r.5.bur.128,mise.hf.ti.r.6.bur.128,mise.hf.ti.r.7.bur.128)),
               mean(c(mise.an.ti.cvh.4.bur.128,mise.an.ti.cvh.5.bur.128,mise.an.ti.cvh.6.bur.128,mise.an.ti.cvh.7.bur.128)),
               mean(c(mise.an.ti.uh.4.bur.128,mise.an.ti.uh.5.bur.128,mise.an.ti.uh.6.bur.128,mise.an.ti.uh.7.bur.128)),
               mean(c(mise.haar.h.4.bur.128,mise.haar.h.5.bur.128,mise.haar.h.6.bur.128,mise.haar.h.7.bur.128)),
               mean(mise.pen.bur.128))

names(mise.bur.128)=c("ash",
                      "BMSM",
                      "BMMIM",
                      "haarfisz_R",
                      "anscombe_cv_h",
                      "anscombe_u_h",
                      "haar_h",
                      "l1_pen")





mise.cb.1=c(mean(mise.ash.cb.1),
            mean(mise.BMSM.cb.1),
            mean(c(mise.BMIE.4.cb.1,mise.BMIE.5.cb.1,mise.BMIE.6.cb.1,mise.BMIE.7.cb.1)),
            mean(c(mise.hf.ti.r.4.cb.1,mise.hf.ti.r.5.cb.1,mise.hf.ti.r.6.cb.1,mise.hf.ti.r.7.cb.1)),
            mean(c(mise.an.ti.cvh.4.cb.1,mise.an.ti.cvh.5.cb.1,mise.an.ti.cvh.6.cb.1,mise.an.ti.cvh.7.cb.1)),
            mean(c(mise.an.ti.uh.4.cb.1,mise.an.ti.uh.5.cb.1,mise.an.ti.uh.6.cb.1,mise.an.ti.uh.7.cb.1)),
            mean(c(mise.haar.h.4.cb.1,mise.haar.h.5.cb.1,mise.haar.h.6.cb.1,mise.haar.h.7.cb.1)),
            mean(mise.pen.cb.1))

names(mise.cb.1)=c("ash",
                   "BMSM",
                   "BMMIM",
                   "haarfisz_R",
                   "anscombe_cv_h",
                   "anscombe_u_h",
                   "haar_h",
                   "l1_pen")


mise.cb.8=c(mean(mise.ash.cb.8),
            mean(mise.BMSM.cb.8),
            mean(c(mise.BMIE.4.cb.8,mise.BMIE.5.cb.8,mise.BMIE.6.cb.8,mise.BMIE.7.cb.8)),
            mean(c(mise.hf.ti.r.4.cb.8,mise.hf.ti.r.5.cb.8,mise.hf.ti.r.6.cb.8,mise.hf.ti.r.7.cb.8)),
            mean(c(mise.an.ti.cvh.4.cb.8,mise.an.ti.cvh.5.cb.8,mise.an.ti.cvh.6.cb.8,mise.an.ti.cvh.7.cb.8)),
            mean(c(mise.an.ti.uh.4.cb.8,mise.an.ti.uh.5.cb.8,mise.an.ti.uh.6.cb.8,mise.an.ti.uh.7.cb.8)),
            mean(c(mise.haar.h.4.cb.8,mise.haar.h.5.cb.8,mise.haar.h.6.cb.8,mise.haar.h.7.cb.8)),
            mean(mise.pen.cb.8))

names(mise.cb.8)=c("ash",
                   "BMSM",
                   "BMMIM",
                   "haarfisz_R",
                   "anscombe_cv_h",
                   "anscombe_u_h",
                   "haar_h",
                   "l1_pen")


mise.cb.128=c(mean(mise.ash.cb.128),
              mean(mise.BMSM.cb.128),
              mean(c(mise.BMIE.4.cb.128,mise.BMIE.5.cb.128,mise.BMIE.6.cb.128,mise.BMIE.7.cb.128)),
              mean(c(mise.hf.ti.r.4.cb.128,mise.hf.ti.r.5.cb.128,mise.hf.ti.r.6.cb.128,mise.hf.ti.r.7.cb.128)),
              mean(c(mise.an.ti.cvh.4.cb.128,mise.an.ti.cvh.5.cb.128,mise.an.ti.cvh.6.cb.128,mise.an.ti.cvh.7.cb.128)),
              mean(c(mise.an.ti.uh.4.cb.128,mise.an.ti.uh.5.cb.128,mise.an.ti.uh.6.cb.128,mise.an.ti.uh.7.cb.128)),
              mean(c(mise.haar.h.4.cb.128,mise.haar.h.5.cb.128,mise.haar.h.6.cb.128,mise.haar.h.7.cb.128)),
              mean(mise.pen.cb.128))

names(mise.cb.128)=c("ash",
                     "BMSM",
                     "BMMIM",
                     "haarfisz_R",
                     "anscombe_cv_h",
                     "anscombe_u_h",
                     "haar_h",
                     "l1_pen")




mise.b.1=c(mean(mise.ash.b.1),
           mean(mise.BMSM.b.1),
           mean(c(mise.BMIE.4.b.1,mise.BMIE.5.b.1,mise.BMIE.6.b.1,mise.BMIE.7.b.1)),
           mean(c(mise.hf.ti.r.4.b.1,mise.hf.ti.r.5.b.1,mise.hf.ti.r.6.b.1,mise.hf.ti.r.7.b.1)),
           mean(c(mise.an.ti.cvh.4.b.1,mise.an.ti.cvh.5.b.1,mise.an.ti.cvh.6.b.1,mise.an.ti.cvh.7.b.1)),
           mean(c(mise.an.ti.uh.4.b.1,mise.an.ti.uh.5.b.1,mise.an.ti.uh.6.b.1,mise.an.ti.uh.7.b.1)),
           mean(c(mise.haar.h.4.b.1,mise.haar.h.5.b.1,mise.haar.h.6.b.1,mise.haar.h.7.b.1)),
           mean(mise.pen.b.1))

names(mise.b.1)=c("ash",
                  "BMSM",
                  "BMMIM",
                  "haarfisz_R",
                  "anscombe_cv_h",
                  "anscombe_u_h",
                  "haar_h",
                  "l1_pen")


mise.b.8=c(mean(mise.ash.b.8),
           mean(mise.BMSM.b.8),
           mean(c(mise.BMIE.4.b.8,mise.BMIE.5.b.8,mise.BMIE.6.b.8,mise.BMIE.7.b.8)),
           mean(c(mise.hf.ti.r.4.b.8,mise.hf.ti.r.5.b.8,mise.hf.ti.r.6.b.8,mise.hf.ti.r.7.b.8)),
           mean(c(mise.an.ti.cvh.4.b.8,mise.an.ti.cvh.5.b.8,mise.an.ti.cvh.6.b.8,mise.an.ti.cvh.7.b.8)),
           mean(c(mise.an.ti.uh.4.b.8,mise.an.ti.uh.5.b.8,mise.an.ti.uh.6.b.8,mise.an.ti.uh.7.b.8)),
           mean(c(mise.haar.h.4.b.8,mise.haar.h.5.b.8,mise.haar.h.6.b.8,mise.haar.h.7.b.8)),
           mean(mise.pen.b.8))

names(mise.b.8)=c("ash",
                  "BMSM",
                  "BMMIM",
                  "haarfisz_R",
                  "anscombe_cv_h",
                  "anscombe_u_h",
                  "haar_h",
                  "l1_pen")


mise.b.128=c(mean(mise.ash.b.128),
             mean(mise.BMSM.b.128),
             mean(c(mise.BMIE.4.b.128,mise.BMIE.5.b.128,mise.BMIE.6.b.128,mise.BMIE.7.b.128)),
             mean(c(mise.hf.ti.r.4.b.128,mise.hf.ti.r.5.b.128,mise.hf.ti.r.6.b.128,mise.hf.ti.r.7.b.128)),
             mean(c(mise.an.ti.cvh.4.b.128,mise.an.ti.cvh.5.b.128,mise.an.ti.cvh.6.b.128,mise.an.ti.cvh.7.b.128)),
             mean(c(mise.an.ti.uh.4.b.128,mise.an.ti.uh.5.b.128,mise.an.ti.uh.6.b.128,mise.an.ti.uh.7.b.128)),
             mean(c(mise.haar.h.4.b.128,mise.haar.h.5.b.128,mise.haar.h.6.b.128,mise.haar.h.7.b.128)),
             mean(mise.pen.b.128))

names(mise.b.128)=c("ash",
                    "BMSM",
                    "BMMIM",
                    "haarfisz_R",
                    "anscombe_cv_h",
                    "anscombe_u_h",
                    "haar_h",
                    "l1_pen")

#save data
save.image("res_pois_full.RData")

rm(list = setdiff(ls(), ls(pattern="^mise")))

save.image("res_pois.RData")

#look at results
sort(mise.s.1)
sort(mise.s.8)
sort(mise.s.128)

sort(mise.ang.1)
sort(mise.ang.8)
sort(mise.ang.128)

sort(mise.hs.1)
sort(mise.hs.8)
sort(mise.hs.128)

sort(mise.bur.1)
sort(mise.bur.8)
sort(mise.bur.128)

sort(mise.cb.1)
sort(mise.cb.8)
sort(mise.cb.128)

sort(mise.b.1)
sort(mise.b.8)
sort(mise.b.128)



