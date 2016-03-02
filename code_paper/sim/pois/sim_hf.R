library(haarfisz)

#define functions
mse=function(x,y) mean((x-y)^2)
l2norm=function(x) sum(x^2)
mise=function(x,y) 10000*mean(apply(x-matrix(rep(y,times=100),nrow=100,byrow=T),1,l2norm)/l2norm(y))


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



hf.la10.est.ti2=function (x, filter.number = 10, family = "DaubLeAsymm", min.level = 2) 
{
  x.w <- wd(x, filter.number, family, type = "station")
  x.w.t <- threshold(x.w, levels = (min.level):(x.w$nlevels - 1), policy = "universal", type = "hard")
  x.w.t.r <- AvBasis(convert(x.w.t))
  return(x.w.t.r)
}


hf.la10.est.ti3=function (x, filter.number = 10, family = "DaubLeAsymm", min.level = 3) 
{
  x.w <- wd(x, filter.number, family, type = "station")
  x.w.t <- threshold(x.w, levels = (min.level):(x.w$nlevels - 1), policy = "universal", type = "hard")
  x.w.t.r <- AvBasis(convert(x.w.t))
  return(x.w.t.r)
}

hf.la10.est.ti4=function (x, filter.number = 10, family = "DaubLeAsymm", min.level = 4) 
{
  x.w <- wd(x, filter.number, family, type = "station")
  x.w.t <- threshold(x.w, levels = (min.level):(x.w$nlevels - 1), policy = "universal", type = "hard")
  x.w.t.r <- AvBasis(convert(x.w.t))
  return(x.w.t.r)
}

hf.la10.est.ti5=function (x, filter.number = 10, family = "DaubLeAsymm", min.level = 5) 
{
  x.w <- wd(x, filter.number, family, type = "station")
  x.w.t <- threshold(x.w, levels = (min.level):(x.w$nlevels - 1), policy = "universal", type = "hard")
  x.w.t.r <- AvBasis(convert(x.w.t))
  return(x.w.t.r)
}

hf.la10.est.ti6=function (x, filter.number = 10, family = "DaubLeAsymm", min.level = 6) 
{
  x.w <- wd(x, filter.number, family, type = "station")
  x.w.t <- threshold(x.w, levels = (min.level):(x.w$nlevels - 1), policy = "universal", type = "hard")
  x.w.t.r <- AvBasis(convert(x.w.t))
  return(x.w.t.r)
}

hf.la10.est.ti7=function (x, filter.number = 10, family = "DaubLeAsymm", min.level = 7) 
{
  x.w <- wd(x, filter.number, family, type = "station")
  x.w.t <- threshold(x.w, levels = (min.level):(x.w$nlevels - 1), policy = "universal", type = "hard")
  x.w.t.r <- AvBasis(convert(x.w.t))
  return(x.w.t.r)
}


#run haarfisz with different options
#some test functions converge with default options

#default options
est.hf.default.r.4.s.8=apply(sim.m.s.8,1,denoise.poisson) 
# est.hf.default.r.5.s.8=apply(sim.m.s.8,1,denoise.poisson)
# est.hf.default.r.6.s.8=apply(sim.m.s.8,1,denoise.poisson)
# est.hf.default.r.7.s.8=apply(sim.m.s.8,1,denoise.poisson)

est.hf.default.r.4.ang.8=apply(sim.m.ang.8,1,denoise.poisson)
# est.hf.default.r.5.ang.8=apply(sim.m.ang.8,1,denoise.poisson)
# est.hf.default.r.6.ang.8=apply(sim.m.ang.8,1,denoise.poisson)
# est.hf.default.r.7.ang.8=apply(sim.m.ang.8,1,denoise.poisson)

est.hf.default.r.4.cb.8=apply(sim.m.cb.8,1,denoise.poisson)
# est.hf.default.r.5.cb.8=apply(sim.m.cb.8,1,denoise.poisson)
# est.hf.default.r.6.cb.8=apply(sim.m.cb.8,1,denoise.poisson)
# est.hf.default.r.7.cb.8=apply(sim.m.cb.8,1,denoise.poisson)

est.hf.default.r.4.b.8=apply(sim.m.b.8,1,denoise.poisson)
# est.hf.default.r.5.b.8=apply(sim.m.b.8,1,denoise.poisson)
# est.hf.default.r.6.b.8=apply(sim.m.b.8,1,denoise.poisson)
# est.hf.default.r.7.b.8=apply(sim.m.b.8,1,denoise.poisson)

est.hf.default.r.4.hs.8=apply(sim.m.hs.8,1,denoise.poisson)
# est.hf.default.r.5.hs.8=apply(sim.m.hs.8,1,denoise.poisson)
# est.hf.default.r.6.hs.8=apply(sim.m.hs.8,1,denoise.poisson)
# est.hf.default.r.7.hs.8=apply(sim.m.hs.8,1,denoise.poisson)

est.hf.default.r.4.bur.8=apply(sim.m.bur.8,1,denoise.poisson)
# est.hf.default.r.5.bur.8=apply(sim.m.bur.8,1,denoise.poisson)
# est.hf.default.r.6.bur.8=apply(sim.m.bur.8,1,denoise.poisson)
# est.hf.default.r.7.bur.8=apply(sim.m.bur.8,1,denoise.poisson)

est.hf.default.r.4.s.128=apply(sim.m.s.128,1,denoise.poisson)
# est.hf.default.r.5.s.128=apply(sim.m.s.128,1,denoise.poisson)
# est.hf.default.r.6.s.128=apply(sim.m.s.128,1,denoise.poisson)
# est.hf.default.r.7.s.128=apply(sim.m.s.128,1,denoise.poisson)

est.hf.default.r.4.ang.128=apply(sim.m.ang.128,1,denoise.poisson)
# est.hf.default.r.5.ang.128=apply(sim.m.ang.128,1,denoise.poisson)
# est.hf.default.r.6.ang.128=apply(sim.m.ang.128,1,denoise.poisson)
# est.hf.default.r.7.ang.128=apply(sim.m.ang.128,1,denoise.poisson)

est.hf.default.r.4.cb.128=apply(sim.m.cb.128,1,denoise.poisson)
# est.hf.default.r.5.cb.128=apply(sim.m.cb.128,1,denoise.poisson)
# est.hf.default.r.6.cb.128=apply(sim.m.cb.128,1,denoise.poisson)
# est.hf.default.r.7.cb.128=apply(sim.m.cb.128,1,denoise.poisson)

est.hf.default.r.4.b.128=apply(sim.m.b.128,1,denoise.poisson)
# est.hf.default.r.5.b.128=apply(sim.m.b.128,1,denoise.poisson)
# est.hf.default.r.6.b.128=apply(sim.m.b.128,1,denoise.poisson)
# est.hf.default.r.7.b.128=apply(sim.m.b.128,1,denoise.poisson)

est.hf.default.r.4.hs.128=apply(sim.m.hs.128,1,denoise.poisson)
# est.hf.default.r.5.hs.128=apply(sim.m.hs.128,1,denoise.poisson)
# est.hf.default.r.6.hs.128=apply(sim.m.hs.128,1,denoise.poisson)
# est.hf.default.r.7.hs.128=apply(sim.m.hs.128,1,denoise.poisson)

est.hf.default.r.4.bur.128=apply(sim.m.bur.128,1,denoise.poisson)
# est.hf.default.r.5.bur.128=apply(sim.m.bur.128,1,denoise.poisson)
# est.hf.default.r.6.bur.128=apply(sim.m.bur.128,1,denoise.poisson)
# est.hf.default.r.7.bur.128=apply(sim.m.bur.128,1,denoise.poisson)


#universal thresholding, with noise estimated
est.hf.u.r.4.s.8=apply(sim.m.s.8,1,denoise.poisson,meth.1=hf.u,cs.1=50,hybrid=F)
# est.hf.u.r.5.s.8=apply(sim.m.s.8,1,denoise.poisson,meth.1=hf.u,cs.1=50,hybrid=F)
# est.hf.u.r.6.s.8=apply(sim.m.s.8,1,denoise.poisson,meth.1=hf.u,cs.1=50,hybrid=F)
# est.hf.u.r.7.s.8=apply(sim.m.s.8,1,denoise.poisson,meth.1=hf.u,cs.1=50,hybrid=F)

est.hf.u.r.4.ang.8=apply(sim.m.ang.8,1,denoise.poisson,meth.1=hf.u,cs.1=50,hybrid=F)
# est.hf.u.r.5.ang.8=apply(sim.m.ang.8,1,denoise.poisson,meth.1=hf.u,cs.1=50,hybrid=F)
# est.hf.u.r.6.ang.8=apply(sim.m.ang.8,1,denoise.poisson,meth.1=hf.u,cs.1=50,hybrid=F)
# est.hf.u.r.7.ang.8=apply(sim.m.ang.8,1,denoise.poisson,meth.1=hf.u,cs.1=50,hybrid=F)

est.hf.u.r.4.cb.8=apply(sim.m.cb.8,1,denoise.poisson,meth.1=hf.u,cs.1=50,hybrid=F)
# est.hf.u.r.5.cb.8=apply(sim.m.cb.8,1,denoise.poisson,meth.1=hf.u,cs.1=50,hybrid=F)
# est.hf.u.r.6.cb.8=apply(sim.m.cb.8,1,denoise.poisson,meth.1=hf.u,cs.1=50,hybrid=F)
# est.hf.u.r.7.cb.8=apply(sim.m.cb.8,1,denoise.poisson,meth.1=hf.u,cs.1=50,hybrid=F)

est.hf.u.r.4.b.8=apply(sim.m.b.8,1,denoise.poisson,meth.1=hf.u,cs.1=50,hybrid=F)
# est.hf.u.r.5.b.8=apply(sim.m.b.8,1,denoise.poisson,meth.1=hf.u,cs.1=50,hybrid=F)
# est.hf.u.r.6.b.8=apply(sim.m.b.8,1,denoise.poisson,meth.1=hf.u,cs.1=50,hybrid=F)
# est.hf.u.r.7.b.8=apply(sim.m.b.8,1,denoise.poisson,meth.1=hf.u,cs.1=50,hybrid=F)

est.hf.u.r.4.hs.8=apply(sim.m.hs.8,1,denoise.poisson,meth.1=hf.u,cs.1=50,hybrid=F)
# est.hf.u.r.5.hs.8=apply(sim.m.hs.8,1,denoise.poisson,meth.1=hf.u,cs.1=50,hybrid=F)
# est.hf.u.r.6.hs.8=apply(sim.m.hs.8,1,denoise.poisson,meth.1=hf.u,cs.1=50,hybrid=F)
# est.hf.u.r.7.hs.8=apply(sim.m.hs.8,1,denoise.poisson,meth.1=hf.u,cs.1=50,hybrid=F)

est.hf.u.r.4.bur.8=apply(sim.m.bur.8,1,denoise.poisson,meth.1=hf.u,cs.1=50,hybrid=F)
# est.hf.u.r.5.bur.8=apply(sim.m.bur.8,1,denoise.poisson,meth.1=hf.u,cs.1=50,hybrid=F)
# est.hf.u.r.6.bur.8=apply(sim.m.bur.8,1,denoise.poisson,meth.1=hf.u,cs.1=50,hybrid=F)
# est.hf.u.r.7.bur.8=apply(sim.m.bur.8,1,denoise.poisson,meth.1=hf.u,cs.1=50,hybrid=F)

est.hf.u.r.4.s.128=apply(sim.m.s.128,1,denoise.poisson,meth.1=hf.u,cs.1=50,hybrid=F)
# est.hf.u.r.5.s.128=apply(sim.m.s.128,1,denoise.poisson,meth.1=hf.u,cs.1=50,hybrid=F)
# est.hf.u.r.6.s.128=apply(sim.m.s.128,1,denoise.poisson,meth.1=hf.u,cs.1=50,hybrid=F)
# est.hf.u.r.7.s.128=apply(sim.m.s.128,1,denoise.poisson,meth.1=hf.u,cs.1=50,hybrid=F)

est.hf.u.r.4.ang.128=apply(sim.m.ang.128,1,denoise.poisson,meth.1=hf.u,cs.1=50,hybrid=F)
# est.hf.u.r.5.ang.128=apply(sim.m.ang.128,1,denoise.poisson,meth.1=hf.u,cs.1=50,hybrid=F)
# est.hf.u.r.6.ang.128=apply(sim.m.ang.128,1,denoise.poisson,meth.1=hf.u,cs.1=50,hybrid=F)
# est.hf.u.r.7.ang.128=apply(sim.m.ang.128,1,denoise.poisson,meth.1=hf.u,cs.1=50,hybrid=F)

est.hf.u.r.4.cb.128=apply(sim.m.cb.128,1,denoise.poisson,meth.1=hf.u,cs.1=50,hybrid=F)
# est.hf.u.r.5.cb.128=apply(sim.m.cb.128,1,denoise.poisson,meth.1=hf.u,cs.1=50,hybrid=F)
# est.hf.u.r.6.cb.128=apply(sim.m.cb.128,1,denoise.poisson,meth.1=hf.u,cs.1=50,hybrid=F)
# est.hf.u.r.7.cb.128=apply(sim.m.cb.128,1,denoise.poisson,meth.1=hf.u,cs.1=50,hybrid=F)

est.hf.u.r.4.b.128=apply(sim.m.b.128,1,denoise.poisson,meth.1=hf.u,cs.1=50,hybrid=F)
# est.hf.u.r.5.b.128=apply(sim.m.b.128,1,denoise.poisson,meth.1=hf.u,cs.1=50,hybrid=F)
# est.hf.u.r.6.b.128=apply(sim.m.b.128,1,denoise.poisson,meth.1=hf.u,cs.1=50,hybrid=F)
# est.hf.u.r.7.b.128=apply(sim.m.b.128,1,denoise.poisson,meth.1=hf.u,cs.1=50,hybrid=F)

est.hf.u.r.4.hs.128=apply(sim.m.hs.128,1,denoise.poisson,meth.1=hf.u,cs.1=50,hybrid=F)
# est.hf.u.r.5.hs.128=apply(sim.m.hs.128,1,denoise.poisson,meth.1=hf.u,cs.1=50,hybrid=F)
# est.hf.u.r.6.hs.128=apply(sim.m.hs.128,1,denoise.poisson,meth.1=hf.u,cs.1=50,hybrid=F)
# est.hf.u.r.7.hs.128=apply(sim.m.hs.128,1,denoise.poisson,meth.1=hf.u,cs.1=50,hybrid=F)

est.hf.u.r.4.bur.128=apply(sim.m.bur.128,1,denoise.poisson,meth.1=hf.u,cs.1=50,hybrid=F)
# est.hf.u.r.5.bur.128=apply(sim.m.bur.128,1,denoise.poisson,meth.1=hf.u,cs.1=50,hybrid=F)
# est.hf.u.r.6.bur.128=apply(sim.m.bur.128,1,denoise.poisson,meth.1=hf.u,cs.1=50,hybrid=F)
# est.hf.u.r.7.bur.128=apply(sim.m.bur.128,1,denoise.poisson,meth.1=hf.u,cs.1=50,hybrid=F)

#tiu, with noise level 1
est.hf.ti.r.4.s.8=apply(sim.m.s.8,1,denoise.poisson,meth.1=hf.la10.ti4,cs.1=50,hybrid=F)
est.hf.ti.r.5.s.8=apply(sim.m.s.8,1,denoise.poisson,meth.1=hf.la10.ti5,cs.1=50,hybrid=F)
est.hf.ti.r.6.s.8=apply(sim.m.s.8,1,denoise.poisson,meth.1=hf.la10.ti6,cs.1=50,hybrid=F)
est.hf.ti.r.7.s.8=apply(sim.m.s.8,1,denoise.poisson,meth.1=hf.la10.ti7,cs.1=50,hybrid=F)

est.hf.ti.r.4.ang.8=apply(sim.m.ang.8,1,denoise.poisson,meth.1=hf.la10.ti4,cs.1=50,hybrid=F)
est.hf.ti.r.5.ang.8=apply(sim.m.ang.8,1,denoise.poisson,meth.1=hf.la10.ti5,cs.1=50,hybrid=F)
est.hf.ti.r.6.ang.8=apply(sim.m.ang.8,1,denoise.poisson,meth.1=hf.la10.ti6,cs.1=50,hybrid=F)
est.hf.ti.r.7.ang.8=apply(sim.m.ang.8,1,denoise.poisson,meth.1=hf.la10.ti7,cs.1=50,hybrid=F)

est.hf.ti.r.4.cb.8=apply(sim.m.cb.8,1,denoise.poisson,meth.1=hf.la10.ti4,cs.1=50,hybrid=F)
est.hf.ti.r.5.cb.8=apply(sim.m.cb.8,1,denoise.poisson,meth.1=hf.la10.ti5,cs.1=50,hybrid=F)
est.hf.ti.r.6.cb.8=apply(sim.m.cb.8,1,denoise.poisson,meth.1=hf.la10.ti6,cs.1=50,hybrid=F)
est.hf.ti.r.7.cb.8=apply(sim.m.cb.8,1,denoise.poisson,meth.1=hf.la10.ti7,cs.1=50,hybrid=F)

est.hf.ti.r.4.b.8=apply(sim.m.b.8,1,denoise.poisson,meth.1=hf.la10.ti4,cs.1=50,hybrid=F)
est.hf.ti.r.5.b.8=apply(sim.m.b.8,1,denoise.poisson,meth.1=hf.la10.ti5,cs.1=50,hybrid=F)
est.hf.ti.r.6.b.8=apply(sim.m.b.8,1,denoise.poisson,meth.1=hf.la10.ti6,cs.1=50,hybrid=F)
est.hf.ti.r.7.b.8=apply(sim.m.b.8,1,denoise.poisson,meth.1=hf.la10.ti7,cs.1=50,hybrid=F)

est.hf.ti.r.4.hs.8=apply(sim.m.hs.8,1,denoise.poisson,meth.1=hf.la10.ti4,cs.1=50,hybrid=F)
est.hf.ti.r.5.hs.8=apply(sim.m.hs.8,1,denoise.poisson,meth.1=hf.la10.ti5,cs.1=50,hybrid=F)
est.hf.ti.r.6.hs.8=apply(sim.m.hs.8,1,denoise.poisson,meth.1=hf.la10.ti6,cs.1=50,hybrid=F)
est.hf.ti.r.7.hs.8=apply(sim.m.hs.8,1,denoise.poisson,meth.1=hf.la10.ti7,cs.1=50,hybrid=F)

est.hf.ti.r.4.bur.8=apply(sim.m.bur.8,1,denoise.poisson,meth.1=hf.la10.ti4,cs.1=50,hybrid=F)
est.hf.ti.r.5.bur.8=apply(sim.m.bur.8,1,denoise.poisson,meth.1=hf.la10.ti5,cs.1=50,hybrid=F)
est.hf.ti.r.6.bur.8=apply(sim.m.bur.8,1,denoise.poisson,meth.1=hf.la10.ti6,cs.1=50,hybrid=F)
est.hf.ti.r.7.bur.8=apply(sim.m.bur.8,1,denoise.poisson,meth.1=hf.la10.ti7,cs.1=50,hybrid=F)

est.hf.ti.r.4.s.128=apply(sim.m.s.128,1,denoise.poisson,meth.1=hf.la10.ti4,cs.1=50,hybrid=F)
est.hf.ti.r.5.s.128=apply(sim.m.s.128,1,denoise.poisson,meth.1=hf.la10.ti5,cs.1=50,hybrid=F)
est.hf.ti.r.6.s.128=apply(sim.m.s.128,1,denoise.poisson,meth.1=hf.la10.ti6,cs.1=50,hybrid=F)
est.hf.ti.r.7.s.128=apply(sim.m.s.128,1,denoise.poisson,meth.1=hf.la10.ti7,cs.1=50,hybrid=F)

est.hf.ti.r.4.ang.128=apply(sim.m.ang.128,1,denoise.poisson,meth.1=hf.la10.ti4,cs.1=50,hybrid=F)
est.hf.ti.r.5.ang.128=apply(sim.m.ang.128,1,denoise.poisson,meth.1=hf.la10.ti5,cs.1=50,hybrid=F)
est.hf.ti.r.6.ang.128=apply(sim.m.ang.128,1,denoise.poisson,meth.1=hf.la10.ti6,cs.1=50,hybrid=F)
est.hf.ti.r.7.ang.128=apply(sim.m.ang.128,1,denoise.poisson,meth.1=hf.la10.ti7,cs.1=50,hybrid=F)

est.hf.ti.r.4.cb.128=apply(sim.m.cb.128,1,denoise.poisson,meth.1=hf.la10.ti4,cs.1=50,hybrid=F)
est.hf.ti.r.5.cb.128=apply(sim.m.cb.128,1,denoise.poisson,meth.1=hf.la10.ti5,cs.1=50,hybrid=F)
est.hf.ti.r.6.cb.128=apply(sim.m.cb.128,1,denoise.poisson,meth.1=hf.la10.ti6,cs.1=50,hybrid=F)
est.hf.ti.r.7.cb.128=apply(sim.m.cb.128,1,denoise.poisson,meth.1=hf.la10.ti7,cs.1=50,hybrid=F)

est.hf.ti.r.4.b.128=apply(sim.m.b.128,1,denoise.poisson,meth.1=hf.la10.ti4,cs.1=50,hybrid=F)
est.hf.ti.r.5.b.128=apply(sim.m.b.128,1,denoise.poisson,meth.1=hf.la10.ti5,cs.1=50,hybrid=F)
est.hf.ti.r.6.b.128=apply(sim.m.b.128,1,denoise.poisson,meth.1=hf.la10.ti6,cs.1=50,hybrid=F)
est.hf.ti.r.7.b.128=apply(sim.m.b.128,1,denoise.poisson,meth.1=hf.la10.ti7,cs.1=50,hybrid=F)

est.hf.ti.r.4.hs.128=apply(sim.m.hs.128,1,denoise.poisson,meth.1=hf.la10.ti4,cs.1=50,hybrid=F)
est.hf.ti.r.5.hs.128=apply(sim.m.hs.128,1,denoise.poisson,meth.1=hf.la10.ti5,cs.1=50,hybrid=F)
est.hf.ti.r.6.hs.128=apply(sim.m.hs.128,1,denoise.poisson,meth.1=hf.la10.ti6,cs.1=50,hybrid=F)
est.hf.ti.r.7.hs.128=apply(sim.m.hs.128,1,denoise.poisson,meth.1=hf.la10.ti7,cs.1=50,hybrid=F)

est.hf.ti.r.4.bur.128=apply(sim.m.bur.128,1,denoise.poisson,meth.1=hf.la10.ti4,cs.1=50,hybrid=F)
est.hf.ti.r.5.bur.128=apply(sim.m.bur.128,1,denoise.poisson,meth.1=hf.la10.ti5,cs.1=50,hybrid=F)
est.hf.ti.r.6.bur.128=apply(sim.m.bur.128,1,denoise.poisson,meth.1=hf.la10.ti6,cs.1=50,hybrid=F)
est.hf.ti.r.7.bur.128=apply(sim.m.bur.128,1,denoise.poisson,meth.1=hf.la10.ti7,cs.1=50,hybrid=F)


#tiu, with noise estimated
est.hf.ti.est.r.4.s.8=apply(sim.m.s.8,1,denoise.poisson,meth.1=hf.la10.est.ti4,cs.1=50,hybrid=F)
est.hf.ti.est.r.5.s.8=apply(sim.m.s.8,1,denoise.poisson,meth.1=hf.la10.est.ti5,cs.1=50,hybrid=F)
est.hf.ti.est.r.6.s.8=apply(sim.m.s.8,1,denoise.poisson,meth.1=hf.la10.est.ti6,cs.1=50,hybrid=F)
est.hf.ti.est.r.7.s.8=apply(sim.m.s.8,1,denoise.poisson,meth.1=hf.la10.est.ti7,cs.1=50,hybrid=F)

est.hf.ti.est.r.4.ang.8=apply(sim.m.ang.8,1,denoise.poisson,meth.1=hf.la10.est.ti4,cs.1=50,hybrid=F)
est.hf.ti.est.r.5.ang.8=apply(sim.m.ang.8,1,denoise.poisson,meth.1=hf.la10.est.ti5,cs.1=50,hybrid=F)
est.hf.ti.est.r.6.ang.8=apply(sim.m.ang.8,1,denoise.poisson,meth.1=hf.la10.est.ti6,cs.1=50,hybrid=F)
est.hf.ti.est.r.7.ang.8=apply(sim.m.ang.8,1,denoise.poisson,meth.1=hf.la10.est.ti7,cs.1=50,hybrid=F)

est.hf.ti.est.r.4.cb.8=apply(sim.m.cb.8,1,denoise.poisson,meth.1=hf.la10.est.ti4,cs.1=50,hybrid=F)
est.hf.ti.est.r.5.cb.8=apply(sim.m.cb.8,1,denoise.poisson,meth.1=hf.la10.est.ti5,cs.1=50,hybrid=F)
est.hf.ti.est.r.6.cb.8=apply(sim.m.cb.8,1,denoise.poisson,meth.1=hf.la10.est.ti6,cs.1=50,hybrid=F)
est.hf.ti.est.r.7.cb.8=apply(sim.m.cb.8,1,denoise.poisson,meth.1=hf.la10.est.ti7,cs.1=50,hybrid=F)

est.hf.ti.est.r.4.b.8=apply(sim.m.b.8,1,denoise.poisson,meth.1=hf.la10.est.ti4,cs.1=50,hybrid=F)
est.hf.ti.est.r.5.b.8=apply(sim.m.b.8,1,denoise.poisson,meth.1=hf.la10.est.ti5,cs.1=50,hybrid=F)
est.hf.ti.est.r.6.b.8=apply(sim.m.b.8,1,denoise.poisson,meth.1=hf.la10.est.ti6,cs.1=50,hybrid=F)
est.hf.ti.est.r.7.b.8=apply(sim.m.b.8,1,denoise.poisson,meth.1=hf.la10.est.ti7,cs.1=50,hybrid=F)

est.hf.ti.est.r.4.hs.8=apply(sim.m.hs.8,1,denoise.poisson,meth.1=hf.la10.est.ti4,cs.1=50,hybrid=F)
est.hf.ti.est.r.5.hs.8=apply(sim.m.hs.8,1,denoise.poisson,meth.1=hf.la10.est.ti5,cs.1=50,hybrid=F)
est.hf.ti.est.r.6.hs.8=apply(sim.m.hs.8,1,denoise.poisson,meth.1=hf.la10.est.ti6,cs.1=50,hybrid=F)
est.hf.ti.est.r.7.hs.8=apply(sim.m.hs.8,1,denoise.poisson,meth.1=hf.la10.est.ti7,cs.1=50,hybrid=F)

est.hf.ti.est.r.4.bur.8=apply(sim.m.bur.8,1,denoise.poisson,meth.1=hf.la10.est.ti4,cs.1=50,hybrid=F)
est.hf.ti.est.r.5.bur.8=apply(sim.m.bur.8,1,denoise.poisson,meth.1=hf.la10.est.ti5,cs.1=50,hybrid=F)
est.hf.ti.est.r.6.bur.8=apply(sim.m.bur.8,1,denoise.poisson,meth.1=hf.la10.est.ti6,cs.1=50,hybrid=F)
est.hf.ti.est.r.7.bur.8=apply(sim.m.bur.8,1,denoise.poisson,meth.1=hf.la10.est.ti7,cs.1=50,hybrid=F)

est.hf.ti.est.r.4.s.128=apply(sim.m.s.128,1,denoise.poisson,meth.1=hf.la10.est.ti4,cs.1=50,hybrid=F)
est.hf.ti.est.r.5.s.128=apply(sim.m.s.128,1,denoise.poisson,meth.1=hf.la10.est.ti5,cs.1=50,hybrid=F)
est.hf.ti.est.r.6.s.128=apply(sim.m.s.128,1,denoise.poisson,meth.1=hf.la10.est.ti6,cs.1=50,hybrid=F)
est.hf.ti.est.r.7.s.128=apply(sim.m.s.128,1,denoise.poisson,meth.1=hf.la10.est.ti7,cs.1=50,hybrid=F)

est.hf.ti.est.r.4.ang.128=apply(sim.m.ang.128,1,denoise.poisson,meth.1=hf.la10.est.ti4,cs.1=50,hybrid=F)
est.hf.ti.est.r.5.ang.128=apply(sim.m.ang.128,1,denoise.poisson,meth.1=hf.la10.est.ti5,cs.1=50,hybrid=F)
est.hf.ti.est.r.6.ang.128=apply(sim.m.ang.128,1,denoise.poisson,meth.1=hf.la10.est.ti6,cs.1=50,hybrid=F)
est.hf.ti.est.r.7.ang.128=apply(sim.m.ang.128,1,denoise.poisson,meth.1=hf.la10.est.ti7,cs.1=50,hybrid=F)

est.hf.ti.est.r.4.cb.128=apply(sim.m.cb.128,1,denoise.poisson,meth.1=hf.la10.est.ti4,cs.1=50,hybrid=F)
est.hf.ti.est.r.5.cb.128=apply(sim.m.cb.128,1,denoise.poisson,meth.1=hf.la10.est.ti5,cs.1=50,hybrid=F)
est.hf.ti.est.r.6.cb.128=apply(sim.m.cb.128,1,denoise.poisson,meth.1=hf.la10.est.ti6,cs.1=50,hybrid=F)
est.hf.ti.est.r.7.cb.128=apply(sim.m.cb.128,1,denoise.poisson,meth.1=hf.la10.est.ti7,cs.1=50,hybrid=F)

est.hf.ti.est.r.4.b.128=apply(sim.m.b.128,1,denoise.poisson,meth.1=hf.la10.est.ti4,cs.1=50,hybrid=F)
est.hf.ti.est.r.5.b.128=apply(sim.m.b.128,1,denoise.poisson,meth.1=hf.la10.est.ti5,cs.1=50,hybrid=F)
est.hf.ti.est.r.6.b.128=apply(sim.m.b.128,1,denoise.poisson,meth.1=hf.la10.est.ti6,cs.1=50,hybrid=F)
est.hf.ti.est.r.7.b.128=apply(sim.m.b.128,1,denoise.poisson,meth.1=hf.la10.est.ti7,cs.1=50,hybrid=F)

est.hf.ti.est.r.4.hs.128=apply(sim.m.hs.128,1,denoise.poisson,meth.1=hf.la10.est.ti4,cs.1=50,hybrid=F)
est.hf.ti.est.r.5.hs.128=apply(sim.m.hs.128,1,denoise.poisson,meth.1=hf.la10.est.ti5,cs.1=50,hybrid=F)
est.hf.ti.est.r.6.hs.128=apply(sim.m.hs.128,1,denoise.poisson,meth.1=hf.la10.est.ti6,cs.1=50,hybrid=F)
est.hf.ti.est.r.7.hs.128=apply(sim.m.hs.128,1,denoise.poisson,meth.1=hf.la10.est.ti7,cs.1=50,hybrid=F)

est.hf.ti.est.r.4.bur.128=apply(sim.m.bur.128,1,denoise.poisson,meth.1=hf.la10.est.ti4,cs.1=50,hybrid=F)
est.hf.ti.est.r.5.bur.128=apply(sim.m.bur.128,1,denoise.poisson,meth.1=hf.la10.est.ti5,cs.1=50,hybrid=F)
est.hf.ti.est.r.6.bur.128=apply(sim.m.bur.128,1,denoise.poisson,meth.1=hf.la10.est.ti6,cs.1=50,hybrid=F)
est.hf.ti.est.r.7.bur.128=apply(sim.m.bur.128,1,denoise.poisson,meth.1=hf.la10.est.ti7,cs.1=50,hybrid=F)



#compute mises
mise.hf.default.r.4.s.8=mise(t(est.hf.default.r.4.s.8),mu.s.8)
# mise.hf.default.r.5.s.8=mise(t(est.hf.default.r.5.s.8),mu.s.8)
# mise.hf.default.r.6.s.8=mise(t(est.hf.default.r.6.s.8),mu.s.8)
# mise.hf.default.r.7.s.8=mise(t(est.hf.default.r.7.s.8),mu.s.8)

mise.hf.default.r.4.ang.8=mise(t(est.hf.default.r.4.ang.8),mu.ang.8)
# mise.hf.default.r.5.ang.8=mise(t(est.hf.default.r.5.ang.8),mu.ang.8)
# mise.hf.default.r.6.ang.8=mise(t(est.hf.default.r.6.ang.8),mu.ang.8)
# mise.hf.default.r.7.ang.8=mise(t(est.hf.default.r.7.ang.8),mu.ang.8)

mise.hf.default.r.4.cb.8=mise(t(est.hf.default.r.4.cb.8),mu.cb.8)
# mise.hf.default.r.5.cb.8=mise(t(est.hf.default.r.5.cb.8),mu.cb.8)
# mise.hf.default.r.6.cb.8=mise(t(est.hf.default.r.6.cb.8),mu.cb.8)
# mise.hf.default.r.7.cb.8=mise(t(est.hf.default.r.7.cb.8),mu.cb.8)

mise.hf.default.r.4.b.8=mise(t(est.hf.default.r.4.b.8),mu.b.8)
# mise.hf.default.r.5.b.8=mise(t(est.hf.default.r.5.b.8),mu.b.8)
# mise.hf.default.r.6.b.8=mise(t(est.hf.default.r.6.b.8),mu.b.8)
# mise.hf.default.r.7.b.8=mise(t(est.hf.default.r.7.b.8),mu.b.8)

mise.hf.default.r.4.hs.8=mise(t(est.hf.default.r.4.hs.8),mu.hs.8)
# mise.hf.default.r.5.hs.8=mise(t(est.hf.default.r.5.hs.8),mu.hs.8)
# mise.hf.default.r.6.hs.8=mise(t(est.hf.default.r.6.hs.8),mu.hs.8)
# mise.hf.default.r.7.hs.8=mise(t(est.hf.default.r.7.hs.8),mu.hs.8)

mise.hf.default.r.4.bur.8=mise(t(est.hf.default.r.4.bur.8),mu.bur.8)
# mise.hf.default.r.5.bur.8=mise(t(est.hf.default.r.5.bur.8),mu.bur.8)
# mise.hf.default.r.6.bur.8=mise(t(est.hf.default.r.6.bur.8),mu.bur.8)
# mise.hf.default.r.7.bur.8=mise(t(est.hf.default.r.7.bur.8),mu.bur.8)

mise.hf.default.r.4.s.128=mise(t(est.hf.default.r.4.s.128),mu.s.128)
# mise.hf.default.r.5.s.128=mise(t(est.hf.default.r.5.s.128),mu.s.128)
# mise.hf.default.r.6.s.128=mise(t(est.hf.default.r.6.s.128),mu.s.128)
# mise.hf.default.r.7.s.128=mise(t(est.hf.default.r.7.s.128),mu.s.128)

mise.hf.default.r.4.ang.128=mise(t(est.hf.default.r.4.ang.128),mu.ang.128)
# mise.hf.default.r.5.ang.128=mise(t(est.hf.default.r.5.ang.128),mu.ang.128)
# mise.hf.default.r.6.ang.128=mise(t(est.hf.default.r.6.ang.128),mu.ang.128)
# mise.hf.default.r.7.ang.128=mise(t(est.hf.default.r.7.ang.128),mu.ang.128)

mise.hf.default.r.4.cb.128=mise(t(est.hf.default.r.4.cb.128),mu.cb.128)
# mise.hf.default.r.5.cb.128=mise(t(est.hf.default.r.5.cb.128),mu.cb.128)
# mise.hf.default.r.6.cb.128=mise(t(est.hf.default.r.6.cb.128),mu.cb.128)
# mise.hf.default.r.7.cb.128=mise(t(est.hf.default.r.7.cb.128),mu.cb.128)

mise.hf.default.r.4.b.128=mise(t(est.hf.default.r.4.b.128),mu.b.128)
# mise.hf.default.r.5.b.128=mise(t(est.hf.default.r.5.b.128),mu.b.128)
# mise.hf.default.r.6.b.128=mise(t(est.hf.default.r.6.b.128),mu.b.128)
# mise.hf.default.r.7.b.128=mise(t(est.hf.default.r.7.b.128),mu.b.128)

mise.hf.default.r.4.hs.128=mise(t(est.hf.default.r.4.hs.128),mu.hs.128)
# mise.hf.default.r.5.hs.128=mise(t(est.hf.default.r.5.hs.128),mu.hs.128)
# mise.hf.default.r.6.hs.128=mise(t(est.hf.default.r.6.hs.128),mu.hs.128)
# mise.hf.default.r.7.hs.128=mise(t(est.hf.default.r.7.hs.128),mu.hs.128)

mise.hf.default.r.4.bur.128=mise(t(est.hf.default.r.4.bur.128),mu.bur.128)
# mise.hf.default.r.5.bur.128=mise(t(est.hf.default.r.5.bur.128),mu.bur.128)
# mise.hf.default.r.6.bur.128=mise(t(est.hf.default.r.6.bur.128),mu.bur.128)
# mise.hf.default.r.7.bur.128=mise(t(est.hf.default.r.7.bur.128),mu.bur.128)


mise.hf.u.r.4.s.8=mise(t(est.hf.u.r.4.s.8),mu.s.8)
# mise.hf.u.r.5.s.8=mise(t(est.hf.u.r.5.s.8),mu.s.8)
# mise.hf.u.r.6.s.8=mise(t(est.hf.u.r.6.s.8),mu.s.8)
# mise.hf.u.r.7.s.8=mise(t(est.hf.u.r.7.s.8),mu.s.8)

mise.hf.u.r.4.ang.8=mise(t(est.hf.u.r.4.ang.8),mu.ang.8)
# mise.hf.u.r.5.ang.8=mise(t(est.hf.u.r.5.ang.8),mu.ang.8)
# mise.hf.u.r.6.ang.8=mise(t(est.hf.u.r.6.ang.8),mu.ang.8)
# mise.hf.u.r.7.ang.8=mise(t(est.hf.u.r.7.ang.8),mu.ang.8)

mise.hf.u.r.4.cb.8=mise(t(est.hf.u.r.4.cb.8),mu.cb.8)
# mise.hf.u.r.5.cb.8=mise(t(est.hf.u.r.5.cb.8),mu.cb.8)
# mise.hf.u.r.6.cb.8=mise(t(est.hf.u.r.6.cb.8),mu.cb.8)
# mise.hf.u.r.7.cb.8=mise(t(est.hf.u.r.7.cb.8),mu.cb.8)

mise.hf.u.r.4.b.8=mise(t(est.hf.u.r.4.b.8),mu.b.8)
# mise.hf.u.r.5.b.8=mise(t(est.hf.u.r.5.b.8),mu.b.8)
# mise.hf.u.r.6.b.8=mise(t(est.hf.u.r.6.b.8),mu.b.8)
# mise.hf.u.r.7.b.8=mise(t(est.hf.u.r.7.b.8),mu.b.8)

mise.hf.u.r.4.hs.8=mise(t(est.hf.u.r.4.hs.8),mu.hs.8)
# mise.hf.u.r.5.hs.8=mise(t(est.hf.u.r.5.hs.8),mu.hs.8)
# mise.hf.u.r.6.hs.8=mise(t(est.hf.u.r.6.hs.8),mu.hs.8)
# mise.hf.u.r.7.hs.8=mise(t(est.hf.u.r.7.hs.8),mu.hs.8)

mise.hf.u.r.4.bur.8=mise(t(est.hf.u.r.4.bur.8),mu.bur.8)
# mise.hf.u.r.5.bur.8=mise(t(est.hf.u.r.5.bur.8),mu.bur.8)
# mise.hf.u.r.6.bur.8=mise(t(est.hf.u.r.6.bur.8),mu.bur.8)
# mise.hf.u.r.7.bur.8=mise(t(est.hf.u.r.7.bur.8),mu.bur.8)

mise.hf.u.r.4.s.128=mise(t(est.hf.u.r.4.s.128),mu.s.128)
# mise.hf.u.r.5.s.128=mise(t(est.hf.u.r.5.s.128),mu.s.128)
# mise.hf.u.r.6.s.128=mise(t(est.hf.u.r.6.s.128),mu.s.128)
# mise.hf.u.r.7.s.128=mise(t(est.hf.u.r.7.s.128),mu.s.128)

mise.hf.u.r.4.ang.128=mise(t(est.hf.u.r.4.ang.128),mu.ang.128)
# mise.hf.u.r.5.ang.128=mise(t(est.hf.u.r.5.ang.128),mu.ang.128)
# mise.hf.u.r.6.ang.128=mise(t(est.hf.u.r.6.ang.128),mu.ang.128)
# mise.hf.u.r.7.ang.128=mise(t(est.hf.u.r.7.ang.128),mu.ang.128)

mise.hf.u.r.4.cb.128=mise(t(est.hf.u.r.4.cb.128),mu.cb.128)
# mise.hf.u.r.5.cb.128=mise(t(est.hf.u.r.5.cb.128),mu.cb.128)
# mise.hf.u.r.6.cb.128=mise(t(est.hf.u.r.6.cb.128),mu.cb.128)
# mise.hf.u.r.7.cb.128=mise(t(est.hf.u.r.7.cb.128),mu.cb.128)

mise.hf.u.r.4.b.128=mise(t(est.hf.u.r.4.b.128),mu.b.128)
# mise.hf.u.r.5.b.128=mise(t(est.hf.u.r.5.b.128),mu.b.128)
# mise.hf.u.r.6.b.128=mise(t(est.hf.u.r.6.b.128),mu.b.128)
# mise.hf.u.r.7.b.128=mise(t(est.hf.u.r.7.b.128),mu.b.128)

mise.hf.u.r.4.hs.128=mise(t(est.hf.u.r.4.hs.128),mu.hs.128)
# mise.hf.u.r.5.hs.128=mise(t(est.hf.u.r.5.hs.128),mu.hs.128)
# mise.hf.u.r.6.hs.128=mise(t(est.hf.u.r.6.hs.128),mu.hs.128)
# mise.hf.u.r.7.hs.128=mise(t(est.hf.u.r.7.hs.128),mu.hs.128)

mise.hf.u.r.4.bur.128=mise(t(est.hf.u.r.4.bur.128),mu.bur.128)
# mise.hf.u.r.5.bur.128=mise(t(est.hf.u.r.5.bur.128),mu.bur.128)
# mise.hf.u.r.6.bur.128=mise(t(est.hf.u.r.6.bur.128),mu.bur.128)
# mise.hf.u.r.7.bur.128=mise(t(est.hf.u.r.7.bur.128),mu.bur.128)


mise.hf.ti.r.4.s.8=mise(t(est.hf.ti.r.4.s.8),mu.s.8)
mise.hf.ti.r.5.s.8=mise(t(est.hf.ti.r.5.s.8),mu.s.8)
mise.hf.ti.r.6.s.8=mise(t(est.hf.ti.r.6.s.8),mu.s.8)
mise.hf.ti.r.7.s.8=mise(t(est.hf.ti.r.7.s.8),mu.s.8)

mise.hf.ti.r.4.ang.8=mise(t(est.hf.ti.r.4.ang.8),mu.ang.8)
mise.hf.ti.r.5.ang.8=mise(t(est.hf.ti.r.5.ang.8),mu.ang.8)
mise.hf.ti.r.6.ang.8=mise(t(est.hf.ti.r.6.ang.8),mu.ang.8)
mise.hf.ti.r.7.ang.8=mise(t(est.hf.ti.r.7.ang.8),mu.ang.8)

mise.hf.ti.r.4.cb.8=mise(t(est.hf.ti.r.4.cb.8),mu.cb.8)
mise.hf.ti.r.5.cb.8=mise(t(est.hf.ti.r.5.cb.8),mu.cb.8)
mise.hf.ti.r.6.cb.8=mise(t(est.hf.ti.r.6.cb.8),mu.cb.8)
mise.hf.ti.r.7.cb.8=mise(t(est.hf.ti.r.7.cb.8),mu.cb.8)

mise.hf.ti.r.4.b.8=mise(t(est.hf.ti.r.4.b.8),mu.b.8)
mise.hf.ti.r.5.b.8=mise(t(est.hf.ti.r.5.b.8),mu.b.8)
mise.hf.ti.r.6.b.8=mise(t(est.hf.ti.r.6.b.8),mu.b.8)
mise.hf.ti.r.7.b.8=mise(t(est.hf.ti.r.7.b.8),mu.b.8)

mise.hf.ti.r.4.hs.8=mise(t(est.hf.ti.r.4.hs.8),mu.hs.8)
mise.hf.ti.r.5.hs.8=mise(t(est.hf.ti.r.5.hs.8),mu.hs.8)
mise.hf.ti.r.6.hs.8=mise(t(est.hf.ti.r.6.hs.8),mu.hs.8)
mise.hf.ti.r.7.hs.8=mise(t(est.hf.ti.r.7.hs.8),mu.hs.8)

mise.hf.ti.r.4.bur.8=mise(t(est.hf.ti.r.4.bur.8),mu.bur.8)
mise.hf.ti.r.5.bur.8=mise(t(est.hf.ti.r.5.bur.8),mu.bur.8)
mise.hf.ti.r.6.bur.8=mise(t(est.hf.ti.r.6.bur.8),mu.bur.8)
mise.hf.ti.r.7.bur.8=mise(t(est.hf.ti.r.7.bur.8),mu.bur.8)

mise.hf.ti.r.4.s.128=mise(t(est.hf.ti.r.4.s.128),mu.s.128)
mise.hf.ti.r.5.s.128=mise(t(est.hf.ti.r.5.s.128),mu.s.128)
mise.hf.ti.r.6.s.128=mise(t(est.hf.ti.r.6.s.128),mu.s.128)
mise.hf.ti.r.7.s.128=mise(t(est.hf.ti.r.7.s.128),mu.s.128)

mise.hf.ti.r.4.ang.128=mise(t(est.hf.ti.r.4.ang.128),mu.ang.128)
mise.hf.ti.r.5.ang.128=mise(t(est.hf.ti.r.5.ang.128),mu.ang.128)
mise.hf.ti.r.6.ang.128=mise(t(est.hf.ti.r.6.ang.128),mu.ang.128)
mise.hf.ti.r.7.ang.128=mise(t(est.hf.ti.r.7.ang.128),mu.ang.128)

mise.hf.ti.r.4.cb.128=mise(t(est.hf.ti.r.4.cb.128),mu.cb.128)
mise.hf.ti.r.5.cb.128=mise(t(est.hf.ti.r.5.cb.128),mu.cb.128)
mise.hf.ti.r.6.cb.128=mise(t(est.hf.ti.r.6.cb.128),mu.cb.128)
mise.hf.ti.r.7.cb.128=mise(t(est.hf.ti.r.7.cb.128),mu.cb.128)

mise.hf.ti.r.4.b.128=mise(t(est.hf.ti.r.4.b.128),mu.b.128)
mise.hf.ti.r.5.b.128=mise(t(est.hf.ti.r.5.b.128),mu.b.128)
mise.hf.ti.r.6.b.128=mise(t(est.hf.ti.r.6.b.128),mu.b.128)
mise.hf.ti.r.7.b.128=mise(t(est.hf.ti.r.7.b.128),mu.b.128)

mise.hf.ti.r.4.hs.128=mise(t(est.hf.ti.r.4.hs.128),mu.hs.128)
mise.hf.ti.r.5.hs.128=mise(t(est.hf.ti.r.5.hs.128),mu.hs.128)
mise.hf.ti.r.6.hs.128=mise(t(est.hf.ti.r.6.hs.128),mu.hs.128)
mise.hf.ti.r.7.hs.128=mise(t(est.hf.ti.r.7.hs.128),mu.hs.128)

mise.hf.ti.r.4.bur.128=mise(t(est.hf.ti.r.4.bur.128),mu.bur.128)
mise.hf.ti.r.5.bur.128=mise(t(est.hf.ti.r.5.bur.128),mu.bur.128)
mise.hf.ti.r.6.bur.128=mise(t(est.hf.ti.r.6.bur.128),mu.bur.128)
mise.hf.ti.r.7.bur.128=mise(t(est.hf.ti.r.7.bur.128),mu.bur.128)


mise.hf.ti.est.r.4.s.8=mise(t(est.hf.ti.est.r.4.s.8),mu.s.8)
mise.hf.ti.est.r.5.s.8=mise(t(est.hf.ti.est.r.5.s.8),mu.s.8)
mise.hf.ti.est.r.6.s.8=mise(t(est.hf.ti.est.r.6.s.8),mu.s.8)
mise.hf.ti.est.r.7.s.8=mise(t(est.hf.ti.est.r.7.s.8),mu.s.8)

mise.hf.ti.est.r.4.ang.8=mise(t(est.hf.ti.est.r.4.ang.8),mu.ang.8)
mise.hf.ti.est.r.5.ang.8=mise(t(est.hf.ti.est.r.5.ang.8),mu.ang.8)
mise.hf.ti.est.r.6.ang.8=mise(t(est.hf.ti.est.r.6.ang.8),mu.ang.8)
mise.hf.ti.est.r.7.ang.8=mise(t(est.hf.ti.est.r.7.ang.8),mu.ang.8)

mise.hf.ti.est.r.4.cb.8=mise(t(est.hf.ti.est.r.4.cb.8),mu.cb.8)
mise.hf.ti.est.r.5.cb.8=mise(t(est.hf.ti.est.r.5.cb.8),mu.cb.8)
mise.hf.ti.est.r.6.cb.8=mise(t(est.hf.ti.est.r.6.cb.8),mu.cb.8)
mise.hf.ti.est.r.7.cb.8=mise(t(est.hf.ti.est.r.7.cb.8),mu.cb.8)

mise.hf.ti.est.r.4.b.8=mise(t(est.hf.ti.est.r.4.b.8),mu.b.8)
mise.hf.ti.est.r.5.b.8=mise(t(est.hf.ti.est.r.5.b.8),mu.b.8)
mise.hf.ti.est.r.6.b.8=mise(t(est.hf.ti.est.r.6.b.8),mu.b.8)
mise.hf.ti.est.r.7.b.8=mise(t(est.hf.ti.est.r.7.b.8),mu.b.8)

mise.hf.ti.est.r.4.hs.8=mise(t(est.hf.ti.est.r.4.hs.8),mu.hs.8)
mise.hf.ti.est.r.5.hs.8=mise(t(est.hf.ti.est.r.5.hs.8),mu.hs.8)
mise.hf.ti.est.r.6.hs.8=mise(t(est.hf.ti.est.r.6.hs.8),mu.hs.8)
mise.hf.ti.est.r.7.hs.8=mise(t(est.hf.ti.est.r.7.hs.8),mu.hs.8)

mise.hf.ti.est.r.4.bur.8=mise(t(est.hf.ti.est.r.4.bur.8),mu.bur.8)
mise.hf.ti.est.r.5.bur.8=mise(t(est.hf.ti.est.r.5.bur.8),mu.bur.8)
mise.hf.ti.est.r.6.bur.8=mise(t(est.hf.ti.est.r.6.bur.8),mu.bur.8)
mise.hf.ti.est.r.7.bur.8=mise(t(est.hf.ti.est.r.7.bur.8),mu.bur.8)

mise.hf.ti.est.r.4.s.128=mise(t(est.hf.ti.est.r.4.s.128),mu.s.128)
mise.hf.ti.est.r.5.s.128=mise(t(est.hf.ti.est.r.5.s.128),mu.s.128)
mise.hf.ti.est.r.6.s.128=mise(t(est.hf.ti.est.r.6.s.128),mu.s.128)
mise.hf.ti.est.r.7.s.128=mise(t(est.hf.ti.est.r.7.s.128),mu.s.128)

mise.hf.ti.est.r.4.ang.128=mise(t(est.hf.ti.est.r.4.ang.128),mu.ang.128)
mise.hf.ti.est.r.5.ang.128=mise(t(est.hf.ti.est.r.5.ang.128),mu.ang.128)
mise.hf.ti.est.r.6.ang.128=mise(t(est.hf.ti.est.r.6.ang.128),mu.ang.128)
mise.hf.ti.est.r.7.ang.128=mise(t(est.hf.ti.est.r.7.ang.128),mu.ang.128)

mise.hf.ti.est.r.4.cb.128=mise(t(est.hf.ti.est.r.4.cb.128),mu.cb.128)
mise.hf.ti.est.r.5.cb.128=mise(t(est.hf.ti.est.r.5.cb.128),mu.cb.128)
mise.hf.ti.est.r.6.cb.128=mise(t(est.hf.ti.est.r.6.cb.128),mu.cb.128)
mise.hf.ti.est.r.7.cb.128=mise(t(est.hf.ti.est.r.7.cb.128),mu.cb.128)

mise.hf.ti.est.r.4.b.128=mise(t(est.hf.ti.est.r.4.b.128),mu.b.128)
mise.hf.ti.est.r.5.b.128=mise(t(est.hf.ti.est.r.5.b.128),mu.b.128)
mise.hf.ti.est.r.6.b.128=mise(t(est.hf.ti.est.r.6.b.128),mu.b.128)
mise.hf.ti.est.r.7.b.128=mise(t(est.hf.ti.est.r.7.b.128),mu.b.128)

mise.hf.ti.est.r.4.hs.128=mise(t(est.hf.ti.est.r.4.hs.128),mu.hs.128)
mise.hf.ti.est.r.5.hs.128=mise(t(est.hf.ti.est.r.5.hs.128),mu.hs.128)
mise.hf.ti.est.r.6.hs.128=mise(t(est.hf.ti.est.r.6.hs.128),mu.hs.128)
mise.hf.ti.est.r.7.hs.128=mise(t(est.hf.ti.est.r.7.hs.128),mu.hs.128)

mise.hf.ti.est.r.4.bur.128=mise(t(est.hf.ti.est.r.4.bur.128),mu.bur.128)
mise.hf.ti.est.r.5.bur.128=mise(t(est.hf.ti.est.r.5.bur.128),mu.bur.128)
mise.hf.ti.est.r.6.bur.128=mise(t(est.hf.ti.est.r.6.bur.128),mu.bur.128)
mise.hf.ti.est.r.7.bur.128=mise(t(est.hf.ti.est.r.7.bur.128),mu.bur.128)

#save data
save.image("res_pois_hf_full.RData")

rm(list = setdiff(ls(), ls(pattern = "^mise")))

save.image("res_pois_hf.RData")

#note: default option for denoise.poisson for intensity (1/8,8) not run for certain test functions
#because of convergence issues

#collect mises for each test function

mise.s.8=c(mean(c(mise.hf.ti.r.4.s.8,mise.hf.ti.r.5.s.8,mise.hf.ti.r.6.s.8,mise.hf.ti.r.7.s.8)),
           mean(c(mise.hf.ti.est.r.4.s.8,mise.hf.ti.est.r.5.s.8,mise.hf.ti.est.r.6.s.8,mise.hf.ti.est.r.7.s.8)),
           mise.hf.u.r.4.s.8)

names(mise.s.8)=c("tiu_1","tiu_null","universal_null")

mise.ang.8=c(mise.hf.default.r.4.ang.8,
             mean(c(mise.hf.ti.r.4.ang.8,mise.hf.ti.r.5.ang.8,mise.hf.ti.r.6.ang.8,mise.hf.ti.r.7.ang.8)),
             mean(c(mise.hf.ti.est.r.4.ang.8,mise.hf.ti.est.r.5.ang.8,mise.hf.ti.est.r.6.ang.8,mise.hf.ti.est.r.7.ang.8)),
             mise.hf.u.r.4.ang.8)

names(mise.ang.8)=c("bt_cv_null","tiu_1","tiu_null","universal_null")

mise.cb.8=c(mise.hf.default.r.4.cb.8,
            mean(c(mise.hf.ti.r.4.cb.8,mise.hf.ti.r.5.cb.8,mise.hf.ti.r.6.cb.8,mise.hf.ti.r.7.cb.8)),
            mean(c(mise.hf.ti.est.r.4.cb.8,mise.hf.ti.est.r.5.cb.8,mise.hf.ti.est.r.6.cb.8,mise.hf.ti.est.r.7.cb.8)),
            mise.hf.u.r.4.cb.8)

names(mise.cb.8)=c("bt_cv_null","tiu_1","tiu_null","universal_null")

mise.b.8=c(mean(c(mise.hf.ti.r.4.b.8,mise.hf.ti.r.5.b.8,mise.hf.ti.r.6.b.8,mise.hf.ti.r.7.b.8)),
           mean(c(mise.hf.ti.est.r.4.b.8,mise.hf.ti.est.r.5.b.8,mise.hf.ti.est.r.6.b.8,mise.hf.ti.est.r.7.b.8)),
           mise.hf.u.r.4.b.8)

names(mise.b.8)=c("tiu_1","tiu_null","universal_null")

mise.hs.8=c(mise.hf.default.r.4.hs.8,
            mean(c(mise.hf.ti.r.4.hs.8,mise.hf.ti.r.5.hs.8,mise.hf.ti.r.6.hs.8,mise.hf.ti.r.7.hs.8)),
            mean(c(mise.hf.ti.est.r.4.hs.8,mise.hf.ti.est.r.5.hs.8,mise.hf.ti.est.r.6.hs.8,mise.hf.ti.est.r.7.hs.8)),
            mise.hf.u.r.4.hs.8)

names(mise.hs.8)=c("bt_cv_null","tiu_1","tiu_null","universal_null")

mise.bur.8=c(mean(c(mise.hf.ti.r.4.bur.8,mise.hf.ti.r.5.bur.8,mise.hf.ti.r.6.bur.8,mise.hf.ti.r.7.bur.8)),
             mean(c(mise.hf.ti.est.r.4.bur.8,mise.hf.ti.est.r.5.bur.8,mise.hf.ti.est.r.6.bur.8,mise.hf.ti.est.r.7.bur.8)),
             mise.hf.u.r.4.bur.8)

names(mise.bur.8)=c("tiu_1","tiu_null","universal_null")

mise.s.128=c(mise.hf.default.r.4.s.128,
             mean(c(mise.hf.ti.r.4.s.128,mise.hf.ti.r.5.s.128,mise.hf.ti.r.6.s.128,mise.hf.ti.r.7.s.128)),
             mean(c(mise.hf.ti.est.r.4.s.128,mise.hf.ti.est.r.5.s.128,mise.hf.ti.est.r.6.s.128,mise.hf.ti.est.r.7.s.128)),
             mise.hf.u.r.4.s.128)

names(mise.s.128)=c("bt_cv_null","tiu_1","tiu_null","universal_null")

mise.ang.128=c(mise.hf.default.r.4.ang.128,
               mean(c(mise.hf.ti.r.4.ang.128,mise.hf.ti.r.5.ang.128,mise.hf.ti.r.6.ang.128,mise.hf.ti.r.7.ang.128)),
               mean(c(mise.hf.ti.est.r.4.ang.128,mise.hf.ti.est.r.5.ang.128,mise.hf.ti.est.r.6.ang.128,mise.hf.ti.est.r.7.ang.128)),
               mise.hf.u.r.4.ang.128)

names(mise.ang.128)=c("bt_cv_null","tiu_1","tiu_null","universal_null")

mise.cb.128=c(mise.hf.default.r.4.cb.128,
              mean(c(mise.hf.ti.r.4.cb.128,mise.hf.ti.r.5.cb.128,mise.hf.ti.r.6.cb.128,mise.hf.ti.r.7.cb.128)),
              mean(c(mise.hf.ti.est.r.4.cb.128,mise.hf.ti.est.r.5.cb.128,mise.hf.ti.est.r.6.cb.128,mise.hf.ti.est.r.7.cb.128)),
              mise.hf.u.r.4.cb.128)

names(mise.cb.128)=c("bt_cv_null","tiu_1","tiu_null","universal_null")

mise.b.128=c(mise.hf.default.r.4.b.128,
             mean(c(mise.hf.ti.r.4.b.128,mise.hf.ti.r.5.b.128,mise.hf.ti.r.6.b.128,mise.hf.ti.r.7.b.128)),
             mean(c(mise.hf.ti.est.r.4.b.128,mise.hf.ti.est.r.5.b.128,mise.hf.ti.est.r.6.b.128,mise.hf.ti.est.r.7.b.128)),
             mise.hf.u.r.4.b.128)

names(mise.b.128)=c("bt_cv_null","tiu_1","tiu_null","universal_null")

mise.hs.128=c(mise.hf.default.r.4.hs.128,
              mean(c(mise.hf.ti.r.4.hs.128,mise.hf.ti.r.5.hs.128,mise.hf.ti.r.6.hs.128,mise.hf.ti.r.7.hs.128)),
              mean(c(mise.hf.ti.est.r.4.hs.128,mise.hf.ti.est.r.5.hs.128,mise.hf.ti.est.r.6.hs.128,mise.hf.ti.est.r.7.hs.128)),
              mise.hf.u.r.4.hs.128)

names(mise.hs.128)=c("bt_cv_null","tiu_1","tiu_null","universal_null")

mise.bur.128=c(mise.hf.default.r.4.bur.128,
               mean(c(mise.hf.ti.r.4.bur.128,mise.hf.ti.r.5.bur.128,mise.hf.ti.r.6.bur.128,mise.hf.ti.r.7.bur.128)),
               mean(c(mise.hf.ti.est.r.4.bur.128,mise.hf.ti.est.r.5.bur.128,mise.hf.ti.est.r.6.bur.128,mise.hf.ti.est.r.7.bur.128)),
               mise.hf.u.r.4.bur.128)

names(mise.bur.128)=c("bt_cv_null","tiu_1","tiu_null","universal_null")


#look at results
mise.s.8
mise.ang.8
mise.cb.8
mise.b.8
mise.hs.8
mise.bur.8

mise.s.128
mise.ang.128
mise.cb.128
mise.b.128
mise.hs.128
mise.bur.128



#collect mises for each test function, with mises for each primary resolution level (4,5,6,7)
mise.s.full.8=matrix(c(mise.hf.ti.r.4.s.8,mise.hf.ti.r.5.s.8,mise.hf.ti.r.6.s.8,mise.hf.ti.r.7.s.8,
                       mise.hf.ti.est.r.4.s.8,mise.hf.ti.est.r.5.s.8,mise.hf.ti.est.r.6.s.8,mise.hf.ti.est.r.7.s.8,
                       rep(mise.hf.u.r.4.s.8,4)),nrow=4)

colnames(mise.s.full.8)=c("tiu_1","tiu_null","universal_null")

mise.ang.full.8=matrix(c(rep(mise.hf.default.r.4.ang.8,4),mise.hf.ti.r.4.ang.8,mise.hf.ti.r.5.ang.8,mise.hf.ti.r.6.ang.8,mise.hf.ti.r.7.ang.8,
                         mise.hf.ti.est.r.4.ang.8,mise.hf.ti.est.r.5.ang.8,mise.hf.ti.est.r.6.ang.8,mise.hf.ti.est.r.7.ang.8,
                         rep(mise.hf.u.r.4.ang.8,4)),nrow=4)

colnames(mise.ang.full.8)=c("bt_cv_null","tiu_1","tiu_null","universal_null")

mise.cb.full.8=matrix(c(rep(mise.hf.default.r.4.cb.8,4),mise.hf.ti.r.4.cb.8,mise.hf.ti.r.5.cb.8,mise.hf.ti.r.6.cb.8,mise.hf.ti.r.7.cb.8,
                        mise.hf.ti.est.r.4.cb.8,mise.hf.ti.est.r.5.cb.8,mise.hf.ti.est.r.6.cb.8,mise.hf.ti.est.r.7.cb.8,
                        rep(mise.hf.u.r.4.cb.8,4)),nrow=4)

colnames(mise.cb.full.8)=c("bt_cv_null","tiu_1","tiu_null","universal_null")

mise.b.full.8=matrix(c(mise.hf.ti.r.4.b.8,mise.hf.ti.r.5.b.8,mise.hf.ti.r.6.b.8,mise.hf.ti.r.7.b.8,
                       mise.hf.ti.est.r.4.b.8,mise.hf.ti.est.r.5.b.8,mise.hf.ti.est.r.6.b.8,mise.hf.ti.est.r.7.b.8,
                       rep(mise.hf.u.r.4.b.8,4)),nrow=4)

colnames(mise.b.full.8)=c("tiu_1","tiu_null","universal_null")

mise.hs.full.8=matrix(c(rep(mise.hf.default.r.4.hs.8,4),mise.hf.ti.r.4.hs.8,mise.hf.ti.r.5.hs.8,mise.hf.ti.r.6.hs.8,mise.hf.ti.r.7.hs.8,
                        mise.hf.ti.est.r.4.hs.8,mise.hf.ti.est.r.5.hs.8,mise.hf.ti.est.r.6.hs.8,mise.hf.ti.est.r.7.hs.8,
                        rep(mise.hf.u.r.4.hs.8,4)),nrow=4)

colnames(mise.hs.full.8)=c("bt_cv_null","tiu_1","tiu_null","universal_null")

mise.bur.full.8=matrix(c(mise.hf.ti.r.4.bur.8,mise.hf.ti.r.5.bur.8,mise.hf.ti.r.6.bur.8,mise.hf.ti.r.7.bur.8,
                         mise.hf.ti.est.r.4.bur.8,mise.hf.ti.est.r.5.bur.8,mise.hf.ti.est.r.6.bur.8,mise.hf.ti.est.r.7.bur.8,
                         rep(mise.hf.u.r.4.bur.8,4)),nrow=4)

colnames(mise.bur.full.8)=c("tiu_1","tiu_null","universal_null")




mise.s.full.128=matrix(c(rep(mise.hf.default.r.4.s.128,4),mise.hf.ti.r.4.s.128,mise.hf.ti.r.5.s.128,mise.hf.ti.r.6.s.128,mise.hf.ti.r.7.s.128,
                         mise.hf.ti.est.r.4.s.128,mise.hf.ti.est.r.5.s.128,mise.hf.ti.est.r.6.s.128,mise.hf.ti.est.r.7.s.128,
                         rep(mise.hf.u.r.4.s.128,4)),nrow=4)

colnames(mise.s.full.128)=c("bt_cv_null","tiu_1","tiu_null","universal_null")

mise.ang.full.128=matrix(c(rep(mise.hf.default.r.4.ang.128,4),mise.hf.ti.r.4.ang.128,mise.hf.ti.r.5.ang.128,mise.hf.ti.r.6.ang.128,mise.hf.ti.r.7.ang.128,
                           mise.hf.ti.est.r.4.ang.128,mise.hf.ti.est.r.5.ang.128,mise.hf.ti.est.r.6.ang.128,mise.hf.ti.est.r.7.ang.128,
                           rep(mise.hf.u.r.4.ang.128,4)),nrow=4)

colnames(mise.ang.full.128)=c("bt_cv_null","tiu_1","tiu_null","universal_null")

mise.cb.full.128=matrix(c(rep(mise.hf.default.r.4.cb.128,4),mise.hf.ti.r.4.cb.128,mise.hf.ti.r.5.cb.128,mise.hf.ti.r.6.cb.128,mise.hf.ti.r.7.cb.128,
                          mise.hf.ti.est.r.4.cb.128,mise.hf.ti.est.r.5.cb.128,mise.hf.ti.est.r.6.cb.128,mise.hf.ti.est.r.7.cb.128,
                          rep(mise.hf.u.r.4.cb.128,4)),nrow=4)

colnames(mise.cb.full.128)=c("bt_cv_null","tiu_1","tiu_null","universal_null")

mise.b.full.128=matrix(c(rep(mise.hf.default.r.4.b.128,4),mise.hf.ti.r.4.b.128,mise.hf.ti.r.5.b.128,mise.hf.ti.r.6.b.128,mise.hf.ti.r.7.b.128,
                         mise.hf.ti.est.r.4.b.128,mise.hf.ti.est.r.5.b.128,mise.hf.ti.est.r.6.b.128,mise.hf.ti.est.r.7.b.128,
                         rep(mise.hf.u.r.4.b.128,4)),nrow=4)

colnames(mise.b.full.128)=c("bt_cv_null","tiu_1","tiu_null","universal_null")

mise.hs.full.128=matrix(c(rep(mise.hf.default.r.4.hs.128,4),mise.hf.ti.r.4.hs.128,mise.hf.ti.r.5.hs.128,mise.hf.ti.r.6.hs.128,mise.hf.ti.r.7.hs.128,
                          mise.hf.ti.est.r.4.hs.128,mise.hf.ti.est.r.5.hs.128,mise.hf.ti.est.r.6.hs.128,mise.hf.ti.est.r.7.hs.128,
                          rep(mise.hf.u.r.4.hs.128,4)),nrow=4)

colnames(mise.hs.full.128)=c("bt_cv_null","tiu_1","tiu_null","universal_null")

mise.bur.full.128=matrix(c(rep(mise.hf.default.r.4.bur.128,4),mise.hf.ti.r.4.bur.128,mise.hf.ti.r.5.bur.128,mise.hf.ti.r.6.bur.128,mise.hf.ti.r.7.bur.128,
                           mise.hf.ti.est.r.4.bur.128,mise.hf.ti.est.r.5.bur.128,mise.hf.ti.est.r.6.bur.128,mise.hf.ti.est.r.7.bur.128,
                           rep(mise.hf.u.r.4.bur.128,4)),nrow=4)

colnames(mise.bur.full.128)=c("bt_cv_null","tiu_1","tiu_null","universal_null")


#look at results
mise.s.full.8
mise.ang.full.8
mise.cb.full.8
mise.b.full.8
mise.hs.full.8
mise.bur.full.8

mise.s.full.128
mise.ang.full.128
mise.cb.full.128
mise.b.full.128
mise.hs.full.128
mise.bur.full.128
