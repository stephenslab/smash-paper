source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/bayesmooth.R"))
source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/bayesmooth_test.R"))

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
mu.t=(1+mu.s)/5
var1=rep(1,n)
var2=(0.01+4*(exp(-550*(t-0.2)^2)+exp(-200*(t-0.5)^2)+exp(-950*(t-0.8)^2)))/1.35
sigma.ini=sqrt(var2)
sigma.t=sigma.ini/mean(sigma.ini)*sd(mu.t)/rsnr^2

X=rnorm(n,mu.t,sigma.t)

mse=function(x,y) mean((x-y)^2)
l2norm=function(x) sum(x^2)
mise=function(x,y) 10000*mean(apply(x-rep(1,100)%o%y,1,l2norm)/l2norm(y))

mu.est0=bayesmooth(X,gridmult=0)
var.est0=bayesmooth(X,v.est=TRUE,gridmult=0)
mu.est1=bayesmooth.test(X,gridmult=0)
var.est1=bayesmooth.test(X,v.est=TRUE,gridmult=0)
mu.est2=bayesmooth.test(X,return.est=FALSE,gridmult=0)
var.est2=bayesmooth.test(X,return.est=FALSE,v.est=TRUE,gridmult=0)

mse(mu.est0,mu.est1)
mse(mu.est0,mu.est2$mu.est)
mse(var.est0,var.est1)
mse(var.est0,var.est2$var.est)

plot(mu.t,type='l',ylim=c(0,1))
lines(mu.est2$mu.est,col=2)
lines(mu.est2$mu.est-2*sqrt(mu.est2$mu.var),col=4)
lines(mu.est2$mu.est+2*sqrt(mu.est2$mu.var),col=4)
lines(mu.est2$mu.est-2*sigma.t,col=3)
lines(mu.est2$mu.est+2*sigma.t,col=3)



plot(sigma.t^2,type='l',ylim=c(-0.015,0.03))
lines(var.est2$var.est,col=2)
lines(var.est2$var.est-2*sqrt(var.est2$var.var),col=4)
lines(var.est2$var.est+2*sqrt(var.est2$var.var),col=4)



vtable=titable(sigma.t^2)$sumtable
y = titable(X.s)
n=1024
J=log2(n)
alpha0=matrix(0,J,n)
alpha1=matrix(0,J,n)
pi0=list(0)
pi1=list(0)
sigvec0=list(0)
sigvec1=list(0)

for(j in 0:(J-1)){
  zdat.ash=ash(y$difftable[j+2,],sqrt(vtable[j+2,]),nullcheck=FALSE, usePointMass = TRUE, localfdr = FALSE, prior="nullbiased")
  alpha0[j+1,]=zdat.ash$PosteriorMean
  gg=zdat.ash$fitted.g
  pi0[[j+1]]=gg$pi
  sigvec0[[j+1]]=gg$sd
  zdat.ash=ash(y$difftable[j+2,],sqrt(2*vtable[j+2,]),nullcheck=FALSE, usePointMass = TRUE, localfdr = FALSE, prior="nullbiased")
  alpha1[j+1,]=zdat.ash$PosteriorMean
  gg=zdat.ash$fitted.g
  pi1[[j+1]]=gg$pi
  sigvec1[[j+1]]=gg$sd
}

pi0[[1]]
pi1[[1]]
sigvec0[[1]]
sigvec1[[1]]
pi0[[2]]
pi1[[2]]
sigvec0[[2]]
sigvec1[[2]]
pi0[[3]]
pi1[[3]]
sigvec0[[3]]
sigvec1[[3]]


plot(as.vector(t(alpha0)),pch=".",cex=2)
points(as.vector(t(alpha1)),pch=".",cex=2,col=2)

plot(as.vector(t(alpha0))-as.vector(t(alpha1)),pch=".",cex=1.5)
plot(reverse.gwave(sum(X.s),alpha0),reverse.gwave(sum(X.s),alpha1),pch=".",cex=1.5)

mu0=cs.smooth(X.s,sigma.t)
mu1=cs.test(X.s,sigma.t)

plot(mu0,type='l')
lines(mu1,col=2)




###variance


n=4096
J=log2(n)
alpha0=matrix(0,J,n)
alpha1=matrix(0,J,n)
pi0=list(0)
pi1=list(0)
sigvec0=list(0)
sigvec1=list(0)

    var.est=(x-lshift(x))^2/2
    var.estl=lshift(var.est)
    var.est.ini=var.est
    vtable=titable(2/9*(var.est+var.estl)^2)$sumtable
    vvtable=1/3*vartable(var.est^2)
    vtable=vtable+vvtable
    vdtable=titable(var.est)$difftable

for(j in 0:(J-1)){
  zdat.ash=ash(vdtable[j+2,],sqrt(vtable[j+2,]),nullcheck=FALSE, usePointMass = TRUE, localfdr = FALSE, prior="nullbiased")
  alpha0[j+1,]=zdat.ash$PosteriorMean
  gg=zdat.ash$fitted.g
  pi0[[j+1]]=gg$pi
  sigvec0[[j+1]]=gg$sd
}

    var.est=x^2
    var.est1=(x-lshift(x))^2/2
    var.est1.ini=var.est1
    vtable1=titable(2/3*var.est^2)$sumtable
    vvtable1=1/3*vartable(var.est^2)
    vtable1=vtable1+vvtable1
    vdtable1=titable(var.est1)$difftable


for(j in 0:(J-1)){
  zdat.ash=ash(vdtable1[j+2,],sqrt(vtable1[j+2,]),nullcheck=FALSE, usePointMass = TRUE, localfdr = FALSE, prior="nullbiased")
  alpha1[j+1,]=zdat.ash$PosteriorMean
  gg=zdat.ash$fitted.g
  pi1[[j+1]]=gg$pi
  sigvec1[[j+1]]=gg$sd
}

pi0[[1]]
pi1[[1]]
sigvec0[[1]]
sigvec1[[1]]
pi0[[2]]
pi1[[2]]
sigvec0[[2]]
sigvec1[[2]]
pi0[[3]]
pi1[[3]]
sigvec0[[3]]
sigvec1[[3]]


plot(as.vector(t(alpha0)),pch=".",cex=2)
points(as.vector(t(alpha1)),pch=".",cex=2,col=2)

