spike.f=function(x) (0.75*exp(-500*(x-0.23)^2)+1.5*exp(-2000*(x-0.33)^2)+3*exp(-8000*(x-0.47)^2)+2.25*exp(-16000*(x-0.69)^2)+0.5*exp(-32000*(x-0.83)^2))
n=1024
t=1:n/n
mu.s=spike.f(1:n/n)
#sigma.t=seq(2,10,length.out=n)
sigma.t=rep(0.17,n)

mu.s0=2
set.seed(425)
X.s=rnorm(n,mu.s0*(mu.s),sigma.t)

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


mse=function(x,y) mean((x-y)^2)
l2norm=function(x) sum(x^2)
mise=function(x,y) 10000*mean(apply(x-rep(1,100)%o%y,1,l2norm)/l2norm(y))




source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/bayesmooth_mc.R"))
source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/bayesmooth.R"))
source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/bayesmooth_new.R"))


#source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/cs.smooth_test.R"))
#source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/var.smooth.R"))


la10.sure0=function (x, min.level = 0) 
{
    x.w <- wd(x)
    x.w.t <- threshold(x.w, levels = (min.level):(x.w$nlevels - 1))
    x.w.t.r <- wr(x.w.t)
    return(x.w.t.r)
}

la10.sure3=function (x) 
{
    x.w <- wd(x)
    x.w.t <- threshold(x.w)
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





#####

hf.la10.ti2=function (x, filter.number = 10, family = "DaubLeAsymm", min.level = 2, noise.level = 1) 
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

sigma=sqrt(2/(3*(n-2))*sum((1/2*X.s[1:(n-2)]-X.s[2:(n-1)]+1/2*X.s[3:n])^2))

system.time(mu.est.0<-cs.smooth(X.s,1))
system.time(mu.est.hf0<-hf.la10.ti2(X.s,noise.level=sigma))
mu.est.hf2=hf.la10.ti2(X.s,noise.level=sigma)
mu.est.hf3=hf.la10.ti3(X.s,noise.level=sigma)
mu.est.hf4=hf.la10.ti4(X.s,noise.level=sigma)
mu.est.hf5=hf.la10.ti5(X.s,noise.level=sigma)
mu.est.hf6=hf.la10.ti6(X.s,noise.level=sigma)
mu.est.hf7=hf.la10.ti7(X.s,noise.level=sigma)

mse(mu.s0*(mu.s),mu.est.0)
mse(mu.s0*(mu.s),mu.est.hf0)
mse(mu.s0*(mu.s),mu.est.hf2)
mse(mu.s0*(mu.s),mu.est.hf3)
mse(mu.s0*(mu.s),mu.est.hf4)
mse(mu.s0*(mu.s),mu.est.hf5)
mse(mu.s0*(mu.s),mu.est.hf6)
mse(mu.s0*(mu.s),mu.est.hf7)

plot(mu.s0*(mu.s),type='l')
lines(mu.est.hf3,col=2)
lines(mu.est.0,col=4)


sigma=sqrt(2/(3*(n-2))*sum((1/2*X.b[1:(n-2)]-X.b[2:(n-1)]+1/2*X.b[3:n])^2))

system.time(mu.est.0<-cs.smooth(X.b,sigma.t))
mu.est.hf0=hf.la10.ti0(X.b,noise.level=sigma)
mu.est.hf2=hf.la10.ti2(X.b,noise.level=sigma)
mu.est.hf3=hf.la10.ti3(X.b,noise.level=sigma)
mu.est.hf4=hf.la10.ti4(X.b,noise.level=sigma)
mu.est.hf5=hf.la10.ti5(X.b,noise.level=sigma)
mu.est.hf6=hf.la10.ti6(X.b,noise.level=sigma)
mu.est.hf7=hf.la10.ti7(X.b,noise.level=sigma)

mse(mu.b0*(mu.b),mu.est.0)
mse(mu.b0*(mu.b),mu.est.hf0)
mse(mu.b0*(mu.b),mu.est.hf2)
mse(mu.b0*(mu.b),mu.est.hf3)
mse(mu.b0*(mu.b),mu.est.hf4)
mse(mu.b0*(mu.b),mu.est.hf5)
mse(mu.b0*(mu.b),mu.est.hf6)
mse(mu.b0*(mu.b),mu.est.hf7)

plot(mu.b0*(mu.b),type='l')
lines(mu.est.hf3,col=2)
lines(mu.est.0,col=4)



n=1024
set.seed(1025)
#mu.t=(1+mu.s)/5
#mu.t=(1+mu.b)/5
#mu.t=(1+mu.blk)/5
rsnr=8
sigma.t=rep(sqrt(l2norm(mu.t)/rsnr^2/n),n)
X.s=matrix(rnorm(100*n,mu.t,sigma.t),nrow=100,byrow=TRUE)
mu.est.ash=matrix(0,nrow=100,ncol=n)
mu.est.ashe=matrix(0,nrow=100,ncol=n)
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

for(i in 1:100){
  sigma=sqrt(2/(3*(n-2))*sum((1/2*X.s[i,1:(n-2)]-X.s[i,2:(n-1)]+1/2*X.s[i,3:n])^2))
  mu.est.ash[i,]=bayesmooth(X.s[i,],sigma)
  mu.est.ashe[i,]=bayesmooth(X.s[i,])
  mu.est.ti0[i,]=la10.ti0(X.s[i,],noise.level=sigma)
  mu.est.ti3[i,]=la10.ti3(X.s[i,],noise.level=sigma)
  mu.est.sure0[i,]=la10.sure0(X.s[i,])
  mu.est.sure3[i,]=la10.sure3(X.s[i,])
  mu.est.ti0u[i,]=la10.ti0u(X.s[i,],noise.level=sigma)
  mu.est.ti3u[i,]=la10.ti3u(X.s[i,],noise.level=sigma)
  mu.est.ti0.mod[i,]=la10.ti0.mod(X.s[i,])
  mu.est.ti3.mod[i,]=la10.ti3.mod(X.s[i,])
  print(i)
}

(mise(mu.est.ash,mu.t))
(mise(mu.est.ashe,mu.t))
(mise(mu.est.ti0u,mu.t))
(mise(mu.est.ti3u,mu.t))
(mise(mu.est.ti0,mu.t))
(mise(mu.est.ti3,mu.t))
(mise(mu.est.ti0.mod,mu.t))
(mise(mu.est.ti3.mod,mu.t))
(mise(mu.est.sure0,mu.t))
(mise(mu.est.sure3,mu.t))


plot(mu.t,type='l')
lines(mu.est.sure0[1,],col=4)

lines(mu.est.ti3[1,],col=4)
lines(mu.est.ashe[1,],col=2)


multimix, 2
1.267606, 27.66753
1.602286, 308.5797
1.881001, 882.6439

1.841075, 27.78215
2.076758, 304.3344
2.626197, 875.991
multimix, var
1.170439, 31.69251
1.775539, 321.9961
2.026061, 909.98

1.547167, 31.54611
2.340342, 318.3228
2.784097, 904.8981
multimix, known
1.030359, 32.56845
1.76683, 342.1839
2.453837, 976.293

1.400861, 32.10729
2.252333, 339.1009
3.156232, 972.5803
2mix, 2
1.230328, 27.34378
1.604813, 311.055
1.99395, 890.7364

1.744702, 27.51118
1.933743, 306.8624
2.676913, 883.6396
2mix, var
0.942461, 31.36283
1.494225, 324.1331
1.573484, 915.9358

1.329028, 30.82365
2.143026, 321.1913
2.359406, 911.5031
2mix, known
0.9895533, 32.12692
1.601807, 343.3145
2.216103, 976.8718
 
1.294851, 31.29577
2.058961, 340.4559
2.933778, 974.3178


tt=0
for(i in 1:10000){
xx=rnorm(1,2,2)
xxx=rnorm(1,2,2)
xxxx=rnorm(1,2,2)

yy=rnorm(1,2,2)
yyy=rnorm(1,2,2)
yyyy=rnorm(1,2,2)

tt[i]=(xx+xxx-yy-yyy)/2/sqrt(((xx-xxx)^2+(xxx-yy)^2+(yy-yyy)^2+(yyy-yyyy)^2)/8)
}
sd(tt)
hist(tt,breaks=20,xlim=c(-10,10),prob=TRUE)
points(seq(-10,10,0.1),dnorm(seq(-10,10,0.1),0,1.2),type='l')



hist((xx-yy)/2/sqrt(((xx-yy)^2+(yy-yyy)^2)/8),breaks=20,xlim=c(-10,10),prob=TRUE)
points(seq(-10,10,0.1),dnorm(seq(-10,10,0.1),0,1.2),type='l')



tt=0
for(i in 1:10000){
n=2^0
xx=rnorm(n,2,2)
yy=rnorm(n,2,2)
xy=c(xx,yy)
tt[i]=(sum(xx)-sum(yy))/sqrt((sum((xy[1:(2*n-1)]-lshift(xy)[1:(2*n-1)])^2)+(yy[n]-rnorm(1,2,2))^2)/2)
}
sd(tt)
hist(tt,breaks=50,prob=TRUE)
points(seq(-10,10,0.1),dnorm(seq(-10,10,0.1),0,1),type='l')


tt=0
for(i in 1:10000){
n=2^0
tsig=10
xx=rnorm(n,2,tsig)
yy=rnorm(n,2,tsig)
sd=sqrt(tsig^2*2*n)
xy=c(xx,yy)
tt[i]=sqrt((sum((xy[1:(2*n-1)]-lshift(xy)[1:(2*n-1)])^2)+(yy[n]-rnorm(1,2,tsig))^2)/2)
}
mean(tt)
sd^2/mean(tt)^2
hist(tt,breaks=50,prob=TRUE)
points(seq(-10,10,0.1),dnorm(seq(-10,10,0.1),0,1),type='l')



###
tt=matrix(0,12,10000)
ts=matrix(0,12,10000)
for(i in 1:10000){
n=2^12
tsig=2
xx=rnorm(n,2,tsig)
vv1=(xx-lshift(xx))^2/2
vt1=titable(vv1)$sumtable
tt[,i]=sqrt(vt1[2:13,1])
ts[,i]=(titable(xx)$difftable)[2:13,1]
if(i%%1000==0){print(i)}
}
sd(ts[1,])^2/mean(tt[1,])^2
sd(ts[2,])^2/mean(tt[2,])^2
sd(ts[3,])^2/mean(tt[3,])^2
sd(ts[4,])^2/mean(tt[4,])^2
sd(ts[5,])^2/mean(tt[5,])^2
sd(ts[6,])^2/mean(tt[6,])^2
sd(ts[7,])^2/mean(tt[7,])^2
sd(ts[8,])^2/mean(tt[8,])^2
sd(ts[9,])^2/mean(tt[9,])^2
sd(ts[10,])^2/mean(tt[10,])^2
sd(ts[11,])^2/mean(tt[11,])^2
sd(ts[12,])^2/mean(tt[12,])^2



###
tt=matrix(0,12,10000)
ts=matrix(0,12,10000)
for(i in 1:10000){
n=2^12
tsig=2
xx=rnorm(n,2,tsig)
vv2=(rshift(xx)-xx)^2/2
vt2=titable(vv2)$sumtable
tt[,i]=sqrt(vt2[2:13,1])
ts[,i]=(titable(xx)$difftable)[2:13,1]
if(i%%1000==0){print(i)}
}
sd(ts[1,])^2/mean(tt[1,])^2
sd(ts[2,])^2/mean(tt[2,])^2
sd(ts[3,])^2/mean(tt[3,])^2
sd(ts[4,])^2/mean(tt[4,])^2
sd(ts[5,])^2/mean(tt[5,])^2
sd(ts[6,])^2/mean(tt[6,])^2
sd(ts[7,])^2/mean(tt[7,])^2
sd(ts[8,])^2/mean(tt[8,])^2
sd(ts[9,])^2/mean(tt[9,])^2
sd(ts[10,])^2/mean(tt[10,])^2
sd(ts[11,])^2/mean(tt[11,])^2
sd(ts[12,])^2/mean(tt[12,])^2

###
J=8
tt=matrix(0,J,10000)
ts=matrix(0,J,10000)
for(i in 1:10000){
n=2^J
tsig=20
xx=rnorm(n,2,tsig)
vv1=(xx-lshift(xx))^2/2
vv2=(rshift(xx)-xx)^2/2
vt1=titable(vv1)$sumtable
vt2=titable(vv2)$sumtable
vt=(vt1+vt2)/2
tt[,i]=sqrt(vt[2:(J+1),1])
ts[,i]=(titable(xx)$difftable)[2:(J+1),1]
if(i%%1000==0){print(i)}
}
sd(ts[1,])^2/mean(tt[1,])^2
sd(ts[2,])^2/mean(tt[2,])^2
sd(ts[3,])^2/mean(tt[3,])^2
sd(ts[4,])^2/mean(tt[4,])^2
sd(ts[5,])^2/mean(tt[5,])^2
sd(ts[6,])^2/mean(tt[6,])^2
sd(ts[7,])^2/mean(tt[7,])^2
sd(ts[8,])^2/mean(tt[8,])^2
sd(ts[9,])^2/mean(tt[9,])^2
sd(ts[10,])^2/mean(tt[10,])^2
sd(ts[11,])^2/mean(tt[11,])^2
sd(ts[12,])^2/mean(tt[12,])^2

sd(ts[1,]/tt[1,])
sd(ts[2,]/tt[2,])
sd(ts[3,]/tt[3,])
sd(ts[4,]/tt[4,])
sd(ts[5,]/tt[5,])
sd(ts[6,]/tt[6,])
sd(ts[7,]/tt[7,])
sd(ts[8,]/tt[8,])
sd(ts[9,]/tt[9,])
sd(ts[10,]/tt[10,])
sd(ts[11,]/tt[11,])
sd(ts[12,]/tt[12,])

hist(ts[1,]/tt[1,],xlim=c(-8,8),breaks=30,prob=TRUE)
lines(seq(-8,8,0.01),dnorm(seq(-8,8,0.01),0,1),col=2)
hist(ts[2,]/tt[2,],xlim=c(-8,8),breaks=30,prob=TRUE)
lines(seq(-8,8,0.01),dnorm(seq(-8,8,0.01),0,1),col=2)
hist(ts[3,]/tt[3,],xlim=c(-8,8),breaks=30,prob=TRUE)
lines(seq(-8,8,0.01),dnorm(seq(-8,8,0.01),0,1),col=2)
hist(ts[4,]/tt[4,],xlim=c(-8,8),breaks=30,prob=TRUE)
lines(seq(-8,8,0.01),dnorm(seq(-8,8,0.01),0,1),col=2)
hist(ts[5,]/tt[5,],xlim=c(-8,8),breaks=30,prob=TRUE)
lines(seq(-8,8,0.01),dnorm(seq(-8,8,0.01),0,1),col=2)
hist(ts[6,]/tt[6,],xlim=c(-8,8),breaks=30,prob=TRUE)
lines(seq(-8,8,0.01),dnorm(seq(-8,8,0.01),0,1),col=2)
hist(ts[7,]/tt[7,],xlim=c(-8,8),breaks=30,prob=TRUE)
lines(seq(-8,8,0.01),dnorm(seq(-8,8,0.01),0,1),col=2)
hist(ts[8,]/tt[8,],xlim=c(-8,8),breaks=30,prob=TRUE)
lines(seq(-8,8,0.01),dnorm(seq(-8,8,0.01),0,1),col=2)
lines(seq(-8,8,0.01),dnorm(seq(-8,8,0.01),0,1),col=2)

#1.3,1.17,1.08,1.04,1.02,1.01,1,1,1,1,1,1




#2
tt=0
ts=0
for(i in 1:10000){
xx=rnorm(2^3,0,2)
xv=(xx-lshift(xx))^2/2
xvl=lshift(xv)
vt=titable(1/4*(xv+xvl)^2)$sumtable
vvt=1/3*vartable(xv^2)
vt=vt+vvt
vd=titable(xv)$difftable
ts[i]=vd[2,1]
tt[i]=sqrt(vt[2,1])
if(i%%1000==0){print(i)}
}
sd(ts)
mean(tt)
hist(tt,breaks=50,xlim=c(-10,10),prob=TRUE)


#4
tt=0
ts=0
for(i in 1:10000){
xx=rnorm(2^4,0,2)
xv=(xx-lshift(xx))^2/2
xvl=lshift(xv)
vt=titable(1/4*(xv+xvl)^2)$sumtable
vvt=1/3*vartable(xv^2)
vt=vt+vvt
vd=titable(xv)$difftable
ts[i]=vd[3,1]
tt[i]=sqrt(vt[3,1])
if(i%%1000==0){print(i)}
}
sd(ts)
mean(tt)
hist(tt,breaks=50,xlim=c(-10,10),prob=TRUE)


#8
tt=0
ts=0
for(i in 1:10000){
xx=rnorm(2^5,0,2)
xv=(xx-lshift(xx))^2/2
xvl=lshift(xv)
vt=titable(1/4*(xv+xvl)^2)$sumtable
vvt=1/3*vartable(xv^2)
vt=vt+vvt
vd=titable(xv)$difftable
ts[i]=vd[4,1]
tt[i]=sqrt(vt[4,1])
if(i%%1000==0){print(i)}
}
sd(ts)
mean(tt)
sd(ts)^2/mean(tt)^2
mean(tt)hist(tt,breaks=50,xlim=c(-10,10),prob=TRUE)




#16
tt=0
ts=0
for(i in 1:10000){
xx=rnorm(2^6,0,2)
xv=(xx-lshift(xx))^2/2
xvl=lshift(xv)
vt=titable(1/4*(xv+xvl)^2)$sumtable
vvt=1/3*vartable(xv^2)
vt=vt+vvt
vd=titable(xv)$difftable
ts[i]=vd[5,1]
tt[i]=sqrt(vt[5,1])
if(i%%1000==0){print(i)}
}
sd(ts)
mean(tt)
sd(ts)^2/mean(tt)^2
hist(tt,breaks=50,xlim=c(-5,5),prob=TRUE)


#32
tt=0
ts=0
for(i in 1:10000){
xx=rnorm(2^7,0,2)
xv=(xx-lshift(xx))^2/2
xvl=lshift(xv)
vt=titable(1/4*(xv+xvl)^2)$sumtable
vvt=1/3*vartable(xv^2)
vt=vt+vvt
vd=titable(xv)$difftable
ts[i]=vd[6,1]
tt[i]=sqrt(vt[6,1])
if(i%%1000==0){print(i)}
}
sd(ts)
mean(tt)
sd(ts)^2/mean(tt)^2
hist(tt,breaks=50,xlim=c(-5,5),prob=TRUE)




###m2
tsig=20
J=8
tt=matrix(0,J,10000)
ts=matrix(0,J,10000)
for(i in 1:10000){
xx=rnorm(2^J,0,tsig)
xv=(xx-lshift(xx))^2/2
vt=titable(2/3*xv^2)$sumtable
vvt=1/3*vartable(xv^2)
vt=vt+vvt
vd=titable(xv)$difftable
ts[,i]=vd[2:(J+1),1]
tt[,i]=sqrt(vt[2:(J+1),1])
if(i%%1000==0){print(i)}
}

sd(ts[1,])^2/mean(tt[1,])^2
sd(ts[2,])^2/mean(tt[2,])^2
sd(ts[3,])^2/mean(tt[3,])^2
sd(ts[4,])^2/mean(tt[4,])^2
sd(ts[5,])^2/mean(tt[5,])^2
sd(ts[6,])^2/mean(tt[6,])^2
sd(ts[7,])^2/mean(tt[7,])^2
sd(ts[8,])^2/mean(tt[8,])^2
sd(ts[9,])^2/mean(tt[9,])^2
sd(ts[10,])^2/mean(tt[10,])^2
sd(ts[11,])^2/mean(tt[11,])^2
sd(ts[12,])^2/mean(tt[12,])^2


sd(ts[1,]/tt[1,])
sd(ts[2,]/tt[2,])
sd(ts[3,]/tt[3,])
sd(ts[4,]/tt[4,])
sd(ts[5,]/tt[5,])
sd(ts[6,]/tt[6,])
sd(ts[7,]/tt[7,])
sd(ts[8,]/tt[8,])
sd(ts[9,]/tt[9,])
sd(ts[10,]/tt[10,])
sd(ts[11,]/tt[11,])
sd(ts[12,]/tt[12,])

hist(ts[1,]/tt[1,],xlim=c(-8,8),breaks=30,prob=TRUE)
lines(seq(-8,8,0.01),dnorm(seq(-8,8,0.01),0,1),col=2)
hist(ts[2,]/tt[2,],xlim=c(-8,8),breaks=30,prob=TRUE)
lines(seq(-8,8,0.01),dnorm(seq(-8,8,0.01),0,1),col=2)
hist(ts[3,]/tt[3,],xlim=c(-8,8),breaks=30,prob=TRUE)
lines(seq(-8,8,0.01),dnorm(seq(-8,8,0.01),0,1),col=2)
hist(ts[4,]/tt[4,],xlim=c(-8,8),breaks=30,prob=TRUE)
lines(seq(-8,8,0.01),dnorm(seq(-8,8,0.01),0,1),col=2)
hist(ts[5,]/tt[5,],xlim=c(-8,8),breaks=30,prob=TRUE)
lines(seq(-8,8,0.01),dnorm(seq(-8,8,0.01),0,1),col=2)
hist(ts[6,]/tt[6,],xlim=c(-8,8),breaks=30,prob=TRUE)
lines(seq(-8,8,0.01),dnorm(seq(-8,8,0.01),0,1),col=2)
hist(ts[7,]/tt[7,],xlim=c(-8,8),breaks=30,prob=TRUE)
lines(seq(-8,8,0.01),dnorm(seq(-8,8,0.01),0,1),col=2)
hist(ts[8,]/tt[8,],xlim=c(-8,8),breaks=30,prob=TRUE)
lines(seq(-8,8,0.01),dnorm(seq(-8,8,0.01),0,1),col=2)
lines(seq(-8,8,0.01),dnorm(seq(-8,8,0.01),0,1),col=2)

#2.17,1.71,1.39,1.21,1.12,1.09,1.05,1.02,1.01,0.99,1.01,0.98



###m3
J=2
tt=matrix(0,J,10000)
ts=matrix(0,J,10000)
for(i in 1:10000){
xx=rnorm(2^J,0,2)
xv=(xx-lshift(xx))^2/2
xvl=lshift(xv)
vt=titable(2/9*(xv+xvl)^2)$sumtable
vvt=1/3*vartable(xv^2)
vt=vt+vvt
vd=titable(xv)$difftable
ts[,i]=vd[2:(J+1),1]
tt[,i]=sqrt(vt[2:(J+1),1])
if(i%%1000==0){print(i)}
}

sd(ts[1,])^2/mean(tt[1,])^2
sd(ts[2,])^2/mean(tt[2,])^2
sd(ts[3,])^2/mean(tt[3,])^2
sd(ts[4,])^2/mean(tt[4,])^2
sd(ts[5,])^2/mean(tt[5,])^2
sd(ts[6,])^2/mean(tt[6,])^2
sd(ts[7,])^2/mean(tt[7,])^2
sd(ts[8,])^2/mean(tt[8,])^2
sd(ts[9,])^2/mean(tt[9,])^2
sd(ts[10,])^2/mean(tt[10,])^2
sd(ts[11,])^2/mean(tt[11,])^2
sd(ts[12,])^2/mean(tt[12,])^2

#1.85,1.68,1.39,1.22,1.12,1.09,1.05,1.02,1.01,0.99,1.00,0.98



###m4
J=12
tt=matrix(0,J,10000)
ts=matrix(0,J,10000)
for(i in 1:10000){
xx=rnorm(2^J,0,2)
xv=xx^2
vt=cxxtitable(2/3*(xv)^2)$sumtable
vd=cxxtitable(xv)$difftable
ts[,i]=vd[2:(J+1),1]
tt[,i]=sqrt(vt[2:(J+1),1])
if(i%%1000==0){print(i)}
}

sd(ts[1,])^2/mean(tt[1,])^2
sd(ts[2,])^2/mean(tt[2,])^2
sd(ts[3,])^2/mean(tt[3,])^2
sd(ts[4,])^2/mean(tt[4,])^2
sd(ts[5,])^2/mean(tt[5,])^2
sd(ts[6,])^2/mean(tt[6,])^2
sd(ts[7,])^2/mean(tt[7,])^2
sd(ts[8,])^2/mean(tt[8,])^2
sd(ts[9,])^2/mean(tt[9,])^2
sd(ts[10,])^2/mean(tt[10,])^2
sd(ts[11,])^2/mean(tt[11,])^2
sd(ts[12,])^2/mean(tt[12,])^2

#2.03,1.50,1.29,1.14,1.10,1.07,1.01,1,0.98,1,0.99,1.00

###m4a
J=12
tt=matrix(0,J,10000)
ts=matrix(0,J,10000)
for(i in 1:10000){
xx=rnorm(2^J,0,2)
xv=xx^2
xvv=2/3*(xv)^2
vt=wd.var(xvv,filter.number=1,family="DaubExPhase",type="station")
vd=wd(xvv,filter.number=1,family="DaubExPhase",type="station")
for(j in 0:(J-1)){
  ts[j+1,i]=accessD(vd,j)[1]
  tt[j+1,i]=accessC(vt,j)[1]
}
if(i%%1000==0){print(i)}
}

###m5
J=2
tt=matrix(0,J,10000)
ts=matrix(0,J,10000)
for(i in 1:10000){
xx=rnorm(2^J,0,2)
xv=xx^2
xv1=(xx-lshift(xx))^2/2
xv2=(rshift(xx)-xx)^2/2
vt=titable(2/3*(xv1)^2)$sumtable
vd=titable(xv)$difftable
ts[,i]=vd[2:(J+1),1]
tt[,i]=sqrt(vt[2:(J+1),1])
if(i%%1000==0){print(i)}
}

sd(ts[1,])^2/mean(tt[1,])^2
sd(ts[2,])^2/mean(tt[2,])^2
sd(ts[3,])^2/mean(tt[3,])^2
sd(ts[4,])^2/mean(tt[4,])^2
sd(ts[5,])^2/mean(tt[5,])^2
sd(ts[6,])^2/mean(tt[6,])^2
sd(ts[7,])^2/mean(tt[7,])^2
sd(ts[8,])^2/mean(tt[8,])^2
sd(ts[9,])^2/mean(tt[9,])^2
sd(ts[10,])^2/mean(tt[10,])^2
sd(ts[11,])^2/mean(tt[11,])^2
sd(ts[12,])^2/mean(tt[12,])^2





###m6
J=2
tt=matrix(0,J,10000)
ts=matrix(0,J,10000)
for(i in 1:10000){
xx=rnorm(2^J,0,2)
xv=xx^2
xv1=(xx-lshift(xx))^2/2
vt=titable(2/3*(xv)^2)$sumtable
vvt=1/3*vartable(xv^2)
vt=vt+vvt
vd=titable(xv1)$difftable
ts[,i]=vd[2:(J+1),1]
tt[,i]=sqrt(vt[2:(J+1),1])
if(i%%1000==0){print(i)}
}

sd(ts[1,])^2/mean(tt[1,])^2
sd(ts[2,])^2/mean(tt[2,])^2
sd(ts[3,])^2/mean(tt[3,])^2
sd(ts[4,])^2/mean(tt[4,])^2
sd(ts[5,])^2/mean(tt[5,])^2
sd(ts[6,])^2/mean(tt[6,])^2
sd(ts[7,])^2/mean(tt[7,])^2
sd(ts[8,])^2/mean(tt[8,])^2
sd(ts[9,])^2/mean(tt[9,])^2
sd(ts[10,])^2/mean(tt[10,])^2
sd(ts[11,])^2/mean(tt[11,])^2
sd(ts[12,])^2/mean(tt[12,])^2




tt=0
for(i in 1:1000){
pp=rpois(2^10,(mu.s+0.01))
yy=ParentTItable(pp)$parent
gg=glm.approx(yy)
tt[i]=gg[1,1]/gg[2,1]
}
sd(tt)


tt=0
vvss=0
vv=0
for(i in 1:10000){
ss=rbinom(1,100,0.9)
pp=(ss)/(100)
yy=logit(pp)
vvss[i]=vss(100,ss,100-ss)
vv[i]=var(logit(rbeta(10000,ss+1,100-ss+1)))
}


sd(logit(rbeta(10000,10,92)))
plot(seq(0,1,0.01),dbeta(seq(0,1,0.01),10,92))