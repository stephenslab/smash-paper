source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/bayesmooth.R"))
source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/multiseq_gaus.R"))
source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/multi_lmm.R"))

n=1024
t=1:n/n

spike.f=function(x) (0.75*exp(-500*(x-0.23)^2)+1.5*exp(-2000*(x-0.33)^2)+3*exp(-8000*(x-0.47)^2)+2.25*exp(-16000*(x-0.69)^2)+0.5*exp(-32000*(x-0.83)^2))
spikem.f=function(x) (0.75*exp(-500*(x-0.23)^2)+2.5*exp(-2000*(x-0.33)^2)+3*exp(-8000*(x-0.47)^2)+1.5*exp(-16000*(x-0.69)^2)+1*exp(-32000*(x-0.83)^2)+0.8*exp(-12000*(x-0.92)^2))

dop.f=function(x) sqrt(x*(1-x))*sin((2*pi*1.05)/(x+0.05))
mu.dop=dop.f(t)
mu.dop=3/(max(mu.dop)-min(mu.dop))*(mu.dop-min(mu.dop))
mu.dop.var=10*dop.f(t)
mu.dop.var=mu.dop.var-min(mu.dop.var)


pos = c(.1, .13, .15, .23, .25, .40, .44, .65, .76, .78, .81)
hgt = 2.97/5*c(4, 5, 3, 4, 5, 4.2, 2.1, 4.3, 3.1, 5.1, 4.2)
wth = c(.005, .005, .006, .01, .01, .03, .01, .01, .005, .008, .005)
mu.b = rep(0,n)
for(j in 1:length(pos)){
  mu.b = mu.b + hgt[j]/(( 1 + (abs(t - pos[j])/wth[j]))^4)
}


pos=c(.1, .13, .15, .23, .25, .40, .44, .65, .76, .78, .81)
hgt = 2.88/5*c(4, (-5), 3, (-4), 5, (-4.2), 2.1, 4.3, (-3.1), 2.1, (-4.2))
mu.blk = rep(0,n)
for(j in 1:length(pos)){
  mu.blk = mu.blk + (1 + sign(t-pos[j]))*(hgt[j]/2)
}

mu.sp=spike.f(t)
mu.spm=spikem.f(t)


var1=rep(1,n)
var2=(0.0001+4*(exp(-550*(t-0.2)^2)+exp(-200*(t-0.5)^2)+exp(-950*(t-0.8)^2)))
var3=(0.00001+2*mu.dop.var)
var4=0.00001+mu.b
#var4=0.00001+1*(mu.blk-min(mu.blk))/max(mu.blk)


sigma.ini.v1=sqrt(var1)
sigma.ini.v2=sqrt(var2)
sigma.ini.v3=sqrt(var3)
sigma.ini.v4=sqrt(var4)


rsnr=sqrt(3)
sigma.0.v1=sigma.ini.v1/mean(sigma.ini.v1)*sd(mu.sp)/rsnr^2
sigma.0.v2=sigma.ini.v2/mean(sigma.ini.v2)*sd(mu.sp)/rsnr^2
sigma.0.v3=sigma.ini.v3/mean(sigma.ini.v3)*sd(mu.sp)/rsnr^2
sigma.0.v4=sigma.ini.v4/mean(sigma.ini.v4)*sd(mu.sp)/rsnr^2

sigma.1.v1=sigma.ini.v1/mean(sigma.ini.v1)*sd(mu.spm)/rsnr^2
sigma.1.v2=sigma.ini.v2/mean(sigma.ini.v2)*sd(mu.spm)/rsnr^2
sigma.1.v3=sigma.ini.v3/mean(sigma.ini.v3)*sd(mu.spm)/rsnr^2
sigma.1.v4=sigma.ini.v4/mean(sigma.ini.v4)*sd(mu.spm)/rsnr^2


set.seed(1112)
#X.1=matrix(rnorm(1*n,mu.sp,rep(sigma.0.v1,1)),nrow=1,byrow=TRUE)
X.1=matrix(rnorm(10*n,mu.sp,c(rep(sigma.0.v1,2),rep(sigma.0.v2,4),rep(sigma.0.v3,4))),nrow=10,byrow=TRUE)
#X.1=matrix(rnorm(10*n,mu.sp,rep(sigma.0.v1,10)),nrow=10,byrow=TRUE)
#X.1=matrix(rnorm(20*n,mu.sp,c(rep(sigma.0.v1,6),rep(sigma.0.v2,7),rep(sigma.0.v3,7))),nrow=20,byrow=TRUE)
#X.1=matrix(rnorm(30*n,mu.sp,c(rep(sigma.0.v1,9),rep(sigma.0.v2,7),rep(sigma.0.v3,14))),nrow=30,byrow=TRUE)
#X.1=matrix(rnorm(50*n,mu.sp,c(rep(sigma.0.v1,16),rep(sigma.0.v2,17),rep(sigma.0.v3,17))),nrow=50,byrow=TRUE)
#X.1=matrix(rnorm(80*n,mu.sp,c(rep(sigma.0.v1,18),rep(sigma.0.v2,33),rep(sigma.0.v3,29))),nrow=80,byrow=TRUE)
#X.1=matrix(rnorm(100*n,mu.sp,c(rep(sigma.0.v1,26),rep(sigma.0.v2,47),rep(sigma.0.v3,27))),nrow=100,byrow=TRUE)
set.seed(1112)
#X.2=matrix(rnorm(1*n,mu.spm,rep(sigma.1.v1,1)),nrow=1,byrow=TRUE)
X.2=matrix(rnorm(10*n,mu.spm,c(rep(sigma.1.v2,3),rep(sigma.1.v3,2),rep(sigma.1.v4,5))),nrow=10,byrow=TRUE)
#X.2=matrix(rnorm(10*n,mu.spm,rep(sigma.1.v1,10)),nrow=10,byrow=TRUE)
#X.2=matrix(rnorm(20*n,mu.spm,c(rep(sigma.1.v1,5),rep(sigma.1.v3,7),rep(sigma.1.v4,8))),nrow=20,byrow=TRUE)
#X.2=matrix(rnorm(30*n,mu.spm,c(rep(sigma.1.v1,6),rep(sigma.1.v3,10),rep(sigma.1.v4,14))),nrow=30,byrow=TRUE)
#X.2=matrix(rnorm(50*n,mu.spm,c(rep(sigma.1.v1,12),rep(sigma.1.v3,18),rep(sigma.1.v4,20))),nrow=50,byrow=TRUE)
#X.2=matrix(rnorm(80*n,mu.spm,c(rep(sigma.1.v1,16),rep(sigma.1.v3,29),rep(sigma.1.v4,35))),nrow=80,byrow=TRUE)
#X.2=matrix(rnorm(100*n,mu.spm,c(rep(sigma.1.v1,18),rep(sigma.1.v3,38),rep(sigma.1.v4,44))),nrow=100,byrow=TRUE)

X=rbind(X.1,X.2)

#g=c(rep(0,1),rep(1,1))
g=c(rep(0,10),rep(1,10))
#g=c(rep(0,20),rep(1,20))
#g=c(rep(0,30),rep(1,30))
#g=c(rep(0,50),rep(1,50))
#g=c(rep(0,80),rep(1,80))
#g=c(rep(0,100),rep(1,100))

mse=function(x,y) mean((x-y)^2)

write.table(X,"../simulation_g/sim_data.txt",col.names=FALSE,row.names=FALSE)

system.time(res1<-multiseq.gaus(X,g,gridmult=2))
system.time(res2<-multiseq.gaus(X,g,gridmult=0))


pdf("data_4096_1_80.pdf")
par(mfrow=c(2,1))
plot(X.1[10,],ylim=c(-0.5,3.5),pch=3,main="data_g0",ylab="y")
plot(X.2[10,],ylim=c(-0.5,3.5),pch=3,main="data_g1",ylab="y")
dev.off()



pdf("baseline_4096_1_80.pdf")
par(mfrow=c(2,1))
plot(mu.sp,ylim=c(-0.5,3.5),type='l',main="baseline",ylab="baseline")
lines(res1$baseline.mean,col=2)
lines(res1$baseline.mean-2*sqrt(res1$baseline.var),col=4)
lines(res1$baseline.mean+2*sqrt(res1$baseline.var),col=4)

plot(mu.sp,ylim=c(-0.5,3.5),type='l',main="baseline",ylab="baseline")
lines(res2$baseline.mean,type='l',col=2)
lines(res2$baseline.mean-2*sqrt(res2$baseline.var),col=4)
lines(res2$baseline.mean+2*sqrt(res2$baseline.var),col=4)
dev.off()


pdf("effect_4096_1_80.pdf")
par(mfrow=c(2,1))
plot(mu.spm-mu.sp,type='l',ylim=c(-1.2,1.6),main="effect",ylab="effect")
lines(res1$effect.mean,col=2)
lines(res1$effect.mean-2*sqrt(res1$effect.var),col=4)
lines(res1$effect.mean+2*sqrt(res1$effect.var),col=4)

plot(mu.spm-mu.sp,type='l',ylim=c(-1.2,1.6),main="effect",ylab="effect")
lines(res2$effect.mean,col=2)
lines(res2$effect.mean-2*sqrt(res2$effect.var),col=4)
lines(res2$effect.mean+2*sqrt(res2$effect.var),col=4)
dev.off()






