source(file.path("~/ashwave/Rcode/bayesmooth.R"))



n=2^12
t=1:n/n


pos=c(.1, .13, .15, .23, .25, .40, .44, .65, .76, .78, .81)
hgt = c(4, (-5), 3, (-4), 5, (-4.2), 2.1, 4.3, (-3.1), 2.1, (-4.2))
mu.blk = rep(0,n)
for(j in 1:length(pos)){
  mu.blk = mu.blk + (1 + sign(t-pos[j]))*(hgt[j]/2)
}


pos = c(.1, .13, .15, .23, .25, .40, .44, .65, .76, .78, .81)
hgt = c(4, 5, 3, 4, 5, 4.2, 2.1, 4.3, 3.1, 5.1, 4.2)
wth = c(.005, .005, .006, .01, .01, .03, .01, .01, .005, .008, .005)
mu.bump = rep(0,n)
for(j in 1:length(pos)){
  mu.bump = mu.bump + hgt[j]/(( 1 + (abs(t - pos[j])/wth[j]))^4)
}

dop.f=function(x) sqrt(x*(1-x))*sin((2*pi*1.05)/(x+0.05))
mu.dop=dop.f(t)

var1.ini=((3-20*t)*(t>=0&t<0.1)+(20*t-1)*(t>=0.1&t<0.25)+(4+(1-4*t)*18/19)*(t>=0.25&t<0.725)+(2.2+10*(t-0.725))*(t>=0.725&t<0.89)+(3.85-85*(t-0.89)/11)*(t>=0.89&t<=1))
var2.ini=(1+4*(exp(-550*(t-0.2)^2)+exp(-200*(t-0.5)^2)+exp(-950*(t-0.8)^2)))

mu.sine=sin(20*t)
mu.cons=rep(0,n)


var1=var1.ini/sqrt(var(var1.ini))
var2=var2.ini/sqrt(var(var2.ini))
var3=mu.bump/sqrt(var(mu.bump))
var4=(mu.dop+2)/sqrt(var(mu.dop))


l2norm=function(x) sum(x^2)
mise=function(x,y) l2norm(x-y)/l2norm(y)
mse=function(x,y) mean((x-y)^2)



set.seed(1121)
sim.m.cons.v1=matrix(rnorm(100*n,mu.cons,sqrt(var1)),nrow=100,byrow=TRUE)
set.seed(1121)
sim.m.cons.v2=matrix(rnorm(100*n,mu.cons,sqrt(var2)),nrow=100,byrow=TRUE)
set.seed(1121)
sim.m.cons.v3=matrix(rnorm(100*n,mu.cons,sqrt(var3)),nrow=100,byrow=TRUE)
set.seed(1121)
sim.m.cons.v4=matrix(rnorm(100*n,mu.cons,sqrt(var4)),nrow=100,byrow=TRUE)

set.seed(1121)
sim.m.sine.v1=matrix(rnorm(100*n,mu.sine,sqrt(var1)),nrow=100,byrow=TRUE)
set.seed(1121)
sim.m.sine.v2=matrix(rnorm(100*n,mu.sine,sqrt(var2)),nrow=100,byrow=TRUE)
set.seed(1121)
sim.m.sine.v3=matrix(rnorm(100*n,mu.sine,sqrt(var3)),nrow=100,byrow=TRUE)
set.seed(1121)
sim.m.sine.v4=matrix(rnorm(100*n,mu.sine,sqrt(var4)),nrow=100,byrow=TRUE)

set.seed(1121)
sim.m.bump.v1=matrix(rnorm(100*n,mu.bump,sqrt(var1)),nrow=100,byrow=TRUE)
set.seed(1121)
sim.m.bump.v2=matrix(rnorm(100*n,mu.bump,sqrt(var2)),nrow=100,byrow=TRUE)
set.seed(1121)
sim.m.bump.v3=matrix(rnorm(100*n,mu.bump,sqrt(var3)),nrow=100,byrow=TRUE)
set.seed(1121)
sim.m.bump.v4=matrix(rnorm(100*n,mu.bump,sqrt(var4)),nrow=100,byrow=TRUE)


set.seed(1121)
sim.m.blk.v1=matrix(rnorm(100*n,mu.blk,sqrt(var1)),nrow=100,byrow=TRUE)
set.seed(1121)
sim.m.blk.v2=matrix(rnorm(100*n,mu.blk,sqrt(var2)),nrow=100,byrow=TRUE)
set.seed(1121)
sim.m.blk.v3=matrix(rnorm(100*n,mu.blk,sqrt(var3)),nrow=100,byrow=TRUE)
set.seed(1121)
sim.m.blk.v4=matrix(rnorm(100*n,mu.blk,sqrt(var4)),nrow=100,byrow=TRUE)


set.seed(1121)
sim.m.dop.v1=matrix(rnorm(100*n,mu.dop,sqrt(var1)),nrow=100,byrow=TRUE)
set.seed(1121)
sim.m.dop.v2=matrix(rnorm(100*n,mu.dop,sqrt(var2)),nrow=100,byrow=TRUE)
set.seed(1121)
sim.m.dop.v3=matrix(rnorm(100*n,mu.dop,sqrt(var3)),nrow=100,byrow=TRUE)
set.seed(1121)
sim.m.dop.v4=matrix(rnorm(100*n,mu.dop,sqrt(var4)),nrow=100,byrow=TRUE)


var.est.ash.cons.v1=matrix(0,100,n)
var.est.ash.cons.v2=matrix(0,100,n)
var.est.ash.cons.v3=matrix(0,100,n)
var.est.ash.cons.v4=matrix(0,100,n)

var.est.ash.sine.v1=matrix(0,100,n)
var.est.ash.sine.v2=matrix(0,100,n)
var.est.ash.sine.v3=matrix(0,100,n)
var.est.ash.sine.v4=matrix(0,100,n)

var.est.ash.bump.v1=matrix(0,100,n)
var.est.ash.bump.v2=matrix(0,100,n)
var.est.ash.bump.v3=matrix(0,100,n)
var.est.ash.bump.v4=matrix(0,100,n)

var.est.ash.blk.v1=matrix(0,100,n)
var.est.ash.blk.v2=matrix(0,100,n)
var.est.ash.blk.v3=matrix(0,100,n)
var.est.ash.blk.v4=matrix(0,100,n)

var.est.ash.dop.v1=matrix(0,100,n)
var.est.ash.dop.v2=matrix(0,100,n)
var.est.ash.dop.v3=matrix(0,100,n)
var.est.ash.dop.v4=matrix(0,100,n)



for(i in 1:100){
var.est.ash.cons.v1[i,]=bayesmooth.jash(sim.m.cons.v1[i,],v.est=TRUE,weight=1)
var.est.ash.cons.v2[i,]=bayesmooth.jash(sim.m.cons.v2[i,],v.est=TRUE,weight=1)
var.est.ash.cons.v3[i,]=bayesmooth.jash(sim.m.cons.v3[i,],v.est=TRUE,weight=1)
var.est.ash.cons.v4[i,]=bayesmooth.jash(sim.m.cons.v4[i,],v.est=TRUE,weight=1)

var.est.ash.sine.v1[i,]=bayesmooth.jash(sim.m.sine.v1[i,],v.est=TRUE,weight=1)
var.est.ash.sine.v2[i,]=bayesmooth.jash(sim.m.sine.v2[i,],v.est=TRUE,weight=1)
var.est.ash.sine.v3[i,]=bayesmooth.jash(sim.m.sine.v3[i,],v.est=TRUE,weight=1)
var.est.ash.sine.v4[i,]=bayesmooth.jash(sim.m.sine.v4[i,],v.est=TRUE,weight=1)

var.est.ash.bump.v1[i,]=bayesmooth.jash(sim.m.bump.v1[i,],v.est=TRUE,weight=1)
var.est.ash.bump.v2[i,]=bayesmooth.jash(sim.m.bump.v2[i,],v.est=TRUE,weight=1)
var.est.ash.bump.v3[i,]=bayesmooth.jash(sim.m.bump.v3[i,],v.est=TRUE,weight=1)
var.est.ash.bump.v4[i,]=bayesmooth.jash(sim.m.bump.v4[i,],v.est=TRUE,weight=1)

var.est.ash.blk.v1[i,]=bayesmooth.jash(sim.m.blk.v1[i,],v.est=TRUE,weight=1)
var.est.ash.blk.v2[i,]=bayesmooth.jash(sim.m.blk.v2[i,],v.est=TRUE,weight=1)
var.est.ash.blk.v3[i,]=bayesmooth.jash(sim.m.blk.v3[i,],v.est=TRUE,weight=1)
var.est.ash.blk.v4[i,]=bayesmooth.jash(sim.m.blk.v4[i,],v.est=TRUE,weight=1)

var.est.ash.dop.v1[i,]=bayesmooth.jash(sim.m.dop.v1[i,],v.est=TRUE,weight=1)
var.est.ash.dop.v2[i,]=bayesmooth.jash(sim.m.dop.v2[i,],v.est=TRUE,weight=1)
var.est.ash.dop.v3[i,]=bayesmooth.jash(sim.m.dop.v3[i,],v.est=TRUE,weight=1)
var.est.ash.dop.v4[i,]=bayesmooth.jash(sim.m.dop.v4[i,],v.est=TRUE,weight=1)

print(i)
}



mse.ash.cons.v1=mean(apply((var.est.ash.cons.v1[1:100,]-rep(1,100)%o%var1)^2,1,mean))
mse.ash.cons.v2=mean(apply((var.est.ash.cons.v2[1:100,]-rep(1,100)%o%var2)^2,1,mean))
mse.ash.cons.v3=mean(apply((var.est.ash.cons.v3[1:100,]-rep(1,100)%o%var3)^2,1,mean))
mse.ash.cons.v4=mean(apply((var.est.ash.cons.v4[1:100,]-rep(1,100)%o%var4)^2,1,mean))

mse.ash.sine.v1=mean(apply((var.est.ash.sine.v1[1:100,]-rep(1,100)%o%var1)^2,1,mean))
mse.ash.sine.v2=mean(apply((var.est.ash.sine.v2[1:100,]-rep(1,100)%o%var2)^2,1,mean))
mse.ash.sine.v3=mean(apply((var.est.ash.sine.v3[1:100,]-rep(1,100)%o%var3)^2,1,mean))
mse.ash.sine.v4=mean(apply((var.est.ash.sine.v4[1:100,]-rep(1,100)%o%var4)^2,1,mean))

mse.ash.bump.v1=mean(apply((var.est.ash.bump.v1[1:100,]-rep(1,100)%o%var1)^2,1,mean))
mse.ash.bump.v2=mean(apply((var.est.ash.bump.v2[1:100,]-rep(1,100)%o%var2)^2,1,mean))
mse.ash.bump.v3=mean(apply((var.est.ash.bump.v3[1:100,]-rep(1,100)%o%var3)^2,1,mean))
mse.ash.bump.v4=mean(apply((var.est.ash.bump.v4[1:100,]-rep(1,100)%o%var4)^2,1,mean))

mse.ash.blk.v1=mean(apply((var.est.ash.blk.v1[1:100,]-rep(1,100)%o%var1)^2,1,mean))
mse.ash.blk.v2=mean(apply((var.est.ash.blk.v2[1:100,]-rep(1,100)%o%var2)^2,1,mean))
mse.ash.blk.v3=mean(apply((var.est.ash.blk.v3[1:100,]-rep(1,100)%o%var3)^2,1,mean))
mse.ash.blk.v4=mean(apply((var.est.ash.blk.v4[1:100,]-rep(1,100)%o%var4)^2,1,mean))

mse.ash.dop.v1=mean(apply((var.est.ash.dop.v1[1:100,]-rep(1,100)%o%var1)^2,1,mean))
mse.ash.dop.v2=mean(apply((var.est.ash.dop.v2[1:100,]-rep(1,100)%o%var2)^2,1,mean))
mse.ash.dop.v3=mean(apply((var.est.ash.dop.v3[1:100,]-rep(1,100)%o%var3)^2,1,mean))
mse.ash.dop.v4=mean(apply((var.est.ash.dop.v4[1:100,]-rep(1,100)%o%var4)^2,1,mean))


save.image("~/ashwave/simulation_1d_g/res_followup_jash_w1.RData")