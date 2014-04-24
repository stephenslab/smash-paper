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


var.est.ash.s8.cons.v1=matrix(0,100,n)
var.est.ash.s8.cons.v2=matrix(0,100,n)
var.est.ash.s8.cons.v3=matrix(0,100,n)
var.est.ash.s8.cons.v4=matrix(0,100,n)

var.est.ash.s8.sine.v1=matrix(0,100,n)
var.est.ash.s8.sine.v2=matrix(0,100,n)
var.est.ash.s8.sine.v3=matrix(0,100,n)
var.est.ash.s8.sine.v4=matrix(0,100,n)

var.est.ash.s8.bump.v1=matrix(0,100,n)
var.est.ash.s8.bump.v2=matrix(0,100,n)
var.est.ash.s8.bump.v3=matrix(0,100,n)
var.est.ash.s8.bump.v4=matrix(0,100,n)

var.est.ash.s8.blk.v1=matrix(0,100,n)
var.est.ash.s8.blk.v2=matrix(0,100,n)
var.est.ash.s8.blk.v3=matrix(0,100,n)
var.est.ash.s8.blk.v4=matrix(0,100,n)

var.est.ash.s8.dop.v1=matrix(0,100,n)
var.est.ash.s8.dop.v2=matrix(0,100,n)
var.est.ash.s8.dop.v3=matrix(0,100,n)
var.est.ash.s8.dop.v4=matrix(0,100,n)


var.est.ash.s8.j.cons.v1=matrix(0,100,n)
var.est.ash.s8.j.cons.v2=matrix(0,100,n)
var.est.ash.s8.j.cons.v3=matrix(0,100,n)
var.est.ash.s8.j.cons.v4=matrix(0,100,n)

var.est.ash.s8.j.sine.v1=matrix(0,100,n)
var.est.ash.s8.j.sine.v2=matrix(0,100,n)
var.est.ash.s8.j.sine.v3=matrix(0,100,n)
var.est.ash.s8.j.sine.v4=matrix(0,100,n)

var.est.ash.s8.j.bump.v1=matrix(0,100,n)
var.est.ash.s8.j.bump.v2=matrix(0,100,n)
var.est.ash.s8.j.bump.v3=matrix(0,100,n)
var.est.ash.s8.j.bump.v4=matrix(0,100,n)

var.est.ash.s8.j.blk.v1=matrix(0,100,n)
var.est.ash.s8.j.blk.v2=matrix(0,100,n)
var.est.ash.s8.j.blk.v3=matrix(0,100,n)
var.est.ash.s8.j.blk.v4=matrix(0,100,n)

var.est.ash.s8.j.dop.v1=matrix(0,100,n)
var.est.ash.s8.j.dop.v2=matrix(0,100,n)
var.est.ash.s8.j.dop.v3=matrix(0,100,n)
var.est.ash.s8.j.dop.v4=matrix(0,100,n)


print("start")

for(i in 1:100){

  var.est.ash.s8.cons.v1[i,]=bayesmooth(sim.m.cons.v1[i,],v.est=TRUE,v.basis=TRUE,filter.number=8,family="DaubLeAsymm",weight=1)
  var.est.ash.s8.cons.v2[i,]=bayesmooth(sim.m.cons.v2[i,],v.est=TRUE,v.basis=TRUE,filter.number=8,family="DaubLeAsymm",weight=1)
  var.est.ash.s8.cons.v3[i,]=bayesmooth(sim.m.cons.v3[i,],v.est=TRUE,v.basis=TRUE,filter.number=8,family="DaubLeAsymm",weight=1)
  var.est.ash.s8.cons.v4[i,]=bayesmooth(sim.m.cons.v4[i,],v.est=TRUE,v.basis=TRUE,filter.number=8,family="DaubLeAsymm",weight=1)

  var.est.ash.s8.sine.v1[i,]=bayesmooth(sim.m.sine.v1[i,],v.est=TRUE,v.basis=TRUE,filter.number=8,family="DaubLeAsymm",weight=1)
  var.est.ash.s8.sine.v2[i,]=bayesmooth(sim.m.sine.v2[i,],v.est=TRUE,v.basis=TRUE,filter.number=8,family="DaubLeAsymm",weight=1)
  var.est.ash.s8.sine.v3[i,]=bayesmooth(sim.m.sine.v3[i,],v.est=TRUE,v.basis=TRUE,filter.number=8,family="DaubLeAsymm",weight=1)
  var.est.ash.s8.sine.v4[i,]=bayesmooth(sim.m.sine.v4[i,],v.est=TRUE,v.basis=TRUE,filter.number=8,family="DaubLeAsymm",weight=1)

  var.est.ash.s8.bump.v1[i,]=bayesmooth(sim.m.bump.v1[i,],v.est=TRUE,v.basis=TRUE,filter.number=8,family="DaubLeAsymm",weight=1)
  var.est.ash.s8.bump.v2[i,]=bayesmooth(sim.m.bump.v2[i,],v.est=TRUE,v.basis=TRUE,filter.number=8,family="DaubLeAsymm",weight=1)
  var.est.ash.s8.bump.v3[i,]=bayesmooth(sim.m.bump.v3[i,],v.est=TRUE,v.basis=TRUE,filter.number=8,family="DaubLeAsymm",weight=1)
  var.est.ash.s8.bump.v4[i,]=bayesmooth(sim.m.bump.v4[i,],v.est=TRUE,v.basis=TRUE,filter.number=8,family="DaubLeAsymm",weight=1)

  var.est.ash.s8.blk.v1[i,]=bayesmooth(sim.m.blk.v1[i,],v.est=TRUE,v.basis=TRUE,filter.number=8,family="DaubLeAsymm",weight=1)
  var.est.ash.s8.blk.v2[i,]=bayesmooth(sim.m.blk.v2[i,],v.est=TRUE,v.basis=TRUE,filter.number=8,family="DaubLeAsymm",weight=1)
  var.est.ash.s8.blk.v3[i,]=bayesmooth(sim.m.blk.v3[i,],v.est=TRUE,v.basis=TRUE,filter.number=8,family="DaubLeAsymm",weight=1)
  var.est.ash.s8.blk.v4[i,]=bayesmooth(sim.m.blk.v4[i,],v.est=TRUE,v.basis=TRUE,filter.number=8,family="DaubLeAsymm",weight=1)

  var.est.ash.s8.dop.v1[i,]=bayesmooth(sim.m.dop.v1[i,],v.est=TRUE,v.basis=TRUE,filter.number=8,family="DaubLeAsymm",weight=1)
  var.est.ash.s8.dop.v2[i,]=bayesmooth(sim.m.dop.v2[i,],v.est=TRUE,v.basis=TRUE,filter.number=8,family="DaubLeAsymm",weight=1)
  var.est.ash.s8.dop.v3[i,]=bayesmooth(sim.m.dop.v3[i,],v.est=TRUE,v.basis=TRUE,filter.number=8,family="DaubLeAsymm",weight=1)
  var.est.ash.s8.dop.v4[i,]=bayesmooth(sim.m.dop.v4[i,],v.est=TRUE,v.basis=TRUE,filter.number=8,family="DaubLeAsymm",weight=1)

  var.est.ash.s8.j.cons.v1[i,]=bayesmooth(sim.m.cons.v1[i,],v.est=TRUE,v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE,weight=1)
  var.est.ash.s8.j.cons.v2[i,]=bayesmooth(sim.m.cons.v2[i,],v.est=TRUE,v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE,weight=1)
  var.est.ash.s8.j.cons.v3[i,]=bayesmooth(sim.m.cons.v3[i,],v.est=TRUE,v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE,weight=1)
  var.est.ash.s8.j.cons.v4[i,]=bayesmooth(sim.m.cons.v4[i,],v.est=TRUE,v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE,weight=1)

  var.est.ash.s8.j.sine.v1[i,]=bayesmooth(sim.m.sine.v1[i,],v.est=TRUE,v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE,weight=1)
  var.est.ash.s8.j.sine.v2[i,]=bayesmooth(sim.m.sine.v2[i,],v.est=TRUE,v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE,weight=1)
  var.est.ash.s8.j.sine.v3[i,]=bayesmooth(sim.m.sine.v3[i,],v.est=TRUE,v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE,weight=1)
  var.est.ash.s8.j.sine.v4[i,]=bayesmooth(sim.m.sine.v4[i,],v.est=TRUE,v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE,weight=1)

  var.est.ash.s8.j.bump.v1[i,]=bayesmooth(sim.m.bump.v1[i,],v.est=TRUE,v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE,weight=1)
  var.est.ash.s8.j.bump.v2[i,]=bayesmooth(sim.m.bump.v2[i,],v.est=TRUE,v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE,weight=1)
  var.est.ash.s8.j.bump.v3[i,]=bayesmooth(sim.m.bump.v3[i,],v.est=TRUE,v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE,weight=1)
  var.est.ash.s8.j.bump.v4[i,]=bayesmooth(sim.m.bump.v4[i,],v.est=TRUE,v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE,weight=1)

  var.est.ash.s8.j.blk.v1[i,]=bayesmooth(sim.m.blk.v1[i,],v.est=TRUE,v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE,weight=1)
  var.est.ash.s8.j.blk.v2[i,]=bayesmooth(sim.m.blk.v2[i,],v.est=TRUE,v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE,weight=1)
  var.est.ash.s8.j.blk.v3[i,]=bayesmooth(sim.m.blk.v3[i,],v.est=TRUE,v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE,weight=1)
  var.est.ash.s8.j.blk.v4[i,]=bayesmooth(sim.m.blk.v4[i,],v.est=TRUE,v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE,weight=1)

  var.est.ash.s8.j.dop.v1[i,]=bayesmooth(sim.m.dop.v1[i,],v.est=TRUE,v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE,weight=1)
  var.est.ash.s8.j.dop.v2[i,]=bayesmooth(sim.m.dop.v2[i,],v.est=TRUE,v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE,weight=1)
  var.est.ash.s8.j.dop.v3[i,]=bayesmooth(sim.m.dop.v3[i,],v.est=TRUE,v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE,weight=1)
  var.est.ash.s8.j.dop.v4[i,]=bayesmooth(sim.m.dop.v4[i,],v.est=TRUE,v.basis=TRUE,filter.number=8,family="DaubLeAsymm",jash=TRUE,weight=1)

  print(i)
}


mse.ash.s8.cons.v1=mean(apply((var.est.ash.s8.cons.v1-rep(1,100)%o%var1)^2,1,mean))
mse.ash.s8.cons.v2=mean(apply((var.est.ash.s8.cons.v2-rep(1,100)%o%var2)^2,1,mean))
mse.ash.s8.cons.v3=mean(apply((var.est.ash.s8.cons.v3-rep(1,100)%o%var3)^2,1,mean))
mse.ash.s8.cons.v4=mean(apply((var.est.ash.s8.cons.v4-rep(1,100)%o%var4)^2,1,mean))

mse.ash.s8.sine.v1=mean(apply((var.est.ash.s8.sine.v1-rep(1,100)%o%var1)^2,1,mean))
mse.ash.s8.sine.v2=mean(apply((var.est.ash.s8.sine.v2-rep(1,100)%o%var2)^2,1,mean))
mse.ash.s8.sine.v3=mean(apply((var.est.ash.s8.sine.v3-rep(1,100)%o%var3)^2,1,mean))
mse.ash.s8.sine.v4=mean(apply((var.est.ash.s8.sine.v4-rep(1,100)%o%var4)^2,1,mean))

mse.ash.s8.bump.v1=mean(apply((var.est.ash.s8.bump.v1-rep(1,100)%o%var1)^2,1,mean))
mse.ash.s8.bump.v2=mean(apply((var.est.ash.s8.bump.v2-rep(1,100)%o%var2)^2,1,mean))
mse.ash.s8.bump.v3=mean(apply((var.est.ash.s8.bump.v3-rep(1,100)%o%var3)^2,1,mean))
mse.ash.s8.bump.v4=mean(apply((var.est.ash.s8.bump.v4-rep(1,100)%o%var4)^2,1,mean))

mse.ash.s8.blk.v1=mean(apply((var.est.ash.s8.blk.v1-rep(1,100)%o%var1)^2,1,mean))
mse.ash.s8.blk.v2=mean(apply((var.est.ash.s8.blk.v2-rep(1,100)%o%var2)^2,1,mean))
mse.ash.s8.blk.v3=mean(apply((var.est.ash.s8.blk.v3-rep(1,100)%o%var3)^2,1,mean))
mse.ash.s8.blk.v4=mean(apply((var.est.ash.s8.blk.v4-rep(1,100)%o%var4)^2,1,mean))

mse.ash.s8.dop.v1=mean(apply((var.est.ash.s8.dop.v1-rep(1,100)%o%var1)^2,1,mean))
mse.ash.s8.dop.v2=mean(apply((var.est.ash.s8.dop.v2-rep(1,100)%o%var2)^2,1,mean))
mse.ash.s8.dop.v3=mean(apply((var.est.ash.s8.dop.v3-rep(1,100)%o%var3)^2,1,mean))
mse.ash.s8.dop.v4=mean(apply((var.est.ash.s8.dop.v4-rep(1,100)%o%var4)^2,1,mean))


mse.ash.s8.j.cons.v1=mean(apply((var.est.ash.s8.j.cons.v1-rep(1,100)%o%var1)^2,1,mean))
mse.ash.s8.j.cons.v2=mean(apply((var.est.ash.s8.j.cons.v2-rep(1,100)%o%var2)^2,1,mean))
mse.ash.s8.j.cons.v3=mean(apply((var.est.ash.s8.j.cons.v3-rep(1,100)%o%var3)^2,1,mean))
mse.ash.s8.j.cons.v4=mean(apply((var.est.ash.s8.j.cons.v4-rep(1,100)%o%var4)^2,1,mean))

mse.ash.s8.j.sine.v1=mean(apply((var.est.ash.s8.j.sine.v1-rep(1,100)%o%var1)^2,1,mean))
mse.ash.s8.j.sine.v2=mean(apply((var.est.ash.s8.j.sine.v2-rep(1,100)%o%var2)^2,1,mean))
mse.ash.s8.j.sine.v3=mean(apply((var.est.ash.s8.j.sine.v3-rep(1,100)%o%var3)^2,1,mean))
mse.ash.s8.j.sine.v4=mean(apply((var.est.ash.s8.j.sine.v4-rep(1,100)%o%var4)^2,1,mean))

mse.ash.s8.j.bump.v1=mean(apply((var.est.ash.s8.j.bump.v1-rep(1,100)%o%var1)^2,1,mean))
mse.ash.s8.j.bump.v2=mean(apply((var.est.ash.s8.j.bump.v2-rep(1,100)%o%var2)^2,1,mean))
mse.ash.s8.j.bump.v3=mean(apply((var.est.ash.s8.j.bump.v3-rep(1,100)%o%var3)^2,1,mean))
mse.ash.s8.j.bump.v4=mean(apply((var.est.ash.s8.j.bump.v4-rep(1,100)%o%var4)^2,1,mean))

mse.ash.s8.j.blk.v1=mean(apply((var.est.ash.s8.j.blk.v1-rep(1,100)%o%var1)^2,1,mean))
mse.ash.s8.j.blk.v2=mean(apply((var.est.ash.s8.j.blk.v2-rep(1,100)%o%var2)^2,1,mean))
mse.ash.s8.j.blk.v3=mean(apply((var.est.ash.s8.j.blk.v3-rep(1,100)%o%var3)^2,1,mean))
mse.ash.s8.j.blk.v4=mean(apply((var.est.ash.s8.j.blk.v4-rep(1,100)%o%var4)^2,1,mean))

mse.ash.s8.j.dop.v1=mean(apply((var.est.ash.s8.j.dop.v1-rep(1,100)%o%var1)^2,1,mean))
mse.ash.s8.j.dop.v2=mean(apply((var.est.ash.s8.j.dop.v2-rep(1,100)%o%var2)^2,1,mean))
mse.ash.s8.j.dop.v3=mean(apply((var.est.ash.s8.j.dop.v3-rep(1,100)%o%var3)^2,1,mean))
mse.ash.s8.j.dop.v4=mean(apply((var.est.ash.s8.j.dop.v4-rep(1,100)%o%var4)^2,1,mean))


save.image("~/ashwave/simulation_1d_g/res_followup_s8_jash.RData")
