source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/bayesmooth.R"))
source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/bayesmooth_new.R"))
source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/bayesmooth_mc.R"))
source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/bayesmooth_test.R"))
source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/bayesmooth_alt.R"))
source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/bayesmooth_alt1.R"))
source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/bayesmooth_alt2.R"))


n=2^12
t=1:n/n


pos=c(.1, .13, .15, .23, .25, .40, .44, .65, .76, .78, .81)
hgt = 2.88/5*c(4, (-5), 3, (-4), 5, (-4.2), 2.1, 4.3, (-3.1), 2.1, (-4.2))
mu.blk = rep(0,n)
for(j in 1:length(pos)){
  mu.blk = mu.blk + (1 + sign(t-pos[j]))*(hgt[j]/2)
}


pos = c(.1, .13, .15, .23, .25, .40, .44, .65, .76, .78, .81)
hgt = c(2, 3.6, 2, 5, 7.5, 6.9, 2, 4.8, 2, 4.1, 2.3)
wth = c(.005, .005, .006, .01, .01, .03, .01, .01, .005, .008, .005)
mu.bump = rep(0,n)
for(j in 1:length(pos)){
  mu.bump = mu.bump + hgt[j]/(( 1 + (abs(t - pos[j])/wth[j]))^4)
}

dop.f=function(x) sqrt(x*(1-x))*sin((2*pi*1.05)/(x+0.05))
mu.dop=dop.f(t)

mu.sine=sin(20*t)


mu.cons=rep(0,n)


var1=1.42*((3-20*t)*(t>=0&t<0.1)+(20*t-1)*(t>=0.1&t<0.25)+(4+(1-4*t)*18/19)*(t>=0.25&t<0.725)+(2.2+10*(t-0.725))*(t>=0.725&t<0.89)+(3.85-85*(t-0.89)/11)*(t>=0.89&t<=1))
var2=(1+4*(exp(-550*(t-0.2)^2)+exp(-200*(t-0.5)^2)+exp(-950*(t-0.8)^2)))/1.35
var3=mu.bump
var4=3.4*(2+mu.dop)


l2norm=function(x) sum(x^2)
mise=function(x,y) l2norm(x-y)/l2norm(y)
mse=function(x,y) mean((x-y)^2)



set.seed(1205)
sim.m.cons.v1=matrix(rnorm(100*n,mu.cons,sqrt(var1)),nrow=100,byrow=TRUE)
set.seed(1205)
sim.m.cons.v2=matrix(rnorm(100*n,mu.cons,sqrt(var2)),nrow=100,byrow=TRUE)
set.seed(1205)
sim.m.cons.v3=matrix(rnorm(100*n,mu.cons,sqrt(var3)),nrow=100,byrow=TRUE)
set.seed(1205)
sim.m.cons.v4=matrix(rnorm(100*n,mu.cons,sqrt(var4)),nrow=100,byrow=TRUE)

set.seed(1205)
sim.m.sine.v1=matrix(rnorm(100*n,mu.sine,sqrt(var1)),nrow=100,byrow=TRUE)
set.seed(1205)
sim.m.sine.v2=matrix(rnorm(100*n,mu.sine,sqrt(var2)),nrow=100,byrow=TRUE)
set.seed(1205)
sim.m.sine.v3=matrix(rnorm(100*n,mu.sine,sqrt(var3)),nrow=100,byrow=TRUE)
set.seed(1205)
sim.m.sine.v4=matrix(rnorm(100*n,mu.sine,sqrt(var4)),nrow=100,byrow=TRUE)

set.seed(1205)
sim.m.bump.v1=matrix(rnorm(100*n,mu.bump,sqrt(var1)),nrow=100,byrow=TRUE)
set.seed(1205)
sim.m.bump.v2=matrix(rnorm(100*n,mu.bump,sqrt(var2)),nrow=100,byrow=TRUE)
set.seed(1205)
sim.m.bump.v3=matrix(rnorm(100*n,mu.bump,sqrt(var3)),nrow=100,byrow=TRUE)
set.seed(1205)
sim.m.bump.v4=matrix(rnorm(100*n,mu.bump,sqrt(var4)),nrow=100,byrow=TRUE)


set.seed(1205)
sim.m.blk.v1=matrix(rnorm(100*n,mu.blk,sqrt(var1)),nrow=100,byrow=TRUE)
set.seed(1205)
sim.m.blk.v2=matrix(rnorm(100*n,mu.blk,sqrt(var2)),nrow=100,byrow=TRUE)
set.seed(1205)
sim.m.blk.v3=matrix(rnorm(100*n,mu.blk,sqrt(var3)),nrow=100,byrow=TRUE)
set.seed(1205)
sim.m.blk.v4=matrix(rnorm(100*n,mu.blk,sqrt(var4)),nrow=100,byrow=TRUE)


set.seed(1205)
sim.m.dop.v1=matrix(rnorm(100*n,mu.dop,sqrt(var1)),nrow=100,byrow=TRUE)
set.seed(1205)
sim.m.dop.v2=matrix(rnorm(100*n,mu.dop,sqrt(var2)),nrow=100,byrow=TRUE)
set.seed(1205)
sim.m.dop.v3=matrix(rnorm(100*n,mu.dop,sqrt(var3)),nrow=100,byrow=TRUE)
set.seed(1205)
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




var.est.asha.cons.v1=matrix(0,100,n)
var.est.asha.cons.v2=matrix(0,100,n)
var.est.asha.cons.v3=matrix(0,100,n)
var.est.asha.cons.v4=matrix(0,100,n)

var.est.asha.sine.v1=matrix(0,100,n)
var.est.asha.sine.v2=matrix(0,100,n)
var.est.asha.sine.v3=matrix(0,100,n)
var.est.asha.sine.v4=matrix(0,100,n)

var.est.asha.bump.v1=matrix(0,100,n)
var.est.asha.bump.v2=matrix(0,100,n)
var.est.asha.bump.v3=matrix(0,100,n)
var.est.asha.bump.v4=matrix(0,100,n)

var.est.asha.blk.v1=matrix(0,100,n)
var.est.asha.blk.v2=matrix(0,100,n)
var.est.asha.blk.v3=matrix(0,100,n)
var.est.asha.blk.v4=matrix(0,100,n)

var.est.asha.dop.v1=matrix(0,100,n)
var.est.asha.dop.v2=matrix(0,100,n)
var.est.asha.dop.v3=matrix(0,100,n)
var.est.asha.dop.v4=matrix(0,100,n)



var.est.asha1.cons.v1=matrix(0,100,n)
var.est.asha1.cons.v2=matrix(0,100,n)
var.est.asha1.cons.v3=matrix(0,100,n)
var.est.asha1.cons.v4=matrix(0,100,n)

var.est.asha1.sine.v1=matrix(0,100,n)
var.est.asha1.sine.v2=matrix(0,100,n)
var.est.asha1.sine.v3=matrix(0,100,n)
var.est.asha1.sine.v4=matrix(0,100,n)

var.est.asha1.bump.v1=matrix(0,100,n)
var.est.asha1.bump.v2=matrix(0,100,n)
var.est.asha1.bump.v3=matrix(0,100,n)
var.est.asha1.bump.v4=matrix(0,100,n)

var.est.asha1.blk.v1=matrix(0,100,n)
var.est.asha1.blk.v2=matrix(0,100,n)
var.est.asha1.blk.v3=matrix(0,100,n)
var.est.asha1.blk.v4=matrix(0,100,n)

var.est.asha1.dop.v1=matrix(0,100,n)
var.est.asha1.dop.v2=matrix(0,100,n)
var.est.asha1.dop.v3=matrix(0,100,n)
var.est.asha1.dop.v4=matrix(0,100,n)



var.est.ashn.cons.v1=matrix(0,100,n)
var.est.ashn.cons.v2=matrix(0,100,n)
var.est.ashn.cons.v3=matrix(0,100,n)
var.est.ashn.cons.v4=matrix(0,100,n)

var.est.ashn.sine.v1=matrix(0,100,n)
var.est.ashn.sine.v2=matrix(0,100,n)
var.est.ashn.sine.v3=matrix(0,100,n)
var.est.ashn.sine.v4=matrix(0,100,n)

var.est.ashn.bump.v1=matrix(0,100,n)
var.est.ashn.bump.v2=matrix(0,100,n)
var.est.ashn.bump.v3=matrix(0,100,n)
var.est.ashn.bump.v4=matrix(0,100,n)

var.est.ashn.blk.v1=matrix(0,100,n)
var.est.ashn.blk.v2=matrix(0,100,n)
var.est.ashn.blk.v3=matrix(0,100,n)
var.est.ashn.blk.v4=matrix(0,100,n)

var.est.ashn.dop.v1=matrix(0,100,n)
var.est.ashn.dop.v2=matrix(0,100,n)
var.est.ashn.dop.v3=matrix(0,100,n)
var.est.ashn.dop.v4=matrix(0,100,n)





var.est.ashm.cons.v1=matrix(0,100,n)
var.est.ashm.cons.v2=matrix(0,100,n)
var.est.ashm.cons.v3=matrix(0,100,n)
var.est.ashm.cons.v4=matrix(0,100,n)

var.est.ashm.sine.v1=matrix(0,100,n)
var.est.ashm.sine.v2=matrix(0,100,n)
var.est.ashm.sine.v3=matrix(0,100,n)
var.est.ashm.sine.v4=matrix(0,100,n)

var.est.ashm.bump.v1=matrix(0,100,n)
var.est.ashm.bump.v2=matrix(0,100,n)
var.est.ashm.bump.v3=matrix(0,100,n)
var.est.ashm.bump.v4=matrix(0,100,n)

var.est.ashm.blk.v1=matrix(0,100,n)
var.est.ashm.blk.v2=matrix(0,100,n)
var.est.ashm.blk.v3=matrix(0,100,n)
var.est.ashm.blk.v4=matrix(0,100,n)

var.est.ashm.dop.v1=matrix(0,100,n)
var.est.ashm.dop.v2=matrix(0,100,n)
var.est.ashm.dop.v3=matrix(0,100,n)
var.est.ashm.dop.v4=matrix(0,100,n)



var.est.asht.cons.v1=matrix(0,100,n)
var.est.asht.cons.v2=matrix(0,100,n)
var.est.asht.cons.v3=matrix(0,100,n)
var.est.asht.cons.v4=matrix(0,100,n)

var.est.asht.sine.v1=matrix(0,100,n)
var.est.asht.sine.v2=matrix(0,100,n)
var.est.asht.sine.v3=matrix(0,100,n)
var.est.asht.sine.v4=matrix(0,100,n)

var.est.asht.bump.v1=matrix(0,100,n)
var.est.asht.bump.v2=matrix(0,100,n)
var.est.asht.bump.v3=matrix(0,100,n)
var.est.asht.bump.v4=matrix(0,100,n)

var.est.asht.blk.v1=matrix(0,100,n)
var.est.asht.blk.v2=matrix(0,100,n)
var.est.asht.blk.v3=matrix(0,100,n)
var.est.asht.blk.v4=matrix(0,100,n)

var.est.asht.dop.v1=matrix(0,100,n)
var.est.asht.dop.v2=matrix(0,100,n)
var.est.asht.dop.v3=matrix(0,100,n)
var.est.asht.dop.v4=matrix(0,100,n)



for(i in 1:100){
var.est.ash.cons.v1[i,]=bayesmooth(sim.m.cons.v1[i,],v.est=TRUE)
var.est.ash.cons.v2[i,]=bayesmooth(sim.m.cons.v2[i,],v.est=TRUE)
var.est.ash.cons.v3[i,]=bayesmooth(sim.m.cons.v3[i,],v.est=TRUE)
var.est.ash.cons.v4[i,]=bayesmooth(sim.m.cons.v4[i,],v.est=TRUE)

var.est.ash.sine.v1[i,]=bayesmooth(sim.m.sine.v1[i,],v.est=TRUE)
var.est.ash.sine.v2[i,]=bayesmooth(sim.m.sine.v2[i,],v.est=TRUE)
var.est.ash.sine.v3[i,]=bayesmooth(sim.m.sine.v3[i,],v.est=TRUE)
var.est.ash.sine.v4[i,]=bayesmooth(sim.m.sine.v4[i,],v.est=TRUE)

var.est.ash.bump.v1[i,]=bayesmooth(sim.m.bump.v1[i,],v.est=TRUE)
var.est.ash.bump.v2[i,]=bayesmooth(sim.m.bump.v2[i,],v.est=TRUE)
var.est.ash.bump.v3[i,]=bayesmooth(sim.m.bump.v3[i,],v.est=TRUE)
var.est.ash.bump.v4[i,]=bayesmooth(sim.m.bump.v4[i,],v.est=TRUE)

var.est.ash.blk.v1[i,]=bayesmooth(sim.m.blk.v1[i,],v.est=TRUE)
var.est.ash.blk.v2[i,]=bayesmooth(sim.m.blk.v2[i,],v.est=TRUE)
var.est.ash.blk.v3[i,]=bayesmooth(sim.m.blk.v3[i,],v.est=TRUE)
var.est.ash.blk.v4[i,]=bayesmooth(sim.m.blk.v4[i,],v.est=TRUE)

var.est.ash.dop.v1[i,]=bayesmooth(sim.m.dop.v1[i,],v.est=TRUE)
var.est.ash.dop.v2[i,]=bayesmooth(sim.m.dop.v2[i,],v.est=TRUE)
var.est.ash.dop.v3[i,]=bayesmooth(sim.m.dop.v3[i,],v.est=TRUE)
var.est.ash.dop.v4[i,]=bayesmooth(sim.m.dop.v4[i,],v.est=TRUE)



var.est.asha.cons.v1[i,]=bayesmooth.alt(sim.m.cons.v1[i,],v.est=TRUE)
var.est.asha.cons.v2[i,]=bayesmooth.alt(sim.m.cons.v2[i,],v.est=TRUE)
var.est.asha.cons.v3[i,]=bayesmooth.alt(sim.m.cons.v3[i,],v.est=TRUE)
var.est.asha.cons.v4[i,]=bayesmooth.alt(sim.m.cons.v4[i,],v.est=TRUE)

var.est.asha.sine.v1[i,]=bayesmooth.alt(sim.m.sine.v1[i,],v.est=TRUE)
var.est.asha.sine.v2[i,]=bayesmooth.alt(sim.m.sine.v2[i,],v.est=TRUE)
var.est.asha.sine.v3[i,]=bayesmooth.alt(sim.m.sine.v3[i,],v.est=TRUE)
var.est.asha.sine.v4[i,]=bayesmooth.alt(sim.m.sine.v4[i,],v.est=TRUE)

var.est.asha.bump.v1[i,]=bayesmooth.alt(sim.m.bump.v1[i,],v.est=TRUE)
var.est.asha.bump.v2[i,]=bayesmooth.alt(sim.m.bump.v2[i,],v.est=TRUE)
var.est.asha.bump.v3[i,]=bayesmooth.alt(sim.m.bump.v3[i,],v.est=TRUE)
var.est.asha.bump.v4[i,]=bayesmooth.alt(sim.m.bump.v4[i,],v.est=TRUE)

var.est.asha.blk.v1[i,]=bayesmooth.alt(sim.m.blk.v1[i,],v.est=TRUE)
var.est.asha.blk.v2[i,]=bayesmooth.alt(sim.m.blk.v2[i,],v.est=TRUE)
var.est.asha.blk.v3[i,]=bayesmooth.alt(sim.m.blk.v3[i,],v.est=TRUE)
var.est.asha.blk.v4[i,]=bayesmooth.alt(sim.m.blk.v4[i,],v.est=TRUE)

var.est.asha.dop.v1[i,]=bayesmooth.alt(sim.m.dop.v1[i,],v.est=TRUE)
var.est.asha.dop.v2[i,]=bayesmooth.alt(sim.m.dop.v2[i,],v.est=TRUE)
var.est.asha.dop.v3[i,]=bayesmooth.alt(sim.m.dop.v3[i,],v.est=TRUE)
var.est.asha.dop.v4[i,]=bayesmooth.alt(sim.m.dop.v4[i,],v.est=TRUE)




var.est.asha1.cons.v1[i,]=bayesmooth.alt1(sim.m.cons.v1[i,],v.est=TRUE)
var.est.asha1.cons.v2[i,]=bayesmooth.alt1(sim.m.cons.v2[i,],v.est=TRUE)
var.est.asha1.cons.v3[i,]=bayesmooth.alt1(sim.m.cons.v3[i,],v.est=TRUE)
var.est.asha1.cons.v4[i,]=bayesmooth.alt1(sim.m.cons.v4[i,],v.est=TRUE)

var.est.asha1.sine.v1[i,]=bayesmooth.alt1(sim.m.sine.v1[i,],v.est=TRUE)
var.est.asha1.sine.v2[i,]=bayesmooth.alt1(sim.m.sine.v2[i,],v.est=TRUE)
var.est.asha1.sine.v3[i,]=bayesmooth.alt1(sim.m.sine.v3[i,],v.est=TRUE)
var.est.asha1.sine.v4[i,]=bayesmooth.alt1(sim.m.sine.v4[i,],v.est=TRUE)

var.est.asha1.bump.v1[i,]=bayesmooth.alt1(sim.m.bump.v1[i,],v.est=TRUE)
var.est.asha1.bump.v2[i,]=bayesmooth.alt1(sim.m.bump.v2[i,],v.est=TRUE)
var.est.asha1.bump.v3[i,]=bayesmooth.alt1(sim.m.bump.v3[i,],v.est=TRUE)
var.est.asha1.bump.v4[i,]=bayesmooth.alt1(sim.m.bump.v4[i,],v.est=TRUE)

var.est.asha1.blk.v1[i,]=bayesmooth.alt1(sim.m.blk.v1[i,],v.est=TRUE)
var.est.asha1.blk.v2[i,]=bayesmooth.alt1(sim.m.blk.v2[i,],v.est=TRUE)
var.est.asha1.blk.v3[i,]=bayesmooth.alt1(sim.m.blk.v3[i,],v.est=TRUE)
var.est.asha1.blk.v4[i,]=bayesmooth.alt1(sim.m.blk.v4[i,],v.est=TRUE)

var.est.asha1.dop.v1[i,]=bayesmooth.alt1(sim.m.dop.v1[i,],v.est=TRUE)
var.est.asha1.dop.v2[i,]=bayesmooth.alt1(sim.m.dop.v2[i,],v.est=TRUE)
var.est.asha1.dop.v3[i,]=bayesmooth.alt1(sim.m.dop.v3[i,],v.est=TRUE)
var.est.asha1.dop.v4[i,]=bayesmooth.alt1(sim.m.dop.v4[i,],v.est=TRUE)




var.est.ashn.cons.v1[i,]=bayesmooth.new(sim.m.cons.v1[i,],v.est=TRUE)
var.est.ashn.cons.v2[i,]=bayesmooth.new(sim.m.cons.v2[i,],v.est=TRUE)
var.est.ashn.cons.v3[i,]=bayesmooth.new(sim.m.cons.v3[i,],v.est=TRUE)
var.est.ashn.cons.v4[i,]=bayesmooth.new(sim.m.cons.v4[i,],v.est=TRUE)

var.est.ashn.sine.v1[i,]=bayesmooth.new(sim.m.sine.v1[i,],v.est=TRUE)
var.est.ashn.sine.v2[i,]=bayesmooth.new(sim.m.sine.v2[i,],v.est=TRUE)
var.est.ashn.sine.v3[i,]=bayesmooth.new(sim.m.sine.v3[i,],v.est=TRUE)
var.est.ashn.sine.v4[i,]=bayesmooth.new(sim.m.sine.v4[i,],v.est=TRUE)

var.est.ashn.bump.v1[i,]=bayesmooth.new(sim.m.bump.v1[i,],v.est=TRUE)
var.est.ashn.bump.v2[i,]=bayesmooth.new(sim.m.bump.v2[i,],v.est=TRUE)
var.est.ashn.bump.v3[i,]=bayesmooth.new(sim.m.bump.v3[i,],v.est=TRUE)
var.est.ashn.bump.v4[i,]=bayesmooth.new(sim.m.bump.v4[i,],v.est=TRUE)

var.est.ashn.blk.v1[i,]=bayesmooth.new(sim.m.blk.v1[i,],v.est=TRUE)
var.est.ashn.blk.v2[i,]=bayesmooth.new(sim.m.blk.v2[i,],v.est=TRUE)
var.est.ashn.blk.v3[i,]=bayesmooth.new(sim.m.blk.v3[i,],v.est=TRUE)
var.est.ashn.blk.v4[i,]=bayesmooth.new(sim.m.blk.v4[i,],v.est=TRUE)

var.est.ashn.dop.v1[i,]=bayesmooth.new(sim.m.dop.v1[i,],v.est=TRUE)
var.est.ashn.dop.v2[i,]=bayesmooth.new(sim.m.dop.v2[i,],v.est=TRUE)
var.est.ashn.dop.v3[i,]=bayesmooth.new(sim.m.dop.v3[i,],v.est=TRUE)
var.est.ashn.dop.v4[i,]=bayesmooth.new(sim.m.dop.v4[i,],v.est=TRUE)


var.est.ashm.cons.v1[i,]=bayesmooth.mc(sim.m.cons.v1[i,],v.est=TRUE)
var.est.ashm.cons.v2[i,]=bayesmooth.mc(sim.m.cons.v2[i,],v.est=TRUE)
var.est.ashm.cons.v3[i,]=bayesmooth.mc(sim.m.cons.v3[i,],v.est=TRUE)
var.est.ashm.cons.v4[i,]=bayesmooth.mc(sim.m.cons.v4[i,],v.est=TRUE)

var.est.ashm.sine.v1[i,]=bayesmooth.mc(sim.m.sine.v1[i,],v.est=TRUE)
var.est.ashm.sine.v2[i,]=bayesmooth.mc(sim.m.sine.v2[i,],v.est=TRUE)
var.est.ashm.sine.v3[i,]=bayesmooth.mc(sim.m.sine.v3[i,],v.est=TRUE)
var.est.ashm.sine.v4[i,]=bayesmooth.mc(sim.m.sine.v4[i,],v.est=TRUE)

var.est.ashm.bump.v1[i,]=bayesmooth.mc(sim.m.bump.v1[i,],v.est=TRUE)
var.est.ashm.bump.v2[i,]=bayesmooth.mc(sim.m.bump.v2[i,],v.est=TRUE)
var.est.ashm.bump.v3[i,]=bayesmooth.mc(sim.m.bump.v3[i,],v.est=TRUE)
var.est.ashm.bump.v4[i,]=bayesmooth.mc(sim.m.bump.v4[i,],v.est=TRUE)

var.est.ashm.blk.v1[i,]=bayesmooth.mc(sim.m.blk.v1[i,],v.est=TRUE)
var.est.ashm.blk.v2[i,]=bayesmooth.mc(sim.m.blk.v2[i,],v.est=TRUE)
var.est.ashm.blk.v3[i,]=bayesmooth.mc(sim.m.blk.v3[i,],v.est=TRUE)
var.est.ashm.blk.v4[i,]=bayesmooth.mc(sim.m.blk.v4[i,],v.est=TRUE)

var.est.ashm.dop.v1[i,]=bayesmooth.mc(sim.m.dop.v1[i,],v.est=TRUE)
var.est.ashm.dop.v2[i,]=bayesmooth.mc(sim.m.dop.v2[i,],v.est=TRUE)
var.est.ashm.dop.v3[i,]=bayesmooth.mc(sim.m.dop.v3[i,],v.est=TRUE)
var.est.ashm.dop.v4[i,]=bayesmooth.mc(sim.m.dop.v4[i,],v.est=TRUE)



var.est.asht.cons.v1[i,]=bayesmooth.test(sim.m.cons.v1[i,],v.est=TRUE)
var.est.asht.cons.v2[i,]=bayesmooth.test(sim.m.cons.v2[i,],v.est=TRUE)
var.est.asht.cons.v3[i,]=bayesmooth.test(sim.m.cons.v3[i,],v.est=TRUE)
var.est.asht.cons.v4[i,]=bayesmooth.test(sim.m.cons.v4[i,],v.est=TRUE)

var.est.asht.sine.v1[i,]=bayesmooth.test(sim.m.sine.v1[i,],v.est=TRUE)
var.est.asht.sine.v2[i,]=bayesmooth.test(sim.m.sine.v2[i,],v.est=TRUE)
var.est.asht.sine.v3[i,]=bayesmooth.test(sim.m.sine.v3[i,],v.est=TRUE)
var.est.asht.sine.v4[i,]=bayesmooth.test(sim.m.sine.v4[i,],v.est=TRUE)

var.est.asht.bump.v1[i,]=bayesmooth.test(sim.m.bump.v1[i,],v.est=TRUE)
var.est.asht.bump.v2[i,]=bayesmooth.test(sim.m.bump.v2[i,],v.est=TRUE)
var.est.asht.bump.v3[i,]=bayesmooth.test(sim.m.bump.v3[i,],v.est=TRUE)
var.est.asht.bump.v4[i,]=bayesmooth.test(sim.m.bump.v4[i,],v.est=TRUE)

var.est.asht.blk.v1[i,]=bayesmooth.test(sim.m.blk.v1[i,],v.est=TRUE)
var.est.asht.blk.v2[i,]=bayesmooth.test(sim.m.blk.v2[i,],v.est=TRUE)
var.est.asht.blk.v3[i,]=bayesmooth.test(sim.m.blk.v3[i,],v.est=TRUE)
var.est.asht.blk.v4[i,]=bayesmooth.test(sim.m.blk.v4[i,],v.est=TRUE)

var.est.asht.dop.v1[i,]=bayesmooth.test(sim.m.dop.v1[i,],v.est=TRUE)
var.est.asht.dop.v2[i,]=bayesmooth.test(sim.m.dop.v2[i,],v.est=TRUE)
var.est.asht.dop.v3[i,]=bayesmooth.test(sim.m.dop.v3[i,],v.est=TRUE)
var.est.asht.dop.v4[i,]=bayesmooth.test(sim.m.dop.v4[i,],v.est=TRUE)

print(i)
}



mse.ashm.cons.v1=median(apply((var.est.ashm.cons.v1-rep(1,100)%o%var1)^2,1,mean))
mse.ashm.cons.v2=median(apply((var.est.ashm.cons.v2-rep(1,100)%o%var2)^2,1,mean))
mse.ashm.cons.v3=median(apply((var.est.ashm.cons.v3-rep(1,100)%o%var3)^2,1,mean))
mse.ashm.cons.v4=median(apply((var.est.ashm.cons.v4-rep(1,100)%o%var4)^2,1,mean))

mse.ashm.sine.v1=median(apply((var.est.ashm.sine.v1-rep(1,100)%o%var1)^2,1,mean))
mse.ashm.sine.v2=median(apply((var.est.ashm.sine.v2-rep(1,100)%o%var2)^2,1,mean))
mse.ashm.sine.v3=median(apply((var.est.ashm.sine.v3-rep(1,100)%o%var3)^2,1,mean))
mse.ashm.sine.v4=median(apply((var.est.ashm.sine.v4-rep(1,100)%o%var4)^2,1,mean))

mse.ashm.bump.v1=median(apply((var.est.ashm.bump.v1-rep(1,100)%o%var1)^2,1,mean))
mse.ashm.bump.v2=median(apply((var.est.ashm.bump.v2-rep(1,100)%o%var2)^2,1,mean))
mse.ashm.bump.v3=median(apply((var.est.ashm.bump.v3-rep(1,100)%o%var3)^2,1,mean))
mse.ashm.bump.v4=median(apply((var.est.ashm.bump.v4-rep(1,100)%o%var4)^2,1,mean))

mse.ashm.blk.v1=median(apply((var.est.ashm.blk.v1-rep(1,100)%o%var1)^2,1,mean))
mse.ashm.blk.v2=median(apply((var.est.ashm.blk.v2-rep(1,100)%o%var2)^2,1,mean))
mse.ashm.blk.v3=median(apply((var.est.ashm.blk.v3-rep(1,100)%o%var3)^2,1,mean))
mse.ashm.blk.v4=median(apply((var.est.ashm.blk.v4-rep(1,100)%o%var4)^2,1,mean))

mse.ashm.dop.v1=median(apply((var.est.ashm.dop.v1-rep(1,100)%o%var1)^2,1,mean))
mse.ashm.dop.v2=median(apply((var.est.ashm.dop.v2-rep(1,100)%o%var2)^2,1,mean))
mse.ashm.dop.v3=median(apply((var.est.ashm.dop.v3-rep(1,100)%o%var3)^2,1,mean))
mse.ashm.dop.v4=median(apply((var.est.ashm.dop.v4-rep(1,100)%o%var4)^2,1,mean))




mse.ashn.cons.v1=median(apply((var.est.ashn.cons.v1-rep(1,100)%o%var1)^2,1,mean))
mse.ashn.cons.v2=median(apply((var.est.ashn.cons.v2-rep(1,100)%o%var2)^2,1,mean))
mse.ashn.cons.v3=median(apply((var.est.ashn.cons.v3-rep(1,100)%o%var3)^2,1,mean))
mse.ashn.cons.v4=median(apply((var.est.ashn.cons.v4-rep(1,100)%o%var4)^2,1,mean))

mse.ashn.sine.v1=median(apply((var.est.ashn.sine.v1-rep(1,100)%o%var1)^2,1,mean))
mse.ashn.sine.v2=median(apply((var.est.ashn.sine.v2-rep(1,100)%o%var2)^2,1,mean))
mse.ashn.sine.v3=median(apply((var.est.ashn.sine.v3-rep(1,100)%o%var3)^2,1,mean))
mse.ashn.sine.v4=median(apply((var.est.ashn.sine.v4-rep(1,100)%o%var4)^2,1,mean))

mse.ashn.bump.v1=median(apply((var.est.ashn.bump.v1-rep(1,100)%o%var1)^2,1,mean))
mse.ashn.bump.v2=median(apply((var.est.ashn.bump.v2-rep(1,100)%o%var2)^2,1,mean))
mse.ashn.bump.v3=median(apply((var.est.ashn.bump.v3-rep(1,100)%o%var3)^2,1,mean))
mse.ashn.bump.v4=median(apply((var.est.ashn.bump.v4-rep(1,100)%o%var4)^2,1,mean))

mse.ashn.blk.v1=median(apply((var.est.ashn.blk.v1-rep(1,100)%o%var1)^2,1,mean))
mse.ashn.blk.v2=median(apply((var.est.ashn.blk.v2-rep(1,100)%o%var2)^2,1,mean))
mse.ashn.blk.v3=median(apply((var.est.ashn.blk.v3-rep(1,100)%o%var3)^2,1,mean))
mse.ashn.blk.v4=median(apply((var.est.ashn.blk.v4-rep(1,100)%o%var4)^2,1,mean))

mse.ashn.dop.v1=median(apply((var.est.ashn.dop.v1-rep(1,100)%o%var1)^2,1,mean))
mse.ashn.dop.v2=median(apply((var.est.ashn.dop.v2-rep(1,100)%o%var2)^2,1,mean))
mse.ashn.dop.v3=median(apply((var.est.ashn.dop.v3-rep(1,100)%o%var3)^2,1,mean))
mse.ashn.dop.v4=median(apply((var.est.ashn.dop.v4-rep(1,100)%o%var4)^2,1,mean))




mse.asht.cons.v1=median(apply((var.est.asht.cons.v1-rep(1,100)%o%var1)^2,1,mean))
mse.asht.cons.v2=median(apply((var.est.asht.cons.v2-rep(1,100)%o%var2)^2,1,mean))
mse.asht.cons.v3=median(apply((var.est.asht.cons.v3-rep(1,100)%o%var3)^2,1,mean))
mse.asht.cons.v4=median(apply((var.est.asht.cons.v4-rep(1,100)%o%var4)^2,1,mean))

mse.asht.sine.v1=median(apply((var.est.asht.sine.v1-rep(1,100)%o%var1)^2,1,mean))
mse.asht.sine.v2=median(apply((var.est.asht.sine.v2-rep(1,100)%o%var2)^2,1,mean))
mse.asht.sine.v3=median(apply((var.est.asht.sine.v3-rep(1,100)%o%var3)^2,1,mean))
mse.asht.sine.v4=median(apply((var.est.asht.sine.v4-rep(1,100)%o%var4)^2,1,mean))

mse.asht.bump.v1=median(apply((var.est.asht.bump.v1-rep(1,100)%o%var1)^2,1,mean))
mse.asht.bump.v2=median(apply((var.est.asht.bump.v2-rep(1,100)%o%var2)^2,1,mean))
mse.asht.bump.v3=median(apply((var.est.asht.bump.v3-rep(1,100)%o%var3)^2,1,mean))
mse.asht.bump.v4=median(apply((var.est.asht.bump.v4-rep(1,100)%o%var4)^2,1,mean))

mse.asht.blk.v1=median(apply((var.est.asht.blk.v1-rep(1,100)%o%var1)^2,1,mean))
mse.asht.blk.v2=median(apply((var.est.asht.blk.v2-rep(1,100)%o%var2)^2,1,mean))
mse.asht.blk.v3=median(apply((var.est.asht.blk.v3-rep(1,100)%o%var3)^2,1,mean))
mse.asht.blk.v4=median(apply((var.est.asht.blk.v4-rep(1,100)%o%var4)^2,1,mean))

mse.asht.dop.v1=median(apply((var.est.asht.dop.v1-rep(1,100)%o%var1)^2,1,mean))
mse.asht.dop.v2=median(apply((var.est.asht.dop.v2-rep(1,100)%o%var2)^2,1,mean))
mse.asht.dop.v3=median(apply((var.est.asht.dop.v3-rep(1,100)%o%var3)^2,1,mean))
mse.asht.dop.v4=median(apply((var.est.asht.dop.v4-rep(1,100)%o%var4)^2,1,mean))




mse.ash.cons.v1=median(apply((var.est.ash.cons.v1-rep(1,100)%o%var1)^2,1,mean))
mse.ash.cons.v2=median(apply((var.est.ash.cons.v2-rep(1,100)%o%var2)^2,1,mean))
mse.ash.cons.v3=median(apply((var.est.ash.cons.v3-rep(1,100)%o%var3)^2,1,mean))
mse.ash.cons.v4=median(apply((var.est.ash.cons.v4-rep(1,100)%o%var4)^2,1,mean))

mse.ash.sine.v1=median(apply((var.est.ash.sine.v1-rep(1,100)%o%var1)^2,1,mean))
mse.ash.sine.v2=median(apply((var.est.ash.sine.v2-rep(1,100)%o%var2)^2,1,mean))
mse.ash.sine.v3=median(apply((var.est.ash.sine.v3-rep(1,100)%o%var3)^2,1,mean))
mse.ash.sine.v4=median(apply((var.est.ash.sine.v4-rep(1,100)%o%var4)^2,1,mean))

mse.ash.bump.v1=median(apply((var.est.ash.bump.v1-rep(1,100)%o%var1)^2,1,mean))
mse.ash.bump.v2=median(apply((var.est.ash.bump.v2-rep(1,100)%o%var2)^2,1,mean))
mse.ash.bump.v3=median(apply((var.est.ash.bump.v3-rep(1,100)%o%var3)^2,1,mean))
mse.ash.bump.v4=median(apply((var.est.ash.bump.v4-rep(1,100)%o%var4)^2,1,mean))

mse.ash.blk.v1=median(apply((var.est.ash.blk.v1-rep(1,100)%o%var1)^2,1,mean))
mse.ash.blk.v2=median(apply((var.est.ash.blk.v2-rep(1,100)%o%var2)^2,1,mean))
mse.ash.blk.v3=median(apply((var.est.ash.blk.v3-rep(1,100)%o%var3)^2,1,mean))
mse.ash.blk.v4=median(apply((var.est.ash.blk.v4-rep(1,100)%o%var4)^2,1,mean))

mse.ash.dop.v1=median(apply((var.est.ash.dop.v1-rep(1,100)%o%var1)^2,1,mean))
mse.ash.dop.v2=median(apply((var.est.ash.dop.v2-rep(1,100)%o%var2)^2,1,mean))
mse.ash.dop.v3=median(apply((var.est.ash.dop.v3-rep(1,100)%o%var3)^2,1,mean))
mse.ash.dop.v4=median(apply((var.est.ash.dop.v4-rep(1,100)%o%var4)^2,1,mean))




mse.asha.cons.v1=median(apply((var.est.asha.cons.v1-rep(1,100)%o%var1)^2,1,mean))
mse.asha.cons.v2=median(apply((var.est.asha.cons.v2-rep(1,100)%o%var2)^2,1,mean))
mse.asha.cons.v3=median(apply((var.est.asha.cons.v3-rep(1,100)%o%var3)^2,1,mean))
mse.asha.cons.v4=median(apply((var.est.asha.cons.v4-rep(1,100)%o%var4)^2,1,mean))

mse.asha.sine.v1=median(apply((var.est.asha.sine.v1-rep(1,100)%o%var1)^2,1,mean))
mse.asha.sine.v2=median(apply((var.est.asha.sine.v2-rep(1,100)%o%var2)^2,1,mean))
mse.asha.sine.v3=median(apply((var.est.asha.sine.v3-rep(1,100)%o%var3)^2,1,mean))
mse.asha.sine.v4=median(apply((var.est.asha.sine.v4-rep(1,100)%o%var4)^2,1,mean))

mse.asha.bump.v1=median(apply((var.est.asha.bump.v1-rep(1,100)%o%var1)^2,1,mean))
mse.asha.bump.v2=median(apply((var.est.asha.bump.v2-rep(1,100)%o%var2)^2,1,mean))
mse.asha.bump.v3=median(apply((var.est.asha.bump.v3-rep(1,100)%o%var3)^2,1,mean))
mse.asha.bump.v4=median(apply((var.est.asha.bump.v4-rep(1,100)%o%var4)^2,1,mean))

mse.asha.blk.v1=median(apply((var.est.asha.blk.v1-rep(1,100)%o%var1)^2,1,mean))
mse.asha.blk.v2=median(apply((var.est.asha.blk.v2-rep(1,100)%o%var2)^2,1,mean))
mse.asha.blk.v3=median(apply((var.est.asha.blk.v3-rep(1,100)%o%var3)^2,1,mean))
mse.asha.blk.v4=median(apply((var.est.asha.blk.v4-rep(1,100)%o%var4)^2,1,mean))

mse.asha.dop.v1=median(apply((var.est.asha.dop.v1-rep(1,100)%o%var1)^2,1,mean))
mse.asha.dop.v2=median(apply((var.est.asha.dop.v2-rep(1,100)%o%var2)^2,1,mean))
mse.asha.dop.v3=median(apply((var.est.asha.dop.v3-rep(1,100)%o%var3)^2,1,mean))
mse.asha.dop.v4=median(apply((var.est.asha.dop.v4-rep(1,100)%o%var4)^2,1,mean))





mse.asha1.cons.v1=median(apply((var.est.asha1.cons.v1-rep(1,100)%o%var1)^2,1,mean))
mse.asha1.cons.v2=median(apply((var.est.asha1.cons.v2-rep(1,100)%o%var2)^2,1,mean))
mse.asha1.cons.v3=median(apply((var.est.asha1.cons.v3-rep(1,100)%o%var3)^2,1,mean))
mse.asha1.cons.v4=median(apply((var.est.asha1.cons.v4-rep(1,100)%o%var4)^2,1,mean))

mse.asha1.sine.v1=median(apply((var.est.asha1.sine.v1-rep(1,100)%o%var1)^2,1,mean))
mse.asha1.sine.v2=median(apply((var.est.asha1.sine.v2-rep(1,100)%o%var2)^2,1,mean))
mse.asha1.sine.v3=median(apply((var.est.asha1.sine.v3-rep(1,100)%o%var3)^2,1,mean))
mse.asha1.sine.v4=median(apply((var.est.asha1.sine.v4-rep(1,100)%o%var4)^2,1,mean))

mse.asha1.bump.v1=median(apply((var.est.asha1.bump.v1-rep(1,100)%o%var1)^2,1,mean))
mse.asha1.bump.v2=median(apply((var.est.asha1.bump.v2-rep(1,100)%o%var2)^2,1,mean))
mse.asha1.bump.v3=median(apply((var.est.asha1.bump.v3-rep(1,100)%o%var3)^2,1,mean))
mse.asha1.bump.v4=median(apply((var.est.asha1.bump.v4-rep(1,100)%o%var4)^2,1,mean))

mse.asha1.blk.v1=median(apply((var.est.asha1.blk.v1-rep(1,100)%o%var1)^2,1,mean))
mse.asha1.blk.v2=median(apply((var.est.asha1.blk.v2-rep(1,100)%o%var2)^2,1,mean))
mse.asha1.blk.v3=median(apply((var.est.asha1.blk.v3-rep(1,100)%o%var3)^2,1,mean))
mse.asha1.blk.v4=median(apply((var.est.asha1.blk.v4-rep(1,100)%o%var4)^2,1,mean))

mse.asha1.dop.v1=median(apply((var.est.asha1.dop.v1-rep(1,100)%o%var1)^2,1,mean))
mse.asha1.dop.v2=median(apply((var.est.asha1.dop.v2-rep(1,100)%o%var2)^2,1,mean))
mse.asha1.dop.v3=median(apply((var.est.asha1.dop.v3-rep(1,100)%o%var3)^2,1,mean))
mse.asha1.dop.v4=median(apply((var.est.asha1.dop.v4-rep(1,100)%o%var4)^2,1,mean))




##############################

#1-no correction for skewness
#2-3 iterations instead of 2
#3-estimate var first
#4-normal version
#5-use RSS for both var est
#6-use RSS for final step






mse.cons.v1=c(mse.ashm.cons.v1,
mse.ashn.cons.v1,
mse.asht.cons.v1,
mse.ash.cons.v1,
mse.asha.cons.v1,
mse.asha1.cons.v1
)

names(mse.cons.v1)=c("no_skew",
"3 iter",
"est var",
"current",
"RSS both",
"RSS final"
)



mse.cons.v2=c(mse.ashm.cons.v2,
mse.ashn.cons.v2,
mse.asht.cons.v2,
mse.ash.cons.v2,
mse.asha.cons.v2,
mse.asha1.cons.v2
)

names(mse.cons.v2)=c("no_skew",
"3 iter",
"est var",
"current",
"RSS both",
"RSS final"
)


mse.cons.v3=c(mse.ashm.cons.v3,
mse.ashn.cons.v3,
mse.asht.cons.v3,
mse.ash.cons.v3,
mse.asha.cons.v3,
mse.asha1.cons.v3
)

names(mse.cons.v3)=c("no_skew",
"3 iter",
"est var",
"current",
"RSS both",
"RSS final"
)



mse.cons.v4=c(mse.ashm.cons.v4,
mse.ashn.cons.v4,
mse.asht.cons.v4,
mse.ash.cons.v4,
mse.asha.cons.v4,
mse.asha1.cons.v4
)

names(mse.cons.v4)=c("no_skew",
"3 iter",
"est var",
"current",
"RSS both",
"RSS final"
)





mse.sine.v1=c(mse.ashm.sine.v1,
mse.ashn.sine.v1,
mse.asht.sine.v1,
mse.ash.sine.v1,
mse.asha.sine.v1,
mse.asha1.sine.v1
)

names(mse.sine.v1)=c("no_skew",
"3 iter",
"est var",
"current",
"RSS both",
"RSS final"
)



mse.sine.v2=c(mse.ashm.sine.v2,
mse.ashn.sine.v2,
mse.asht.sine.v2,
mse.ash.sine.v2,
mse.asha.sine.v2,
mse.asha1.sine.v2
)

names(mse.sine.v2)=c("no_skew",
"3 iter",
"est var",
"current",
"RSS both",
"RSS final"
)


mse.sine.v3=c(mse.ashm.sine.v3,
mse.ashn.sine.v3,
mse.asht.sine.v3,
mse.ash.sine.v3,
mse.asha.sine.v3,
mse.asha1.sine.v3
)

names(mse.sine.v3)=c("no_skew",
"3 iter",
"est var",
"current",
"RSS both",
"RSS final"
)



mse.sine.v4=c(mse.ashm.sine.v4,
mse.ashn.sine.v4,
mse.asht.sine.v4,
mse.ash.sine.v4,
mse.asha.sine.v4,
mse.asha1.sine.v4
)

names(mse.sine.v4)=c("no_skew",
"3 iter",
"est var",
"current",
"RSS both",
"RSS final"
)





mse.bump.v1=c(mse.ashm.bump.v1,
mse.ashn.bump.v1,
mse.asht.bump.v1,
mse.ash.bump.v1,
mse.asha.bump.v1,
mse.asha1.bump.v1
)

names(mse.bump.v1)=c("no_skew",
"3 iter",
"est var",
"current",
"RSS both",
"RSS final"
)



mse.bump.v2=c(mse.ashm.bump.v2,
mse.ashn.bump.v2,
mse.asht.bump.v2,
mse.ash.bump.v2,
mse.asha.bump.v2,
mse.asha1.bump.v2
)

names(mse.bump.v2)=c("no_skew",
"3 iter",
"est var",
"current",
"RSS both",
"RSS final"
)


mse.bump.v3=c(mse.ashm.bump.v3,
mse.ashn.bump.v3,
mse.asht.bump.v3,
mse.ash.bump.v3,
mse.asha.bump.v3,
mse.asha1.bump.v3
)

names(mse.bump.v3)=c("no_skew",
"3 iter",
"est var",
"current",
"RSS both",
"RSS final"
)



mse.bump.v4=c(mse.ashm.bump.v4,
mse.ashn.bump.v4,
mse.asht.bump.v4,
mse.ash.bump.v4,
mse.asha.bump.v4,
mse.asha1.bump.v4
)

names(mse.bump.v4)=c("no_skew",
"3 iter",
"est var",
"current",
"RSS both",
"RSS final"
)





mse.blk.v1=c(mse.ashm.blk.v1,
mse.ashn.blk.v1,
mse.asht.blk.v1,
mse.ash.blk.v1,
mse.asha.blk.v1,
mse.asha1.blk.v1
)

names(mse.blk.v1)=c("no_skew",
"3 iter",
"est var",
"current",
"RSS both",
"RSS final"
)



mse.blk.v2=c(mse.ashm.blk.v2,
mse.ashn.blk.v2,
mse.asht.blk.v2,
mse.ash.blk.v2,
mse.asha.blk.v2,
mse.asha1.blk.v2
)

names(mse.blk.v2)=c("no_skew",
"3 iter",
"est var",
"current",
"RSS both",
"RSS final"
)


mse.blk.v3=c(mse.ashm.blk.v3,
mse.ashn.blk.v3,
mse.asht.blk.v3,
mse.ash.blk.v3,
mse.asha.blk.v3,
mse.asha1.blk.v3
)

names(mse.blk.v3)=c("no_skew",
"3 iter",
"est var",
"current",
"RSS both",
"RSS final"
)



mse.blk.v4=c(mse.ashm.blk.v4,
mse.ashn.blk.v4,
mse.asht.blk.v4,
mse.ash.blk.v4,
mse.asha.blk.v4,
mse.asha1.blk.v4
)

names(mse.blk.v4)=c("no_skew",
"3 iter",
"est var",
"current",
"RSS both",
"RSS final"
)





mse.dop.v1=c(mse.ashm.dop.v1,
mse.ashn.dop.v1,
mse.asht.dop.v1,
mse.ash.dop.v1,
mse.asha.dop.v1,
mse.asha1.dop.v1
)

names(mse.dop.v1)=c("no_skew",
"3 iter",
"est var",
"current",
"RSS both",
"RSS final"
)



mse.dop.v2=c(mse.ashm.dop.v2,
mse.ashn.dop.v2,
mse.asht.dop.v2,
mse.ash.dop.v2,
mse.asha.dop.v2,
mse.asha1.dop.v2
)

names(mse.dop.v2)=c("no_skew",
"3 iter",
"est var",
"current",
"RSS both",
"RSS final"
)


mse.dop.v3=c(mse.ashm.dop.v3,
mse.ashn.dop.v3,
mse.asht.dop.v3,
mse.ash.dop.v3,
mse.asha.dop.v3,
mse.asha1.dop.v3
)

names(mse.dop.v3)=c("no_skew",
"3 iter",
"est var",
"current",
"RSS both",
"RSS final"
)



mse.dop.v4=c(mse.ashm.dop.v4,
mse.ashn.dop.v4,
mse.asht.dop.v4,
mse.ash.dop.v4,
mse.asha.dop.v4,
mse.asha1.dop.v4
)

names(mse.dop.v4)=c("no_skew",
"3 iter",
"est var",
"current",
"RSS both",
"RSS final"
)


sort(mse.cons.v1)
sort(mse.cons.v2)
sort(mse.cons.v3)
sort(mse.cons.v4)

sort(mse.sine.v1)
sort(mse.sine.v2)
sort(mse.sine.v3)
sort(mse.sine.v4)

sort(mse.bump.v1)
sort(mse.bump.v2)
sort(mse.bump.v3)
sort(mse.bump.v4)

sort(mse.blk.v1)
sort(mse.blk.v2)
sort(mse.blk.v3)
sort(mse.blk.v4)

sort(mse.dop.v1)
sort(mse.dop.v2)
sort(mse.dop.v3)
sort(mse.dop.v4)

