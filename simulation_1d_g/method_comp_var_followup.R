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


var.est.asha2.cons.v1=matrix(0,100,n)
var.est.asha2.cons.v2=matrix(0,100,n)
var.est.asha2.cons.v3=matrix(0,100,n)
var.est.asha2.cons.v4=matrix(0,100,n)

var.est.asha2.sine.v1=matrix(0,100,n)
var.est.asha2.sine.v2=matrix(0,100,n)
var.est.asha2.sine.v3=matrix(0,100,n)
var.est.asha2.sine.v4=matrix(0,100,n)

var.est.asha2.bump.v1=matrix(0,100,n)
var.est.asha2.bump.v2=matrix(0,100,n)
var.est.asha2.bump.v3=matrix(0,100,n)
var.est.asha2.bump.v4=matrix(0,100,n)

var.est.asha2.blk.v1=matrix(0,100,n)
var.est.asha2.blk.v2=matrix(0,100,n)
var.est.asha2.blk.v3=matrix(0,100,n)
var.est.asha2.blk.v4=matrix(0,100,n)

var.est.asha2.dop.v1=matrix(0,100,n)
var.est.asha2.dop.v2=matrix(0,100,n)
var.est.asha2.dop.v3=matrix(0,100,n)
var.est.asha2.dop.v4=matrix(0,100,n)



for(i in 1:100){
var.est.ash.cons.v1[i,]=bayesmooth(sim.m.cons.v1[i,],v.est=TRUE,weight=0)
var.est.ash.cons.v2[i,]=bayesmooth(sim.m.cons.v2[i,],v.est=TRUE,weight=0)
var.est.ash.cons.v3[i,]=bayesmooth(sim.m.cons.v3[i,],v.est=TRUE,weight=0)
var.est.ash.cons.v4[i,]=bayesmooth(sim.m.cons.v4[i,],v.est=TRUE,weight=0)

var.est.ash.sine.v1[i,]=bayesmooth(sim.m.sine.v1[i,],v.est=TRUE,weight=0)
var.est.ash.sine.v2[i,]=bayesmooth(sim.m.sine.v2[i,],v.est=TRUE,weight=0)
var.est.ash.sine.v3[i,]=bayesmooth(sim.m.sine.v3[i,],v.est=TRUE,weight=0)
var.est.ash.sine.v4[i,]=bayesmooth(sim.m.sine.v4[i,],v.est=TRUE,weight=0)

var.est.ash.bump.v1[i,]=bayesmooth(sim.m.bump.v1[i,],v.est=TRUE,weight=0)
var.est.ash.bump.v2[i,]=bayesmooth(sim.m.bump.v2[i,],v.est=TRUE,weight=0)
var.est.ash.bump.v3[i,]=bayesmooth(sim.m.bump.v3[i,],v.est=TRUE,weight=0)
var.est.ash.bump.v4[i,]=bayesmooth(sim.m.bump.v4[i,],v.est=TRUE,weight=0)

var.est.ash.blk.v1[i,]=bayesmooth(sim.m.blk.v1[i,],v.est=TRUE,weight=0)
var.est.ash.blk.v2[i,]=bayesmooth(sim.m.blk.v2[i,],v.est=TRUE,weight=0)
var.est.ash.blk.v3[i,]=bayesmooth(sim.m.blk.v3[i,],v.est=TRUE,weight=0)
var.est.ash.blk.v4[i,]=bayesmooth(sim.m.blk.v4[i,],v.est=TRUE,weight=0)

var.est.ash.dop.v1[i,]=bayesmooth(sim.m.dop.v1[i,],v.est=TRUE,weight=0)
var.est.ash.dop.v2[i,]=bayesmooth(sim.m.dop.v2[i,],v.est=TRUE,weight=0)
var.est.ash.dop.v3[i,]=bayesmooth(sim.m.dop.v3[i,],v.est=TRUE,weight=0)
var.est.ash.dop.v4[i,]=bayesmooth(sim.m.dop.v4[i,],v.est=TRUE,weight=0)



var.est.asha.cons.v1[i,]=bayesmooth(sim.m.cons.v1[i,],v.est=TRUE,weight=1)
var.est.asha.cons.v2[i,]=bayesmooth(sim.m.cons.v2[i,],v.est=TRUE,weight=1)
var.est.asha.cons.v3[i,]=bayesmooth(sim.m.cons.v3[i,],v.est=TRUE,weight=1)
var.est.asha.cons.v4[i,]=bayesmooth(sim.m.cons.v4[i,],v.est=TRUE,weight=1)

var.est.asha.sine.v1[i,]=bayesmooth(sim.m.sine.v1[i,],v.est=TRUE,weight=1)
var.est.asha.sine.v2[i,]=bayesmooth(sim.m.sine.v2[i,],v.est=TRUE,weight=1)
var.est.asha.sine.v3[i,]=bayesmooth(sim.m.sine.v3[i,],v.est=TRUE,weight=1)
var.est.asha.sine.v4[i,]=bayesmooth(sim.m.sine.v4[i,],v.est=TRUE,weight=1)

var.est.asha.bump.v1[i,]=bayesmooth(sim.m.bump.v1[i,],v.est=TRUE,weight=1)
var.est.asha.bump.v2[i,]=bayesmooth(sim.m.bump.v2[i,],v.est=TRUE,weight=1)
var.est.asha.bump.v3[i,]=bayesmooth(sim.m.bump.v3[i,],v.est=TRUE,weight=1)
var.est.asha.bump.v4[i,]=bayesmooth(sim.m.bump.v4[i,],v.est=TRUE,weight=1)

var.est.asha.blk.v1[i,]=bayesmooth(sim.m.blk.v1[i,],v.est=TRUE,weight=1)
var.est.asha.blk.v2[i,]=bayesmooth(sim.m.blk.v2[i,],v.est=TRUE,weight=1)
var.est.asha.blk.v3[i,]=bayesmooth(sim.m.blk.v3[i,],v.est=TRUE,weight=1)
var.est.asha.blk.v4[i,]=bayesmooth(sim.m.blk.v4[i,],v.est=TRUE,weight=1)

var.est.asha.dop.v1[i,]=bayesmooth(sim.m.dop.v1[i,],v.est=TRUE,weight=1)
var.est.asha.dop.v2[i,]=bayesmooth(sim.m.dop.v2[i,],v.est=TRUE,weight=1)
var.est.asha.dop.v3[i,]=bayesmooth(sim.m.dop.v3[i,],v.est=TRUE,weight=1)
var.est.asha.dop.v4[i,]=bayesmooth(sim.m.dop.v4[i,],v.est=TRUE,weight=1)




var.est.asha1.cons.v1[i,]=bayesmooth(sim.m.cons.v1[i,],v.est=TRUE,weight=0)
var.est.asha1.cons.v2[i,]=bayesmooth(sim.m.cons.v2[i,],v.est=TRUE,weight=0)
var.est.asha1.cons.v3[i,]=bayesmooth(sim.m.cons.v3[i,],v.est=TRUE,weight=0)
var.est.asha1.cons.v4[i,]=bayesmooth(sim.m.cons.v4[i,],v.est=TRUE,weight=0)

var.est.asha1.sine.v1[i,]=bayesmooth(sim.m.sine.v1[i,],v.est=TRUE,weight=0)
var.est.asha1.sine.v2[i,]=bayesmooth(sim.m.sine.v2[i,],v.est=TRUE,weight=0)
var.est.asha1.sine.v3[i,]=bayesmooth(sim.m.sine.v3[i,],v.est=TRUE,weight=0)
var.est.asha1.sine.v4[i,]=bayesmooth(sim.m.sine.v4[i,],v.est=TRUE,weight=0)

var.est.asha1.bump.v1[i,]=bayesmooth(sim.m.bump.v1[i,],v.est=TRUE,weight=0)
var.est.asha1.bump.v2[i,]=bayesmooth(sim.m.bump.v2[i,],v.est=TRUE,weight=0)
var.est.asha1.bump.v3[i,]=bayesmooth(sim.m.bump.v3[i,],v.est=TRUE,weight=0)
var.est.asha1.bump.v4[i,]=bayesmooth(sim.m.bump.v4[i,],v.est=TRUE,weight=0)

var.est.asha1.blk.v1[i,]=bayesmooth(sim.m.blk.v1[i,],v.est=TRUE,weight=0)
var.est.asha1.blk.v2[i,]=bayesmooth(sim.m.blk.v2[i,],v.est=TRUE,weight=0)
var.est.asha1.blk.v3[i,]=bayesmooth(sim.m.blk.v3[i,],v.est=TRUE,weight=0)
var.est.asha1.blk.v4[i,]=bayesmooth(sim.m.blk.v4[i,],v.est=TRUE,weight=0)

var.est.asha1.dop.v1[i,]=bayesmooth(sim.m.dop.v1[i,],v.est=TRUE,weight=0)
var.est.asha1.dop.v2[i,]=bayesmooth(sim.m.dop.v2[i,],v.est=TRUE,weight=0)
var.est.asha1.dop.v3[i,]=bayesmooth(sim.m.dop.v3[i,],v.est=TRUE,weight=0)
var.est.asha1.dop.v4[i,]=bayesmooth(sim.m.dop.v4[i,],v.est=TRUE,weight=0)




var.est.asha2.cons.v1[i,]=bayesmooth(sim.m.cons.v1[i,],v.est=TRUE,weight=0.5)
var.est.asha2.cons.v2[i,]=bayesmooth(sim.m.cons.v2[i,],v.est=TRUE,weight=0.5)
var.est.asha2.cons.v3[i,]=bayesmooth(sim.m.cons.v3[i,],v.est=TRUE,weight=0.5)
var.est.asha2.cons.v4[i,]=bayesmooth(sim.m.cons.v4[i,],v.est=TRUE,weight=0.5)

var.est.asha2.sine.v1[i,]=bayesmooth(sim.m.sine.v1[i,],v.est=TRUE,weight=0.5)
var.est.asha2.sine.v2[i,]=bayesmooth(sim.m.sine.v2[i,],v.est=TRUE,weight=0.5)
var.est.asha2.sine.v3[i,]=bayesmooth(sim.m.sine.v3[i,],v.est=TRUE,weight=0.5)
var.est.asha2.sine.v4[i,]=bayesmooth(sim.m.sine.v4[i,],v.est=TRUE,weight=0.5)

var.est.asha2.bump.v1[i,]=bayesmooth(sim.m.bump.v1[i,],v.est=TRUE,weight=0.5)
var.est.asha2.bump.v2[i,]=bayesmooth(sim.m.bump.v2[i,],v.est=TRUE,weight=0.5)
var.est.asha2.bump.v3[i,]=bayesmooth(sim.m.bump.v3[i,],v.est=TRUE,weight=0.5)
var.est.asha2.bump.v4[i,]=bayesmooth(sim.m.bump.v4[i,],v.est=TRUE,weight=0.5)

var.est.asha2.blk.v1[i,]=bayesmooth(sim.m.blk.v1[i,],v.est=TRUE,weight=0.5)
var.est.asha2.blk.v2[i,]=bayesmooth(sim.m.blk.v2[i,],v.est=TRUE,weight=0.5)
var.est.asha2.blk.v3[i,]=bayesmooth(sim.m.blk.v3[i,],v.est=TRUE,weight=0.5)
var.est.asha2.blk.v4[i,]=bayesmooth(sim.m.blk.v4[i,],v.est=TRUE,weight=0.5)

var.est.asha2.dop.v1[i,]=bayesmooth(sim.m.dop.v1[i,],v.est=TRUE,weight=0.5)
var.est.asha2.dop.v2[i,]=bayesmooth(sim.m.dop.v2[i,],v.est=TRUE,weight=0.5)
var.est.asha2.dop.v3[i,]=bayesmooth(sim.m.dop.v3[i,],v.est=TRUE,weight=0.5)
var.est.asha2.dop.v4[i,]=bayesmooth(sim.m.dop.v4[i,],v.est=TRUE,weight=0.5)

print(i)
}



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


mse.asha2.cons.v1=median(apply((var.est.asha2.cons.v1-rep(1,100)%o%var1)^2,1,mean))
mse.asha2.cons.v2=median(apply((var.est.asha2.cons.v2-rep(1,100)%o%var2)^2,1,mean))
mse.asha2.cons.v3=median(apply((var.est.asha2.cons.v3-rep(1,100)%o%var3)^2,1,mean))
mse.asha2.cons.v4=median(apply((var.est.asha2.cons.v4-rep(1,100)%o%var4)^2,1,mean))

mse.asha2.sine.v1=median(apply((var.est.asha2.sine.v1-rep(1,100)%o%var1)^2,1,mean))
mse.asha2.sine.v2=median(apply((var.est.asha2.sine.v2-rep(1,100)%o%var2)^2,1,mean))
mse.asha2.sine.v3=median(apply((var.est.asha2.sine.v3-rep(1,100)%o%var3)^2,1,mean))
mse.asha2.sine.v4=median(apply((var.est.asha2.sine.v4-rep(1,100)%o%var4)^2,1,mean))

mse.asha2.bump.v1=median(apply((var.est.asha2.bump.v1-rep(1,100)%o%var1)^2,1,mean))
mse.asha2.bump.v2=median(apply((var.est.asha2.bump.v2-rep(1,100)%o%var2)^2,1,mean))
mse.asha2.bump.v3=median(apply((var.est.asha2.bump.v3-rep(1,100)%o%var3)^2,1,mean))
mse.asha2.bump.v4=median(apply((var.est.asha2.bump.v4-rep(1,100)%o%var4)^2,1,mean))

mse.asha2.blk.v1=median(apply((var.est.asha2.blk.v1-rep(1,100)%o%var1)^2,1,mean))
mse.asha2.blk.v2=median(apply((var.est.asha2.blk.v2-rep(1,100)%o%var2)^2,1,mean))
mse.asha2.blk.v3=median(apply((var.est.asha2.blk.v3-rep(1,100)%o%var3)^2,1,mean))
mse.asha2.blk.v4=median(apply((var.est.asha2.blk.v4-rep(1,100)%o%var4)^2,1,mean))

mse.asha2.dop.v1=median(apply((var.est.asha2.dop.v1-rep(1,100)%o%var1)^2,1,mean))
mse.asha2.dop.v2=median(apply((var.est.asha2.dop.v2-rep(1,100)%o%var2)^2,1,mean))
mse.asha2.dop.v3=median(apply((var.est.asha2.dop.v3-rep(1,100)%o%var3)^2,1,mean))
mse.asha2.dop.v4=median(apply((var.est.asha2.dop.v4-rep(1,100)%o%var4)^2,1,mean))



mse.cons.v1=c(
mse.ash.cons.v1,
mse.asha.cons.v1,
mse.asha1.cons.v1,
mse.asha2.cons.v1
)

names(mse.cons.v1)=c(
"current",
"RSS both",
"RSS final",
"RSS weighted"
)



mse.cons.v2=c(
mse.ash.cons.v2,
mse.asha.cons.v2,
mse.asha1.cons.v2,
mse.asha2.cons.v2
)

names(mse.cons.v2)=c(
"current",
"RSS both",
"RSS final",
"RSS weighted"
)


mse.cons.v3=c(
mse.ash.cons.v3,
mse.asha.cons.v3,
mse.asha1.cons.v3,
mse.asha2.cons.v3
)

names(mse.cons.v3)=c(
"current",
"RSS both",
"RSS final",
"RSS weighted"
)



mse.cons.v4=c(
mse.ash.cons.v4,
mse.asha.cons.v4,
mse.asha1.cons.v4,
mse.asha2.cons.v4
)

names(mse.cons.v4)=c(
"current",
"RSS both",
"RSS final",
"RSS weighted"
)



mse.sine.v1=c(
mse.ash.sine.v1,
mse.asha.sine.v1,
mse.asha1.sine.v1,
mse.asha2.sine.v1
)

names(mse.sine.v1)=c(
"current",
"RSS both",
"RSS final",
"RSS weighted"
)



mse.sine.v2=c(
mse.ash.sine.v2,
mse.asha.sine.v2,
mse.asha1.sine.v2,
mse.asha2.sine.v2
)

names(mse.sine.v2)=c(
"current",
"RSS both",
"RSS final",
"RSS weighted"
)


mse.sine.v3=c(
mse.ash.sine.v3,
mse.asha.sine.v3,
mse.asha1.sine.v3,
mse.asha2.sine.v3
)

names(mse.sine.v3)=c(
"current",
"RSS both",
"RSS final",
"RSS weighted"
)



mse.sine.v4=c(
mse.ash.sine.v4,
mse.asha.sine.v4,
mse.asha1.sine.v4,
mse.asha2.sine.v4
)

names(mse.sine.v4)=c(
"current",
"RSS both",
"RSS final",
"RSS weighted"
)





mse.bump.v1=c(
mse.ash.bump.v1,
mse.asha.bump.v1,
mse.asha1.bump.v1,
mse.asha2.bump.v1
)

names(mse.bump.v1)=c(
"current",
"RSS both",
"RSS final",
"RSS weighted"
)



mse.bump.v2=c(
mse.ash.bump.v2,
mse.asha.bump.v2,
mse.asha1.bump.v2,
mse.asha2.bump.v2
)

names(mse.bump.v2)=c(
"current",
"RSS both",
"RSS final",
"RSS weighted"
)


mse.bump.v3=c(
mse.ash.bump.v3,
mse.asha.bump.v3,
mse.asha1.bump.v3,
mse.asha2.bump.v3
)

names(mse.bump.v3)=c(
"current",
"RSS both",
"RSS final",
"RSS weighted"
)



mse.bump.v4=c(
mse.ash.bump.v4,
mse.asha.bump.v4,
mse.asha1.bump.v4,
mse.asha2.bump.v4
)

names(mse.bump.v4)=c(
"current",
"RSS both",
"RSS final",
"RSS weighted"
)






mse.blk.v1=c(
mse.ash.blk.v1,
mse.asha.blk.v1,
mse.asha1.blk.v1,
mse.asha2.blk.v1
)


names(mse.blk.v1)=c(
"current",
"RSS both",
"RSS final",
"RSS weighted"
)



mse.blk.v2=c(
mse.ash.blk.v2,
mse.asha.blk.v2,
mse.asha1.blk.v2,
mse.asha2.blk.v2
)

names(mse.blk.v2)=c(
"current",
"RSS both",
"RSS final",
"RSS weighted"
)


mse.blk.v3=c(
mse.ash.blk.v3,
mse.asha.blk.v3,
mse.asha1.blk.v3,
mse.asha2.blk.v3
)

names(mse.blk.v3)=c(
"current",
"RSS both",
"RSS final",
"RSS weighted"
)



mse.blk.v4=c(
mse.ash.blk.v4,
mse.asha.blk.v4,
mse.asha1.blk.v4,
mse.asha2.blk.v4
)

names(mse.blk.v4)=c(
"current",
"RSS both",
"RSS final",
"RSS weighted"
)






mse.dop.v1=c(
mse.ash.dop.v1,
mse.asha.dop.v1,
mse.asha1.dop.v1,
mse.asha2.dop.v1
)

names(mse.dop.v1)=c(
"current",
"RSS both",
"RSS final",
"RSS weighted"
)



mse.dop.v2=c(
mse.ash.dop.v2,
mse.asha.dop.v2,
mse.asha1.dop.v2,
mse.asha2.dop.v2
)

names(mse.dop.v2)=c(
"current",
"RSS both",
"RSS final",
"RSS weighted"
)


mse.dop.v3=c(
mse.ash.dop.v3,
mse.asha.dop.v3,
mse.asha1.dop.v3,
mse.asha2.dop.v3
)

names(mse.dop.v3)=c(
"current",
"RSS both",
"RSS final",
"RSS weighted"
)



mse.dop.v4=c(
mse.ash.dop.v4,
mse.asha.dop.v4,
mse.asha1.dop.v4,
mse.asha2.dop.v4
)

names(mse.dop.v4)=c(
"current",
"RSS both",
"RSS final",
"RSS weighted"
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


