source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/bayesmooth.R"))
source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/bayesmooth_alt.R"))
source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/bayesmooth_alt1.R"))
source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/bayesmooth_est_var.R"))

mu.mad=0
var.mad=0
mu.mad.alt=0
var.mad.alt=0
mu.mad.alt1=0
var.mad.alt1=0
mu.mad.ev=0
var.mad.ev=0

for(i in 1:400){
x=sort(runif(200,-2,2))
err=rnorm(200)
sigma.t=0.4*exp(-2*x^2)+0.2
mu.t=x+2*exp(-16*x^2)
y=mu.t+sigma.t*err

xgrid=seq(-1.8,1.8,length.out=101)
mu.t.inter=xgrid+2*exp(-16*xgrid^2)
var.t.inter=(0.4*exp(-2*xgrid^2)+0.2)^2

y.exp=c(y,y[200:145])
y.data=c(y.exp,y.exp[256:1])

mu.est=bayesmooth(y.data)
mu.est=mu.est[1:200]
var.est=bayesmooth(y.data,v.est=TRUE)
var.est=var.est[1:200]

mu.est.inter=approx(x,mu.est,xgrid,'linear')$y
var.est.inter=approx(x,var.est,xgrid,'linear')$y

mu.mad[i]=1/101*sum(abs(mu.est.inter-mu.t.inter))
var.mad[i]=1/101*sum(abs(var.est.inter-var.t.inter))

mu.est=bayesmooth.alt(y.data)
mu.est=mu.est[1:200]
var.est=bayesmooth.alt(y.data,v.est=TRUE)
var.est=var.est[1:200]

mu.est.inter=approx(x,mu.est,xgrid,'linear')$y
var.est.inter=approx(x,var.est,xgrid,'linear')$y

mu.mad.alt[i]=1/101*sum(abs(mu.est.inter-mu.t.inter))
var.mad.alt[i]=1/101*sum(abs(var.est.inter-var.t.inter))

mu.est=bayesmooth.alt1(y.data)
mu.est=mu.est[1:200]
var.est=bayesmooth.alt1(y.data,v.est=TRUE)
var.est=var.est[1:200]

mu.est.inter=approx(x,mu.est,xgrid,'linear')$y
var.est.inter=approx(x,var.est,xgrid,'linear')$y

mu.mad.alt1[i]=1/101*sum(abs(mu.est.inter-mu.t.inter))
var.mad.alt1[i]=1/101*sum(abs(var.est.inter-var.t.inter))


mu.est=bayesmooth.est.var(y.data)
mu.est=mu.est[1:200]
var.est=bayesmooth.est.var(y.data,v.est=TRUE)
var.est=var.est[1:200]

mu.est.inter=approx(x,mu.est,xgrid,'linear')$y
var.est.inter=approx(x,var.est,xgrid,'linear')$y

mu.mad.ev[i]=1/101*sum(abs(mu.est.inter-mu.t.inter))
var.mad.ev[i]=1/101*sum(abs(var.est.inter-var.t.inter))

print(i)
}




mu.mad.equal=0
var.mad.equal=0
mu.mad.equal.alt=0
var.mad.equal.alt=0
mu.mad.equal.alt1=0
var.mad.equal.alt1=0
mu.mad.equal.ev=0
var.mad.equal.ev=0

for(i in 1:400){
x=seq(-2,2,length.out=256)
err=rnorm(256)
sigma.t=0.4*exp(-2*x^2)+0.2
mu.t=x+2*exp(-16*x^2)
y.data=mu.t+sigma.t*err

xgrid=seq(-1.8,1.8,length.out=101)
mu.t.inter=xgrid+2*exp(-16*xgrid^2)
var.t.inter=(0.4*exp(-2*xgrid^2)+0.2)^2


mu.est=bayesmooth(y.data)
var.est=bayesmooth(y.data,v.est=TRUE)

mu.est.inter=approx(x,mu.est,xgrid,'linear')$y
var.est.inter=approx(x,var.est,xgrid,'linear')$y

mu.mad.equal[i]=1/101*sum(abs(mu.est.inter-mu.t.inter))
var.mad.equal[i]=1/101*sum(abs(var.est.inter-var.t.inter))

mu.est=bayesmooth.alt(y.data)
var.est=bayesmooth.alt(y.data,v.est=TRUE)

mu.est.inter=approx(x,mu.est,xgrid,'linear')$y
var.est.inter=approx(x,var.est,xgrid,'linear')$y

mu.mad.equal.alt[i]=1/101*sum(abs(mu.est.inter-mu.t.inter))
var.mad.equal.alt[i]=1/101*sum(abs(var.est.inter-var.t.inter))

mu.est=bayesmooth.alt1(y.data)
var.est=bayesmooth.alt1(y.data,v.est=TRUE)

mu.est.inter=approx(x,mu.est,xgrid,'linear')$y
var.est.inter=approx(x,var.est,xgrid,'linear')$y

mu.mad.equal.alt1[i]=1/101*sum(abs(mu.est.inter-mu.t.inter))
var.mad.equal.alt1[i]=1/101*sum(abs(var.est.inter-var.t.inter))

mu.est=bayesmooth.est.var(y.data)
var.est=bayesmooth.est.var(y.data,v.est=TRUE)

mu.est.inter=approx(x,mu.est,xgrid,'linear')$y
var.est.inter=approx(x,var.est,xgrid,'linear')$y

mu.mad.equal.ev[i]=1/101*sum(abs(mu.est.inter-mu.t.inter))
var.mad.equal.ev[i]=1/101*sum(abs(var.est.inter-var.t.inter))

print(i)

}


plot(x,mu.t,type='l')
lines(x,mu.est[1:200],col=2)

plot(x,sigma.t^2,type='l')
lines(x,var.est[1:200],col=2)
