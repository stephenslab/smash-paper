source(file.path("D:/Grad School/Spring 2013/multiscale_ash/ash/bayesmooth.R"))
mse=function(x,y) mean((x-y)^2)

n=200
mse.hsm=0
for(i in 1:500){
xt=sort(rnorm(n,0.5,0.2))
xt=(xt-min(xt))/(max(xt)-min(xt))

mu.hsm=4 * sin(4 * pi * xt) - 2*sign(xt - 0.3) - 2*sign(0.72 - xt)

mu.t=mu.hsm
var1=rep(1,n)
var2=(xt<=0.5)+2*(xt>0.5)
sigma.ini=sqrt(var2)
rsnr=sqrt(4)
sigma.t=sigma.ini/mean(sigma.ini)*sd(mu.t)/rsnr^2

y=rnorm(n,mu.t,sigma.t)
y.exp=c(y,y[200:145])
y.final=c(y.exp,y.exp[256:1])

mu.est<-bayesmooth(y.final)
mu.est=mu.est[1:n]
mse.hsm[i]=mse(mu.est,mu.t)

print(i)
}

sqrt(median(mse.hsm))


mu.cor=function(t) 623.87*t^3*(1-2*t)*(t>=0&t<=0.5)+187.161*(0.125-t^3)*t^4*(t>0.5&t<=0.8)+3708.470441*(t-1)^3*(t>0.8&t<=1)



plot(xt,mu.t,type='l')
lines(xt,mu.est,col=2)