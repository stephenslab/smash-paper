library(Rcpp)
Rcpp.package.skeleton("smash", path=".", code_files=c("ashsmooth.R", "deltamethod.R", "glm_approx.R", "jash.R", "wd_var.R", "documentation.R"), cpp_files="ashsmooth.cpp", example_code = FALSE, attributes = TRUE)


require(roxygen2)
require(devtools)
roxygenize()
document()












install.packages("D:/Grad School/Spring 2013/ashwave/package/smash_1.0.tar.gz",repos=NULL,type="source")
library(smash)
?smash





















spike.f=function(x) (0.75*exp(-500*(x-0.23)^2)+1.5*exp(-2000*(x-0.33)^2)+3*exp(-8000*(x-0.47)^2)+2.25*exp(-16000*(x-0.69)^2)+0.5*exp(-32000*(x-0.83)^2))
n=2^8
t=1:n/n
mu.s=spike.f(t)

pos = c(.1, .13, .15, .23, .25, .40, .44, .65, .76, .78, .81)
hgt = 2.97/5*c(4, 5, 3, 4, 5, 4.2, 2.1, 4.3, 3.1, 5.1, 4.2)
wth = c(.005, .005, .006, .01, .01, .03, .01, .01, .005, .008, .005)
mu.b = rep(0,n)
for(j in 1:length(pos)){
  mu.b = mu.b + hgt[j]/(( 1 + (abs(t - pos[j])/wth[j]))^4)
}


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

mu.cblk=mu.blk
mu.cblk[mu.cblk<0]=0



mu.cor=623.87*t^3*(1-2*t)*(t>=0&t<=0.5)+187.161*(0.125-t^3)*t^4*(t>0.5&t<=0.8)+3708.470441*(t-1)^3*(t>0.8&t<=1)
mu.cor=(0.6/(max(mu.cor)-min(mu.cor)))*mu.cor
mu.cor=mu.cor-min(mu.cor)+0.2


mse=function(x,y) mean((x-y)^2)
l2norm=function(x) sum(x^2)
mise=function(x,y) 10000*mean(apply(x-rep(1,100)%o%y,1,l2norm)/l2norm(y))






mu.t=(1+mu.s)/5
#mu.t=(1+mu.b)/5
#mu.t=(1+mu.blk)/5
#mu.t=(1+mu.ang)/5
#mu.t=(1+mu.dop)/5

#mu.t=mu.blip
#mu.t=mu.cor

rsnr=sqrt(1)

var1=rep(1,n)
var2=(0.0001+4*(exp(-550*(t-0.2)^2)+exp(-200*(t-0.5)^2)+exp(-950*(t-0.8)^2)))/1.35
var4=3.4*(2+mu.dop)
var5=0.00001+1*(mu.cblk-min(mu.cblk))/max(mu.cblk)
sigma.ini=sqrt(var5)
sigma.t=sigma.ini/mean(sigma.ini)*sd(mu.t)/rsnr^2

#sigma.t=sd(mu.t)/rsnr^2

set.seed(377)
X.s=rnorm(n,mu.t,sigma.t)

system.time(mu.est.haar<-ashsmooth.gaus(X.s,filter.number=1,family="DaubExPhase"))


