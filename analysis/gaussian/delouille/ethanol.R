library(smash)
library(lattice)

x.ini = sort(ethanol$E)
y.ini = ethanol$NOx[order(ethanol$E)]

x = unique(x.ini)
y = 0
for(i in 1:length(x)){
  y[i] = median(y.ini[x.ini == x[i]])
}


y.exp=c(y,y[length(y):(2*length(y)-128+1)])
y.final = c(y.exp, y.exp[length(y.exp):1])


y.var.est=ashsmooth.gaus(y.final,v.est=TRUE)

y.var.est=y.var.est[1:length(y)]

plot(x,y,ylim=c(-1,5))
lines(x,y.est,col=2)
plot(x,y.var.est,type='l')



