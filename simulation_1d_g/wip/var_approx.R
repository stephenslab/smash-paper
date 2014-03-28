
yseq=seq(1,20,1)
lik=0
for(i in 1:length(yseq)){
  xx=(rnorm(100000,0,10)^2)
  yy=(rnorm(100000,0,yseq[i])^2)
  res=hist(yy-xx,breaks=100,plot=FALSE)
  midd=res$mid
  dens=res$density
  yfind=min(which(midd>=10))
  lik[i]=dens[yfind]
}
plot(yseq^2-10^2,lik,type='l')
dnorm(yseq^2-10^2,10,)


  xx=(rnorm(100000,0,sqrt(1))^2)
  yy=(rnorm(100000,0,sqrt(3))^2)

skewness(yy-xx)

hist(yy-xx,breaks=100,freq=FALSE)
lines(seq(-20,60,0.01),dnorm(seq(-20,60,0.01),2,sqrt((2/3*3^2+2/3*1^2))))


  xx=(rnorm(100000,0,sqrt(1))^2)
  xxx=(rnorm(100000,0,sqrt(1))^2)
  xxxx=(rnorm(100000,0,sqrt(1))^2)
  xxxxx=(rnorm(100000,0,sqrt(1))^2)
  xxxxxx=(rnorm(100000,0,sqrt(1))^2)
  xxxxxxx=(rnorm(100000,0,sqrt(1))^2)
  xxxxxxxx=(rnorm(100000,0,sqrt(1))^2)
  xxxxxxxxx=(rnorm(100000,0,sqrt(1))^2)
  xx1=(rnorm(100000,0,sqrt(1))^2)
  xxx1=(rnorm(100000,0,sqrt(1))^2)
  xxxx1=(rnorm(100000,0,sqrt(1))^2)
  xxxxx1=(rnorm(100000,0,sqrt(1))^2)
  xxxxxx1=(rnorm(100000,0,sqrt(1))^2)
  xxxxxxx1=(rnorm(100000,0,sqrt(1))^2)
  xxxxxxxx1=(rnorm(100000,0,sqrt(1))^2)
  xxxxxxxxx1=(rnorm(100000,0,sqrt(1))^2)
  yy=(rnorm(100000,0,sqrt(10))^2)
  yyy=(rnorm(100000,0,sqrt(10))^2)
  yyyy=(rnorm(100000,0,sqrt(10))^2)
  yyyyy=(rnorm(100000,0,sqrt(10))^2)
  yyyyyy=(rnorm(100000,0,sqrt(10))^2)
  yyyyyyy=(rnorm(100000,0,sqrt(10))^2)
  yyyyyyyy=(rnorm(100000,0,sqrt(10))^2)
  yyyyyyyyy=(rnorm(100000,0,sqrt(10))^2)
  yy1=(rnorm(100000,0,sqrt(10))^2)
  yyy1=(rnorm(100000,0,sqrt(10))^2)
  yyyy1=(rnorm(100000,0,sqrt(10))^2)
  yyyyy1=(rnorm(100000,0,sqrt(10))^2)
  yyyyyy1=(rnorm(100000,0,sqrt(10))^2)
  yyyyyyy1=(rnorm(100000,0,sqrt(10))^2)
  yyyyyyyy1=(rnorm(100000,0,sqrt(10))^2)
  yyyyyyyyy1=(rnorm(100000,0,sqrt(10))^2)


var(xx)
mean(2/3*xx^2)
var(yy)
mean(2/3*yy^2)
var(yy-xx)

hist(2/3*yy^2+2/3*xx^2,breaks=1000,xlim=c(0,50))



mean(2/3*yy^2+2/3*xx^2)/median(2/3*yy^2+2/3*xx^2)

mean(2/3*yy^2+2/3*xx^2+2/3*yyy^2+2/3*xxx^2)/median(2/3*yy^2+2/3*xx^2+2/3*yyy^2+2/3*xxx^2)

mean(2/3*yy^2+2/3*xx^2+2/3*yyy^2+2/3*xxx^2+2/3*yyyy^2+2/3*xxxx^2+2/3*yyyyy^2+2/3*xxxxx^2)/median(2/3*yy^2+2/3*xx^2+2/3*yyy^2+2/3*xxx^2+2/3*yyyy^2+2/3*xxxx^2+2/3*yyyyy^2+2/3*xxxxx^2)

mean(2/3*yy^2+2/3*xx^2+2/3*yyy^2+2/3*xxx^2+2/3*yyyy^2+2/3*xxxx^2+2/3*yyyyy^2+2/3*xxxxx^2+
2/3*yyyyyy^2+2/3*xxxxxx^2+2/3*yyyyyyy^2+2/3*xxxxxxx^2+2/3*yyyyyyyy^2+2/3*xxxxxxxx^2+2/3*yyyyyyyyy^2+2/3*xxxxxxxxx^2)/
median(2/3*yy^2+2/3*xx^2+2/3*yyy^2+2/3*xxx^2+2/3*yyyy^2+2/3*xxxx^2+2/3*yyyyy^2+2/3*xxxxx^2+
2/3*yyyyyy^2+2/3*xxxxxx^2+2/3*yyyyyyy^2+2/3*xxxxxxx^2+2/3*yyyyyyyy^2+2/3*xxxxxxxx^2+2/3*yyyyyyyyy^2+2/3*xxxxxxxxx^2)

mean(2/3*yy^2+2/3*xx^2+2/3*yyy^2+2/3*xxx^2+2/3*yyyy^2+2/3*xxxx^2+2/3*yyyyy^2+2/3*xxxxx^2+
2/3*yyyyyy^2+2/3*xxxxxx^2+2/3*yyyyyyy^2+2/3*xxxxxxx^2+2/3*yyyyyyyy^2+2/3*xxxxxxxx^2+2/3*yyyyyyyyy^2+2/3*xxxxxxxxx^2+
2/3*yy1^2+2/3*xx1^2+2/3*yyy1^2+2/3*xxx1^2+2/3*yyyy1^2+2/3*xxxx1^2+2/3*yyyyy1^2+2/3*xxxxx1^2+
2/3*yyyyyy1^2+2/3*xxxxxx1^2+2/3*yyyyyyy1^2+2/3*xxxxxxx1^2+2/3*yyyyyyyy1^2+2/3*xxxxxxxx1^2+2/3*yyyyyyyyy1^2+2/3*xxxxxxxxx1^2)/
median(2/3*yy^2+2/3*xx^2+2/3*yyy^2+2/3*xxx^2+2/3*yyyy^2+2/3*xxxx^2+2/3*yyyyy^2+2/3*xxxxx^2+
2/3*yyyyyy^2+2/3*xxxxxx^2+2/3*yyyyyyy^2+2/3*xxxxxxx^2+2/3*yyyyyyyy^2+2/3*xxxxxxxx^2+2/3*yyyyyyyyy^2+2/3*xxxxxxxxx^2+
2/3*yy1^2+2/3*xx1^2+2/3*yyy1^2+2/3*xxx1^2+2/3*yyyy1^2+2/3*xxxx1^2+2/3*yyyyy1^2+2/3*xxxxx1^2+
2/3*yyyyyy1^2+2/3*xxxxxx1^2+2/3*yyyyyyy1^2+2/3*xxxxxxx1^2+2/3*yyyyyyyy1^2+2/3*xxxxxxxx1^2+2/3*yyyyyyyyy1^2+2/3*xxxxxxxxx1^2)



