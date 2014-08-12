library(mixtools)

vv=psigamma(0.5,1)
J=12

lambda.mix=list(0)
sigma.mix=0
scale.mix=list(0)
ncomp=c(rep(3,2),rep(2,6))

size=2^(0:(J-1))
for(k in 1:J){
  K=size[k]
  x=matrix(0,K,1000000)
  y=matrix(0,K,1000000)
  for(i in 1:K){
    x[i,]=rchisq(1000000,1)
    y[i,]=rchisq(1000000,1) 
  }
  normmix<-normalmixEM(apply(log(x)-log(y),2,sum),mu=0,k=ncomp[k],maxit=10000,eps=1e-7)
  lambda.mix[[k]]=normmix$lambda
  sigma.mix[k]=normmix$sigma
  scale.mix[[k]]=normmix$scale
}


save(lambda.mix,sigma.mix,scale.mix,file="norm_mix.Robj")

