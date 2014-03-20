multiseq.gaus=function(x,g,center=FALSE,repara=FALSE,rand.eff=TRUE,disp="add",pointmass=TRUE,prior="nullbiased",nullcheck=TRUE,gridmult=2,mixsd=NULL,VB=FALSE){
  if(!is.null(g)&is.factor(g)&center==TRUE){    #considers different cases for covariate g
    g.num = as.numeric(levels(g))[g]
    w = unique(sort(g.num-mean(g.num)))
  }else if(!is.null(g)&is.factor(g)&center==FALSE){
    w = c(0,1)
  }else if(!is.null(g)&!is.factor(g)&center==TRUE){  
    w = unique(sort(g-mean(g)))
  }else if(!is.null(g)&!is.factor(g)&center==FALSE){
    #weights for quantitative covariate
    w = c(0,1)
  }else{
    w = NA
  }
  
  n=dim(x)[2]
  J=log2(n)
  nsig=dim(x)[1]
  var.est=t(apply(x,1,bayesmooth,v.est=TRUE,pointmass=pointmass,prior=prior,nullcheck=nullcheck,gridmult=gridmult,VB=VB))
  var.est[var.est<=0]=1e-6
  wc=matrix(0,nrow=nsig,ncol=J*n)
  wc.var=matrix(0,nrow=nsig,ncol=J*n)
  for(i in 1:nsig){
    wc[i,]=as.vector(t((cxxtitable(x[i,])$difftable)[2:(J+1),]))/2
    wc.var[i,]=as.vector(t((cxxtitable(var.est[i,])$sumtable)[2:(J+1),]))/4
  }

    x.sum = matrix(rowSums(x),ncol=1)
    v.sum = matrix(rowSums(var.est),ncol=1)
    zdat.rate = as.vector(multi.lmm(x.sum,v.sum,g=g,center=center,repara=repara,rand.eff=rand.eff,disp=disp))
    alpha.m = zdat.rate[1]
    alpha.v = zdat.rate[2]^2
    beta.m = zdat.rate[3]
    beta.v = zdat.rate[4]^2
    if(repara==TRUE){
      mbvar = zdat.rate[5]
    }else{
      mbvar=0
    }
    b.m = alpha.m+(w[1]+mbvar)*beta.m
    #b.v = alpha.v+(w[1]+mbvar)^2*beta.v
    b.v = 0
    e.m = beta.m
    #e.v = beta.v
    e.v = 0

  bmean=matrix(0,nrow=J,ncol=n)
  bvar=matrix(0,nrow=J,ncol=n)
  emean=matrix(0,nrow=J,ncol=n)
  evar=matrix(0,nrow=J,ncol=n)
  zdat=multi.lmm(wc,wc.var,g,center=center,repara=repara,rand.eff=rand.eff,disp=disp)
  for(j in 1:J){
    ind = ((j-1)*n+1):(j*n)  
    #apply ash to vector of intercept estimates and SEs
    zdat.ash=fast.ash(zdat[1,ind],zdat[2,ind],prior=prior,pointmass=pointmass,nullcheck=nullcheck,gridmult=gridmult,mixsd=NULL,VB=VB)
    alpha.mv=list(mean=zdat.ash$PosteriorMean,var=zdat.ash$PosteriorSD^2) #find mean and variance of alpha
    zdat.ash=fast.ash(zdat[3,ind],zdat[4,ind],prior=prior,pointmass=pointmass,nullcheck=nullcheck,gridmult=gridmult,mixsd=NULL,VB=VB)
    beta.mv=list(mean=zdat.ash$PosteriorMean,var=zdat.ash$PosteriorSD^2) #find mean and variance of effect beta
    if(repara==TRUE){  #if reparametrization is used then we want gamma returned as well
      mbvar.ind=is.na(zdat[5,ind])
      mbvar=zdat[5,ind]
      mbvar[mbvar.ind]=0
    }else{
      mbvar=0
    }
    gamma.mv=list(mean=alpha.mv$mean+(w[1]+mbvar)*beta.mv$mean,var=alpha.mv$var+(w[1]+mbvar)^2*beta.mv$var)  #define baseline estimate by accounting for weights and gamma
    bmean[j,]=gamma.mv$mean
    bvar[j,]=gamma.mv$var

    emean[j,]=beta.mv$mean
    evar[j,]=beta.mv$var
  }
  brmean=-bmean
  brvar=bvar
  ermean=-emean
  ervar=evar

  baseline.mean=cxxreverse.gwave(b.m,bmean,brmean)
  baseline.var=cxxreverse.gvwave(b.v,bvar,brvar)
  effect.mean=cxxreverse.gwave(e.m,emean,ermean)
  effect.var=cxxreverse.gvwave(e.v,evar,ervar)

  return(list(baseline.mean=baseline.mean, baseline.var=baseline.var, effect.mean=effect.mean, effect.var=effect.var))  
}