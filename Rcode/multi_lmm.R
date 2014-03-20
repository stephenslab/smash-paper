multi.lmm=function(x,vm,g,center=FALSE,repara=FALSE,rand.eff=TRUE,disp="add"){
  if(is.factor(g)){g=as.numeric(levels(g))[g]}   #if g is a 2-level factor convert to numeric
  n=ncol(x)
  ng=length(g)
  df=ng-2
  w1=1/vm
  if(center==TRUE){g=g-mean(g)}
  w1s=colSums(w1)                                #define sum of weights for each of the n models (to be used later)
  gwmean=colSums(w1*g)/w1s                       #define weighted center of g for each of the n models (to be used later)
  ggwmeanm=g%o%rep(1,n)-rep(1,ng)%o%gwmean       #define weighted difference between each g and its weighted center for each of the n models (to be used later)
  betahat=colSums(w1*x*ggwmeanm)/colSums(w1*ggwmeanm^2)         #compute betahat using formula in documentation
  g.betahat=g%o%betahat                          #compute betahat*g	for each of the n models
  muhat=colSums(w1*(x-g.betahat))/w1s           #compute muhat using formula in documentation
  wrse=sqrt(colSums((x-rep(1,ng)%o%muhat-g.betahat)^2*w1)/df)   #compute residual standard error
  wrse[is.na(wrse)]=1
  if(rand.eff==FALSE|ng==length(unique(g))){     #force dispersion to be absent (also in the case with only 2 observations)
    wrse=1
  }else{           
    wrse[(wrse==Inf)|(wrse<1)]=1          #do not allow for "underdispersion"
  }
  wrse2=wrse^2
  wgg=colSums(w1*ggwmeanm^2)                     #compute sum_j w_j^2*(g_j-gwmean)^2 (to be used later)
  wgg.ind=wgg<1e-6                               
  wgg[wgg.ind]=0                                 #change small values of wgg to 0 for numeric reasons
  sebetahat=sqrt(wrse2/wgg)                      #compute se(betahat) using formula in documentation
  sebetahat[wgg.ind]=NA 
  semuhat=sqrt((1/w1s+gwmean^2/wgg)*wrse2)       #compute se(muhat) using formula in documentation
  semuhat[wgg.ind]=NA
  covmubeta=colSums(w1*ggwmeanm)/w1s/wgg*wrse2-gwmean*sebetahat^2   #compute covariance between muhat and betahat
  if(disp=="mult"){                       #return estimates if multiplicative dispersion is assumed
    coef=c(muhat,betahat)
    se=c(semuhat,sebetahat)
    if(repara==TRUE){                            #return reparametrized muhat and behat as well as their SEs, together with gamma as defined in documentation
      mbvar=covmubeta/sebetahat^2                
      coef[1:n]=muhat-betahat*mbvar              
      se[1:n]=sqrt(semuhat^2+mbvar^2*sebetahat^2-2*mbvar*covmubeta)
    }
  }else{
    vv=pmax((wrse2-1)*colMeans(vm),0)     #computes crude estimate of sigma_u^2 as in documentation
    ww=1/(vm+rep(1,ng)%o%vv)                    #recompute weights with estimate of sigma_u^2
    wws=colSums(ww)  
    gwmean=colSums(ww*g)/wws                     #define weighted center of g for each of the n models (to be used later)
    ggwmeanm=g%o%rep(1,n)-rep(1,ng)%o%gwmean     #define weighted difference between each g and its weighted center for each of the n models (to be used later)
    betahat=colSums(ww*x*ggwmeanm)/colSums(ww*ggwmeanm^2)         #compute betahat using formula in documentation
    g.betahat=g%o%betahat                        #compute betahat*g	for each of the n models
    muhat=colSums(ww*(x-g.betahat))/wws         #compute muhat using formula in documentation
    wrse=sqrt(colSums((x-rep(1,ng)%o%muhat-g.betahat)^2*ww)/df)   #compute residual standard error    
    if(rand.eff==FALSE|ng==length(unique(g))){   #force dispersion to be absent (also in the case with only 1 observation in each group)
      wrse=1
    }else{           
      wrse[(wrse==Inf)|(vv==0)]=1                #do not allow for "underdispersion"
    }
    wrse2=wrse^2
    wgg=colSums(ww*ggwmeanm^2)                   #compute sum_j w_j^2*(g_j-gwmean)^2 (to be used later)
    wgg.ind=wgg<1e-6      
    wgg[wgg.ind]=0                               #change small values of wgg to 0 for numeric reasons
    sebetahat=sqrt(wrse2/wgg)                    #compute se(betahat) using formula in documentation
    sebetahat[wgg.ind]=NA
    semuhat=sqrt((1/wws+gwmean^2/wgg)*wrse2)
    semuhat[wgg.ind]=NA
    covmubeta=colSums(ww*ggwmeanm)/wws/wgg*wrse2-gwmean*sebetahat^2       #compute covariance between muhat and betahat
    coef=c(muhat,betahat)
    se=c(semuhat,sebetahat)
    if(repara==TRUE){                            #return reparametrized muhat and behat as well as their SEs, together with gamma as defined in documentation
      mbvar=covmubeta/sebetahat^2
      coef[1:n]=muhat-betahat*mbvar
      se[1:n]=sqrt(semuhat^2+mbvar^2*sebetahat^2-2*mbvar*covmubeta)
    }
  }
  na.ind=is.na(coef[1:n])|is.na(se[1:n])|is.na(coef[(n+1):(2*n)])|is.na(se[(n+1):(2*n)])
  na.ind2=rep(na.ind,times=2)         
  coef[na.ind2]=NA
  se[na.ind2]=NA
  if(repara==TRUE){mbvar[na.ind]=NA}
  out=array(rbind(coef,se),dim=c(2,n,2))
  if(repara==FALSE){
    return(apply(out,2,rbind))
    #return(list(muhat=muhat,betahat=betahat,semuhat=semuhat,sebetahat=sebetahat))
  }else{
    return(matrix(rbind(apply(out,2,rbind),mbvar),ncol=n))
    #return(mbvar)
  }
}