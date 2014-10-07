require(ashr)
require(Rcpp)
require(inline)
source(file.path("~/ashwave/Rcode/glm_approx.R"))
source(file.path("~/ashwave/Rcode/deltamethod.R"))

#interleave two vectors
interleave=function(x,y){
  return(as.vector(rbind(x,y)))
}

#shift a vector right and left respectively
rshift = function(x){L=length(x); return(c(x[L],x[-L]))}
lshift = function(x){return(c(x[-1],x[1]))}



#The following produces both a TItable,
#and a "parent" table whose pairwise comparisons would be used to create a TI table
#So, for example, in the ith row, elements 1, 2 would be the parents of the
#first element in the (i+1)the row of the TI table
# INPUT: sig, an n vector of Poisson counts at n locations
# OUTPUT: a list, with elements
# TItable - the usual TI table
# parent - the parent values used to make the TI table
ParentTItable=function(sig){
  n = length(sig)
  J = log2(n)
  
# Create decomposition table of signal, using pairwise sums,
# keeping just the values that are *not* redundant under the
# shift-invariant scheme.  This is very similar to TI-tables
# in Donoho and Coifman's TI-denoising framework.
  dmat = matrix(0, nrow=J+1, ncol=n)
  dmat[1,] = sig
  #dmat[1,] = as.matrix(sig)  
  dmat2 = matrix(0, nrow=J, ncol=2*n) #the parent table
  
  for(D in 0:(J-1)){
    nD = 2^(J-D); 
    nDo2 = nD/2;
    twonD = 2*nD;
    for(l in 0:(2^D-1)){
      ind = (l*nD+1):((l+1)*nD)
      ind2 = (l*twonD+1):((l+1)*twonD)
      x = dmat[D+1,ind]
      lsumx = x[seq(from=1,to=nD-1, by=2)] + x[seq(from=2,to=nD,by=2)]
      rx = rshift(x);
      rsumx = rx[seq(from=1,to=nD-1, by=2)] + rx[seq(from=2,to=nD,by=2)]
      dmat[D+2,ind] = c(lsumx,rsumx)
      dmat2[D+1,ind2] = c(x,rx)
    }
  }
  return(list(TItable=dmat,parent=dmat2))
}




#cxxParentTItable
#cxxParentTItable implements ParentTItable in C++
#note that input to cxxParentTItable is an nsig by n matrix sig
#while input to ParentTItable is a 1 by n vector instead
#src is a string containing the C++ code
src1 <- '
        NumericVector signal=sig; 
        int n=(int) signal.size();
        int J=(int) log2((double)n);

        NumericVector parent(2*J*n);
        NumericMatrix TItable(J+1,n);
        TItable(0,_) = signal;
        for (int D=0; D<J; D++){
           int nD=(int) pow(2., (int) (J-D)), pD=(int) pow(2.,(int) D);
           for (int l=0; l<pD; l++){
              int a=l*nD+1, b=2*l*nD+1, c=2*D*n+b, d;
              for (int i=0; i<nD-1; i++){
                 d=TItable(D,a+i-1);
                 parent(c+i-1)=d;
                 parent(c+i+nD)=d;
              }
              //i=nD-1
              d=TItable(D,a+nD-2);
              parent(c+nD-2)=d;
              parent(c+nD-1)=d;

              for (int i=0; i<nD; i++)
                TItable(D+1,a+i-1)=parent(c+2*i-1)+parent(c+2*i);
          }
        }
        return(parent);
        '
cxxSParentTItable <- cxxfunction(signature(sig="numeric"),
                                body=src1,
                                plugin="Rcpp",
                                inc="#include <cmath>")

# src2 <- '
#         NumericVector signal=sig; 
#         int n=(int) signal.size();
#         int J=(int) log2((double)n);
#         
#         NumericVector parent(2*J*n);
#         NumericMatrix TItable(J+1,n);
#         TItable(0,_) = signal;
#         for (int D=0; D<J; D++){
#           int nD=(int) pow(2., (int) (J-D)), pD=(int) pow(2.,(int) D);
#           for (int l=0; l<pD; l++){
#             int a=l*nD+1, b=2*l*nD+1, c=2*D*n+b, d;
#             for (int i=0; i<nD-1; i++){
#               d=TItable(D,a+i-1);
#               parent(c+i-1)=d;
#               parent(c+i+nD)=d;
#             }
#             //i=nD-1
#             d=TItable(D,a+nD-2);
#             parent(c+nD-2)=d;
#             parent(c+nD-1)=d;
#             
#             
#             for (int i=0; i<nD; i++)
#               TItable(D+1,a+i-1)=parent(c+2*i-1)+parent(c+2*i);
#           }
#         }
#         return(parent);
#         '
# cxxSParentTItable <- cxxfunction(signature(sig="numeric"),
#                                 body=src2,
#                                 plugin="Rcpp",
#                                 inc="#include <cmath>")




# function to turn a parent TItable into something that one can apply glm to
# INPUT: dmat, a k by n matrix
# OUTPUT: a 2 by nk/2 matrix
simplify=function(dmat){
  matrix(t(dmat),nrow=2)
}

#reverse wavelet transform a set of probabilities in TItable format
#lp is a set of estimated wavelet proportions for Poisson data, on the log scale.
#That is, lp is an estimate of log(p) 
#First row of lp gives the high frequency proportions (so 1/(1+2), 3/(3+4) etc)
#second row gives next level of resolution etc
#lq is an estimate of log(1-p). If lq is not given, set lq = log(1-exp(lp))
#INPUT: est, an n-vector. Usually a constant vector with each element equal to the estimated log(total rate)
# lp: a J by n matrix of estimated log(p) 
# lq: a J by n matrix of estimated log(1-p)
reverse.pwave=function(est,lp,lq=NULL){
  if(is.null(lq)){
    lq = log(1-exp(lp))
  }
  if(length(est)==1){est = rep(est,ncol(lp))}
  
  J=nrow(lp)

  for(D in J:1){
    #print(exp(est))
    #readline("press a key")
    nD = 2^(J-D+1) 
    nDo2 = nD/2
    for(l in 0:(2^(D-1)-1)){
      # Set indexing so as to pick off blocks of size 2^(J-D+1)
      # when shrinking estimates at depth D+1 down to finer
      # scale at depth D.
      ind = (l*nD+1):((l+1)*nD)

      estvec = est[ind]
      lpvec = lp[D,ind]
      lqvec = lq[D,ind]
      # In the first half of the vector of D+1-depth estimates,
      # we can shrink using the D-depth counts in the order
      # in which they appear.
      estl = estvec[1:nDo2]
      lpl=lpvec[1:nDo2]
      lql=lqvec[1:nDo2]
      nestl = interleave(estl+lpl,estl+lql) #interleaves the two
    
      # In the second half of the vector of D+1-depth counts,
      # the shrunken values are for the right shifted vector so these
      # have to be left shifted before averaging
    
      estr = estvec[(nDo2+1):nD]
      lpr = lpvec[(nDo2+1):nD]
      lqr = lqvec[(nDo2+1):nD]
      nestr = interleave(estr+lpr,estr+lqr) #interleaves the two
      nestr = lshift(nestr)
    
      # Combine the estimates from both halves of the D+1-depth
      # counts, and store.
      est[ind] = 0.5*( nestl + nestr )
    }
  }
  return(est)
}


#reverse.pwave, as implemented in C++
src2 <- '
        NumericMatrix pp=pmat;
        NumericMatrix qq=qmat;
        NumericVector est1=estimate;
        int np=(int) pp.ncol();
        int J=(int) pp.nrow();
        NumericVector est(np,est1(0));
        //for(int D=J; D-->0;){
        for(int D=0; D<J; D++){
          //int nD=pow(2., (int) (J-D)), pD=pow(2., (int) D);
          int nD=pow(2., (int) (D+1)), pD=pow(2., (int) (J-1-D));
          int nDo2=nD/2;
          NumericVector tempvecl(nD), tempvecr(nD);
          for(int l=0; l<pD; l++){
            int a=l*nD+1; 
            double dep, deq, dp, dq;
            for (int i=0; i<nDo2; i++){
              dep=est(a+i-1);
              //dp=pp(D,a+i-1);
              dp=pp(J-1-D,a+i-1);
              //dq=qq(D,a+i-1);
              dq=qq(J-1-D,a+i-1);
              tempvecl(2*i)=dep+dp;
              tempvecl(2*i+1)=dep+dq;
            }
            for (int i=nDo2; i<nD-1; i++){
              dep=est(a+i);
              //dp=pp(D,a+i);
              dp=pp(J-1-D,a+i);
              deq=est(a+i-1);
              //dq=qq(D,a+i-1);
              dq=qq(J-1-D,a+i-1);
              tempvecr(2*(i-nDo2))=deq+dq;
              tempvecr(2*(i-nDo2)+1)=dep+dp;
            }         
            //i=nD-1
            dep=est(a+nDo2-1);
            //dp=pp(D,a+nDo2-1);
            dp=pp(J-1-D,a+nDo2-1);
            deq=est(a+nD-2); 
            //dq=qq(D,a+nD-2);
            dq=qq(J-1-D,a+nD-2);
            tempvecr(nD-2)=deq+dq;
            tempvecr(nD-1)=dep+dp;
            for(int i=0; i<nD; i++)
               est(a+i-1)=0.5*(tempvecl(i)+tempvecr(i));
          }
        }
        return(est);
        '
cxxreverse_pwave <- cxxfunction(signature(estimate="numeric",pmat="numeric",qmat="numeric"),
                                body=src2,
                                plugin="Rcpp",
                                inc="#include <cmath>")


#take a univariate inhomogeneous Poisson process x and smooths it
#length of x must be power of 2
#INPUT: x a vector of nsig by n=2^J Poisson counts
#reflect: boolean, indicates if the signals should be reflected; otherwise periodicity is assumed
#lev: integer from 0 to J-1, coarsest level of resolution
#log: boolean, determines if smoothed signal is returned on log scale or not
#pseudocounts, all: variables for passing to glm calculations
#pseudocounts: boolean, a number to be added to counts
#all is a parameter indicating if pseudocounts should be added too all entries or only cases when either number of successes or number of failures (but not both) is 0
#prior, nullcheck, usePointMass, mixsd, VB: variables for passing to ash
#prior: "nullbiased" or "uniform", type of prior used as part of EM algorithm in ash
#nullcheck: boolean, check that any fitted model exceeds the "null" likelihood
#usePointMass: bool, indicating whether to use a point mass at zero as one of components for a mixture distribution
#mixsd: vector of sigma components to be specified for mixture model; defaults to NULL, in which case an automatic procedure is used
#VB: boolean, indicates if a variational bayes alternative should be used to EM
#return.est: boolean, indicates if just estimate (without variance) should be returned


binnorm.smooth = function(x,reflect=FALSE,lev=0,log=FALSE,pseudocounts=0.5,all=FALSE,pointmass=TRUE,prior="nullbiased",gridmult=2,nullcheck=TRUE,mixsd=NULL,VB=FALSE,return.est=TRUE){
  if(is.matrix(x)){
    if(nrow(x)==1){  #change matrix x to vector
      x=as.vector(x)
    }else{
      stop("x cannot have multiple rows")
    }
  }
  
  n = length(x)
  J = log2(n)
  if(!isTRUE(all.equal(J,trunc(J)))){stop("Error: length of x must be a power of 2")}
  if(reflect==TRUE){
    x.mir=x[n:1]
    x=c(x,x.mir)
    sig.ind=1:n
  }
  n = length(x)
  J = log2(n)
  
  #create the parent TI table for each signal, and put into rows of matrix y
  ls=sum(x)
  
  #y = as.vector(t(ParentTItable(x)$parent))
  y = cxxSParentTItable(x)
  
  zdat=glm.approx(y,g=NULL,minobs=1,pseudocounts=pseudocounts,center=FALSE,all=all,forcebin=TRUE,repara=TRUE,lm.approx=FALSE,disp="add")
  #define empty matrices for posterior means and variances of log(p) and log(q)
  lp.mean = matrix(0,ncol=n,nrow=J) # estimated E(log(p)) where p is the prob of going left
  lp.var = matrix(0,ncol=n,nrow=J) # estimated Var(log(p)) 
  lq.mean = matrix(0,ncol=n,nrow=J) # estimated E(log(q)) where q=1-p is the prob of going right
  lq.var = matrix(0,ncol=n,nrow=J) # estimated Var(log(q))
  
  # loop through resolutions, smoothing each resolution separately
  for(j in 1:(J-lev)){
    ind = ((j-1)*n+1):(j*n)
    #apply ash to vector of intercept estimates and SEs
    zdat.ash=ash(zdat[1,ind],zdat[2,ind], prior=prior, multiseqoutput=TRUE, pointmass=pointmass, nullcheck=nullcheck, gridmult=gridmult, mixsd=mixsd, VB=VB)
    alpha.mv=list(mean=zdat.ash$PosteriorMean,var=zdat.ash$PosteriorSD^2) #find mean and variance of alpha
    #if we want estimates on the log scale.
    if(log==TRUE){
      lp.mv = ff.moments(alpha.mv$mean, alpha.mv$var)  #find mean and variance of log(p)
      lp.mean[j,] = lp.mv$mean
      lp.var[j,]= lp.mv$var
      
      lq.mv = ff.moments(-alpha.mv$mean, alpha.mv$var)  #find mean and variance of log(q)  
      lq.mean[j,] = lq.mv$mean
      lq.var[j,]= lq.mv$var 
    }else{  #returns estimates on the original scale (non-log).
      lp.mv = ff.moments_exp(alpha.mv$mean, alpha.mv$var)  #find mean and meansq of p
      lp.mean[j,] = lp.mv$mean
      lp.var[j,]= lp.mv$meansq
      
      lq.mv = ff.moments_exp(-alpha.mv$mean, alpha.mv$var)  #find mean and meansq of q  
      lq.mean[j,] = lq.mv$mean
      lq.var[j,]= lq.mv$meansq
    }
  }
  #Do not smooth for coarser levels, with everything the same as above but using the estimate and its variance as the posterior mean and variance ie flat prior   
  if(lev!=0){
    for(j in (J-lev+1):J){
      ind = ((j-1)*n+1):(j*n)
      alpha.mv=list(mean=fill.nas(zdat[1,ind]),var=fill.nas(zdat[2,ind])^2) #find mean and variance of alpha
      if(log==TRUE){
        lp.mv = ff.moments(alpha.mv$mean, alpha.mv$var)  #find mean and variance of log(p)
        lp.mean[j,] = lp.mv$mean
        lp.var[j,]= lp.mv$var
        
        lq.mv = ff.moments(-alpha.mv$mean, alpha.mv$var)  #find mean and variance of log(q)  
        lq.mean[j,] = lq.mv$mean
        lq.var[j,]= lq.mv$var
      }else{
        lp.mv = ff.moments_exp(alpha.mv$mean, alpha.mv$var)  #find mean and meansq of p
        lp.mean[j,] = lp.mv$mean
        lp.var[j,]= lp.mv$meansq
        
        lq.mv = ff.moments_exp(-alpha.mv$mean, alpha.mv$var)  #find mean and meansq of q  
        lq.mean[j,] = lq.mv$mean
        lq.var[j,]= lq.mv$meansq
      }
    }
  }
  if(log==TRUE){ #reconstructs estimate from the "wavelet" space on the log level
    est.mean=cxxreverse_pwave(log(ls),lp.mean,lq.mean)
    est.var=cxxreverse_pwave(0,lp.var,lq.var)     
  }else{         #reconstruction on non-log level
    est.mean=exp(cxxreverse_pwave(log(ls),log(lp.mean),log(lq.mean)))
    est.ms=exp(cxxreverse_pwave(2*log(ls),log(lp.var),log(lq.var)))
    est.var=pmax(est.ms-est.mean^2,0)
  }
  if(reflect==TRUE){
    est.mean=est.mean[sig.ind]
    est.var=est.var[sig.ind]
  }
  if(return.est==TRUE){  #if return.est=TRUE then only estimate (and not variance) is returned
    return(est=est.mean)
  }else{
    return(list(est=est.mean,var=est.var,lp.mean=lp.mean,lq.mean=lq.mean))  
  }
}




