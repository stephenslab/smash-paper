source("D:/Grad School/Spring 2013/multiscale_ash/simulation_1d_g/wd_var.R")
library(wavethresh)
require(ashr)
require(Rcpp)
require(inline)


#interleave two vectors
interleave=function(x,y){
  return(as.vector(rbind(x,y)))
}


#shift a vector right and left respectively
rshift = function(x){L=length(x); return(c(x[L],x[-L]))}
lshift = function(x){return(c(x[-1],x[1]))}
rshiftx = function(x,sh){
  L=length(x)
  if(sh==0){
     return(x)
  }else{
    return(c(x[(L-sh+1):L],x[-((L-sh+1):L)]))
  }
}
lshiftmat = function(x,sh){
  if(sh==0){
    return(x)
  }else{
    return(cbind(x[,-(1:sh)],x[,1:sh]))
  }
}
rshiftmat = function(x,sh){
  if(sh==0){
    return(x)
  }else{
    return(cbind(x[,(L-sh+1):L],x[,-((L-sh+1):L)]))
  }
}


#The following produces both a TItable,
#and a "parent" table whose pairwise comparisons would be used to create a TI table
#So, for example, in the ith row, elements 1, 2 would be the parents of the
#first element in the (i+1)the row of the TI table
# INPUT: sig, an n vector of Poisson counts at n locations
# OUTPUT: a list, with elements
# TItable - the usual TI table
# parent - the parent values used to make the TI table
titable=function(sig){
  n = length(sig)
  J = log2(n)
  
# Create decomposition table of signal, using pairwise sums,
# keeping just the values that are *not* redundant under the
# shift-invariant scheme.  This is very similar to TI-tables
# in Donoho and Coifman's TI-denoising framework.
  dmat = matrix(0, nrow=J+1, ncol=n)
  ddmat = matrix(0, nrow=J+1, ncol=n)
  dmat[1,] = sig
  ddmat[1,] = sig

  #dmat[1,] = as.matrix(sig)  
  #dmat2 = matrix(0, nrow=J, ncol=2*n) #the parent table
  
  for(D in 0:(J-1)){
    nD = 2^(J-D); 
    nDo2 = nD/2;
    twonD = 2*nD;
    for(l in 0:(2^D-1)){
      ind = (l*nD+1):((l+1)*nD)
      #ind2 = (l*twonD+1):((l+1)*twonD)
      x = dmat[D+1,ind]
      lsumx = x[seq(from=1,to=nD-1, by=2)] + x[seq(from=2,to=nD,by=2)]
      ldiffx = x[seq(from=1,to=nD-1, by=2)] - x[seq(from=2,to=nD,by=2)]
      rx = rshift(x);
      rsumx = rx[seq(from=1,to=nD-1, by=2)] + rx[seq(from=2,to=nD,by=2)]
      rdiffx = rx[seq(from=1,to=nD-1, by=2)] - rx[seq(from=2,to=nD,by=2)]
      dmat[D+2,ind] = c(lsumx,rsumx)
      ddmat[D+2,ind] = c(ldiffx,rdiffx)
    }
  }
  return(list(sumtable=dmat,difftable=ddmat))
}


#titable implemented in C++
src1 <- '
        NumericVector signal=sig; 
        int n=(int) signal.size();
        int J=(int) log2((double)n);

        NumericVector tempvec(2*n);
        NumericMatrix sumtable(J+1,n);
        NumericMatrix difftable(J+1,n);
        sumtable(0,_) = signal;
        difftable(0,_) = signal;

        for (int D=0; D<J; D++){
           int nD=(int) pow(2., (int) (J-D)), pD=(int) pow(2.,(int) D);
           for (int l=0; l<pD; l++){
              int a=l*nD+1;
              double d;
              for (int i=0; i<nD-1; i++){
                 d=sumtable(D,a+i-1);
                 tempvec(i)=d;
                 tempvec(i+nD+1)=d;
              }
              //i=nD-1
              d=sumtable(D,a+nD-2);
              tempvec(nD-1)=d;
              tempvec(nD)=d;

              for (int i=0; i<nD; i++){
                sumtable(D+1,a+i-1)=tempvec(2*i)+tempvec(2*i+1);
                difftable(D+1,a+i-1)=tempvec(2*i)-tempvec(2*i+1);
              }
           }
        }
        return(List::create(Named("sumtable")=sumtable, Named("difftable")=difftable));
        '
cxxtitable <- cxxfunction(signature(sig="numeric"),
                                body=src1,
                                plugin="Rcpp",
                                inc="#include <cmath>")






tirtable=function(sig){
  n = length(sig)
  J = log2(n)
  
  dmat = matrix(0, nrow=J+1, ncol=n)
  ddmat = matrix(0, nrow=J+1, ncol=n)
  dmat[1,] = sig
  ddmat[1,] = sig

  #dmat[1,] = as.matrix(sig)  
  #dmat2 = matrix(0, nrow=J, ncol=2*n) #the parent table
  
  for(D in 0:(J-1)){
    nD = 2^(J-D); 
    nDo2 = nD/2;
    twonD = 2*nD;
    for(l in 0:(2^D-1)){
      ind = (l*nD+1):((l+1)*nD)
      #ind2 = (l*twonD+1):((l+1)*twonD)
      x = dmat[D+1,ind]
      lsumx = x[seq(from=1,to=nD-1, by=2)] + x[seq(from=2,to=nD,by=2)]
      ldiffx = log(x[seq(from=1,to=nD-1, by=2)]) - log(x[seq(from=2,to=nD,by=2)])
      rx = rshift(x);
      rsumx = rx[seq(from=1,to=nD-1, by=2)] + rx[seq(from=2,to=nD,by=2)]
      rdiffx = log(rx[seq(from=1,to=nD-1, by=2)]) - log(rx[seq(from=2,to=nD,by=2)])
      dmat[D+2,ind] = c(lsumx,rsumx)
      ddmat[D+2,ind] = c(ldiffx,rdiffx)
    }
  }
  return(ddmat)
}


#tirtable implemented in C++
src4 <- '
        NumericVector signal=sig; 
        int n=(int) signal.size();
        int J=(int) log2((double)n);

        NumericVector tempvec(2*n);
        NumericMatrix sumtable(J+1,n);
        NumericMatrix vrtable(J+1,n);
        sumtable(0,_) = signal;
        vrtable(0,_) = signal;

        for (int D=0; D<J; D++){
           int nD=(int) pow(2., (int) (J-D)), pD=(int) pow(2.,(int) D);
           for (int l=0; l<pD; l++){
              int a=l*nD+1;
              double d;
              for (int i=0; i<nD-1; i++){
                 d=sumtable(D,a+i-1);
                 tempvec(i)=d;
                 tempvec(i+nD+1)=d;
              }
              //i=nD-1
              d=sumtable(D,a+nD-2);
              tempvec(nD-1)=d;
              tempvec(nD)=d;

              for (int i=0; i<nD; i++){
                sumtable(D+1,a+i-1)=tempvec(2*i)+tempvec(2*i+1);
                vrtable(D+1,a+i-1)=log((double)tempvec(2*i))-log((double)tempvec(2*i+1));
              }
           }
        }
        return(vrtable);
        '
cxxtirtable <- cxxfunction(signature(sig="numeric"),
                                body=src4,
                                plugin="Rcpp",
                                inc="#include <cmath>")


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
reverse.gwave=function(est,lp,lq=NULL){
  if(is.null(lq)){
    lq = -lp
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

      estvec = est[ind]/2
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


#reverse.gwave, as implemented in C++
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
              dep=est(a+i-1)/2;
              //dp=pp(D,a+i-1);
              dp=pp(J-1-D,a+i-1);
              //dq=qq(D,a+i-1);
              dq=qq(J-1-D,a+i-1);
              tempvecl(2*i)=dep+dp;
              tempvecl(2*i+1)=dep+dq;
            }
            for (int i=nDo2; i<nD-1; i++){
              dep=est(a+i)/2;
              //dp=pp(D,a+i);
              dp=pp(J-1-D,a+i);
              deq=est(a+i-1)/2;
              //dq=qq(D,a+i-1);
              dq=qq(J-1-D,a+i-1);
              tempvecr(2*(i-nDo2))=deq+dq;
              tempvecr(2*(i-nDo2)+1)=dep+dp;
            }         
            //i=nD-1
            dep=est(a+nDo2-1)/2;
            //dp=pp(D,a+nDo2-1);
            dp=pp(J-1-D,a+nDo2-1);
            deq=est(a+nD-2)/2; 
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
cxxreverse.gwave <- cxxfunction(signature(estimate="numeric",pmat="numeric",qmat="numeric"),
                                body=src2,
                                plugin="Rcpp",
                                inc="#include <cmath>")



reverse.gvwave=function(est,lp,lq=NULL){
  if(is.null(lq)){
    lq = -lp
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

      estvec = est[ind]/4
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





src3 <- '
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
              dep=est(a+i-1)/4;
              //dp=pp(D,a+i-1);
              dp=pp(J-1-D,a+i-1);
              //dq=qq(D,a+i-1);
              dq=qq(J-1-D,a+i-1);
              tempvecl(2*i)=dep+dp;
              tempvecl(2*i+1)=dep+dq;
            }
            for (int i=nDo2; i<nD-1; i++){
              dep=est(a+i)/4;
              //dp=pp(D,a+i);
              dp=pp(J-1-D,a+i);
              deq=est(a+i-1)/4;
              //dq=qq(D,a+i-1);
              dq=qq(J-1-D,a+i-1);
              tempvecr(2*(i-nDo2))=deq+dq;
              tempvecr(2*(i-nDo2)+1)=dep+dp;
            }         
            //i=nD-1
            dep=est(a+nDo2-1)/4;
            //dp=pp(D,a+nDo2-1);
            dp=pp(J-1-D,a+nDo2-1);
            deq=est(a+nD-2)/4; 
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
cxxreverse.gvwave <- cxxfunction(signature(estimate="numeric",pmat="numeric",qmat="numeric"),
                                body=src3,
                                plugin="Rcpp",
                                inc="#include <cmath>")



vartable=function(sig){
  n = length(sig)
  J = log2(n)

  dmat = matrix(0, nrow=J+1, ncol=n)
  ddmat = matrix(0, nrow=J+1, ncol=n)
  dmat[1,] = sig
  shift = list(0)
  
  for(D in 0:(J-1)){
    shift[[D+2]]=0
    nD = 2^(J-D); 
    nDo2 = nD/2;
    if(D==0){
      pos = (1:n)%%2==0
    }else{
      pos = (1:n)%%2^(D+1)==(1+2^D)
    }
    for(l in 0:(2^D-1)){
      ind = (l*nD+1):((l+1)*nD)
      sh = shift[[D+1]][l+1]
      x.sig = rshiftx(sig,sh)
      x = dmat[D+1,ind]     
      x.ddmat = ddmat[D+1,ind]
      xv.sel = x.sig[pos]
      x.temp = x.ddmat[seq(from=1,to=nD-1, by=2)] + x.ddmat[seq(from=2,to=nD,by=2)]
      xv.dmat = x.temp-xv.sel
      xv.ddmat = x.temp+xv.sel
      rx = rshift(x);
      rsh = sh+2^D
      rx.sig = rshiftx(sig,rsh)
      rx.ddmat = rshift(x.ddmat)
      rxv.sel = rx.sig[pos]
      rx.temp = rx.ddmat[seq(from=1,to=nD-1, by=2)] + rx.ddmat[seq(from=2,to=nD,by=2)]
      rxv.dmat = rx.temp-rxv.sel
      rxv.ddmat = rx.temp+rxv.sel
      dmat[D+2,ind] = c(xv.dmat,rxv.dmat)
      ddmat[D+2,ind] = c(xv.ddmat,rxv.ddmat)
      shift[[D+2]][(2*(l+1)-1):(2*(l+1))]=c(sh,rsh)
    }
  }
  return(vartable=dmat)
}



bayesmooth.est.var = function(x,sigma=NULL,v.est=FALSE,return.est=TRUE,basis="haar",prior="nullbiased",pointmass=TRUE,nullcheck=TRUE,gridmult=0,mixsd=NULL,VB=FALSE){
  n = length(x)
  J = log2(n)
  if(!isTRUE(all.equal(J,trunc(J)))){stop("Error: number of columns of x must be power of 2")}
  if(length(sigma)==1){sigma=rep(sigma,n)}


  tsum = sum(x)
  y = cxxtitable(x)$difftable
  if(basis=="symm8"){x.w.d = wd(x, filter.number=8, family="DaubLeAsymm", type = "station")}

  wmean = matrix(0,J,n) 
  wvar = matrix(0,J,n) 
  wmean1 = matrix(0,J,n) 
  wmean2 = matrix(0,J,n) 

  if(is.null(sigma)){
    var.est1.ini=(x-lshift(x))^2/2
    var.est1l.ini=lshift(var.est1.ini)
    var.est2.ini=(rshift(x)-x)^2/2
    var.est2l.ini=lshift(var.est2.ini)
    vtable1=cxxtitable(2/9*(var.est1.ini+var.est1l.ini)^2)$sumtable
    vtable2=cxxtitable(2/9*(var.est2.ini+var.est2l.ini)^2)$sumtable
    vvtable1=1/3*vartable(var.est1.ini^2)
    vvtable2=1/3*vartable(var.est2.ini^2)
    vtable1=vtable1+vvtable1
    vtable2=vtable2+vvtable2
    vdtable1=titable(var.est1.ini)$difftable
    vdtable2=titable(var.est2.ini)$difftable
    fac=c(1.85,1.7,1.4,1.2,1.12,1.1,1.05,rep(1,J-7))
    for(j in 0:(J-1)){
      zdat.ash=ash(vdtable1[j+2,],sqrt(fac[j+1]*(vtable1[j+2,])),prior=prior,pointmass=pointmass,nullcheck=nullcheck,VB=VB,mixsd=mixsd,gridmult=gridmult)
      wmean1[j+1,] = zdat.ash$PosteriorMean/2
      zdat.ash=ash(vdtable2[j+2,],sqrt(fac[j+1]*(vtable2[j+2,])),prior=prior,pointmass=pointmass,nullcheck=nullcheck,VB=VB,mixsd=mixsd,gridmult=gridmult)
      wmean2[j+1,] = zdat.ash$PosteriorMean/2
    }
    wmean=(wmean1+wmean2)/2
    wwmean=-wmean
    var.est=cxxreverse.gwave(sum((var.est1.ini+var.est2.ini)/2),wmean,wwmean)
    var.est[var.est<=0]=1e-6
    sigma=sqrt(var.est)
    if(basis=="haar"){
      vtable=cxxtitable(sigma^2)$sumtable
      for(j in 0:(J-1)){
        fac=c(1.3,1.17,1.08,1.04,rep(1,J-4))
        zdat.ash=fast.ash(y[j+2,],sqrt(vtable[j+2,]),prior=prior,pointmass=pointmass,nullcheck=nullcheck,VB=VB,mixsd=mixsd,gridmult=gridmult)
        wmean[j+1,] = zdat.ash$PosteriorMean/2
      }
      wwmean=-wmean
      mu.est=cxxreverse.gwave(tsum,wmean,wwmean)
    }else if(basis=="symm8"){
      x.w=x.w.d
      x.w.v1 = wd.var(var.est1.ini, type = "station")
      x.w.v2 = wd.var(var.est2.ini, type = "station")
      for(j in 0:(J-1)){
        zdat.ash=fast.ash(accessD(x.w,j),sqrt((accessC(x.w.v1,j)+accessC(x.w.v2,j))/2),prior=prior,pointmass=pointmass,nullcheck=nullcheck,VB=VB,mixsd=mixsd,gridmult=gridmult)
        x.pm = zdat.ash$PosteriorMean
        x.w = putD(x.w,j,x.pm)
      }
      mu.est=AvBasis(convert(x.w))
    }
    var.est=(x-mu.est)^2
    var.var.est=2/3*var.est^2
    vtable=cxxtitable(var.var.est)$sumtable
    vdtable=cxxtitable(var.est)$difftable
    vrtable=cxxtirtable(var.est)
    virtable=3*(abs(vrtable)>=log(25))+2*(abs(vrtable)>=log(3)&abs(vrtable)<log(25))+1*(abs(vrtable)<log(3))
    fac=c(2,1.5,1.3,1.15,1.1,1.07,rep(1,J-6))
    for(j in 0:(J-1)){
      zdat.ash=fast.ash(vdtable[j+2,],sqrt((fac[j+1]*vtable[j+2,]/virtable[j+2,])),prior=prior,pointmass=pointmass,nullcheck=nullcheck,VB=VB,mixsd=mixsd,gridmult=gridmult)
      wmean[j+1,] = zdat.ash$PosteriorMean/2
    }
    #for(j in (J-3):(J-1)){
    #  wmean[j+1,] = vdtable[j+2,]/2
    #}
    wwmean=-wmean
    var.est=cxxreverse.gwave(sum((var.est1.ini+var.est2.ini)/2),wmean,wwmean)
    var.est[var.est<=0]=1e-6
    sigma=sqrt(var.est)
  }
  if(basis=="haar"){
    vtable=cxxtitable(sigma^2)$sumtable
    for(j in 0:(J-1)){
      zdat.ash=fast.ash(y[j+2,],sqrt(vtable[j+2,]),prior=prior,pointmass=pointmass,nullcheck=nullcheck,VB=VB,mixsd=mixsd,gridmult=gridmult)
      wmean[j+1,] = zdat.ash$PosteriorMean/2
      wvar[j+1,] = zdat.ash$PosteriorSD^2/4
    }
    wwmean=-wmean
    wwvar=wvar
    mu.est=cxxreverse.gwave(tsum,wmean,wwmean)
    mu.est.var=cxxreverse.gvwave(0,wvar,wwvar)
  }else if(basis=="symm8"){
    x.w=x.w.d
    x.w.v = wd.var(sigma^2, type = "station")
    for(j in 0:(J-1)){
      zdat.ash=fast.ash(accessD(x.w,j),sqrt(accessC(x.w.v,j)),prior=prior,pointmass=pointmass,nullcheck=nullcheck,VB=VB,mixsd=mixsd,gridmult=gridmult)
      x.pm = zdat.ash$PosteriorMean
      x.w = putD(x.w,j,x.pm)
    }
    mu.est=AvBasis(convert(x.w))
  }
  if(v.est==FALSE){
    if(return.est==TRUE){
      return(mu.est)
    }else{
      if(basis=="symm8"){stop("Error: Posterior variance only returned with haar basis")}
      return(list(mu.est=mu.est,mu.var=mu.est.var))
    }
  }else{
    var.est=(x-mu.est)^2
    vtable=cxxtitable(2/3*var.est^2)$sumtable
    vdtable=cxxtitable(var.est)$difftable
    vrtable=cxxtirtable(var.est)
    virtable=3*(abs(vrtable)>=log(25))+2*(abs(vrtable)>=log(3)&abs(vrtable)<log(25))+1*(abs(vrtable)<log(3))
    fac=c(2,1.5,1.3,1.15,1.1,1.07,rep(1,J-6))
    for(j in 0:(J-1)){
      zdat.ash=fast.ash(vdtable[j+2,],sqrt((fac[j+1]*vtable[j+2,]/virtable[j+2,])),prior=prior,pointmass=pointmass,nullcheck=nullcheck,VB=VB,mixsd=mixsd,gridmult=gridmult)
      wmean[j+1,] = zdat.ash$PosteriorMean/2
      wvar[j+1,] = zdat.ash$PosteriorSD^2/4
    }
    #for(j in (J-4):(J-1)){
    #  wmean[j+1,] = vdtable[j+2,]/2
    #}
    wwmean=-wmean
    wwvar=wvar
    var.est=cxxreverse.gwave(sum((var.est1.ini+var.est2.ini)/2),wmean,wwmean)
    var.est.var=cxxreverse.gvwave(0,wvar,wwvar)
    if(return.est==TRUE){
      return(var.est)
    }else{
      return(list(var.est=var.est,var.var=var.est.var))
    }
  }
}
 