source("../Rcode/wd_var.R")
source("../Rcode/jash.R")

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




bayesmooth.jash = function(x,mu.t,prior="nullbiased",pointmass=TRUE,nullcheck=TRUE,gridmult=0,mixsd=NULL,VB=FALSE,weight=0.5){
  n = length(x)
  J = log2(n)
  if(!isTRUE(all.equal(J,trunc(J)))){stop("Error: number of columns of x must be power of 2")}


  tsum = sum(x)
  y = cxxtitable(x)$difftable

  wmean = matrix(0,J,n) 

    var.est=(x-mu.t)^2
    vtable=cxxtitable(2/3*var.est^2)$sumtable
    vdtable=cxxtitable(var.est)$difftable
    for(j in 0:(J-1)){
      zdat.ash=jasha(vdtable[j+2,],sqrt((vtable[j+2,])),df=50)
      wmean[j+1,] = zdat.ash$PosteriorMean/2
    }
    wwmean=-wmean
    var.est=cxxreverse.gwave(sum(var.est),wmean,wwmean)
    return(var.est)

} 