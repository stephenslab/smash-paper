library(wavethresh)
require(ashr)
require(Rcpp)
require(inline)
library(Matrix)
source("~/ashwave/Rcode/jash.r")



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





wd.D=function (data, filter.number = 10, family = "DaubLeAsymm", type = "wavelet", 
    bc = "periodic", verbose = FALSE, min.scale = 0, precond = TRUE) 
{
    if (verbose == TRUE) 
        cat("wd: Argument checking...")
    if (!is.atomic(data)) 
        stop("Data is not atomic")
    DataLength <- length(data)
    nlevels <- nlevelsWT(data)
    if (is.na(nlevels)) 
        stop("Data length is not power of two")
    if (type != "wavelet" && type != "station") 
        stop("Unknown type of wavelet decomposition")
    if (type == "station" && bc != "periodic") 
        stop("Can only do periodic boundary conditions with station")
    if (verbose == TRUE) 
        cat("...done\nFilter...")
    if (bc != "interval") 
        filter <- filter.select(filter.number = filter.number, 
            family = family)
    if (verbose == TRUE) 
        cat("...selected\nFirst/last database...")
    fl.dbase <- first.last(LengthH = length(filter$H), DataLength = DataLength, 
        type = type, bc = bc)
    if (bc == "interval") {
        ans <- wd.int(data = data, preferred.filter.number = filter.number, 
            min.scale = min.scale, precond = precond)
        fl.dbase <- first.last(LengthH = length(filter$H), DataLength = DataLength, 
            type = type, bc = bc, current.scale = min.scale)
        filter <- list(name = paste("CDV", filter.number, sep = ""), 
            family = "CDV", filter.number = filter.number)
        l <- list(transformed.vector = ans$transformed.vector, 
            current.scale = ans$current.scale, filters.used = ans$filters.used, 
            preconditioned = ans$preconditioned, date = ans$date, 
            nlevels = IsPowerOfTwo(length(ans$transformed.vector)), 
            fl.dbase = fl.dbase, type = type, bc = bc, filter = filter)
        class(l) <- "wd"
        return(l)
    }
    dtsp <- tsp(data)
    C <- rep(0, fl.dbase$ntotal)
    C[1:DataLength] <- data
    if (verbose == TRUE) 
        error <- 1
    else error <- 0
    if (verbose == TRUE) 
        cat("built\n")
    if (verbose == TRUE) 
        cat("Decomposing...\n")
    nbc <- switch(bc, periodic = 1, symmetric = 2)
    if (is.null(nbc)) 
        stop("Unknown boundary condition")
    ntype <- switch(type, wavelet = 1, station = 2)
    if (is.null(filter$G)) {
        wavelet.decomposition <- .C("wavedecomp", C = as.double(C), 
            D = as.double(rep(0, fl.dbase$ntotal.d)), H = as.double(filter$H), 
            LengthH = as.integer(length(filter$H)), nlevels = as.integer(nlevels), 
            firstC = as.integer(fl.dbase$first.last.c[, 1]), 
            lastC = as.integer(fl.dbase$first.last.c[, 2]), offsetC = as.integer(fl.dbase$first.last.c[, 
                3]), firstD = as.integer(fl.dbase$first.last.d[, 
                1]), lastD = as.integer(fl.dbase$first.last.d[, 
                2]), offsetD = as.integer(fl.dbase$first.last.d[, 
                3]), ntype = as.integer(ntype), nbc = as.integer(nbc), 
            error = as.integer(error), PACKAGE = "wavethresh")
    }
    else {
        wavelet.decomposition <- .C("comwd", CR = as.double(Re(C)), 
            CI = as.double(Im(C)), LengthC = as.integer(fl.dbase$ntotal), 
            DR = as.double(rep(0, fl.dbase$ntotal.d)), DI = as.double(rep(0, 
                fl.dbase$ntotal.d)), LengthD = as.integer(fl.dbase$ntotal.d), 
            HR = as.double(Re(filter$H)), HI = as.double(-Im(filter$H)), 
            GR = as.double(Re(filter$G)), GI = as.double(-Im(filter$G)), 
            LengthH = as.integer(length(filter$H)), nlevels = as.integer(nlevels), 
            firstC = as.integer(fl.dbase$first.last.c[, 1]), 
            lastC = as.integer(fl.dbase$first.last.c[, 2]), offsetC = as.integer(fl.dbase$first.last.c[, 
                3]), firstD = as.integer(fl.dbase$first.last.d[, 
                1]), lastD = as.integer(fl.dbase$first.last.d[, 
                2]), offsetD = as.integer(fl.dbase$first.last.d[, 
                3]), ntype = as.integer(ntype), nbc = as.integer(nbc), 
            error = as.integer(error), PACKAGE = "wavethresh")
    }
    if (verbose == TRUE) 
        cat("done\n")
    error <- wavelet.decomposition$error
    if (error != 0) {
        cat("Error ", error, " occured in wavedecomp\n")
        stop("Error")
    }
    if (is.null(filter$G)) {
        l <- list(C = wavelet.decomposition$C, D = wavelet.decomposition$D, 
            nlevels = nlevelsWT(wavelet.decomposition), fl.dbase = fl.dbase, 
            filter = filter, type = type, bc = bc, date = date())
    }
    else {
        l <- list(C = complex(real = wavelet.decomposition$CR, 
            imaginary = wavelet.decomposition$CI), D = complex(real = wavelet.decomposition$DR, 
            imaginary = wavelet.decomposition$DI), nlevels = nlevelsWT(wavelet.decomposition), 
            fl.dbase = fl.dbase, filter = filter, type = type, 
            bc = bc, date = date())
    }
    class(l) <- "wd"
    if (!is.null(dtsp)) 
        tsp(l) <- dtsp
    l$D
}





ndwt.mat=function(n,filter.number,family){
  J=log2(n)
  X=diag(rep(1,n))
  W=matrix(0,n*J,n)
  #for(i in 1:n){
  #  W[,i]=wd(X[i,],filter.number,family,type="station")$D
  #}
  W=apply(X,1,wd.D,filter.number=filter.number,family=family,type="station")
  #W=Matrix(W,sparse=TRUE) 
  return(W^2)
}



bayesmooth = function(x,sigma=NULL,v.est=FALSE,v.basis=FALSE,return.est=TRUE,filter.number=1,family="DaubExPhase",prior="nullbiased",pointmass=TRUE,nullcheck=TRUE,gridmult=0,mixsd=NULL,VB=FALSE,jash=FALSE,weight=0.5){
  n = length(x)
  J = log2(n)
  if(!isTRUE(all.equal(J,trunc(J)))){stop("Error: number of columns of x must be power of 2")}
  if(length(sigma)==1){sigma=rep(sigma,n)}
  if(filter.number==1&family=="DaubExPhase"){
    basis="haar"
  }else{
    basis=paste0("family",filter.number)
  }


  tsum = sum(x)
  y = cxxtitable(x)$difftable
  if(basis!="haar"){
    x.w.d = wd(x, filter.number=filter.number, family=family, type = "station")
    W2=ndwt.mat(n,filter.number=filter.number,family=family)
  }

  wmean = matrix(0,J,n) 
  wvar = matrix(0,J,n) 

  if(is.null(sigma)){
    var.est1.ini=(x-lshift(x))^2/2
    var.est2.ini=(rshift(x)-x)^2/2
    if(basis=="haar"){
      vtable1=cxxtitable(var.est1.ini)$sumtable
      vtable2=cxxtitable(var.est2.ini)$sumtable
      vtable=(vtable1+vtable2)/2
      for(j in 0:(J-1)){
        ind.nnull=(y[j+2,]!=0)|(vtable[j+2,]!=0)
        zdat.ash=fast.ash(y[j+2,ind.nnull],sqrt(vtable[j+2,ind.nnull]),prior=prior,pointmass=pointmass,nullcheck=nullcheck,VB=VB,mixsd=mixsd,gridmult=gridmult)
        wmean[j+1,ind.nnull]=zdat.ash$PosteriorMean/2
        wmean[j+1,!ind.nnull]=0
      }
      wwmean=-wmean
      mu.est=cxxreverse.gwave(tsum,wmean,wwmean)
    }else{
      x.w=x.w.d
      x.w.v=apply((rep(1,n*J)%o%((var.est1.ini+var.est2.ini)/2))*W2,1,sum)  #diagonal of W*V*W'   
      x.pm=rep(0,n)  
      for(j in 0:(J-1)){
        index=(((J-1)-j)*n+1):((J-j)*n)
        x.w.j=accessD(x.w,j)
        x.w.v.j=x.w.v[index]
        ind.nnull=(x.w.j!=0)|(x.w.v.j!=0)
        zdat.ash=fast.ash(x.w.j[ind.nnull],sqrt(x.w.v.j[ind.nnull]),prior=prior,pointmass=pointmass,nullcheck=nullcheck,VB=VB,mixsd=mixsd,gridmult=gridmult)
        x.pm[ind.nnull] = zdat.ash$PosteriorMean
        x.pm[!ind.nnull] = 0
        x.w = putD(x.w,j,x.pm)
      }
      mu.est=AvBasis(convert(x.w))
    }
    var.est=(x-mu.est)^2
    var.var.est=2/3*var.est^2
    if(basis=="haar"|v.basis==FALSE){
      vtable=cxxtitable(var.var.est)$sumtable
      vdtable=cxxtitable(var.est)$difftable
      vrtable=cxxtirtable(var.est)
      fac1=c(3,2,1.5)
      fac2=c(2,1.5,1)
      #fac1=c(1,1,1)
      #fac2=c(1,1,1)
      for(j in 0:(J-1)){
        if(j>=0&j<=1){
          virtable=vrtable[j+2,]
          virtable=fac1[1]*(abs(virtable)>=log(4))+fac1[2]*(abs(virtable)>=log(2)&abs(virtable)<log(4))+fac1[3]*(abs(virtable)<log(2))
        }else if(j>=2&j<=5){
          virtable=vrtable[j+2,]
          virtable=fac2[1]*(abs(virtable)>=log(4))+fac2[2]*(abs(virtable)>=log(2)&abs(virtable)<log(4))+fac2[3]*(abs(virtable)<log(2))
        }else{
          virtable=1
        }
        if(jash==TRUE){
          zdat.ash=jasha(vdtable[j+2,],sqrt((vtable[j+2,])),df=50)
        }else{
          zdat.ash=fast.ash(vdtable[j+2,],sqrt((vtable[j+2,]/virtable)),prior=prior,pointmass=pointmass,nullcheck=nullcheck,VB=VB,mixsd=mixsd,gridmult=gridmult)
        }
        wmean[j+1,] = zdat.ash$PosteriorMean/2
      }
      wwmean=-wmean
      var.est=cxxreverse.gwave(weight*sum(var.est)+(1-weight)*sum((var.est1.ini+var.est2.ini)/2),wmean,wwmean)
    }else{
      x.w=wd(var.est, filter.number=filter.number, family=family, type = "station")
      x.w.v=apply((rep(1,n*J)%o%var.var.est)*W2,1,sum)  #diagonal of W*V*W'   
      for(j in 0:(J-1)){
        index=(((J-1)-j)*n+1):((J-j)*n)
        x.w.j=accessD(x.w,j)
        x.w.v.j=x.w.v[index]
        if(jash==TRUE){
          zdat.ash=jasha(x.w.j,sqrt(x.w.v.j),df=50)
        }else{
          zdat.ash=fast.ash(x.w.j,sqrt(x.w.v.j),prior=prior,pointmass=pointmass,nullcheck=nullcheck,VB=VB,mixsd=mixsd,gridmult=gridmult)
        }
        x.pm = zdat.ash$PosteriorMean
        x.w = putD(x.w,j,x.pm)
      }
      var.est=AvBasis(convert(x.w))
    }
    var.est[var.est<=0]=1e-8
    sigma=sqrt(var.est)
  }
  if(basis=="haar"|v.basis==FALSE){
    vtable=cxxtitable(sigma^2)$sumtable
    for(j in 0:(J-1)){
      ind.nnull=(y[j+2,]!=0)|(vtable[j+2,]!=0)
      zdat.ash=fast.ash(y[j+2,ind.nnull],sqrt(vtable[j+2,ind.nnull]),prior=prior,pointmass=pointmass,nullcheck=nullcheck,VB=VB,mixsd=mixsd,gridmult=gridmult)
      wmean[j+1,ind.nnull]=zdat.ash$PosteriorMean/2
      wmean[j+1,!ind.nnull]=0
      wvar[j+1,ind.nnull]=zdat.ash$PosteriorSD^2/4
      wvar[j+1,!ind.nnull]=0
    }
    wwmean=-wmean
    wwvar=wvar
    mu.est=cxxreverse.gwave(tsum,wmean,wwmean)
    mu.est.var=cxxreverse.gvwave(0,wvar,wwvar)
  }else{
    x.w=x.w.d
    x.w.v=apply((rep(1,n*J)%o%(sigma^2))*W2,1,sum)  #diagonal of W*V*W'   
    x.pm=rep(0,n)
    for(j in 0:(J-1)){
      index=(((J-1)-j)*n+1):((J-j)*n)
      x.w.j=accessD(x.w,j)
      x.w.v.j=x.w.v[index]
      ind.nnull=(x.w.j!=0)|(x.w.v.j!=0)
      zdat.ash=fast.ash(x.w.j[ind.nnull],sqrt(x.w.v.j[ind.nnull]),prior=prior,pointmass=pointmass,nullcheck=nullcheck,VB=VB,mixsd=mixsd,gridmult=gridmult)
      x.pm[ind.nnull] = zdat.ash$PosteriorMean
      x.pm[!ind.nnull] = 0
      x.w = putD(x.w,j,x.pm)
    }
    mu.est=AvBasis(convert(x.w))
  }
  if(v.est==FALSE){
    if(return.est==TRUE){
      return(mu.est)
    }else{
      if(basis!="haar"){stop("Error: Posterior variance only returned with haar basis")}
      return(list(mu.est=mu.est,mu.var=mu.est.var))
    }
  }else{
    var.est=(x-mu.est)^2
    var.var.est=2/3*var.est^2
    if(basis=="haar"){
      vtable=cxxtitable(var.var.est)$sumtable
      vdtable=cxxtitable(var.est)$difftable
      vrtable=cxxtirtable(var.est)
      fac1=c(3,2,1.5)
      fac2=c(2,1.5,1)
      #fac1=c(1,1,1)
      #fac2=c(1,1,1)
      for(j in 0:(J-1)){
        if(j>=0&j<=1){
          virtable=vrtable[j+2,]
          virtable=fac1[1]*(abs(virtable)>=log(4))+fac1[2]*(abs(virtable)>=log(2)&abs(virtable)<log(4))+fac1[3]*(abs(virtable)<log(2))
        }else if(j>=2&j<=5){
          virtable=vrtable[j+2,]
          virtable=fac2[1]*(abs(virtable)>=log(4))+fac2[2]*(abs(virtable)>=log(2)&abs(virtable)<log(4))+fac2[3]*(abs(virtable)<log(2))
        }else{
          virtable=1
        }
        if(jash==TRUE){
          zdat.ash=jasha(vdtable[j+2,],sqrt((vtable[j+2,])),df=50)
        }else{
          zdat.ash=fast.ash(vdtable[j+2,],sqrt((vtable[j+2,]/virtable)),prior=prior,pointmass=pointmass,nullcheck=nullcheck,VB=VB,mixsd=mixsd,gridmult=gridmult)
        }
        wmean[j+1,] = zdat.ash$PosteriorMean/2
      }
      wwmean=-wmean
      wwvar=wvar
      var.est=cxxreverse.gwave(sum(var.est),wmean,wwmean)
      var.est.var=cxxreverse.gvwave(0,wvar,wwvar)
    }else{
      x.w=wd(var.est, filter.number=filter.number, family=family, type = "station")
      x.w.v=apply((rep(1,n*J)%o%var.var.est)*W2,1,sum)  #diagonal of W*V*W'   
      for(j in 0:(J-1)){
        index=(((J-1)-j)*n+1):((J-j)*n)
        x.w.j=accessD(x.w,j)
        x.w.v.j=x.w.v[index]
        if(jash==TRUE){
          zdat.ash=jasha(x.w.j,sqrt(x.w.v.j),df=50)
        }else{
          zdat.ash=fast.ash(x.w.j,sqrt(x.w.v.j),prior=prior,pointmass=pointmass,nullcheck=nullcheck,VB=VB,mixsd=mixsd,gridmult=gridmult)
        }
        x.pm = zdat.ash$PosteriorMean
        x.w = putD(x.w,j,x.pm)
      }
      var.est=AvBasis(convert(x.w))
    }
    if(return.est==TRUE){
      return(var.est)
    }else{
      return(list(var.est=var.est,var.var=var.est.var))
    }
  }
}
 