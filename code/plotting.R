#This script contains some necessary code for plotting figures


#' Interleave two vectors
#' @param x,y: two vectors of the same length
#' @return a vector of length twice that of x (or y)
interleave=function(x,y){
  return(as.vector(rbind(x,y)))
}


#' Shift a vector one unit to the right
#' @param x: a vector
#' @return a vector of the same length as that of x 
rshift = function(x){L=length(x); return(c(x[L],x[-L]))}
lshift = function(x){return(c(x[-1],x[1]))}


#' Produces two TI tables. One table contains the difference between adjacent pairs of data in the same resolution, and the other table contains the sum.
#' @param sig: a signal of length a power of 2
#' @return a list of two TI tables in the form of matrices
titable=function(sig){
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