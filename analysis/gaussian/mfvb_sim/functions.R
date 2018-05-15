# The source code for the ZOSull and meanVarMFVB functions was
# generously shared by Marianne Menictas and Matt Wand from the
# University of Technology, Sydney.

# Specify true nonparametric mean and variance functions.
fTrue    <- function(x) return(sin(3*pi*x^2))
gTrue    <- function(x) return(exp(0.1 + cos(4*pi*x)))
loggTrue <- function(x) return(0.1 + cos(4*pi*x))

########## R function: ZOSull ##########

# For creation of O'Sullivan-type Z matrices.

# Last changed: 16 DEC 2009

ZOSull <- function(x,range.x,intKnots,drv=0)
{
   # Check legality of range.x.
  
   if (!missing(range.x))
   {
      if (length(range.x)!=2) stop("range.x must be of length 2.")
      if (range.x[1]>range.x[2]) stop("range.x[1] exceeds range.x[1].")
      if (range.x[1]>min(x)) stop("range.x[1] must be <= than min(x).")
      if (range.x[2]<max(x)) stop("range.x[2] must be >= than max(x).")

   }
   
   if (drv>2) stop("splines not smooth enough for more than 2 derivatives")

   library(splines)

   # Set defaults for `range.x' and `intKnots'

   if (missing(range.x))
      range.x <- c(1.05*min(x)-0.05*max(x),1.05*max(x)-0.05*min(x))
   
   if (missing(intKnots))
   {
      numIntKnots <- min(length(unique(x)),35)
      intKnots <- quantile(unique(x),seq(0,1,length=
                  (numIntKnots+2))[-c(1,(numIntKnots+2))])
   }
   numIntKnots <- length(intKnots) 

   # Obtain the penalty matrix.

   allKnots <- c(rep(range.x[1],4),intKnots,rep(range.x[2],4)) 
   K <- length(intKnots) ; L <- 3*(K+8)
   xtilde <- (rep(allKnots,each=3)[-c(1,(L-1),L)]+ 
              rep(allKnots,each=3)[-c(1,2,L)])/2
   wts <- rep(diff(allKnots),each=3)*rep(c(1,4,1)/6,K+7)
   Bdd <- spline.des(allKnots,xtilde,derivs=rep(2,length(xtilde)),
                     outer.ok=TRUE)$design  
   Omega     <- t(Bdd*wts)%*%Bdd     

   # Use the spectral decomposition of Omega to obtain Z.

   svdOmega <- svd(Omega) 
   indsZ <- 1:(numIntKnots+2)
   UZ <- svdOmega$u[,indsZ] 
   LZ <- t(t(UZ)/sqrt(svdOmega$d[indsZ]))

   # Perform stability check.   

   indsX <- (numIntKnots+3):(numIntKnots+4)
   UX <- svdOmega$u[,indsX]   
   L <- cbind(UX,LZ)
   stabCheck <- t(crossprod(L,t(crossprod(L,Omega))))          
   if (sum(stabCheck^2) > 1.0001*(numIntKnots+2))
       print("WARNING: NUMERICAL INSTABILITY ARISING\\
              FROM SPECTRAL DECOMPOSITION")

   # Obtain B and post-multiply by LZ matrix to get Z.

   B <- spline.des(allKnots,x,derivs=rep(drv,length(x)),
                     outer.ok=TRUE)$design  
   
   Z <- B%*%LZ

   # Add the `range.x' and 'intKnots' as attributes
   # of the return object.

   attr(Z,"range.x") <- range.x
   attr(Z,"intKnots") <- intKnots

   # Return Z matrix with 2 attributes.

   return(Z)
}

########## End of ZOSull ##########

########## R function: meanVarMFVB ##########

# For fitting using mean field variational Bayes 
# for simultaneous mean and variance function
# estimation for univariate nonparametric regression.

# Last changed: 10 JUL 2014

meanVarMFVB <- function(y,X,ncZu,ncZv,Au.hyp,Av.hyp,sigsq.gamma,
                        sigsq.beta,maxIter=1000,tolerance=1e-7)
{
   logML <- function(n,y,Cumat,Cvmat,ncX,ncZu,ncZv,Au.hyp,Av.hyp,
                     sigsq.beta,sigsq.gamma,mu.q.nu,mu.q.omega,
                     Sigma.q.nu,Sigma.q.omega,B.q.sigsqu,B.q.sigsqv,
		     mu.q.recip.sigsqu,mu.q.recip.sigsqv,
                     mu.q.recip.au,mu.q.recip.av)
   {
      ans <- - 0.5*n*log(2*pi) - 0.5*oneVec%*%Cvmat%*%mu.q.omega
      ans <- (-0.5*t(E.qv.r)%*%exp(Cvmat%*%mu.q.omega 
              + 0.5*diag(Cvmat%*%Sigma.q.omega%*%t(Cvmat))))
      ans <- - 0.5*ncX*log(sigsq.beta) 
      ans <- (- (1/(2*sigsq.beta))*(sum(mu.q.nu[1:ncX]^2)
              + sum(diag(Sigma.q.nu[1:ncX,1:ncX]))))
      ans <- + 0.5*determinant(Sigma.q.nu)$modulus 
      ans <- (+ 0.5*sum(diag(Sigma.q.nu*(sum(mu.q.nu^2) 
              + sum(diag(Sigma.q.nu))))))
      ans <- - t(mu.q.nu)%*%t(Cumat)%*%diag(psi.q.omega)%*%y
      ans <- - 0.5*ncX*log(sigsq.gamma) 
      ans <- (-(1/(2*sigsq.gamma))*(sum(mu.q.omega[1:ncX]^2) 
              + sum(diag(Sigma.q.omega[1:ncX,1:ncX])))) 
      ans <- (+ 0.5*determinant(Sigma.q.omega)$modulus 
              + 0.5*(ncZv + ncX) - 2*log(pi))
      ans <- + lgamma((ncZu+1)/2) + lgamma((ncZv+1)/2)
      ans <- -0.5*(ncZu+1)*log(B.q.sigsqu) - 0.5*(ncZv+1)*log(B.q.sigsqv)
      ans <- -log(Au.hyp) -log(Av.hyp) + mu.q.recip.sigsqu*mu.q.recip.au
      ans <- + mu.q.recip.sigsqv*mu.q.recip.av 
      ans <- - log(mu.q.recip.sigsqu + (1/Au.hyp^2))
      ans <- - log(mu.q.recip.sigsqv + (1/Av.hyp^2))
   
      return(ans)
   }
   
   # Set iteration parameters and do initialisations:

   Sigma.q.omega.Inv <- diag(ncCv)
   mu.q.nu <- rep(1,ncCu)
   mu.q.omega <- rep(0,ncCv)
   mu.q.recip.sigsqu <- 1
   mu.q.recip.sigsqv <- 1
   mu.q.recip.au <- 1
   mu.q.recip.av <- 1
   
   oneVec <- rep(1,n)

   # Set up convergence criteria:

   logMLcurr <- -Inf
   converged <- FALSE ; itNum <- 0 
   logMLgrid <- rep(NA,maxIter)

   while (!converged)
   {
      itNum <- itNum + 1
	  
      # Set current log(ML) value to previous one:

      logMLprev <- logMLcurr

      # Perform updates for (beta,u) parameters:

      psi.q.omega <- (as.vector(exp(-Cvmat%*%mu.q.omega 
                     + 0.5*diag(Cvmat%*%solve(Sigma.q.omega.Inv,t(Cvmat))))))			   
				   
      ridgeVecu <- c(rep((1/sigsq.beta),ncX),rep(mu.q.recip.sigsqu,ncZu))

      CuTpsiCu <- crossprod(Cumat,(psi.q.omega*Cumat))
      Sigma.q.nu.Inv <- CuTpsiCu + diag(ridgeVecu)
      
      mu.q.nu <- solve(Sigma.q.nu.Inv)%*%crossprod((Cumat),(diag(psi.q.omega)%*%y))
	  
      E.qv.r <- (as.vector(diag(tcrossprod(y-Cumat%*%mu.q.nu) 
                 + tcrossprod(Cumat%*%solve(Sigma.q.nu.Inv),Cumat))))
	  
      ridgeVecv <- c(rep((1/sigsq.gamma),ncX),rep(mu.q.recip.sigsqv,ncZv))
	  
      CvTpsiCv <- crossprod(Cvmat,diag(E.qv.r*psi.q.omega)%*%Cvmat)
      Sigma.q.omega.Inv <- CvTpsiCv + diag(ridgeVecv)
	  
      mu.q.omega <- (mu.q.omega + solve(Sigma.q.omega.Inv,
                    (crossprod(Cvmat,(E.qv.r*psi.q.omega - oneVec)) 
					 - ridgeVecv*mu.q.omega)))

      # Perform updates for sigsq parameters:

      Sigma.q.nu <- solve(Sigma.q.nu.Inv)
      B.q.sigsqu <- (0.5*sum(mu.q.nu[-(1:ncX)]^2)  
                   + 0.5*sum(diag(Sigma.q.nu[-(1:ncX),-(1:ncX)])) 
                   + mu.q.recip.au)
      mu.q.recip.sigsqu <- 0.5*(ncZu + 1)/B.q.sigsqu 
	  
      Sigma.q.omega <- solve(Sigma.q.omega.Inv)
      B.q.sigsqv <- (0.5*sum(mu.q.omega[-(1:ncX)]^2)  
                   + 0.5*sum(diag(Sigma.q.omega[-(1:ncX),-(1:ncX)])) 
                   + mu.q.recip.av)
      mu.q.recip.sigsqv <- 0.5*(ncZv + 1)/B.q.sigsqv 

      # Perform updatea for "a" parameters:

      mu.q.recip.au <- 1/(mu.q.recip.sigsqu + (1/Au.hyp^2))
      mu.q.recip.av <- 1/(mu.q.recip.sigsqv + (1/Av.hyp^2))

      # Obtain current log(ML):

      logMLcurr <- logML(n,y,Cumat,Cvmat,ncX,ncZu,ncZv,Au.hyp,Av.hyp,
                     sigsq.beta,sigsq.gamma,mu.q.nu,mu.q.omega,
                     Sigma.q.nu,Sigma.q.omega,B.q.sigsqu,B.q.sigsqv,
					 mu.q.recip.sigsqu,mu.q.recip.sigsqv,
                     mu.q.recip.au,mu.q.recip.av)
					 
      logMLgrid[itNum] <- logMLcurr

      # Compute relative error:

      relErr <- abs((logMLcurr/logMLprev)-1)

      # Check `converged' conditions:
                 
      if (itNum>=maxIter) 
      {
         converged <- TRUE
         print("WARNING: maximum number of iterations exceeded.")
      }

      if (relErr<tolerance) converged <- TRUE 
   }

   mu.q.nu <- as.vector(mu.q.nu)
   mu.q.omega <- as.vector(mu.q.omega)

   logMLgrid <- logMLgrid[3:itNum]
   
   return(list(mu.q.omega=mu.q.omega,mu.q.nu=mu.q.nu,
               Sigma.q.omega=Sigma.q.omega,Sigma.q.nu=Sigma.q.nu,
			   Sigma.q.omega.Inv=Sigma.q.omega.Inv,
			   Sigma.q.nu.Inv=Sigma.q.nu.Inv,
			   mu.q.recip.sigsqu=mu.q.recip.sigsqu,
               mu.q.recip.sigsqv=mu.q.recip.sigsqv,
               mu.q.recip.au=mu.q.recip.au,
               mu.q.recip.av=mu.q.recip.av,
			   mu.q.recip.sigsqv=mu.q.recip.sigsqv,
			   mu.q.recip.sigsqu=mu.q.recip.sigsqu,
			   E.qv.r=E.qv.r,psi.q.omega=psi.q.omega)) 	   
}

############# End of meanVarMFVB ############

