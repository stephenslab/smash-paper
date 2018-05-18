# TO DO: Explain here what this script is for.

# SET UP ENVIRONMENT
# ------------------
# TO DO: Explain here what these lines of code are doing.
library(smashr)
source("functions.R")

# SCRIPT PARAMETERS
# -----------------
# Set hyperparameters.
Au.hyp      <- 1e5
Av.hyp      <- 1e5
sigsq.gamma <- 1e10
sigsq.beta  <- 1e10

# Set colours.
mainCol <- "darkslateblue"
ptCol   <- "paleturquoise3"
lineCol <- "skyblue"
axisCol <- "black"

# Set plotting size parameters.
cex.pt      <- 0.75
cex.mainVal <- 1.7
cex.labVal  <- 1.3
xlabVal     <- "x"

# Source MFVB files here.
mse.mu.uneven.mfvb=0
mse.mu.uneven.haar=0
mse.mu.uneven.s8=0

mse.sd.uneven.mfvb=0
mse.sd.uneven.haar=0
mse.sd.uneven.s8=0
mse.sd.uneven.j=0
mse.sd.uneven.s8.j=0

lmse.sd.uneven.mfvb=0
lmse.sd.uneven.haar=0
lmse.sd.uneven.s8=0
lmse.sd.uneven.j=0
lmse.sd.uneven.s8.j=0

for (j in 1:100) {

# GENERATE DATA SET
# -----------------
n <- 500
set.seed(3*j)
xOrig <- runif(n)
set.seed(3*j)
yOrig <- fTrue(xOrig) + sqrt(exp(loggTrue(xOrig)))*rnorm(n)

aOrig <- min(xOrig) ; bOrig <- max(xOrig)
mean.x <- mean(xOrig) ; sd.x <- sd(xOrig)
mean.y <- mean(yOrig) ; sd.y <- sd(yOrig)

a <- (aOrig - mean.x)/sd.x ; b <- (bOrig - mean.x)/sd.x
x <- (xOrig - mean.x)/sd.x ; y <- (yOrig - mean.y)/sd.y

numIntKnotsU <- 17
intKnotsU <- quantile(x,seq(0,1,length=numIntKnotsU+2)[-c(1,numIntKnotsU+2)])
Zu <- ZOSull(x,intKnots=intKnotsU,range.x=c(a,b))
numKnotsU <- ncol(Zu)

numIntKnotsV <- numIntKnotsU 
intKnotsV <- quantile(x,seq(0,1,length=numIntKnotsV+2)[-c(1,numIntKnotsV+2)])
Zv <- ZOSull(x,intKnots=intKnotsV,range.x=c(a,b))
numKnotsV <- ncol(Zv) 

# RUN MEAN FIELD VARIATIONAL BAYES
# --------------------------------
X <- cbind(rep(1,n),x)
Cumat <- cbind(X,Zu)
Cvmat <- cbind(X,Zv)
ncX <- ncol(X)
ncZu <- ncol(Zu)
ncZv <- ncol(Zv)
ncCu <- ncol(Cumat)
ncCv <- ncol(Cvmat)

MFVBfit <- meanVarMFVB(y,X,ncZu,ncZv,Au.hyp,Av.hyp,
                       sigsq.gamma,sigsq.beta)

ng <- 201
xgOrig <- seq(aOrig,bOrig,length=ng)
xg <- (xgOrig - mean.x)/sd.x
Xg <- cbind(rep(1,ng),xg)
Zug <- ZOSull(xg,intKnots=intKnotsU,range.x=c(a,b))
Cug <- cbind(Xg,Zug)
Zvg <- ZOSull(xg,intKnots=intKnotsV,range.x=c(a,b))
Cvg <- cbind(Xg,Zvg)

mu.q.nu <- MFVBfit$mu.q.nu
mu.q.omega <- MFVBfit$mu.q.omega
Sigma.q.nu <- MFVBfit$Sigma.q.nu
Sigma.q.omega <- MFVBfit$Sigma.q.omega

fhatMFVBg <- Cug%*%mu.q.nu

fhatMFVBgOrig <- fhatMFVBg*sd.y + mean.y

logghatMFVBg <- Cvg%*%mu.q.omega 
logghatMFVBgOrig <- logghatMFVBg + 2*log(sd.y)

sdloggMFVBgOrig <- sqrt(diag(Cvg%*%Sigma.q.omega%*%t(Cvg))) 
credLowloggMFVBgOrig <- logghatMFVBgOrig - qnorm(0.975)*sdloggMFVBgOrig
credUpploggMFVBgOrig <- logghatMFVBgOrig + qnorm(0.975)*sdloggMFVBgOrig

sqrtghatMFVBg <- exp(0.5*Cvg%*%mu.q.omega 
                   + 0.125*diag(Cvg%*%Sigma.q.omega%*%t(Cvg)))

sqrtghatMFVBgOrig <- sqrtghatMFVBg*sd.y

x.mod=unique(sort(xOrig))
y.mod=0
for(i in 1:length(x.mod)){
  y.mod[i]=median(yOrig[xOrig==x.mod[i]])
}

y.exp=c(y.mod,y.mod[length(y.mod):(2*length(y.mod)-2^9+1)])
y.final=c(y.exp,y.exp[length(y.exp):1])

mu.est=smash.gaus(y.final,filter.number=1,family="DaubExPhase")
var.est=smash.gaus(y.final,v.est=TRUE)
mu.est=mu.est[1:500]
var.est=var.est[1:500]

mu.est.inter=approx(x.mod,mu.est,xgOrig,'linear')$y
var.est.inter=approx(x.mod,var.est,xgOrig,'linear')$y

mse.mu.uneven.mfvb[j]=mean((fhatMFVBgOrig-fTrue(xgOrig))^2)
mse.mu.uneven.haar[j]=mean((mu.est.inter-fTrue(xgOrig))^2)

mse.sd.uneven.mfvb[j]=mean((sqrtghatMFVBgOrig-exp((loggTrue(xgOrig))/2))^2)
mse.sd.uneven.haar[j]=mean((sqrt(var.est.inter)-exp((loggTrue(xgOrig))/2))^2)
lmse.sd.uneven.mfvb[j]=mean((log(sqrtghatMFVBgOrig)-(loggTrue(xgOrig))/2)^2)
lmse.sd.uneven.haar[j]=mean((log(sqrt(var.est.inter))-(loggTrue(xgOrig))/2)^2)

mu.est=smash.gaus(y.final,filter.number=8,family="DaubLeAsymm")
var.est=smash.gaus(y.final,v.est=TRUE,v.basis=TRUE,filter.number=8,
                   family="DaubLeAsymm")
var.est.s8.j=smash.gaus(y.final,v.est=TRUE,v.basis=TRUE,jash=TRUE,
                        filter.number=8,family="DaubLeAsymm")
var.est.j=smash.gaus(y.final,v.est=TRUE,jash=TRUE)
mu.est=mu.est[1:500]
var.est=var.est[1:500]
var.est.j=var.est.j[1:500]
var.est.s8.j=var.est.s8.j[1:500]

mu.est.inter=approx(x.mod,mu.est,xgOrig,'linear')$y
var.est.inter=approx(x.mod,var.est,xgOrig,'linear')$y
var.est.inter.j=approx(x.mod,var.est.j,xgOrig,'linear')$y
var.est.inter.s8.j=approx(x.mod,var.est.s8.j,xgOrig,'linear')$y

mse.mu.uneven.s8[j]=mean((mu.est.inter-fTrue(xgOrig))^2)
mse.sd.uneven.s8[j]=mean((sqrt(var.est.inter)-exp((loggTrue(xgOrig))/2))^2)
mse.sd.uneven.j[j]=mean((sqrt(var.est.inter.j)-exp((loggTrue(xgOrig))/2))^2)
mse.sd.uneven.s8.j[j] =
  mean((sqrt(var.est.inter.s8.j) - exp((loggTrue(xgOrig))/2))^2)
lmse.sd.uneven.s8[j]=mean((log(sqrt(var.est.inter))-(loggTrue(xgOrig))/2)^2)
lmse.sd.uneven.j[j]=mean((log(sqrt(var.est.inter.j))-(loggTrue(xgOrig))/2)^2)
lmse.sd.uneven.s8.j[j] =
  mean((log(sqrt(var.est.inter.s8.j))-(loggTrue(xgOrig))/2)^2)

print(j)
}

mean(mse.mu.uneven.mfvb)
mean(mse.mu.uneven.haar)
mean(mse.mu.uneven.s8)
mean(mse.sd.uneven.mfvb)
mean(mse.sd.uneven.haar)
mean(mse.sd.uneven.s8)
mean(mse.sd.uneven.j)
mean(mse.sd.uneven.s8.j)

stop()

############################################################################

mse.mu.even.mfvb=0
mse.mu.even.haar=0
mse.mu.even.s8=0

mse.sd.even.mfvb=0
mse.sd.even.haar=0
mse.sd.even.s8=0
mse.sd.even.j=0
mse.sd.even.s8.j=0
lmse.sd.even.mfvb=0
lmse.sd.even.haar=0
lmse.sd.even.s8=0
lmse.sd.even.j=0
lmse.sd.even.s8.j=0

for(j in 1:100){

# Set hyperparameters:

Au.hyp <- 1e5
Av.hyp <- 1e5
sigsq.gamma <- 1e10
sigsq.beta <- 1e10

# Set colours:

mainCol <- "darkslateblue"
ptCol <- "paleturquoise3"
lineCol <- "skyblue"
axisCol <- "black"

# Set plotting size parameters:

cex.pt <- 0.75
cex.mainVal <- 1.7
cex.labVal <- 1.3

xlabVal <- "x"

# Generate data:

n <- 2^10
xOrig <- (1:n)/n
set.seed(30*j)
yOrig <- fTrue(xOrig) + sqrt(exp(loggTrue(xOrig)))*rnorm(n)

aOrig <- min(xOrig) ; bOrig <- max(xOrig)
mean.x <- mean(xOrig) ; sd.x <- sd(xOrig)
mean.y <- mean(yOrig) ; sd.y <- sd(yOrig)

a <- (aOrig - mean.x)/sd.x ; b <- (bOrig - mean.x)/sd.x
x <- (xOrig - mean.x)/sd.x ; y <- (yOrig - mean.y)/sd.y



numIntKnotsU <- 17
intKnotsU <- quantile(x,seq(0,1,length=numIntKnotsU+2)[-c(1,numIntKnotsU+2)])
Zu <- ZOSull(x,intKnots=intKnotsU,range.x=c(a,b))
numKnotsU <- ncol(Zu)

numIntKnotsV <- numIntKnotsU 
intKnotsV <- quantile(x,seq(0,1,length=numIntKnotsV+2)[-c(1,numIntKnotsV+2)])
Zv <- ZOSull(x,intKnots=intKnotsV,range.x=c(a,b))
numKnotsV <- ncol(Zv) 

# Do mean field variational Bayes:

X <- cbind(rep(1,n),x)
Cumat <- cbind(X,Zu)
Cvmat <- cbind(X,Zv)
ncX <- ncol(X)
ncZu <- ncol(Zu)
ncZv <- ncol(Zv)
ncCu <- ncol(Cumat)
ncCv <- ncol(Cvmat)

MFVBfit <- meanVarMFVB(y,X,ncZu,ncZv,Au.hyp,Av.hyp,
                       sigsq.gamma,sigsq.beta)


ng <- 2^10
xgOrig <- seq(aOrig,bOrig,length=ng)
xg <- (xgOrig - mean.x)/sd.x
Xg <- cbind(rep(1,ng),xg)
Zug <- ZOSull(xg,intKnots=intKnotsU,range.x=c(a,b))
Cug <- cbind(Xg,Zug)
Zvg <- ZOSull(xg,intKnots=intKnotsV,range.x=c(a,b))
Cvg <- cbind(Xg,Zvg)

mu.q.nu <- MFVBfit$mu.q.nu
mu.q.omega <- MFVBfit$mu.q.omega
Sigma.q.nu <- MFVBfit$Sigma.q.nu
Sigma.q.omega <- MFVBfit$Sigma.q.omega

# Plot the mean function estimate.

fhatMFVBg <- Cug%*%mu.q.nu

fhatMFVBgOrig <- fhatMFVBg*sd.y + mean.y


logghatMFVBg <- Cvg%*%mu.q.omega 
logghatMFVBgOrig <- logghatMFVBg + 2*log(sd.y)

sdloggMFVBgOrig <- sqrt(diag(Cvg%*%Sigma.q.omega%*%t(Cvg))) 
credLowloggMFVBgOrig <- logghatMFVBgOrig - qnorm(0.975)*sdloggMFVBgOrig
credUpploggMFVBgOrig <- logghatMFVBgOrig + qnorm(0.975)*sdloggMFVBgOrig

sqrtghatMFVBg <- exp(0.5*Cvg%*%mu.q.omega 
                   + 0.125*diag(Cvg%*%Sigma.q.omega%*%t(Cvg)))

sqrtghatMFVBgOrig <- sqrtghatMFVBg*sd.y



mu.est=smash.gaus(yOrig,filter.number=1,family="DaubExPhase")
var.est=smash.gaus(yOrig,v.est=TRUE)


mse.mu.even.mfvb[j]=mean((fhatMFVBgOrig-fTrue(xgOrig))^2)
mse.mu.even.haar[j]=mean((mu.est-fTrue(xgOrig))^2)

mse.sd.even.mfvb[j]=mean((sqrtghatMFVBgOrig-exp((loggTrue(xgOrig))/2))^2)
mse.sd.even.haar[j]=mean((sqrt(var.est)-exp((loggTrue(xgOrig))/2))^2)
lmse.sd.even.mfvb[j]=mean((log(sqrtghatMFVBgOrig)-(loggTrue(xgOrig))/2)^2)
lmse.sd.even.haar[j]=mean((log(sqrt(var.est))-(loggTrue(xgOrig))/2)^2)

mu.est=smash.gaus(yOrig,filter.number=8,family="DaubLeAsymm")
var.est=smash.gaus(yOrig,v.est=TRUE,v.basis=TRUE,filter.number=8,family="DaubLeAsymm")
var.est.s8.j=smash.gaus(yOrig,v.est=TRUE,v.basis=TRUE,jash=TRUE,filter.number=8,family="DaubLeAsymm")
var.est.j=smash.gaus(yOrig,v.est=TRUE,jash=TRUE)

mse.mu.even.s8[j]=mean((mu.est-fTrue(xgOrig))^2)
mse.sd.even.s8[j]=mean((sqrt(var.est)-exp((loggTrue(xgOrig))/2))^2)
mse.sd.even.j[j]=mean((sqrt(var.est.j)-exp((loggTrue(xgOrig))/2))^2)
mse.sd.even.s8.j[j]=mean((sqrt(var.est.s8.j)-exp((loggTrue(xgOrig))/2))^2)
lmse.sd.even.s8[j]=mean((log(sqrt(var.est))-(loggTrue(xgOrig))/2)^2)
lmse.sd.even.j[j]=mean((log(sqrt(var.est.j))-(loggTrue(xgOrig))/2)^2)
lmse.sd.even.s8.j[j]=mean((log(sqrt(var.est.s8.j))-(loggTrue(xgOrig))/2)^2)

print(j)
}

mean(mse.mu.even.mfvb)
mean(mse.mu.even.haar)
mean(mse.mu.even.s8)
mean(mse.sd.even.mfvb)
mean(mse.sd.even.haar)
mean(mse.sd.even.s8)
mean(mse.sd.even.j)
mean(mse.sd.even.s8.j)

save.image("comparison.RData")

#################

xgrid = (0:10000)/10000

pdf("paper/mfvb_eg.pdf", width = 8, height = 5)
plot(xgrid, fTrue(xgrid), type = "l", ylim = c(-5, 5), ylab = "y(X)",
     xlab = "X")
lines(xgrid, fTrue(xgrid) + 2 * sqrt(gTrue(xgrid)), col = 2, lty = 2)
lines(xgrid, fTrue(xgrid) - 2 * sqrt(gTrue(xgrid)), col = 2, lty = 2)
dev.off()
