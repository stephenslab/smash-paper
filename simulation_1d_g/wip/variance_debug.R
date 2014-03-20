J=log2(n)
wmean=matrix(0,nrow=J,ncol=n)
wvmean=matrix(0,nrow=J,ncol=n)  
    x=X.s
    var.est=log((x)^2)-digamma(0.5)-log(2)
    var.var.est=rep(4.94,n)
    vtable=cxxtitable(var.var.est)$sumtable
    vdtable=cxxtitable(var.est)$difftable
    fac=c(2,1.5,1.3,1.15,1.1,1.07,rep(1,J-6))
    for(j in 0:(J-1)){
      zdat.ash=ash(vdtable[j+2,],sqrt((1.8*vtable[j+2,])), gridmult=0, pointmass=FALSE, prior="nullbiased")
      wmean[j+1,] = zdat.ash$PosteriorMean/2
      wvmean[j+1,] = zdat.ash$PosteriorSD^2/4
    }
    #for(j in (J-3):(J-1)){
    #  wmean[j+1,] = vdtable[j+2,]/2
    #}
    wwmean=-wmean
    wwvmean=wvmean
    #var.est=exp(cxxreverse.gwave(sum(log((var.est1.ini+var.est2.ini)/2)),wmean,wwmean))
    var.est=exp(cxxreverse.gwave(sum(log(sigma.t^2)),wmean,wwmean))
    var.var=exp(cxxreverse.gvwave(sum(var.var.est),wvmean,wwvmean)/2)
    var.est=var.est*var.var
plot((sigma.t^2),type='l')
lines((var.est),col=2)


    x=rnorm(4096,var2,1)
    var.est=((x))
    var.var.est=rep(1,4096)
    vtable=cxxtitable(var.var.est)$sumtable
    vdtable=cxxtitable(var.est)$difftable
    fac=c(2,1.5,1.3,1.15,1.1,1.07,rep(1,J-6))
    for(j in 0:(J-1)){
      zdat.ash=ash(vdtable[j+2,],sqrt((1*vtable[j+2,])), gridmult=0, pointmass=FALSE, prior="nullbiased")
      wmean[j+1,] = zdat.ash$PosteriorMean/2
    }
    #for(j in (J-3):(J-1)){
    #  wmean[j+1,] = vdtable[j+2,]/2
    #}
    wwmean=-wmean
    #var.est=exp(cxxreverse.gwave(sum(log((var.est1.ini+var.est2.ini)/2)),wmean,wwmean))
    var.est=(cxxreverse.gwave(sum(var2),wmean,wwmean))
plot(var2,type='l')
lines(var.est,col=2)



    var.est=x^2
    var.var.est=2/3*var.est^2
    vtable=cxxtitable(var.var.est)$sumtable
    vdtable=cxxtitable(var.est)$difftable
    vrtable=tirtable(var.est)
    virtable=3*(abs(vrtable)>=log(25))+2*(abs(vrtable)>=log(3)&abs(vrtable)<log(25))+1*(abs(vrtable)<log(3))
    fac=c(2,1.5,1.3,1.15,1.1,1.07,rep(1,J-6))
    for(j in 0:(J-1)){
      zdat.ash=ash(vdtable[j+2,],sqrt((fac[j+1]*vtable[j+2,]/virtable[j+2,])) gridmult=0,, pointmass=FALSE, prior="nullbiased")
      wmean[j+1,] = zdat.ash$PosteriorMean/2
    }
    #for(j in (J-3):(J-1)){
    #  wmean[j+1,] = vdtable[j+2,]/2
    #}
    wwmean=-wmean
    #var.est=exp(cxxreverse.gwave(sum(log((var.est1.ini+var.est2.ini)/2)),wmean,wwmean))
    var.est=(cxxreverse.gwave(sum(var2),wmean,wwmean))
plot(var2,type='l')
lines(var.est,col=2)


mse(var.est,var2)


var.t=log(sigma.t^2)
tvdtable=cxxtitable(var.t)$difftable

plot(as.vector(t(vdtable)),pch='.',cex=1)
points(as.vector(t(tvdtable[2:13,])),pch='.',cex=1,col=2)

j=0
      zdat.ash=ash(vdtable[j+2,],sqrt((1*vtable[j+2,])), pointmass=FALSE, localfdr = FALSE,nullcheck=TRUE, prior="nullbiased")
zdat.ash$fitted.g$pi