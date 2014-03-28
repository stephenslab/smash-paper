x=X.s

  n = length(x)
  J = log2(n)

  wmean = matrix(0,J,n) 

    var.est=(x-mu.t)^2

    vtable=cxxtitable(2/3*var.est^2)$sumtable
    vdtable=cxxtitable(var.est)$difftable
    vrtable=cxxtirtable(var.est)
    fac1=c(3,2,1.5)
    fac2=c(2,1.5,1)
    #fac1=rep(1,3)
    #fac2=rep(1,3)
j=10

      if(j>=0&j<=1){
        virtable=vrtable[j+2,]
        virtable=fac1[1]*(abs(virtable)>=log(4))+fac1[2]*(abs(virtable)>=log(2)&abs(virtable)<log(4))+fac1[3]*(abs(virtable)<log(2))
      }else if(j>=2&j<=5){
        virtable=vrtable[j+2,]
        virtable=fac2[1]*(abs(virtable)>=log(4))+fac2[2]*(abs(virtable)>=log(2)&abs(virtable)<log(4))+fac2[3]*(abs(virtable)<log(2))
      }else{
        virtable=1
      }
      vt=vash(vtable[j+2,],10)$PosteriorMean


plot(vt,vtable[j+2,]/virtable)
abline(0,1,lty=2)

plot(vt,vtable[j+2,],xlim=c(6000,18000))
abline(0,1,lty=2)
