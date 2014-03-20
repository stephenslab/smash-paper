inter.samp=function(x,y,nn){
  n=length(x)
  xmax=max(x)
  xmin=min(x)
  xgrid.old=c(0,(x-xmin)*(n-1)/(xmax-xmin)+1,n+1)
  mgrid=(xgrid.old[1:(n+1)]+xgrid.old[2:(n+2)])/2
  mgrid.l=mgrid[1:n]
  mgrid.u=mgrid[2:(n+1)]
  delta=n/nn
  xgrid.new=0.5+(1:nn-0.5)*delta
  y.new=0
  for(i in 1:nn){
    y.new[i]=y[xgrid.new[i]>=mgrid.l&xgrid.new[i]<mgrid.u]
  }
  return(list(x=xgrid.new,y=y.new))
}