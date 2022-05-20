makeplot3D <- function(xdata, dm=c(1,2,3), nsegs=30, quantile.increment= 0.005, 
                       nprojs=2000, reltol=0.001){
  
  pp <- seq(from=0.5, to=0.95, by=quantile.increment)
  npp <- length(pp)
  
  maxprob <- max(pp)
  
  
  for(i in 1:npp)	{
    muqie3D(xdata, dm=dm, nsegs=nsegs, plot.it=TRUE, probs=pp[i], nprojs=nprojs, reltol=reltol)
  }
  
}
