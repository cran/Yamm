makeplot <- function(xdata, dm=c(1,2), nsegs=20, quantile.increment= 0.001, nprojs=2000, reltol=0.001){
  
  pp <- seq(from=0.5, to=0.95, by=quantile.increment)
  npp <- length(pp)
  
  
  for(i in 1:npp)	{
    muqie(xdata, dm=dm, nsegs=nsegs, plot.it=TRUE, probs=pp[i], nprojs=nprojs, reltol=reltol)
    text(10, 5, paste("alpha=",pp[i]))
  }
}
