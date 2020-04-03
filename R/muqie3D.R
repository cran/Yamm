muqie3D <-function(xdata, dm=c(1,2,3), probs=0.5, nsegs=30, nprojs=2000, reltol=0.001, 
                   plot.it=FALSE, full.return=FALSE){
  
  if (ncol(xdata) < 3 )
    stop("Data matrix xdata should have at least three columns")
  
  xd.dnames <- dimnames(xdata)
  

  # Reduce xdata just down to the two variables we want
  xdata <- xdata[, dm]
  xdata <- as.matrix(xdata)
  
  my.yamm <- yamm(xdata, nprojs=nprojs, reltol=reltol)
  
  cdata <- xdata - matrix(my.yamm, nrow=nrow(xdata), ncol=3, byrow=TRUE)
  
  cd.dnames <- list(xd.dnames[[1]], xd.dnames[[2]][dm])
  
  dimnames(cdata) <- cd.dnames
  
  thetavec <- seq(from=0, to=2*pi, length=nsegs+1)[1:nsegs]
  phivec <- seq(from=0, to=pi/2, length=nsegs+1)[1:nsegs]
  
  pmat <- matrix(0, nrow=3, ncol=nsegs*nsegs)
  
  count <- 1
  
  for(i in 1:nsegs)
    for(j in 1:nsegs)	{
      pmat[, count] <- c(sin(thetavec[i])*sin(phivec[j]),
                          cos(thetavec[i])*sin(phivec[j]),
                          cos(phivec[j]))
      count <- count+1
    }
  
  
  uvd <- cdata %*% pmat
  
  
  ans <- rbind(pmat, apply(uvd, 2, quantile, probs=probs))
  
  dimnames(ans)[[1]][1:4] <- c("x", "y", "z", paste(as.character(100*probs), "%", sep=""))
  
  xx <- ans[4,]*ans[1,]+my.yamm[1]
  yy <- ans[4,]*ans[2,]+my.yamm[2]
  zz <- ans[4,]*ans[3,]+my.yamm[3]
  
  zlim <- c(my.yamm[3], max(zz))
  
  if (plot.it==TRUE){
	
    persp(interp(xx, yy, zz, duplicate="strip"), shade=0.3, col=3, xlab="", axes=FALSE, box=FALSE, zlim=zlim)
    title(xlab=paste("Quantile: ", as.character(probs)))
  }

  
  if (full.return==TRUE){
    ll <- list(ans=ans, uvd=uvd, cdata=cdata, yamm=my.yamm)
    return(ll)
  }

  else
    return(ans)
}
