muqie <- function(xdata, dm=c(1,2), probs=0.5, nsegs=20, nprojs=2000, reltol=0.001,
                  plot.it=FALSE, full.return=FALSE, xlab=NULL, ylab=NULL){
  
  if (ncol(xdata) < 2 )
    stop("Data matrix xdata should have at least two columns")
  
  xd.dnames <- dimnames(xdata)
  
  if (is.null(xlab))	{
    if (is.null(xd.dnames))
      xlab<-paste("Variable", dm[1])
    else
      xlab <- xd.dnames[[2]][dm[1]]
  }
  if (is.null(ylab))	{
    if (is.null(xd.dnames))
      ylab<-paste("Variable", dm[2])
    else
      ylab <- xd.dnames[[2]][dm[2]]
  }
  #
  # Now reduce xdata just down to the two variables we want
  #
  xdata <- xdata[, dm]
  xdata <- as.matrix(xdata)
  my.yamm <- yamm(xdata, nprojs=nprojs, reltol=reltol)
  
  cdata <- xdata - matrix(my.yamm, nrow=nrow(xdata), ncol=2, byrow=TRUE)
  
  cd.dnames <- list(xd.dnames[[1]], xd.dnames[[2]][dm])
  
  dimnames(cdata) <- cd.dnames
  
  thetavec <- seq(from=0, to=2*pi, length=nsegs+1)[1:nsegs]
  
  pmat <- rbind( cos(thetavec), sin(thetavec))
  
  uvd <- cdata %*% pmat
  
  
  ans <- rbind(thetavec, pmat, apply(uvd, 2, quantile, probs=probs))
  
  dimnames(ans)[[1]][1:4] <- c("theta", "x", "y", paste(as.character(100*probs), "%", sep=""))
  
  if (plot.it==TRUE){
    
    plot(xdata[,1], xdata[,2], type="n", xlab=xlab, ylab=ylab)
    points(xdata[,1], xdata[,2], col="gray")
    points(my.yamm[1], my.yamm[2], col="blue", pch=19, cex=1.5)
    
    
    
    polygon( my.yamm[1]+ans[4,]*ans[2,], my.yamm[2] + ans[4,]*ans[3,])
    
    x1 = my.yamm[1] + ans[4,]*ans[2,]
    y1 = my.yamm[2] + ans[4,]*ans[3,]
    
    segments(x0=my.yamm[1], y0=my.yamm[2], x1=x1, y1=y1, lty=2, col="gray")
    
    if (probs>=0.5)
      opp <- which (ans[4,]< 0.0)
    else
      opp <- which (ans[4,]>= 0.0)
    
    if (length(opp)>0){
      
      segments(x0=my.yamm[1], y0=my.yamm[2], x1=x1[opp], y1=y1[opp], lty=2,col="red")
    }
    
  }
  
  
  if (full.return==TRUE){
    ll <- list(ans=ans, uvd=uvd, cdata=cdata, yamm=my.yamm)
    return(ll)
  }
  else
    return(ans)
}
