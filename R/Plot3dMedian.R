Plot3dMedian <-
function(data,xvec,yvec,zvec,yamm.nprojs=2000,PmedMCInt.nprojs=20000,no.subinterval=c(18,36),opt.method="BFGS",xlab="Component1",ylab="Component2",zlab="Component3"){

  if (dim(data)[[2]]!=3)
    stop("The dimensionality of the data should be 3.")
  if (length(no.subinterval)!=2)
    stop("The number of subintervals should be a vector of length 2.")
  
 # med.Oja <- ojaMedian(data)
  med.Spatial <- L1median(data)$estimate
  med.CWmed <- med(data, "CWmed")$median
  med.Tukey <- med(data, "Tukey",approx = TRUE)$median
  mean <- apply(data,2,mean)
  med.Pmed_MCInt <- PmedMCInt(x=data, nprojs=PmedMCInt.nprojs)
  med.Pmed_Trapz<-PmedTrapz(data,no.subinterval)
  med.Yamm <- yamm(x=data, nprojs=yamm.nprojs, opt.method=opt.method, reltol = 1e-3)
  
  
  ## plot
  layout(matrix(c(1,2,3,4,5,5), ncol=2, byrow=TRUE), heights=c(0.45, 0.45, 0.1))
  
  oldpar <- par(mai=rep(0.52, 4))
  on.exit(par(oldpar))
  
  plot(data[,1], data[,2], pch=20,col="grey",xlab=xlab, ylab=ylab, xlim=c(min(xvec),max(xvec)), ylim=c(min(yvec),max(yvec)))
  #points(med.Oja[1], med.Oja[2], col="hot pink", pch=1,cex=1.5,lwd=2)
  points(med.Spatial[1], med.Spatial[2], col="lime green", pch=2,cex=1.5,lwd=2)
  points(med.CWmed[1], med.CWmed[2], col="blue", pch=3,cex=1.5,lwd=2)
  points(med.Tukey[1], med.Tukey[2], col="purple", pch=4,cex=1.5,lwd=2)
  points(mean[1], mean[2], col="orange", pch=17,cex=1.5)
  points(med.Pmed_MCInt[1], med.Pmed_MCInt[2], col="turquoise1", pch=16,cex=1.5)
  points(med.Pmed_Trapz[1], med.Pmed_Trapz[2], col="aquamarine4", pch=15, cex=1.5)
  points(med.Yamm[1], med.Yamm[2], col="black", pch=8, bg=5,cex=1.5)
  
  plot(data[,3], data[,2], pch=20,col="grey",xlab=zlab, ylab=ylab, xlim=c(min(zvec),max(zvec)), ylim=c(min(yvec),max(yvec)))
  #points(med.Oja[3], med.Oja[2], col="hot pink", pch=1,cex=1.5,lwd=2)
  points(med.Spatial[3], med.Spatial[2], col="lime green", pch=2,cex=1.5,lwd=2)
  points(med.CWmed[3], med.CWmed[2], col="blue", pch=3,cex=1.5,lwd=2)
  points(med.Tukey[3], med.Tukey[2], col="purple", pch=4,cex=1.5,lwd=2)
  points(mean[3], mean[2], col="orange", pch=17,cex=1.5)
  points(med.Pmed_MCInt[3], med.Pmed_MCInt[2], col="turquoise1", pch=16,cex=1.5)
  points(med.Pmed_Trapz[3], med.Pmed_Trapz[2], col="aquamarine4", pch=15, cex=1.5)
  points(med.Yamm[3], med.Yamm[2], col="black", pch=8, bg=5,cex=1.5)
  
  plot(data[,1], data[,3], pch=20,col="grey",xlab=xlab, ylab=zlab, xlim=c(min(xvec),max(xvec)), ylim=c(min(zvec),max(zvec)))
  #points(med.Oja[1], med.Oja[3], col="hot pink", pch=1,cex=1.5,lwd=2)
  points(med.Spatial[1], med.Spatial[3], col="lime green", pch=2,cex=1.5,lwd=2)
  points(med.CWmed[1], med.CWmed[3], col="blue", pch=3,cex=1.5,lwd=2)
  points(med.Tukey[1], med.Tukey[3], col="purple", pch=4,cex=1.5,lwd=2)
  points(mean[1], mean[3], col="orange", pch=17,cex=1.5)
  points(med.Pmed_MCInt[1], med.Pmed_MCInt[3], col="turquoise1", pch=16,cex=1.5)
  points(med.Pmed_Trapz[1], med.Pmed_Trapz[3], col="aquamarine4", pch=15, cex=1.5)
  points(med.Yamm[1], med.Yamm[3], col="black", pch=8, bg=5,cex=1.5)
  
  
  plot.new()
  
  par(mai=c(0,0,0,0))
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  col1 <- c("hot pink","lime green","blue","purple","orange","black","turquoise1","aquamarine4")
  pch1 <- 1:4
  legend(x = "center", ncol=4, inset = 0, x.intersp = 0.6, y.intersp=0.8, bty="n", cex=1.4, legend = c("Oja","Spatial","CWmed","Tukey","Mean","Yamm","PmedMCInt","PmedTrapz"), col=col1, pch = c(pch1,17,8,16,15))
  
}
