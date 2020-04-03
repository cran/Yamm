Plot2dMedian <-
function(data,xvec,yvec,yamm.nprojs=2000,PmedMCInt.nprojs=20000,no.subinterval=36,opt.method="BFGS",xlab="Component1",ylab="Component2"){
  
  if (dim(data)[[2]]!=2)
    stop("The dimensionality of the data should be 2.")
  if (length(no.subinterval)!=1)
    stop("The length of the number of subintervals should be 1.")
  
  layout(matrix(c(1,2), ncol=1), heights=c(4.5,1))
  oldpar <- par(mai=c(0.8,0.8,0.3,0.5))
  on.exit(par(oldpar))
  plot(data[,1], data[,2], pch=20,col="grey",xlab=xlab,ylab=ylab, xlim=c(min(xvec),max(xvec)), ylim=c(min(yvec),max(yvec)))
  
  med.Liu <- med(data, "Liu")$median
  points(med.Liu[1], med.Liu[2], col="red", pch=0,cex=1.5,lwd=2)
  
  med.Oja <- ojaMedian(data)
  points(med.Oja[1], med.Oja[2], col="hot pink", pch=1,cex=1.5,lwd=2)
  
  med.Spatial <- l1median(data)
  points(med.Spatial[1], med.Spatial[2], col="lime green", pch=2,cex=1.5,lwd=2)
  
  med.CWmed <- med(data, "CWmed")$median
  points(med.CWmed[1], med.CWmed[2], col="blue", pch=3,cex=1.5,lwd=2)
  
  med.Tukey <- med(data, "Tukey", approx=TRUE)$median
  points(med.Tukey[1], med.Tukey[2], col="purple", pch=4,cex=1.5,lwd=2)
  
  
  mn1 = mean(data[,1])
  mn2 = mean(data[,2])
  points(mn1, mn2, col="orange", pch=17,cex=1.5)
  
  
  med.Pmed_MCInt <- PmedMCInt(x=data,nprojs = PmedMCInt.nprojs)
  points(med.Pmed_MCInt[1], med.Pmed_MCInt[2], col="turquoise1", pch=16, cex=1.5)
  
  med.Pmed_Trapz <- PmedTrapz(x=data,no.subinterval = no.subinterval)
  points(med.Pmed_Trapz[1], med.Pmed_Trapz[2], col="aquamarine4", pch=15, cex=1.5)
  
  
  med.Yamm <- yamm(x=data, nprojs=yamm.nprojs, opt.method=opt.method, reltol = 1e-3)
  points(med.Yamm[1], med.Yamm[2], col="black", pch=8,cex=1.5)
  
  par(mai=c(0,0,0,0))
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  
  col1 <- c("red","hot pink","lime green","blue","purple","orange","turquoise1","aquamarine4","black")
  pch1 <- 0:4
  legend("center", ncol=3, inset = 0, x.intersp = 0.6, y.intersp=0.7, bty="n", cex=1, legend = c("Liu","Oja","Spatial","CWmed","Tukey","Mean","PmedMCInt","PmedTrapz","Yamm"), col=col1, pch = c(pch1,17,16,15,8))
  
}
