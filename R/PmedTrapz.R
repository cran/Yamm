PmedTrapz <-
function(x,no.subinterval){
  dx <- dim(x)
  nor <- dx[[1]]
  noc <- dx[[2]]
  length.sub <- length(no.subinterval)
  if (noc > 3 || noc < 2)
    stop("The number of columns of data matrix has to be 2 or 3.")
  if ((noc-1)!= length.sub)
    stop("The length of the division vector should equal to the number of columns of data matrix minus one.")
  
  x <- t(x)
  if (noc==2){
    out<-rep(0.0,2)
    ans <- .C("CPmedTrapz2D",
              as.double(x),
              result=as.double(out),
              as.integer(no.subinterval),
              as.integer(nor)
    )
  }
  if (noc==3){
    out<-rep(0.0,3)
    ans <- .C("CPmedTrapz3D",
              as.double(x),
              result=as.double(out),
              as.integer(no.subinterval[1]),
              as.integer(no.subinterval[2]),
              as.integer(nor)
    )
  }
  ans[["result"]]
}
