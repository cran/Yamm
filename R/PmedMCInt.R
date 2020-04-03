PmedMCInt <-
function(x,nprojs=20000){
  nor<-dim(x)[1]
  noc<-dim(x)[2]
  x<-t(x)
  out<-rep(0.0,noc)
  ans <- .C("CPmedMCInt",
            as.double(x),
            result=as.double(out),
            as.integer(nprojs),
            as.integer(nor),
            as.integer(noc)
  )
  ans[["result"]]
}
