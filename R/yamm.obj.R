yamm.obj <- function (x, mu, nprojs = 2000, doabs = 0)
{
  nor <- dim(x)[[1]]
  noc <- dim(x)[[2]]
  
  y <- rep(0, nor)
  a <- rep(0, noc)
  x <-as.matrix(x)
  xs <- x	# Gets overwritten
  
  if (noc != length(mu))
    stop("Dimensionality of mu and x not the same.")
  x <- t(x)
  ans <- .C("Cyammobj", x=as.double(x), mu=as.double(mu),
            result = as.double(0), nprojs = as.integer(nprojs),
            nor=as.integer(nor), noc=as.integer(noc),
            doabs=as.integer(doabs), y=as.double(y), a=as.double(a),
            xs=as.double(xs))
  
  return(ans[["result"]])
}
