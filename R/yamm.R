yamm <- function (x, nprojs = 2000, reltol = 1e-06, abstol=-Inf,  xstart = l1median(x), 
                  opt.method = "BFGS", doabs = 0, full.results = FALSE) 
{
  #	cat("Reltol = ", reltol, ". Abstol = ", abstol, "\n")
  if (length(xstart) != ncol(x)) 
    stop("Length of xstart has to be same as number of columns of x")

  ans <- optim(par = xstart, fn = yamm.obj, x = x, nprojs = nprojs, 
               method = opt.method, control = list(reltol = reltol, abstol=abstol), 
               doabs = doabs)
  if (full.results == "FALSE") {
    ans <- ans$par
  }
  return(ans)
}