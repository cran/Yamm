\name{yamm}
\alias{yamm}

\title{
Yet Another Multivariate Median
}
\description{
Another method for computing the projection median for any dimensional dataset. Basically, it minimises the objective function \code{\link{yamm.obj}} over a unit hypersphere and finds the optimal shift vector \code{mu} in \code{yamm.obj}. \code{\link{optim}} in the \pkg{stats} package is used in this function to minimise \code{yamm.obj}.  
}
\usage{
yamm(x, nprojs = 2000, reltol = 1e-6, abstol=-Inf, xstart = l1median(x),
     opt.method = "BFGS", doabs = 0, full.results=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The data as a matrix or data frame, with each row being viewed as one multivariate observation.
}
  \item{nprojs}{
The number of projections for the shifted data matrix while using the Monte Carlo method to approximate the integration. The default value is 2000, more projections may be required for complicated data to ensure accuracy, which, however, increases the running time.
}
  \item{reltol}{
The tolerance of the optimisation process gets supplied to control arguments of \code{optim}. The default value is \eqn{1e-6}. Loosening tolerance will make the running process faster. Generally, \eqn{1e-3} is enough to obtain a good approximation for a short running time.
}
  \item{abstol}{
The absolute convergence tolerance of the optimisation process gets supplied to control arguments of \code{optim}. The default value is negative infinity.
}
  \item{xstart}{
The starting value for the optimiser. The default value is Spatial median of the data using function \code{l1median}. Other multivariate medians or mean values can also be used. Note, you should be aware of the outliers when using the mean values as a starting point, which may slow down the optimisation process or result in a less accurate median.
}
  \item{opt.method}{
The method chosen for the optimiser when computing the \code{yamm}, with default function \dQuote{\code{BFGS}}. Apart from \dQuote{\code{BFGS}}, other functions in \code{optim} like \dQuote{\code{Nelder-Mead}}, \dQuote{\code{CG}}, \dQuote{\code{L-BFGS-B}}, and \dQuote{\code{SANN}} can also be used.

}
  \item{doabs}{
If 0 (default), the function \code{yamm.obj} integrates the square of the univariate median of the projection to the shifted data set over a unit hypersphere; if 1, \code{yamm.obj} integrates the absolute value of the univariate median instead.
}

 \item{full.results}{
Logical. If \code{FALSE} (default), the function \code{yamm} only returns the best set of yamm location estimator found; if \code{TRUE}, a list of full reults from the function \code{optim} is displayed.
}
}

\value{
If \code{full.results} = \code{FALSE}, it returns the best set of yamm location estimator found, otherwise, it returns a list comprising of
\item{par }{The best set of parameters found, which is the yamm location estimator.}
\item{value }{The value of objective function \code{yamm.obj} corresponding to \code{par}.}
\item{counts }{A two-element integer vector giving the number of calls to the objective function and gradient of the function respectively. This excludes those calls needed to compute the Hessian, if requested, and any calls to the objective function to compute a finite-difference approximation to the gradient.}
\item{convergence }{An integer code. 0 indicates successful completion (which is always the case for method \dQuote{\code{SANN}} and \dQuote{\code{Brent}}). Possible error codes are

1
indicates that the iteration limit had been reached.

10
indicates degeneracy of the \code{Nelder-Mead} simplex.

51
indicates a warning from the \dQuote{\code{L-BFGS-B}} method; see component \code{message} for further details.

52
indicates an error from the \dQuote{\code{L-BFGS-B}} method; see component \code{message} for further details.}
\item{message }{A character string giving any additional information returned by the optimiser, or \code{NULL}}
}


\seealso{
\code{\link{yamm.obj}},
\code{\link{optim}}.
}

\examples{
data(beetle)
#
# Set seed for reproduction.
set.seed(5)
#
# Yamm approximated using 1000 projections.
yamm(beetle,nprojs = 1000,reltol=1e-3,doabs=0,full.results=TRUE)
#
# $par
# [1] 180.30601 124.23781  50.16349 135.53947  13.45252  95.64742
#
# $value
# [1] 5.704375
#
# $counts
# function gradient 
#      69        4 
#
# $convergence
# [1] 0
#
# $message
# NULL
}
\references{
Chen, F. and Nason, Guy P. (2020) A new method for computing the projection medi
an, its influence curve and techniques for the production of projected quantile
plots. \emph{PLOS One}, (to appear)
}

\keyword{ yamm }
