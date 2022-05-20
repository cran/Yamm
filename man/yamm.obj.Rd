\name{yamm.obj}
\alias{yamm.obj}

\title{
Objective Function for Yamm
}
\description{
The objective function when computing \code{\link{yamm}}, which is the integral of the squared or absolute value of the univariate median of the projection of the shifted data set over a unit hypersphere. It is implemented internally using C code \code{Cyammobj} and hence is much faster than coding with R only.
}
\usage{
yamm.obj(x, mu, nprojs = 2000, doabs = 0)
}

\arguments{
  \item{x}{
The data as a matrix or data frame, with each row being viewed as one multivariate observation.
}
  \item{mu}{
A shift vector with length \eqn{n}, where \eqn{n} should equal to the number of columns (variables) of the data matrix. Each row of the data matrix \code{x} is shifted by \code{mu} to obtain the shifted data matrix.
}
  \item{nprojs}{
The number of projections for the shifted data matrix while using the Monte Carlo method to approximate the integration. The default value is 2000.
}
  \item{doabs}{
If 0 (default), function \code{yamm.obj} integrates square of the univariate median of the projection to the shifted data set over a unit hypersphere; if 1, \code{yamm.obj} integrates absolute value of the univariate median instead.
}
}

\value{
A univariate integral of the squared or absolute value of the median of the projection of the shifted data set over a unit hypersphere is returned from the \code{.C} calling function
}


\seealso{
\code{\link{yamm}}
}
\examples{
data(beetle)
#
# Set seed for reproduction.
set.seed(5)
#
# Objective function for yamm with a chosen shift vector.
#
yamm.obj(beetle, mu=rep(10,6), nprojs=5000, doabs=1)
# [1] 88.38346
}

\references{
Chen, F. and Nason, Guy P. (2020) A new method for computing the projection medi
an, its influence curve and techniques for the production of projected quantile
plots. \emph{PLOS One}, \doi{10.1371/journal.pone.0229845}
}

\keyword{ yamm }
