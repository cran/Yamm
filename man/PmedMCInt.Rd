\name{PmedMCInt}
\alias{PmedMCInt}

\title{
Projection Median Approximated by Monte Carlo Integration
}
\description{
This function approximates the projection median using Monte Carlo integration, which can be used for any dimensions. \code{PmedMCInt} is implemented internally using C code \code{CPmedMCInt} and hence is much faster than coding with R only.
}
\usage{
PmedMCInt(x, nprojs = 20000)
}

\arguments{
  \item{x}{
The data as a matrix or data frame, with each row being viewed as one multivariate observation.
}
  \item{nprojs}{
The number of projections when using the Monte Carlo method to approximate the integration. The default value is 20000, since \code{PmedMCInt} requires large a number of projections to ensure the accuracy. More projections may increase the accuracy, as well as the running time.
}
}
\details{
The projection median was introduced by Durocher and Kirkpatrick (2009) and generalised by Basu, Bhattacharya and Talukdar (2012). \code{PmedMCInt} produces the projection median using Monte Carlo approximation, which is valid in any multi-dimensional data. However, a large number of projections is sometimes required to ensure accuracy, which will also increase the running time. In this case, \code{PmedTrapz} is preferred for the two- or three-dimensional data, which is fast and accurate in general. In higher dimensions, \code{\link{yamm}} is another alternative for computing the projection median.
}
\value{
A vector of the projection median for \eqn{n}-dimensional data.
}
\references{
Durocher, S. and Kirkpatrick, D. (2009), The projection median of a set of points, \emph{Computational Geometry},\bold{42}, 364-375.

Basu, R., Bhattacharya, B.B., and Talukdar, T. (2012) The projection median of a set of points in Rd \emph{CCCG.}, \bold{47}, 329-346. \doi{10.1007/s00454-011-9380-6}
}


\seealso{
\code{\link{PmedTrapz}},
\code{\link{yamm}}
}
\examples{
# Load a 2-dimensional data set
data(clusters2d)
#
# Set seed for reproduction.
set.seed(5)
#
# Projection median approximated by Monte Carlo Integration.
PmedMCInt(clusters2d, nprojs = 50000)
# [1]  4.3246488 -0.1535201
#
#
# Load a 6-dimensional data set
data(beetle)
#
set.seed(5)
PmedMCInt(beetle,nprojs = 150000)
# [1] 179.92439 125.16939  50.01176 136.55460  13.22277  95.04224
}

\keyword{ yamm }

