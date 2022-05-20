\name{PmedTrapz}
\alias{PmedTrapz}

\title{
Projection Median Approximated by Trapezoidal Rule
}
\description{
This function approximates the projection median using trapezoidal rule, which is only valid for the two- and three-dimensional cases. \code{PmedTrapz} is implemented internally using C code \code{CPmedTrapz2D} and \code{CPmedTrapz3D} and hence is much faster than coding with R only.
}
\usage{
PmedTrapz(x, no.subinterval)
}

\arguments{
  \item{x}{
The data as a matrix or data frame, with each row being viewed as one multivariate observation.
}
  \item{no.subinterval}{
A vector of subintervals chosen for implementing the trapezoidal rule. It is a number in the two-dimensional case, and has a length of two for the three-dimensional data, since the trapezoidal rule is only required once in 2D and needs to be applied twice for the double integral in 3D. Small values (e.g. less than 10) for each entry of \code{no.subinterval} should not be used, to safeguard the accuracy.
}
}
\details{
The projection median was introduced by Durocher and Kirkpatrick (2009) and generalised by Basu, Bhattacharya and Talukdar (2012). \code{PmedTrapz} produces the projection median directly from the definition using the trapezoidal rule, but current function is only valid in the two-dimensional and three-dimensional case. For more general dimensionalities, you can refer to function \code{\link{PmedMCInt}} and \code{\link{yamm}}. 
}
\value{
A vector of the projection median in the two or three dimensions.
}
\references{
Chen, F. and Nason, Guy P. (2020) A new method for computing the projection medi
an, its influence curve and techniques for the production of projected quantile
plots. \emph{PLOS One}, \doi{10.1371/journal.pone.0229845}

Durocher, S. and Kirkpatrick, D. (2009), The projection median of a set of points, \emph{Computational Geometry},\bold{42}, 364-375.

Basu, R., Bhattacharya, B.B., and Talukdar, T. (2012) The projection median of a set of points in Rd \emph{CCCG.}, \bold{47}, 329-346. \doi{10.1007/s00454-011-9380-6}
}


\seealso{
\code{\link{PmedMCInt}}
}
\examples{
# Load a 2-dimensional dataset
data(clusters2d)
#
# Projection median approximated by the trapezoidal rule.
PmedTrapz(clusters2d,no.subinterval=180)
# [1]  4.1556553 -0.3566614
#
# Load a 3-dimensional dataset
data(clusters3d)
#
PmedTrapz(clusters3d,c(180,360))
# [1] -0.906680  1.584866  2.695584
}


\keyword{ yamm }

