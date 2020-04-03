\name{Yamm-package}
\alias{Yamm-package}
\alias{Yamm}
\docType{package}
\title{
\packageTitle{Yamm}
}
\description{
This package provides functions for computing the projection median. \code{PmedTrapz} approximates the projection median by the trapezoidal rule, which is only valid for the two- and three-dimensional cases, while \code{PmedMCInt} use Monte Carlo approximation, and it is valid for any multivariate median. \code{yamm} provides another method to compute the projection median based on an optimiser technique. This package also provides functions for plotting different multivariate medians, such as the Spatial, Component-wise, Tukey's, Oja's median, etc., for randomly generated data sets in both the two-dimensional and three-dimensional cases. In addition, this package also allows users to produce the two-dimensional and three-dimensional quantile plots with function \code{muqie} and \code{muqie3D} respectively.
}

\details{
The DESCRIPTION file:
\packageDESCRIPTION{Yamm}
\packageIndices{Yamm}
}

\author{
\packageAuthor{Yamm}

Maintainer: \packageMaintainer{Yamm}
}
\references{
Basu, R., Bhattacharya, B.B., and Talukdar, T. (2012) The projection median of a set of points in Rd \emph{CCCG.}, \bold{47}, 329-346. \doi{10.1007/s00454-011-9380-6}

Chen, F. and Nason, Guy P. (2020) A new method for computing the projection median, its influence curve and techniques for the production of projected quantile plots. \emph{PLOS One}, (to appear)

Croux, C., Filzmoser, P., and Oliveira, M., (2007). Algorithms for Projection-Pursuit Robust Principal Component Analysis,  \emph{Chemometrics and Intelligent Laboratory Systems}, \bold{87}, 218-225.

Durocher, S. and Kirkpatrick, D. (2009), The projection median of a set of points, \emph{Computational Geometry}, \bold{42}, 364-375.

Fischer,D., Mosler, K., Mottonen, J.K., Nordhausen, K., Pokotylo, O., and Vogel, D. (2016) Computing the Oja Median in R: The Package OjaNP, \emph{ArXiv:1606.07620}

Rousseeuw, P.J. and Ruts, I. (1996), Algorithm AS 307: Bivariate location depth, \emph{Appl. Stat.-J. Roy. St. C}, \bold{45}, 516-526.

Rousseeuw, P.J. and Ruts, I. (1998), Constructing the bivariate Tukey median, \emph{Stat. Sinica}, \bold{8}, 828-839.

Rousseeuw, P.J., Ruts, I., and Tukey, J.W. (1999), The Bagplot: A Bivariate Boxplot, \emph{The Am. Stat.}, \bold{53}, 382-387.

Struyf, A. and Rousseeuw, P.J. (2000), High-dimensional computation of the deepest location, \emph{Comput. Statist. Data Anal.}, \bold{34}, 415-436.
}

\keyword{ package }
\seealso{
\code{\link{yamm}},
\code{\link{PmedTrapz}},
\code{\link{PmedMCInt}},
}
\examples{
# Load a 2-dimensional data set.
data(clusters2d)
#
# Set seed for reproduction.
set.seed(5)
#
# Projection median approximated by Monte Carlo Integration.
PmedMCInt(clusters2d, nprojs = 30000)
# [1]  4.3369501 -0.1578591
#
#
# Projection median approximated by the trapezoidal rule.
PmedTrapz(clusters2d,no.subinterval=180)
# [1]  4.1556553 -0.3566614
#
#
# Yamm.
set.seed(5)
yamm(clusters2d,nprojs = 2500,reltol=1e-3,doabs=1,full.results=FALSE)
# [1]  4.3871582 -0.1070497
#
#
# Plot 2-D medians
# Remove the outliers of the dataset.
cluster_without_outlier <- clusters2d[c(1:101),]
myxvec <- c(min(cluster_without_outlier[,1]), 
          max(cluster_without_outlier[,1]))
myyvec <- c(min(cluster_without_outlier[,2]), 
          max(cluster_without_outlier[,2]))
#
# Plot the figure.
set.seed(5)
Plot2dMedian(clusters2d, myxvec, myyvec, yamm.nprojs = 2000, 
             PmedMCInt.nprojs = 20000, no.subinterval = 36, 
             opt.method = "BFGS", xlab = "Component1",
             ylab = "Component2")

}
