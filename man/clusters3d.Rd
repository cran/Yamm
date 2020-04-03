\name{clusters3d}
\alias{clusters3d}
\docType{data}
\title{
Four Clusters of 3-dimensional Data
}
\description{
This dataset with 105 observations contains four clusters, which are generated from different Laplace distributions randomly, and five outliers (located in the last five rows).
}
\usage{data("clusters3d")}
\format{
The four clusters are generated from different multivariate Laplace distributions. The first cluster has 20 observations, where the mean values \eqn{\mu} of the Laplace distribution are equal to \eqn{(-8, -8, -8)} and the covariance matrix \eqn{\Sigma} is the product of two times identity matrix. The second cluster has 35 observations, where \eqn{\mu = (-5, 5, 5)} and \eqn{\Sigma} is the identity matrix. The third cluster has 30 observations, where \eqn{\mu = (12, -12, 12)} and \eqn{\Sigma} is the identity matrix. The fourth cluster has 30 observations, where \eqn{\mu = (18, 18, -18)} and \eqn{\Sigma} is the identity matrix. The five outliers are from the \eqn{\mu = (100, 100, -100)} and \eqn{\Sigma} is the product of ten times identity matrix.
}

\references{
Chen, F. and Nason, Guy P. (2020) A new method for computing the projection medi
an, its influence curve and techniques for the production of projected quantile
plots. \emph{PLOS One}, (to appear)
}

\examples{
data(clusters3d)

}
\keyword{datasets}
