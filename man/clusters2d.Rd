\name{clusters2d}
\alias{clusters2d}
\docType{data}
\title{
Three Clusters of 2-dimensional Data
}
\description{
This dataset with 103 observations contains three clusters, which are generated from different independent normal distributions randomly, and two outliers (located in the last two rows).
}
\usage{data("clusters2d")}
\format{
The first cluster has 26 observations, and the two variables are generated from \eqn{N(3, 1)} and \eqn{N(4, 1)} respectively. The second cluster has 36 observations, and the two variables are generated from \eqn{N(10, 1.5)} and \eqn{N(-2, 1.5)} respectively. The third cluster has 39 observations, and the two variables are generated from \eqn{N(2, 0.5)} and \eqn{N(-2, 0.5)} respectively. The two outliers are \eqn{c(100.3, 99.1)} and \eqn{c(97.5, 98.4)}.
}
\references{
Chen, F. and Nason, Guy P. (2020) A new method for computing the projection medi
an, its influence curve and techniques for the production of projected quantile
plots. \emph{PLOS One}, \doi{10.1371/journal.pone.0229845}
}

\examples{
data(clusters2d)
}
\keyword{datasets}
