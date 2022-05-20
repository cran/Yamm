\name{beetle}
\alias{beetle}
\docType{data}
\title{
Six Measurements of Beetles
}
\description{
Multivariate dataset containing six measurements on each of three species of flea-beetles: concinna, heptapotamica, and heikertingeri.
The original data set contains one column identifying the species of the observations, which is irrelevant
and has been deleted here.
}
\usage{data("beetle")}
\format{
  A data frame with 74 observations on the following 6 variables.
  \describe{
    \item{\code{tars1}}{Width of the first joint of the first tarsus in microns (the sum of measurements for both tarsi).}
    \item{\code{tars2}}{The same for the second joint.}
    \item{\code{head}}{The maximal width of the head between the external edges of the eyes in 0.01 mm.}
    \item{\code{aede1}}{The maximal width of the aedeagus in the fore-part in microns.}
    \item{\code{aede2}}{The front angle of the aedeagus (1 unit = 7.5 degrees).}
    \item{\code{aede3}}{The aedeagus width from the side in microns.}
  }
}

\source{Lubischew, A.A.(1962) On the Use of Discriminant Functions in Taxonomy, \emph{Biometrics},\bold{18}, 455-477.
}
\references{
Cook, D.H. and Swayne, D.F. (2007). Interactive and Dynamic Graphics for Data Analysis: With Examples Using R and GGobi. \url{http://ggobi.org/book/data/flea.xml}
}
\examples{
data(beetle)
}
\keyword{datasets}
