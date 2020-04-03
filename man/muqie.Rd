\name{muqie}
\alias{muqie}

\title{
Two-dimensional Quantile
}
\description{
This function plots the collection of all MUltivariate QuantIlE points in two dimensions (\code{muqie}) over all unit-length direction vectors \code{u}, which projects the \code{\link{yamm}}-centred multivariate data onto the chosen vector \code{u} to obtain a univariate set. The \code{muqie} point is merely the vector \code{u} rescaled to have length equal to the quantile of the univariate set.
}
\usage{
muqie(xdata, dm=c(1,2), probs=0.5, nsegs=20,
      nprojs=2000, reltol=0.001, plot.it=FALSE, 
      full.return=FALSE, xlab=NULL, ylab=NULL)
}

\arguments{
  \item{xdata}{
The data as a matrix or dataframe with the number of columns greater than or equal to two, with each row being viewed as one multivariate observation.
}
  \item{dm}{
A numeric vector with two entries representing the selected columns of the data considered. The default value is \eqn{c(1,2)}, which means the first two columns of data are chosen if the dimension of data is more than two.
}
  \item{probs}{
The quantile of the data after projected to obtain a univariate set.
}
  \item{nsegs}{
The number of the two-dimensional unit-length direction vectors \code{u}, which is computed by dividing a unit circle into \code{nsegs} equal sectors.
}
 \item{nprojs}{
The number of projections for the dataset when computing \code{yamm}. The default value is 2000.
}
 \item{reltol}{
The tolerance of the optimisation process in the function \code{yamm}. The default value is 0.001.
}
\item{plot.it}{
Logical. If \code{TRUE}, the function \code{muqie} will produce a two-dimensional quantile plot.
}
\item{full.return}{
Logical. If \code{TRUE}, the function \code{muqie} will return a list of full results. See \dQuote{Value}.
}

\item{xlab}{
x-axis label for the quantile plot.
}
\item{ylab}{
y-axis label for the quantile plot.
}
}

\value{
If \code{full.results} = \code{TRUE}, it returns a list comprising of
\item{ans }{ A data matrix with four rows. The first row represents the angle of the unit-length projection vector \code{u} to the positive x-axis, while the second and third row are the x- and y-coordinates of the projection vector respectively. The last row is univariate quantile of the projected data matrix. 
}
\item{uvd }{ A data matrix after projecting the \code{yamm}-centred multivariate data onto a set of projection vectors \code{u}.}
\item{cdata }{The \code{yamm}-centred multivariate data matrix.}
\item{yamm }{The yamm value of the multivariate data. See \code{\link{yamm}} for more details.}
If \code{full.results} = \code{FALSE} (default), it will only return \code{ans}.
}


\seealso{
\code{\link{yamm}}
}
\examples{
data(beetle)
#
# Compute the 0.7-quantile for the first two columns of the beetle data. 
muqie(beetle,dm=c(1,4), probs=0.7)
}
\references{
Chen, F. and Nason, Guy P. (2020) A new method for computing the projection medi
an, its influence curve and techniques for the production of projected quantile
plots. \emph{PLOS One}, (to appear)
}

\keyword{yamm}
