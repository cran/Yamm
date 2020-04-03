\name{makeplot}
\alias{makeplot}

\title{
Plot Two-dimensional Quantile
}
\description{
This function calls \code{\link{muqie}} for multiple values of quantiles from 0.5 to 0.95 and then produces a set of plots with these quantiles for producing an animated GIF using package \pkg{animation}.
}
\usage{
makeplot(xdata, dm=c(1,2), nsegs=20,
         quantile.increment= 0.001, 
         nprojs=2000, reltol=0.001)
}

\arguments{
  \item{xdata}{
The data as a matrix or dataframe with the number of columns greater than or equal to two, with each row being viewed as one multivariate observation.
}
  \item{dm}{
A numeric vector with two entries representing the selected columns of the data considered. The default value is \eqn{c(1,2)}, which means the first two columns of data are chosen if the dimension of data is greater than two.
}
  \item{nsegs}{
The number of the unit-length direction vectors \code{u}, which is computed by dividing a unit circle into \code{nsegs} equal sectors.
}
  \item{quantile.increment}{
A numeric value specifies the increment of the set of different quantiles.
}
  \item{nprojs}{
The number of projections for the dataset when computing \code{yamm}. The default value is 2000.
}
 \item{reltol}{
The tolerance of the optimisation process in the function \code{yamm}. The default value is 0.001.
}

}

\value{
This function returns a set of plots with various specified quantiles.
}


\seealso{
\code{\link{yamm}}
\code{\link{muqie}}
}
\examples{

# Load a data frame with 103 rows and 2 columns.
# The last two rows of the data are the outliers.
data(clusters2d)
#
# Remove the outliers of the dataset.
cluster_without_outlier <- clusters2d[c(1:101),]
#
# Produce an animation of a set of multivariate quantile plots.

if (requireNamespace("animation"))	{

library("animation")
# Generate temporary file
f <- tempfile(fileext=".gif")
#
# Now generate movie into the temporary file.
# Here nprojs=40: for a real example, for production quality you should increase
# it to 1000, 2000 or even higher
#
# Here quantile.increment=0.1, for production quality this should be reduced to
# e.g. 0.01, of even smaller
#
saveGIF(makeplot(cluster_without_outlier, nprojs=40, quantile.increment=0.1), 
diff.col=3, interval=0.1,width=500, height=500, movie.name=f)
cat("Movie saved to: ", f, "\n")
}
}

\references{
Chen, F. and Nason, Guy P. (2020) A new method for computing the projection medi
an, its influence curve and techniques for the production of projected quantile
plots. \emph{PLOS One}, (to appear)
}

\keyword{ yamm }
