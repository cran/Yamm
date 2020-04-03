\name{Plot2dMedian}
\alias{Plot2dMedian}

\title{
Plot Two-dimensional Medians
}
\description{
This function plots various multivariate medians in the two-dimensional case. The grey dots presented in the figure are the data points and the Spatial, Component-wise (CWmed), Tukey's, Oja's, Liu's, Projection median as well as the mean value of the data set are plotted in the figure.
}
\usage{
Plot2dMedian(data, xvec, yvec, yamm.nprojs = 2000,
             PmedMCInt.nprojs = 20000, 
             no.subinterval = 36, opt.method = "BFGS", 
             xlab = "Component1", ylab = "Component2")
}

\arguments{
  \item{data}{
The data as a matrix or data frame, with each row being viewed as one multivariate observation.
}
  \item{xvec}{
A numeric vector containing the maximum and minimum values you desire for the x-axis. 
}
  \item{yvec}{
A numeric vector containing the maximum and minimum values you desire for the y-axis.
}
  \item{yamm.nprojs}{
The number of projections for the dataset when computing \code{\link{yamm}}. The default value is 2000.
}

 \item{PmedMCInt.nprojs}{
The number of projections for the dataset when computing \code{\link{PmedMCInt}}. The default value is 20000, since \code{PmedMCInt} requires large number of projections while doing the Monte Carlo integration to ensure accuracy.
}

  \item{no.subinterval}{
The number of subintervals while using the trapezoidal rule to approximate the projection median with \code{\link{PmedTrapz}} function. The default value is 36, and small values (e.g. less than 10) of \code{no.subinterval} should not be used, to safeguard accuracy.
}
  \item{opt.method}{
The method chosen for the optimiser when computing the \code{\link{yamm}}, with default function \dQuote{\code{BFGS}}. \code{\link{optim}} is used to minimise the objective function \code{\link{yamm.obj}}. Apart from \dQuote{\code{BFGS}}, other functions in \code{optim} like \dQuote{\code{Nelder-Mead}}, \dQuote{\code{CG}}, \dQuote{\code{L-BFGS-B}}, and \dQuote{\code{SANN}} can also be used.
}
  \item{xlab}{
Title for x-axis. Must be a character string.
}
  \item{ylab}{
Title for y-axis. Must be a character string.
}
}
\details{
The Spatial median is obtained using \code{\link{l1median}} in the Rpackage \pkg{pacPP}. The Component-wise (CWmed), Liu's and Tukey's median are produced using function \code{\link{med}} in the Rpackage \pkg{depth}. Oja's median is produced using function \code{\link{ojaMedian}} in the Rpackage \code{OjaNP}. When computing the projection median, three approximations are implemented and diplayed in the plot, where \code{\link{PmedMCInt}} uses Monte Carlo method, \code{\link{PmedTrapz}} is computed by the trapezoidal rule, and \code{\link{yamm}} uses an optimiser.

The argument \code{xvec} and \code{yvec} are useful when there are outliers in the data set, which are not expected to be shown in the figure in some cases. Determining the x-axis and y-axis allows you to zoom in the plot and see the difference between multivariate medians and mean value.
}



\seealso{
\code{\link{PmedTrapz}},
\code{\link{PmedMCInt}},
\code{\link{yamm}},
\code{\link{yamm.obj}}, 
\code{\link{optim}}.
}
\references{
Chen, F. and Nason, Guy P. (2020) A new method for computing the projection medi
an, its influence curve and techniques for the production of projected quantile
plots. \emph{PLOS One}, (to appear)
}
\examples{
# Load a data frame with 103 rows and 2 columns.
# The last two rows of the data are the outliers.
data(clusters2d)
#
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
\keyword{yamm}

