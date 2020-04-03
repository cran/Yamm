\name{Plot3dMedian}
\alias{Plot3dMedian}

\title{
Plot Three-dimensional Medians
}
\description{
This function plots multivariate medians in the three-dimensional case. The grey dots presented in the figure are the data points and the Spatial, Component-wise (CWmed), Tukey's, Oja's, Liu's, Projection medians as well as the mean value of the data set are plotted in the figure.
}
\usage{
Plot3dMedian(data, xvec, yvec, zvec, yamm.nprojs = 2000,
             PmedMCInt.nprojs = 20000, no.subinterval = c(18,36),
             opt.method = "BFGS", xlab = "Component1", 
             ylab = "Component2", zlab = "Component3")
}
%- maybe also 'usage' for other objects documented here.
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
  \item{zvec}{
A numeric vector containing the maximum and minimum values you desire for the z-axis.
}
   \item{yamm.nprojs}{
The number of projections for the dataset when computing \code{\link{yamm}}. The default value is 2000.
}

 \item{PmedMCInt.nprojs}{
The number of projections for the dataset when computing \code{\link{PmedMCInt}}. The default value is 20000, since \code{PmedMCInt} requires large number of projections while doing the Monte Carlo integration to ensure accuracy.
}
  \item{no.subinterval}{
A numeric vector of two entries which represents the number of subintervals chosen while using the trapezoidal rule to approximate the projection median with \code{\link{PmedTrapz}} function. The default vector is \eqn{c(36,36)}. Note small values (e.g. less than 10) for each entry of \code{no.subinterval} should not be used, to safeguard accuracy.
}

 \item{opt.method}{
The method chosen for the optimiser when computing the \code{\link{yamm}}, with default function \dQuote{\code{BFGS}}. \code{\link{optim}} is used to minimise the objective function \code{\link{yamm.obj}}. Apart from \code{BFGS}, other functions in \code{optim} like \dQuote{\code{Nelder-Mead}}, \dQuote{\code{CG}}, \dQuote{\code{L-BFGS-B}}, and \dQuote{\code{SANN}} can also be used.
}
  \item{xlab}{
Title for x-axis. Must be a character string.
}
  \item{ylab}{
Title for y-axis. Must be a character string.
}
  \item{zlab}{
Title for z-axis. Must be a character string.
}
}
\details{
The Spatial median is obtained using \code{\link{l1median}} in the Rpackage \pkg{pacPP}. The Component-wise (CWmed), and Tukey's median are produced using function \code{\link{med}} in the Rpackage \pkg{depth}. Oja's median is produced using function \code{\link{ojaMedian}} in the Rpackage \code{OjaNP}. Liu's median is not available in higher dimensions (\eqn{>2}), so it is not shown here. When computing the projection median, three approximations are implemented and diplayed in the plot, where \code{\link{PmedMCInt}} uses Monte Carlo method, \code{\link{PmedTrapz}} is computed by the trapezoidal rule, and \code{\link{yamm}} uses an optimiser.

The argument \code{xvec}, \code{yvec} and \code{zvec} are useful when there are outliers in the dataset, which are not expected to be shown in the figure in some cases. Determining the x-axis y-axis, and z-axis allows you to zoom in the plot and see the difference between multivariate medians and mean value.
}


\seealso{
\code{\link{PmedMCInt}},
\code{\link{PmedTrapz}}
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
# Load a data frame with 105 rows and 3 columns.
# The last five rows of the data are the outliers.
data(clusters3d)
#
# Remove the outliers of the dataset.
cluster_without_outlier <- clusters3d[c(1:100),]
myxvec <- c(min(cluster_without_outlier[,1]), 
          max(cluster_without_outlier[,1]))
myyvec <- c(min(cluster_without_outlier[,2]), 
          max(cluster_without_outlier[,2]))
myzvec <- c(min(cluster_without_outlier[,3]), 
          max(cluster_without_outlier[,3]))
#
# Plot the figure.
set.seed(15)
Plot3dMedian(cluster_without_outlier, myxvec, myyvec, myzvec, 
             yamm.nprojs = 2000, PmedMCInt.nprojs = 15000, 
             no.subinterval = c(18,36),opt.method = "BFGS",
             xlab = "Component1",ylab = "Component2", 
             zlab = "Component3")
}
\keyword{ yamm }
