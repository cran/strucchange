\name{plot.mefp}
\alias{plot.mefp}
\alias{lines.mefp}
\title{Plot Methods for mefp Objects}
\description{This is a method of the generic \code{\link{plot}} function for
  for \code{"mefp"} objects as returned by \code{\link{mefp}} or
  \code{\link{monitor}}. It plots the emprical fluctuation process (or a
  functional therof) as a time series plot, and includes boundaries
  corresponding to the significance level of the monitoring procedure.
}
\usage{
\method{plot}{mefp}(x, boundary=TRUE, functional="max", main=NULL,
    ylab="empirical fluctuation process", ylim=NULL, ...)
}
\arguments{
  \item{x}{an object of class \code{"mefp"}.}
  \item{boundary}{if \code{FALSE}, plotting of boundaries is suppressed.}
  \item{functional}{indicates which functional should be applied to a
    multivariate empirical process. If set to \code{NULL} all dimensions
    of the process (one process per coefficient in the linear model) are
    plotted. }
  \item{main, ylab, ylim, ...}{high-level \code{\link{plot}} function parameters.}
}
\seealso{\code{\link{mefp}}}
\examples{
df1 <- data.frame(y=rnorm(300))
df1[150:300,"y"] <- df1[150:300,"y"]+1
me1 <- mefp(y~1, data=df1[1:50,,drop=FALSE], type="ME", h=1,
              alpha=0.05)
me2 <- monitor(me1, data=df1)

plot(me2)
}
\keyword{hplot}
