\name{boundary.efp}
\alias{boundary.efp}
\title{Boundary for Empirical Fluctuation Processes}
\description{Computes boundary for an object of class \code{"efp"}}
\usage{
\method{boundary}{efp}(x, alpha = 0.05, alt.boundary = FALSE,
   functional = "max", ...)
}

\arguments{
  \item{x}{an object of class \code{"efp"}.}
  \item{alpha}{numeric from interval (0,1) indicating the confidence level for
     which the boundary of the corresponding test will be computed.}
  \item{alt.boundary}{logical. If set to \code{TRUE} alternative boundaries
     (instead of the standard linear boundaries) will be computed (for Brownian
     bridge type processes only).}
  \item{functional}{indicates which functional should be applied to the
     empirical fluctuation process. See also \code{\link{plot.efp}}.}
  \item{\dots}{currently not used.}
}

\value{an object of class \code{"ts"} with the same time properties as
the process in \code{x}}

\seealso{\code{\link{efp}}, \code{\link{plot.efp}}}

\examples{
## Load dataset "nhtemp" with average yearly temperatures in New Haven
data(nhtemp)
## plot the data
plot(nhtemp)

## test the model null hypothesis that the average temperature remains constant
## over the years
## compute OLS-CUSUM fluctuation process
temp.cus <- efp(nhtemp ~ 1, type = "OLS-CUSUM")
## plot the process without boundaries
plot(temp.cus, alpha = 0.01, boundary = FALSE)
## add the boundaries in another colour
bound <- boundary(temp.cus, alpha = 0.01)
lines(bound, col=4)
lines(-bound, col=4)
}
\keyword{regression}
