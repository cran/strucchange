\name{boundary.efp}
\alias{boundary}
\alias{boundary.efp}
\title{Boundary Function for Structural Change Tests}
\description{A generic function computing boundaries for structural change
tests}
\usage{
boundary(x, ...)
boundary(x, alpha = 0.05, alt.boundary = FALSE)}

\arguments{
  \item{x}{an object of class \code{"efp"}.}
  \item{alpha}{numeric from interval (0,1) indicating the confidence level for
     which the boundary of the corresponding test will be computed.}
  \item{alt.boundary}{logical. If set to \code{TRUE} alternative boundaries
     (instead of the standard linear boundaries) will be computed (for CUSUM
     processes only).}
}

\value{an object of class \code{"ts"} with the same time properties as
the time series in \code{x}}

\author{Achim Zeileis \email{zeileis@ci.tuwien.ac.at}}

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
lines(bound, col=2)
lines(-bound, col=2)
}
\keyword{regression}
