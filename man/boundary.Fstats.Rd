\name{boundary.Fstats}
\alias{boundary.Fstats}
\title{Boundary for F Statistics}
\description{Computes boundary for an object of class \code{"Fstats"}}
\usage{
\method{boundary}{Fstats}(x, alpha = 0.05, pval = FALSE, aveF = FALSE,
    asymptotic = FALSE, ...)}
\arguments{
 \item{x}{an object of class \code{"Fstats"}.}
 \item{alpha}{numeric from interval (0,1) indicating the confidence level for
     which the boundary of the supF test will be computed.}
 \item{pval}{logical. If set to \code{TRUE} a boundary for the corresponding p
     values will be computed.}
 \item{aveF}{logical. If set to \code{TRUE} the boundary of the aveF (instead of
    the supF) test will be computed. The resulting boundary then is a boundary
    for the mean of the F statistics rather than for the F statistics
    themselves.}
 \item{asymptotic}{logical. If set to \code{TRUE} the asymptotic (chi-square)
     distribution instead of the exact (F) distribution will be used to compute
     the p values (only if \code{pval} is \code{TRUE}).}
 \item{...}{currently not used.}
}

\value{an object of class \code{"ts"} with the same time properties as
the time series in \code{x}}

\seealso{\code{\link{Fstats}}, \code{\link{plot.Fstats}}}

\examples{
## Load dataset "nhtemp" with average yearly temperatures in New Haven
data(nhtemp)
## plot the data
plot(nhtemp)

## test the model null hypothesis that the average temperature remains constant
## over the years for potential break points between 1941 (corresponds to
## from = 0.5) and 1962 (corresponds to to = 0.85)
## compute F statistics
fs <- Fstats(nhtemp ~ 1, from = 0.5, to = 0.85)
## plot the p values without boundary
plot(fs, pval = TRUE, alpha = 0.01)
## add the boundary in another colour
lines(boundary(fs, pval = TRUE, alpha = 0.01), col = 2)
}
\keyword{regression}
