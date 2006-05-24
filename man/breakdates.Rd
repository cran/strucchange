\name{breakdates}
\alias{breakdates}
\alias{breakdates.breakpoints}
\alias{breakdates.confint.breakpoints}
\title{Breakdates Corresponding to Breakpoints}
\description{
 A generic function for computing the breakdates corresponding
 to breakpoints (and their confidence intervals).
}
\usage{
breakdates(obj, format.times = FALSE, ...)
}
\arguments{
  \item{obj}{An object of class \code{"breakpoints"}, \code{"breakpointsfull"} or their
             confidence intervals as returned by \code{\link{confint}}.}
  \item{format.times}{logical. If set to \code{TRUE} a vector of
    strings with the formatted breakdates. See details for more
    information.}
  \item{\dots}{currently not used.}
}
\details{
  Breakpoints are the number of observations that are the last in one
  segment and breakdates are the corresponding points on the underlying
  time scale. The breakdates can be formatted which enhances readability
  in particular for quarterly or monthly time series. For example the
  breakdate \code{2002.75} of a monthly time series will be formatted to
  \code{"2002(10)"}.
}

\value{
  A vector or matrix containing the breakdates.
}

\seealso{\code{\link{breakpoints}}, \code{\link{confint}}}

\examples{
## Nile data with one breakpoint: the annual flows drop in 1898
## because the first Ashwan dam was built
data("Nile")
plot(Nile)

bp.nile <- breakpoints(Nile ~ 1)
summary(bp.nile)
plot(bp.nile)

## compute breakdates corresponding to the
## breakpoints of minimum BIC segmentation
breakdates(bp.nile)

## confidence intervals
ci.nile <- confint(bp.nile)
breakdates(ci.nile)
ci.nile

plot(Nile)
lines(ci.nile)
}

\keyword{regression}
