\name{breakfactor}
\alias{breakfactor}
\title{Factor Coding of Segmentations}
\description{
 Generates a factor encoding the segmentation given by
 a set of breakpoints.
}
\usage{
breakfactor(obj, breaks = NULL, labels = NULL, ...)
}
\arguments{
  \item{obj}{An object of class \code{"breakpoints"} or
    \code{"breakpointsfull"} respectively.}
  \item{breaks}{an integer specifying the number of breaks
    to extract (only if \code{obj} is of class \code{"breakpointsfull"}),
    by default the minimum BIC partition is used.}
  \item{labels}{a vector of labels for the returned factor,
    by default the segments are numbered starting from
    \code{"segment1"}.}
  \item{\dots}{further arguments passed to \code{factor}.}
}

\value{
  A factor encoding the segmentation.
}

\seealso{\code{\link{breakpoints}}}

\examples{
if(! "package:stats" \%in\% search()) library(ts)

## Nile data with one breakpoint: the annual flows drop in 1898
## because the first Ashwan dam was built
data(Nile)
plot(Nile)

## compute breakpoints
bp.nile <- breakpoints(Nile ~ 1)

## fit and visualize segmented and unsegmented model
fm0 <- lm(Nile ~ 1)
fm1 <- lm(Nile ~ breakfactor(bp.nile, breaks = 1))

lines(fitted(fm0), col = 3)
lines(fitted(fm1), col = 4)
lines(bp.nile, breaks = 1)
}

\keyword{regression}
