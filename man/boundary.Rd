\name{boundary}
\alias{boundary}
\title{Boundary Function for Structural Change Tests}
\description{A generic function computing boundaries for structural change
tests}
\usage{
boundary(x, ...)}

\arguments{
  \item{x}{an object. Use \code{\link{methods}} to see which
    \code{\link{class}} has a method for boundary.}
  \item{...}{additional arguments affecting the boundary.}
}

\value{an object of class \code{"ts"} with the same time properties as
the time series in \code{x}}

\seealso{\code{\link{boundary.efp}}, \code{\link{boundary.mefp}},
  \code{\link{boundary.Fstats}}}

\keyword{regression}
