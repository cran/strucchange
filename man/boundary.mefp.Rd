\name{boundary.mefp}
\alias{boundary.mefp}
\title{Boundary Function for Monitoring of Structural Changes}
\description{Computes boundary for an object of class \code{"mefp"}}
\usage{
\method{boundary}{mefp}(x, ...)}
\arguments{
  \item{x}{an object of class \code{"mefp"}.}
  \item{...}{currently not used.}
}
\value{an object of class \code{"ts"} with the same time properties as
the monitored process}
\author{Friedrich Leisch}
\seealso{\code{\link{mefp}}, \code{\link{plot.mefp}}}
\examples{
df1 <- data.frame(y=rnorm(300))
df1[150:300,"y"] <- df1[150:300,"y"]+1
me1 <- mefp(y~1, data=df1[1:50,,drop=FALSE], type="ME", h=1,
              alpha=0.05)
me2 <- monitor(me1, data=df1)

plot(me2, boundary=FALSE)
lines(boundary(me2), col="green", lty="44")
}
\keyword{regression}

