\name{supLM}
\alias{supLM}
\alias{maxMOSUM}

\title{Generators for efpFunctionals along Continuous Variables}

\description{
Generators for \code{efpFunctional} objects suitable for aggregating
empirical fluctuation processes to test statistics along continuous
variables (i.e., along time in time series applications).
}

\usage{
supLM(from = 0.15, to = NULL) 

maxMOSUM(width = 0.15)
}

\arguments{
  \item{from, to}{numeric from interval (0, 1) specifying start and end
    of trimmed sample period. By default, \code{to} is \code{1 - from}, i.e.,
    with the default \code{from = 0.15} the first and last 15 percent of
    observations are trimmed.}
  \item{width}{a numeric from interval (0,1) specifying the bandwidth.
    Determines the size of the moving data window relative to sample size.}
}

\details{
  \code{supLM} and \code{maxMOSUM} generate \code{\link{efpFunctional}}
  objects for Andrews' supLM test and a (maximum) MOSUM test, respectively,
  with the specified optional parameters (\code{from} and \code{to},
  and \code{width}, respectively). The resulting objects can be used in
  combination with empirical fluctuation processes of class \code{\link{gefp}}
  for significance testing and visualization. The corresponding statistics
  are useful for carrying out structural change tests along a continuous
  variable (i.e., along time in time series applications). Further typical
  \code{\link{efpFunctional}}s for this setting are the double-maximum
  functional \code{\link{maxBB}} and the Cramer-von Mises functional
  \code{\link{meanL2BB}}.
}

\value{
  An object of class \code{efpFunctional}.
}

\references{
Merkle E.C., Zeileis A. (2013), Tests of Measurement Invariance without Subgroups:
A Generalization of Classical Methods. \emph{Psychometrika}, \bold{78}(1), 59--82.
doi:10.1007/S11336-012-9302-4

Zeileis A. (2005), A Unified Approach to Structural Change Tests Based on
ML Scores, F Statistics, and OLS Residuals. \emph{Econometric Reviews}, \bold{24},
445--466. doi:10.1080/07474930500406053.

Zeileis A. (2006), Implementing a Class of Structural Change Tests: An
Econometric Computing Approach. \emph{Computational Statistics & Data Analysis}, 
\bold{50}, 2987--3008. doi:10.1016/j.csda.2005.07.001.

Zeileis A., Hornik K. (2007), Generalized M-Fluctuation Tests for Parameter
Instability, \emph{Statistica Neerlandica}, \bold{61}, 488--508.
doi:10.1111/j.1467-9574.2007.00371.x.
}

\seealso{\code{\link{efpFunctional}}, \code{\link{gefp}}}

\examples{
## seatbelt data
data("UKDriverDeaths")
seatbelt <- log10(UKDriverDeaths)
seatbelt <- cbind(seatbelt, lag(seatbelt, k = -1), lag(seatbelt, k = -12))
colnames(seatbelt) <- c("y", "ylag1", "ylag12")
seatbelt <- window(seatbelt, start = c(1970, 1), end = c(1984,12))

## empirical fluctuation process
scus.seat <- gefp(y ~ ylag1 + ylag12, data = seatbelt)

## supLM test
plot(scus.seat, functional = supLM(0.1))
## MOSUM test
plot(scus.seat, functional = maxMOSUM(0.25))
## double maximum test
plot(scus.seat)
## range test
plot(scus.seat, functional = rangeBB)
## Cramer-von Mises statistic (Nyblom-Hansen test)
plot(scus.seat, functional = meanL2BB)
}

\keyword{regression}
