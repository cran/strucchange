\name{logLik.breakpoints}
\alias{logLik.breakpoints}
\alias{logLik.breakpointsfull}
\alias{AIC.breakpointsfull}
\title{Log Likelihood and Information Criteria for Breakpoints}

\description{
  Computation of log likelihood and AIC type information criteria
  for partitions given by breakpoints.
}

\usage{
\method{logLik}{breakpointsfull}(object, breaks = NULL, ...)
\method{AIC}{breakpointsfull}(object, breaks = NULL, ..., k = 2)
}

\arguments{
  \item{object}{an object of class \code{"breakpoints"} or \code{"breakpointsfull"}.}
  \item{breaks}{if \code{object} is of class \code{"breakpointsfull"} the
    number of breaks can be specified.}
  \item{\dots}{\emph{currently not used}.}
  \item{k}{the penalty parameter to be used, the default \code{k = 2}
    is the classical AIC, \code{k = log(n)} gives the BIC, if \code{n}
    is the number of observations.}
}

\details{
  As for linear models the log likelihood is computed on a normal model and
  the degrees of freedom are the number of regression coefficients multiplied
  by the number of segements plus the number of estimated breakpoints plus
  1 for the error variance.

  If \code{AIC} is applied to an object of class \code{"breakpointsfull"}
  \code{breaks} can be a vector of integers and the AIC for each corresponding
  partition will be returned. By default the maximal number of breaks stored
  in the \code{object} is used. See below for an example.
}

\value{
  An object of class \code{"logLik"} or a simple vector containing
  the AIC respectively.
}

\seealso{\code{\link{breakpoints}}}

\examples{
if(! "package:stats" \%in\% search()) library(ts)

## Nile data with one breakpoint: the annual flows drop in 1898
## because the first Ashwan dam was built
data(Nile)
plot(Nile)

bp.nile <- breakpoints(Nile ~ 1)
summary(bp.nile)
plot(bp.nile)

## BIC of partitions with0 to 5 breakpoints
plot(0:5, AIC(bp.nile, k = log(bp.nile$nobs)), type = "b")
## AIC
plot(0:5, AIC(bp.nile), type = "b")

## BIC, AIC, log likelihood of a single partition
bp.nile1 <- breakpoints(bp.nile, breaks = 1)
AIC(bp.nile1, k = log(bp.nile1$nobs))
AIC(bp.nile1)
logLik(bp.nile1)
}

\keyword{regression}
