\name{confint.breakpointsfull}
\alias{confint.breakpointsfull}
\alias{lines.confint.breakpoints}
\alias{print.confint.breakpoints}
\title{Confidence Intervals for Breakpoints}
\description{
  Computes confidence intervals for breakpoints.
}
\usage{
\method{confint}{breakpointsfull}(object, parm = NULL, level = 0.95,
    breaks = NULL, het.reg = TRUE, het.err = TRUE, vcov = NULL, sandwich = TRUE, ...)
\method{lines}{confint.breakpoints}(x, col = 2, angle = 90, length = 0.05,
    code = 3, at = NULL, breakpoints = TRUE, ...)
}
\arguments{
  \item{object}{an object of class \code{"breakpointsfull"} as computed by
    \code{\link{breakpoints}} from a \code{formula}.}
  \item{parm}{the same as \code{breaks}, only one of the two should be
    specified.}
  \item{level}{the confidence level required.}
  \item{breaks}{an integer specifying the number of breaks to be used.
    By default the breaks of the minimum BIC partition are used.}
  \item{het.reg}{logical. Should heterogenous regressors be assumed? If set
    to \code{FALSE} the distribution of the regressors is assumed to be
    homogenous over the segments.}
  \item{het.err}{logical. Should heterogenous errors be assumed? If set
    to \code{FALSE} the distribution of the errors is assumed to be
    homogenous over the segments.}
  \item{vcov}{a function to extract the covariance matrix
     for the coefficients of a fitted model of class \code{"lm"}.}
  \item{sandwich}{logical. Is the function \code{vcov} the sandwich
     estimator or only the middle part?}
  \item{x}{an object of class \code{"confint.breakpoints"} as returned by
    \code{confint}.}
  \item{col, angle, length, code}{arguments passed to \code{\link{arrows}}.}
  \item{at}{position on the y axis, where the confidence arrows should be
    drawn. By default they are drawn at the bottom of the plot.}
  \item{breakpoints}{logical. If \code{TRUE} vertical lines for the breakpoints
    are drawn.}
  \item{\dots}{\emph{currently not used}.}
}

\details{
  As the breakpoints are integers (observation numbers) the corresponding
  confidence intervals are also rounded to integers.

  The distribution function used for the computation of confidence
  intervals of breakpoints is given in Bai (1997). The procedure, in
  particular the usage of heterogenous regressors and/or errors, is
  described in more detail in Bai & Perron (2003).

  The breakpoints should be computed from a formula with \code{breakpoints},
  then the confidence intervals for the breakpoints can be derived by
  \code{confint} and these can be visualized by \code{lines}. For an
  example see below.
}

\value{
  A matrix containing the breakpoints and their lower and upper
  confidence boundary for the given level.
}

\references{
Bai J. (1997), Estimation of a Change Point in Multiple Regression Models,
  \emph{Review of Economics and Statistics}, \bold{79}, 551-563.

Bai J., Perron P. (2003), Computation and Analysis of Multiple Structural Change
  Models, \emph{Journal of Applied Econometrics}, \bold{18}, 1-22.
}

\seealso{\code{\link{breakpoints}}}

\examples{
if(! "package:stats" \%in\% search()) library(ts)

## Nile data with one breakpoint: the annual flows drop in 1898
## because the first Ashwan dam was built
data(Nile)
plot(Nile)

## dating breaks
bp.nile <- breakpoints(Nile ~ 1)
ci.nile <- confint(bp.nile, breaks = 1)
lines(ci.nile)
}

\keyword{regression}
