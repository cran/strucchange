\name{breakpoints}
\alias{breakpoints}
\alias{breakpoints.formula}
\alias{breakpoints.breakpointsfull}
\alias{breakpoints.Fstats}
\alias{summary.breakpoints}
\alias{summary.breakpointsfull}
\alias{plot.breakpointsfull}
\alias{plot.summary.breakpointsfull}
\alias{print.breakpoints}
\alias{print.summary.breakpointsfull}
\alias{lines.breakpoints}
\alias{coef.breakpointsfull}
\alias{vcov.breakpointsfull}
\alias{fitted.breakpointsfull}
\alias{residuals.breakpointsfull}
\encoding{latin1}
\title{Dating Breaks}

\description{
  Computation of breakpoints in regression relationships. Given a number
  of breaks the function computes the optimal breakpoints.
}

\usage{
\method{breakpoints}{formula}(formula, h = 0.15, breaks = NULL,
    data = list(), \dots)
\method{breakpoints}{breakpointsfull}(obj, breaks = NULL, \dots)
\method{summary}{breakpointsfull}(object, breaks = NULL, sort = TRUE,
    format.times = NULL, \dots)
\method{lines}{breakpoints}(x, breaks = NULL, lty = 2, \dots)

\method{coef}{breakpointsfull}(object, breaks = NULL, names = NULL, \dots)
\method{fitted}{breakpointsfull}(object, breaks = NULL, \dots)
\method{residuals}{breakpointsfull}(object, breaks = NULL, \dots)
\method{vcov}{breakpointsfull}(object, breaks = NULL, names = NULL,
    het.reg = TRUE, het.err = TRUE, vcov = NULL, sandwich = TRUE, \dots)
}

\arguments{
  \item{formula}{a symbolic description for the model in which breakpoints
    will be estimated.}
  \item{h}{minimal segment size either given as fraction relative to the
    sample size or as an integer giving the minimal number of observations
    in each segment.}
  \item{breaks}{integer specifying the maximal number of breaks to be calculated.
    By default the maximal number allowed by \code{h} is used.}
  \item{data}{an optional data frame containing the variables in the model. By
     default the variables are taken from the environment which \code{breakpoints} is
     called from.}
  \item{\dots}{currently not used.}
  \item{obj, object}{an object of class \code{"breakpointsfull"}.}
  \item{sort}{logical. If set to \code{TRUE} \code{summary} tries to match
    the breakpoints from partitions with different numbers of breaks.}
  \item{format.times}{logical. If set to \code{TRUE} a vector of
    strings with the formatted breakdates is printed. See \code{\link{breakdates}}
    for more information.}
  \item{x}{an object of class \code{"breakpoints"}.}
  \item{lty}{line type.}
  \item{names}{a character vector giving the names of the segments. If of length
    1 it is taken to be a generic prefix, e.g. \code{"segment"}.}
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
}

\details{
  All procedures in this package are concerned with testing or assessing
  deviations from stability in the classical linear regression model

  \deqn{y_i = x_i^\top \beta + u_i}{y_i = x_i' b + u_i}

  In many applications it is reasonable to assume
  that there are \eqn{m} breakpoints, where the coefficients shift from
  one stable regression relationship to a different one. Thus,
  there are \eqn{m+1} segments in which the regression coefficients are
  constant, and the model can be rewritten as

  \deqn{y_i = x_i^\top \beta_j + u_i
  \qquad (i = i_{j-1} + 1, \dots, i_j, \quad j = 1, \dots, m+1)}{y_i =
    x_i' b_j + u_i   (i = i_{j-1} + 1, \dots, i_j,   j = 1, \dots, m+1)}

  where \eqn{j} denotes the segment index. In practice the breakpoints \eqn{i_j}
  are rarely given exogenously, but have to be estimated.
  \code{breakpoints} estimates these breakpoints by minimizing the residual sum of
  squares (RSS) of the equation above.

  The foundation for estimating breaks in time series regression models
  was given by Bai (1994) and was extended to multiple breaks by Bai (1997ab)
  and Bai & Perron (1998). \code{breakpoints} implements the algorithm
  described in Bai & Perron (2003) for simultanous estimation of
  multiple breakpoints. The distribution function used for the confidence
  intervals for the breakpoints is given in Bai (1997b). The ideas behind
  this implementation are described in Zeileis et al. (2003).

  The algorithm for computing the optimal breakpoints given the number
  of breaks is based on a dynamic programming approach. The underlying
  idea is that of the Bellman principle. The main computational effort
  is to compute a triangular RSS matrix, which gives the residual
  sum of squares for a segment starting at observation \eqn{i} and
  ending at \eqn{i'} with \eqn{i} < \eqn{i'}.

  Given a \code{formula} as the first argument, \code{breakpoints} computes
  an object of class \code{"breakpointsfull"} which inherits from \code{"breakpoints"}.
  This contains in particular the triangular RSS
  matrix and functions to extract an optimal segmentation. A \code{summary}
  of this object will give the breakpoints (and associated) breakdates
  for all segmentations up to the maximal number of breaks together
  with the associated RSS and BIC. These will be plotted if \code{plot}
  is applied and thus visualize the minimum BIC estimator of the number
  of breakpoints. From an object of class \code{"breakpointsfull"} an
  arbitrary number of \code{breaks} (admissable by the minimum segment
  size \code{h}) can be extracted by another application of
  \code{breakpoints}, returning an object of class \code{"breakpoints"}.
  This contains only the breakpoints for the specified number of breaks
  and some model properties (number of observations, regressors, time
  series properties and the associated RSS) but not the triangular RSS
  matrix and related extractor functions. The set of breakpoints which
  is associated by default with a \code{"breakpointsfull"} object is
  the minimum BIC partition.

  Breakpoints are the number of observations that are the last in one
  segment, it is also possible to compute the corresponding \code{breakdates}
  which are the breakpoints on the underlying time scale. The breakdates
  can be formatted which enhances readability in particular for quarterly
  or monthly time series. For example the breakdate \code{2002.75} of a monthly
  time series will be formatted to \code{"2002(10)"}. See \code{\link{breakdates}}
  for more details.

  From a \code{"breakpointsfull"} object confidence intervals for the breakpoints
  can be computed using the method of \code{\link{confint}}.
  The breakdates corresponding to the breakpoints can again be computed
  by \code{\link{breakdates}}. The breakpoints and their confidence
  intervals can be visualized by \code{lines}. Convenience functions are
  provided for extracting the coefficients and covariance matrix, fitted 
  values and residuals of segmented models.

  The log likelihood as well as some information criteria can be computed
  using the methods for the \code{\link{logLik}} and \code{\link{AIC}}. As
  for linear models the log likelihood is computed on a normal model and
  the degrees of freedom are the number of regression coefficients multiplied
  by the number of segements plus the number of estimated breakpoints plus
  1 for the error variance. More details can be found on the help page of
  the method \code{\link{logLik.breakpoints}}.

  As the maximum of a sequence of F statistics is equivalent to the minimum
  OLS estimator of the breakpoint in a 2-segment partition it can be
  extracted by \code{breakpoints} from an object of class \code{"Fstats"}
  as computed by \code{\link{Fstats}}. However, this cannot be used to extract
  a larger number of breakpoints.

  For illustration see the commented examples below and Zeileis et al. (2003).
}
\section{value}{
  An object of class \code{"breakpoints"} is a list with the following
  elements:
  \describe{
   \item{breakpoints}{the breakpoints of the optimal partition with the
     number of breaks specified,}
   \item{RSS}{the associated RSS,}
   \item{nobs}{the number of observations,}
   \item{nreg}{the number of regressors,}
   \item{call}{the function call,}
   \item{datatsp}{the time series properties \code{tsp} of the data,
     if any, \code{c(0, 1, nobs)} otherwise.}
  }

  If applied to a \code{formula} as first argument, \code{breakpoints} returns an object of class
  \code{"breakpointsfull"} (which inherits from \code{"breakpoints"}), that
  contains some additional (or slightly different) elements such as:
  \describe{
   \item{breakpoints}{the breakpoints of the minimum BIC partition,}
   \item{RSS}{a function which takes two arguments \code{i,j} and computes
     the residual sum of squares for a segment starting at observation \code{i} and
     ending at \code{j} by looking up the corresponding element in the triangular
     RSS matrix \code{RSS.triang},}
   \item{RSS.triang}{a list encoding the triangular RSS matrix.}
  }
}

\references{
Bai J. (1994), Least Squares Estimation of a Shift in Linear Processes,
  \emph{Journal of Time Series Analysis}, \bold{15}, 453-472.

Bai J. (1997a), Estimating Multiple Breaks One at a Time,
  \emph{Econometric Theory}, \bold{13}, 315-352.

Bai J. (1997b), Estimation of a Change Point in Multiple Regression Models,
  \emph{Review of Economics and Statistics}, \bold{79}, 551-563.

Bai J., Perron P. (1998), Estimating and Testing Linear Models With Multiple Structural
   Changes, \emph{Econometrica}, \bold{66}, 47-78.

Bai J., Perron P. (2003), Computation and Analysis of Multiple Structural Change
  Models, \emph{Journal of Applied Econometrics}, \bold{18}, 1-22.

Zeileis A., Kleiber C., Krämer W., Hornik K. (2003), Testing and Dating of
  Structural Changes in Practice, \emph{Computational Statistics and Data Analysis},
  \bold{44}, 109-123.
}

\examples{
if(! "package:stats" \%in\% search()) library(ts)

## Nile data with one breakpoint: the annual flows drop in 1898
## because the first Ashwan dam was built
data(Nile)
plot(Nile)

## F statistics indicate one breakpoint
fs.nile <- Fstats(Nile ~ 1)
plot(fs.nile)
breakpoints(fs.nile)
lines(breakpoints(fs.nile))

## or
bp.nile <- breakpoints(Nile ~ 1)
summary(bp.nile)

## the BIC also chooses one breakpoint
plot(bp.nile)
breakpoints(bp.nile)

## fit null hypothesis model and model with 1 breakpoint
fm0 <- lm(Nile ~ 1)
fm1 <- lm(Nile ~ breakfactor(bp.nile, breaks = 1))
plot(Nile)
lines(ts(fitted(fm0), start = 1871), col = 3)
lines(ts(fitted(fm1), start = 1871), col = 4)
lines(bp.nile)

## confidence interval
ci.nile <- confint(bp.nile)
ci.nile
lines(ci.nile)


## UK Seatbelt data: a SARIMA(1,0,0)(1,0,0)_12 model
## (fitted by OLS) is used and reveals (at least) two
## breakpoints - one in 1973 associated with the oil crisis and
## one in 1983 due to the introduction of compulsory
## wearing of seatbelts in the UK.
data(UKDriverDeaths)
seatbelt <- log10(UKDriverDeaths)
seatbelt <- cbind(seatbelt, lag(seatbelt, k = -1), lag(seatbelt, k = -12))
colnames(seatbelt) <- c("y", "ylag1", "ylag12")
seatbelt <- window(seatbelt, start = c(1970, 1), end = c(1984,12))
plot(seatbelt[,"y"], ylab = expression(log[10](casualties)))

## testing
re.seat <- efp(y ~ ylag1 + ylag12, data = seatbelt, type = "RE")
plot(re.seat)

## dating
bp.seat <- breakpoints(y ~ ylag1 + ylag12, data = seatbelt, h = 0.1)
summary(bp.seat)
lines(bp.seat, breaks = 2)

## minimum BIC partition
plot(bp.seat)
breakpoints(bp.seat)
## the BIC would choose 0 breakpoints although the RE and supF test
## clearly reject the hypothesis of structural stability. Bai &
## Perron (2003) report that the BIC has problems in dynamic regressions.
## due to the shape of the RE process of the F statistics choose two
## breakpoints and fit corresponding models
bp.seat2 <- breakpoints(bp.seat, breaks = 2)
fm0 <- lm(y ~ ylag1 + ylag12, data = seatbelt)
fm1 <- lm(y ~ breakfactor(bp.seat2)/(ylag1 + ylag12) - 1, data = seatbelt)

## plot
plot(seatbelt[,"y"], ylab = expression(log[10](casualties)))
time.seat <- as.vector(time(seatbelt))
lines(time.seat, fitted(fm0), col = 3)
lines(time.seat, fitted(fm1), col = 4)
lines(bp.seat2)

## confidence intervals
ci.seat2 <- confint(bp.seat, breaks = 2)
ci.seat2
lines(ci.seat2)
}

\concept{breakpoint estimation}
\concept{changepoint estimation}
\concept{segmented regression}
\keyword{regression}
