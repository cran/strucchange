\name{Fstats}
\alias{Fstats}
\alias{print.Fstats}
\title{F Statistics}
\description{Computes a series of F statistics for a specified data window.}
\usage{
Fstats(formula, from = 0.15, to = NULL, data,
    cov.type = c("const", "HC", "HC1"))}
\arguments{
 \item{formula}{a symbolic description for the model to be tested}
 \item{from, to}{numeric. If \code{from} is smaller than 1 they are
     interpreted as percentages of data and by default \code{to} is taken to be
     1 - \code{from}. F statistics will be calculated for the observations
     \code{(n*from):(n*to)}, when \code{n} is the number of observations in the
     model. If \code{from} is greater than 1 it is interpreted to be the index
     and \code{to} defaults to \code{n - from}. If \code{from} is a vector with
     two elements, then \code{from} and \code{to} are interpreted as time
     specifications like in \code{\link{ts}}, see also the examples.}
 \item{data}{an optional data frame containing the variables in the model. By
     default the variables are taken from the environment which \code{Fstats} is
     called from.}
 \item{cov.type}{a string indicating which type of covariance matrix
     estimator should be used. Constant homoskedastic variances are assumed
     if set to \code{"const"} and White's heteroskedasticity consistent
     estimator is used if set to \code{"HC"}. And \code{"HC1"} stands for a
     standardized estimator of \code{"HC"}, see also \code{\link{covHC}}.}
}
\details{For every potential change point in \code{from:to} a F statistic (Chow
test statistic) is computed. For this an OLS model is fitted for the
observations before and after the potential change point, i.e. \code{2k}
parameters have to be estimated, and the error sum of squares is computed (ESS).
Another OLS model for all obervations with a restricted sum of squares (RSS) is
computed, hence \code{k} parameters have to be estimated here. If \code{n} is 
the number of observations and \code{k} the number of regressors in the model,
the formula is:

\deqn{F = \frac{(RSS - ESS)}{ESS/(n - 2 k)}}{F = (RSS-ESS)/ESS * (n-2*k)}
}

\value{\code{Fstats} returns an object of class \code{"Fstats"}, which contains
mainly a time series of F statistics. The function \code{\link{plot}} has a
method to plot the F statistics or the
corresponding p values; with \code{sctest} a
supF-, aveF- or expF-test on structural change can be performed.}

\references{
Andrews D.W.K. (1993), Tests for parameter instability and structural
change with unknown change point, \emph{Econometrica}, \bold{61}, 821-856.

Hansen B. (1992), Tests for parameter instability in regressions with I(1)
processes, \emph{Journal of Business & Economic Statistics}, \bold{10}, 321-335.

Hansen B. (1997), Approximate asymptotic p values for structural-change
tests, \emph{Journal of Business & Economic Statistics}, \bold{15}, 60-67. }

\seealso{\code{\link{plot.Fstats}}, \code{\link{sctest.Fstats}},
\code{\link{boundary.Fstats}}}

\examples{
if(! "package:stats" \%in\% search()) library(ts)

## Nile data with one breakpoint: the annual flows drop in 1898
## because the first Ashwan dam was built
data(Nile)
plot(Nile)

## test the null hypothesis that the annual flow remains constant
## over the years
fs.nile <- Fstats(Nile ~ 1)
plot(fs.nile)
sctest(fs.nile)
## visualize the breakpoint implied by the argmax of the F statistics
plot(Nile)
lines(breakpoints(fs.nile))

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

## compute F statistics for potential breakpoints between
## 1971(6) (corresponds to from = 0.1) and 1983(6) (corresponds to
## to = 0.9 = 1 - from, the default)
## compute F statistics
fs <- Fstats(y ~ ylag1 + ylag12, data = seatbelt, from = 0.1)
## this gives the same result
fs <- Fstats(y ~ ylag1 + ylag12, data = seatbelt, from = c(1971, 6),
             to = c(1983, 6))
## plot the F statistics
plot(fs, alpha = 0.01)
## plot F statistics with aveF boundary
plot(fs, aveF = TRUE)
## perform the expF test
sctest(fs, type = "expF")
}
\keyword{regression}
