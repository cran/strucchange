\name{Fstats}
\alias{Fstats}
\alias{print.Fstats}
\title{F Statistics}
\description{Computes a series of F statistics for a specified data window.}
\usage{
Fstats(formula, from = 0.15, to = NULL, data)}
\arguments{
 \item{formula}{a symbolic description for the model to be tested}
 \item{from, to}{numeric. If \code{from} is smaller than 1 they are
     interpreted as percentages of data and by default \code{to} is taken to be
     1 - \code{from}. F statistics will be calculated for the observations
     \code{(n*from):(n*to)}, when \code{n} is the number of observations in the
     model. If \code{from} is greater than 1 it is interpreted to be the index
     and \code{to} defaults to \code{n - from}. If \code{from} is a vector with
     two elements, then \code{from} and \{to} are interpreted as time
     specifications like in \code{\link{ts}}, see also the examples.}
 \item{data}{an optional data frame containing the variables in the model. By
     default the variables are taken from the environment which \code{Fstats} is
     called from.} }
\details{For every potential change point in \code{from:to} a F statistic (Chow
test statistic) is computed. For this an OLS model is fitted for the
observations before and after the potential change point, i.e. \code{2k}
parameters have to be estimated, and the error sum of squares is computed (ESS).
Another OLS model for all obervations with a restricted sum of squares (RSS) is
computed, hence \code{k} parameters have to be estimated here. If \code{n} is 
the number of observations and \code{k} the number of regressors in the model,
the formula is:

\deqn{F = \frac{(RSS - ESS)/k}{ESS/(n - 2 k)}}{F = (RSS-ESS)/ESS * (n-2*k)/k}
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

\author{Achim Zeileis \email{zeileis@ci.tuwien.ac.at}}

\seealso{\code{\link{plot.Fstats}}, \code{\link{sctest.Fstats}},
\code{\link{boundary.Fstats}}}

\examples{
## Load dataset "nhtemp" with average yearly temperatures in New Haven
data(nhtemp)
## plot the data
plot(nhtemp)

## test the model null hypothesis that the average temperature remains constant
## over the years for potential break points between 1941 (corresponds to from =
## 0.5) and 1962 (corresponds to to = 0.85)
## compute F statistics
fs <- Fstats(nhtemp ~ 1, from = 0.5, to = 0.85)
## this gives the same result
fs <- Fstats(nhtemp ~ 1, from = c(1941,1), to = c(1962,1))
## plot the F statistics
plot(fs, alpha = 0.01)
## and the corresponding p values
plot(fs, pval = TRUE, alpha = 0.01)
## perform the aveF test
sctest(fs, type = "aveF")
}
\keyword{regression}
