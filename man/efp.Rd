\name{efp}
\alias{efp}
\alias{print.efp}
\encoding{latin1}
\title{Empirical Fluctuation Processes}
\description{Computes an empirical fluctuation process according
to a specified method from the generalized fluctuation test
framework, which includes CUSUM and MOSUM tests based on recursive
or OLS residuals, parameter estimates or ML scores (OLS first order
conditions).}
\usage{
efp(formula, data, type = , h = 0.15,
    dynamic = FALSE, rescale = TRUE)}
\arguments{
  \item{formula}{a symbolic description for the model to be tested.}
  \item{data}{an optional data frame containing the variables in the model. By
     default the variables are taken from the environment which \code{efp} is
     called from.}
  \item{type}{specifies which type of fluctuation process will be
     computed, the default is \code{"Rec-CUSUM"}. For details see below.}
  \item{h}{a numeric from interval (0,1) specifying the bandwidth. determines the
     size of the data window relative to sample size (for MOSUM and ME processes
     only).}
  \item{dynamic}{logical. If \code{TRUE} the lagged observations are included as
     a regressor.}
  \item{rescale}{logical. If \code{TRUE} the estimates will be standardized by
    the regressor matrix of the corresponding subsample according to Kuan & Chen
    (1994); if \code{FALSE} the whole regressor matrix will be used.
    (only if \code{type} is either \code{"RE"} or \code{"ME"})}
}

\details{If \code{type} is one of \code{"Rec-CUSUM"}, \code{"OLS-CUSUM"},
\code{"Rec-MOSUM"} or \code{"OLS-MOSUM"} the function \code{efp} will return a
one-dimensional empirical process of sums of residuals. Either it will be based
on recursive residuals or on OLS residuals and the process will contain
CUmulative SUMs or MOving SUMs of residuals in a certain data window.
For the MOSUM and ME processes all estimations are done for the
observations in a moving data window, whose size is determined by \code{h} and
which is shifted over the whole sample.

If \code{type} is either \code{"RE"} or \code{"ME"} a
\emph{k}-dimensional process will be returned, if \emph{k} is the number of
regressors in the model, as it is based on recursive OLS estimates of the
regression coefficients or moving OLS estimates respectively. The recursive
estimates test is also called fluctuation test, therefore setting \code{type}
to \code{"fluctuation"} was used to specify it in earlier versions of
strucchange. It still can be used now, but will be forced to \code{"RE"}.

If \code{type} is \code{"Score-CUSUM"} or \code{"Score-MOSUM"} a \emph{k+1}-dimensional
process will be returned, one for each score of the regression coefficients and one for
the scores of the variance. The process gives the decorrelated cumulative sums of the ML
scores (in a Gaussian model) or first order conditions respectively (in an OLS framework).

If there is a single structural change point \eqn{t^*}, the recursive CUSUM path
starts to depart from its mean 0 at \eqn{t^*}. The Brownian bridge type paths
will have their respective peaks around \eqn{t^*}.
The Brownian bridge increments type paths should have a strong change at \eqn{t^*}.

The function \code{\link{plot}}
has a method to plot the empirical fluctuation process; with
\code{sctest} the corresponding test on structural change can be
performed.
}

\value{
  \code{efp} returns a list of class \code{"efp"} with components including:
  \item{process}{the fitted empirical fluctuation process of class
    \code{"ts"} or \code{"mts"} respectively,}
  \item{type}{a string with the \code{type} of the process fitted,}
  \item{nreg}{the number of regressors,}
  \item{nobs}{the number of observations,}
  \item{par}{the bandwidth \code{h} used.}
}

\references{

Brown R.L., Durbin J., Evans J.M. (1975), Techniques for
testing constancy of regression relationships over time, \emph{Journal of the
Royal Statistical Society}, B, \bold{37}, 149-163.

Chu C.-S., Hornik K., Kuan C.-M. (1995), MOSUM tests for parameter
constancy, \emph{Biometrika}, \bold{82}, 603-617.

Chu C.-S., Hornik K., Kuan C.-M. (1995), The moving-estimates test for
parameter stability, \emph{Econometric Theory}, \bold{11}, 669-720.

Hansen B. (1992), Testing for Parameter Instability in Linear Models,
\emph{Journal of Policy Modeling}, \bold{14}, 517-533.

Hjort N.L., Koning A. (2002), Tests for Constancy of Model Parameters
Over Time, \emph{Nonparametric Statistics}, \bold{14}, 113-132.

Krämer W., Ploberger W., Alt R. (1988), Testing for structural change in
dynamic models, \emph{Econometrica}, \bold{56}, 1355-1369.

Kuan C.-M., Hornik K. (1995), The generalized fluctuation test: A
unifying view, \emph{Econometric Reviews}, \bold{14}, 135 - 161.

Kuan C.-M., Chen (1994), Implementing the fluctuation and moving estimates
tests in dynamic econometric models, \emph{Economics Letters}, \bold{44},
235-239.

Ploberger W., Krämer W. (1992), The CUSUM test with OLS residuals,
\emph{Econometrica}, \bold{60}, 271-285.

Zeileis A., Leisch F., Hornik K., Kleiber C. (2002), \code{strucchange}:
An R Package for Testing for Structural Change in Linear Regression Models,
\emph{Journal of Statistical Software}, \bold{7}(2), 1-38.
URL \url{http://www.jstatsoft.org/v07/i02/}.

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

\seealso{\code{\link{gefp}}, \code{\link{plot.efp}}, \code{\link{print.efp}},
\code{\link{sctest.efp}}, \code{\link{boundary.efp}}}

\examples{
## Nile data with one breakpoint: the annual flows drop in 1898
## because the first Ashwan dam was built
data("Nile")
plot(Nile)

## test the null hypothesis that the annual flow remains constant
## over the years
## compute OLS-based CUSUM process and plot
## with standard and alternative boundaries
ocus.nile <- efp(Nile ~ 1, type = "OLS-CUSUM")
plot(ocus.nile)
plot(ocus.nile, alpha = 0.01, alt.boundary = TRUE)
## calculate corresponding test statistic
sctest(ocus.nile)

## UK Seatbelt data: a SARIMA(1,0,0)(1,0,0)_12 model
## (fitted by OLS) is used and reveals (at least) two
## breakpoints - one in 1973 associated with the oil crisis and
## one in 1983 due to the introduction of compulsory
## wearing of seatbelts in the UK.
data("UKDriverDeaths")
seatbelt <- log10(UKDriverDeaths)
seatbelt <- cbind(seatbelt, lag(seatbelt, k = -1), lag(seatbelt, k = -12))
colnames(seatbelt) <- c("y", "ylag1", "ylag12")
seatbelt <- window(seatbelt, start = c(1970, 1), end = c(1984,12))
plot(seatbelt[,"y"], ylab = expression(log[10](casualties)))

## use RE process
re.seat <- efp(y ~ ylag1 + ylag12, data = seatbelt, type = "RE")
plot(re.seat)
plot(re.seat, functional = NULL)
sctest(re.seat)
}

\concept{CUSUM}
\concept{MOSUM}
\concept{recursive estimates}
\concept{moving estimates}
\concept{fluctuation test}
\concept{maximum likelihood scores}
\concept{structural change}
\concept{CUSUM}
\concept{MOSUM}
\concept{recursive estimates}
\concept{moving estimates}
\concept{fluctuation test}
\concept{maximum likelihood scores}
\concept{structural change}
\keyword{regression}
