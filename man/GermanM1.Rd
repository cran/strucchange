\name{GermanM1}
\alias{GermanM1}
\alias{historyM1}
\alias{monitorM1}
\title{German M1 Money Demand}
\usage{data(GermanM1)}
\description{
German M1 money demand.
}
\format{
\code{GermanM1} is a data frame containing 12 quarterly time series
from 1961(1) to 1995(4) and two further variables. \code{historyM1}
is the subset of \code{GermanM1} up to 1990(2), i.e., the data before
the German monetary unification on 1990-06-01. \code{monitorM1}
is the complement of \code{historyM1}, i.e., the data after
the unification. All three data frames contain the variables
 \describe{
    \item{m}{time series. Logarithm of real M1 per capita,}
    \item{p}{time series. Logarithm of a price index,}
    \item{y}{time series. Logarithm of real per capita gross
             national product,}
    \item{R}{time series. Long-run interest rate,}
    \item{dm}{time series. First differences of \code{m},}
    \item{dy2}{time series. First differences of lag 2 of \code{y},}
    \item{dR}{time series. First differences of \code{R},}
    \item{dR1}{time series. First differences of lag 1 of \code{R},}
    \item{dp}{time series. First differences of \code{p},}
    \item{m1}{time series. Lag 1 of \code{m},}
    \item{y1}{time series. Lag 1 of \code{y},}
    \item{R1}{time series. Lag 1 of \code{R},}
    \item{season}{factor coding the seasonality,}
    \item{ecm.res}{vector containing the OLS residuals of
                   the Lütkepohl et al. (1999) model fitted in the history
		   period.}
}
}

\details{
Lütkepohl et al. (1999) investigate the linearity and
stability of German M1 money demand: they find a stable regression relation
for the time before the monetary union on 1990-06-01 but a clear structural
instability afterwards.

Zeileis et al. (2002) use a model with
\code{ecm.res} instead of \code{m1}, \code{y1} and \code{R1}, which
leads to equivalent results in the history period but slightly
different results in the monitoring period. The reason for the
replacement is that stationary regressors are needed for the
structural change tests. See references and the examples below for
more details.
}

\source{The data is provided by the German central bank and is
available online in the data archive of the Journal of Applied
Econometrics
\url{http://qed.econ.queensu.ca/jae/1999-v14.5/lutkepohl-terasvirta-wolters/}.}

\references{
Lütkepohl H., Teräsvirta T., Wolters J. (1999), Investigating
Stability and Linearity of a German M1 Money Demand Function,
\emph{Journal of Applied Econometrics}, \bold{14}, 511-525.

Zeileis A., Leisch F., Kleiber C., Hornik K. (2002), Monitoring
Structural Change in Dynamic Econometric Models,
Report 64, SFB "Adaptive Information Systems and Modelling in Economics
and Management Science", Vienna University of Economics,
\url{http://www.wu-wien.ac.at/am/reports.htm#78}.}

\examples{
data(GermanM1)
## Lütkepohl et al. (1999) use the following model
LTW.model <- dm ~ dy2 + dR + dR1 + dp + m1 + y1 + R1 + season
## Zeileis et al. (2002) use
M1.model <- dm ~ dy2 + dR + dR1 + dp + ecm.res + season


## historical tests
ols <- efp(LTW.model, data = GermanM1, type = "OLS-CUSUM")
plot(ols)
re <- efp(LTW.model, data = GermanM1, type = "fluctuation")
plot(re)
fs <- Fstats(LTW.model, data = GermanM1, from = 0.1)
plot(fs)

## monitoring
M1 <- historyM1
ols.efp <- efp(M1.model, type = "OLS-CUSUM", data = M1)
newborder <- function(k) 1.5778*k/118
ols.mefp <- mefp(ols.efp, period = 2)
ols.mefp2 <- mefp(ols.efp, border = newborder)
M1 <- GermanM1
ols.mon <- monitor(ols.mefp)
ols.mon2 <- monitor(ols.mefp2)
plot(ols.mon)
lines(boundary(ols.mon2), col = 2)

## dating
bp <- breakpoints(LTW.model, data = GermanM1)
summary(bp)
plot(bp)

plot(fs)
lines(confint(bp))
}

\keyword{datasets}
