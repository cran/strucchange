\name{PhillipsCurve}
\alias{PhillipsCurve}
\title{UK Phillips Curve Equation Data}
\usage{data("PhillipsCurve")}
\description{
Macroeconomic time series from the United Kingdom with variables
for estimating the Phillips curve equation.
}
\format{
A multivariate annual time series from 1857 to 1987 with the columns
 \describe{
    \item{p}{Logarithm of the consumer price index,}
    \item{w}{Logarithm of nominal wages,}
    \item{u}{Unemployment rate,}
    \item{dp}{First differences of \code{p},}
    \item{dw}{First differences of \code{w},}
    \item{du}{First differences of \code{u}}
    \item{u1}{Lag 1 of \code{u},}
    \item{dp1}{Lag 1 of \code{dp}.}
}
}

\source{The data is available online in the data archive of the
Journal of Applied Econometrics
\url{http://qed.econ.queensu.ca/jae/2003-v18.1/bai-perron/}.}

\references{
Alogoskoufis G.S., Smith R. (1991), The Phillips Curve, the Persistence of
  Inflation, and the Lucas Critique: Evidence from Exchange Rate Regimes,
  \emph{American Economic Review}, \bold{81}, 1254-1275.

Bai J., Perron P. (2003), Computation and Analysis of Multiple Structural Change
  Models, \emph{Journal of Applied Econometrics}, \bold{18}, 1-22.
}

\examples{
## load and plot data
data("PhillipsCurve")
uk <- window(PhillipsCurve, start = 1948)
plot(uk[, "dp"])

## AR(1) inflation model
## estimate breakpoints
bp.inf <- breakpoints(dp ~ dp1, data = uk, h = 8)
plot(bp.inf)
summary(bp.inf)

## fit segmented model with three breaks
fac.inf <- breakfactor(bp.inf, breaks = 2, label = "seg")
fm.inf <- lm(dp ~ 0 + fac.inf/dp1, data = uk)
summary(fm.inf)

## Results from Table 2 in Bai & Perron (2003):
## coefficient estimates
coef(bp.inf, breaks = 2)
## corresponding standard errors
sqrt(sapply(vcov(bp.inf, breaks = 2), diag))
## breakpoints and confidence intervals
confint(bp.inf, breaks = 2)

## Phillips curve equation
## estimate breakpoints
bp.pc <- breakpoints(dw ~ dp1 + du + u1, data = uk, h = 5, breaks = 5)
## look at RSS and BIC
plot(bp.pc)
summary(bp.pc)

## fit segmented model with three breaks
fac.pc <- breakfactor(bp.pc, breaks = 2, label = "seg")
fm.pc <- lm(dw ~ 0 + fac.pc/dp1 + du + u1, data = uk)
summary(fm.pc)

## Results from Table 3 in Bai & Perron (2003):
## coefficient estimates
coef(fm.pc)
## corresponding standard errors
sqrt(diag(vcov(fm.pc)))
## breakpoints and confidence intervals
confint(bp.pc, breaks = 2, het.err = FALSE)
}

\keyword{datasets}
