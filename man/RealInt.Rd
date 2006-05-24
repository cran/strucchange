\name{RealInt}
\alias{RealInt}
\title{US Ex-post Real Interest Rate}
\description{
US ex-post real interest rate: the three-month treasury bill
deflated by the CPI inflation rate.
}
\usage{data("RealInt")}
\format{
A quarterly time series from 1961(1) to 1986(3).
}

\source{The data is available online in the data archive of the
Journal of Applied Econometrics
\url{http://qed.econ.queensu.ca/jae/2003-v18.1/bai-perron/}.}

\references{
Bai J., Perron P. (2003), Computation and Analysis of Multiple Structural Change
  Models, \emph{Journal of Applied Econometrics}, \bold{18}, 1-22.
}

\examples{
## load and plot data
data("RealInt")
plot(RealInt)

## estimate breakpoints
bp.ri <- breakpoints(RealInt ~ 1, h = 15)
plot(bp.ri)
summary(bp.ri)

## fit segmented model with three breaks
fac.ri <- breakfactor(bp.ri, breaks = 3, label = "seg")
fm.ri <- lm(RealInt ~ 0 + fac.ri)
summary(fm.ri)

## setup kernel HAC estimator
vcov.ri <- function(x, ...) kernHAC(x, kernel = "Quadratic Spectral",
  prewhite = 1, approx = "AR(1)", ...)

## Results from Table 1 in Bai & Perron (2003):
## coefficient estimates
coef(bp.ri, breaks = 3)
## corresponding standard errors
sapply(vcov(bp.ri, breaks = 3, vcov = vcov.ri), sqrt)
## breakpoints and confidence intervals
confint(bp.ri, breaks = 3, vcov = vcov.ri)

## Visualization
plot(RealInt)
lines(as.vector(time(RealInt)), fitted(fm.ri), col = 4)
lines(confint(bp.ri, breaks = 3, vcov = vcov.ri))
}

\keyword{datasets}
