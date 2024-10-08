\name{durab}
\alias{durab}
\title{US Labor Productivity}
\usage{data("durab")}
\description{
US labor productivity in the manufacturing/durables sector.
}
\format{
\code{durab} is a multivariate monthly time series from 1947(3)
to 2001(4) with variables
 \describe{
    \item{y}{growth rate of the Industrial Production Index to
             average weekly labor hours in the manufacturing/durables sector,}
    \item{lag}{lag 1 of the series \code{y},}
}
}

\source{The data set is available from Bruce Hansen's homepage
\url{https://www.ssc.wisc.edu/~bhansen/}. For more information
see Hansen (2001).}

\references{
Hansen B. (2001), The New Econometrics of Structural Change:
Dating Breaks in U.S. Labor Productivity,
\emph{Journal of Economic Perspectives}, \bold{15}, 117--128.

Zeileis A., Leisch F., Kleiber C., Hornik K. (2005), Monitoring
Structural Change in Dynamic Econometric Models,
\emph{Journal of Applied Econometrics}, \bold{20}, 99--121.
}

\examples{
data("durab")
## use AR(1) model as in Hansen (2001) and Zeileis et al. (2005)
durab.model <- y ~ lag

## historical tests
## OLS-based CUSUM process
ols <- efp(durab.model, data = durab, type = "OLS-CUSUM")
plot(ols)
## F statistics
fs <- Fstats(durab.model, data = durab, from = 0.1)
plot(fs)

\donttest{
## F statistics based on heteroskadisticy-consistent covariance matrix
fsHC <- Fstats(durab.model, data = durab, from = 0.1,
	       vcov = function(x, ...) vcovHC(x, type = "HC", ...))
plot(fsHC)
}

## monitoring
Durab <- window(durab, start=1964, end = c(1979, 12))
ols.efp <- efp(durab.model, type = "OLS-CUSUM", data = Durab)
newborder <- function(k) 1.723 * k/192
ols.mefp <- mefp(ols.efp, period=2)
ols.mefp2 <- mefp(ols.efp, border=newborder)
Durab <- window(durab, start=1964)
ols.mon <- monitor(ols.mefp)
ols.mon2 <- monitor(ols.mefp2)
plot(ols.mon)
lines(boundary(ols.mon2), col = 2)
## Note: critical value for linear boundary taken from Table III
## in Zeileis et al. 2005: (1.568 + 1.896)/2 = 1.732 is a linear
## interpolation between the values for T = 2 and T = 3 at
## alpha = 0.05. A typo switched 1.732 to 1.723.

\donttest{
## dating
bp <- breakpoints(durab.model, data = durab)
summary(bp)
plot(summary(bp))

plot(ols)
lines(breakpoints(bp, breaks = 1), col = 3)
lines(breakpoints(bp, breaks = 2), col = 4)
plot(fs)
lines(breakpoints(bp, breaks = 1), col = 3)
lines(breakpoints(bp, breaks = 2), col = 4)
}
}
\keyword{datasets}
