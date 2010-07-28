\name{SP2001}
\alias{SP2001}
\title{S\&P 500 Stock Prices}
\description{
  A multivariate series of all S\&P 500 stock prices in the
  second half of the year 2001, i.e., before and after the
  terrorist attacks of 2001-09-11. 
}
\usage{data("SP2001")}
\format{
 A multivariate daily \code{"zoo"} series with \code{"Date"} index 
 from 2001-07-31 to 2001-12-31 (103 observations) of all 500 S\&P
 stock prices. 
}
\source{Yahoo! Finance: \url{http://finance.yahoo.com/}.}

\references{
Zeileis A., Leisch F., Kleiber C., Hornik K. (2005), Monitoring
Structural Change in Dynamic Econometric Models,
\emph{Journal of Applied Econometrics}, \bold{20}, 99--121.
}

\seealso{\code{\link[tseries]{get.hist.quote}}}

\examples{
## load and transform data
## (DAL: Delta Air Lines, LU: Lucent Technologies)
data("SP2001")
stock.prices <- SP2001[, c("DAL", "LU")]
stock.returns <- diff(log(stock.prices))

## price and return series
plot(stock.prices, ylab = c("Delta Air Lines", "Lucent Technologies"), main = "")
plot(stock.returns, ylab = c("Delta Air Lines", "Lucent Technologies"), main = "")

## monitoring of DAL series
myborder <- function(k) 1.939*k/28
x <- as.vector(stock.returns[, "DAL"][1:28])
dal.cusum <- mefp(x ~ 1, type = "OLS-CUSUM", border = myborder)
dal.mosum <- mefp(x ~ 1, type = "OLS-MOSUM", h = 0.5, period = 4)
x <- as.vector(stock.returns[, "DAL"])
dal.cusum <- monitor(dal.cusum)
dal.mosum <- monitor(dal.mosum)

## monitoring of LU series
x <- as.vector(stock.returns[, "LU"][1:28])
lu.cusum <- mefp(x ~ 1, type = "OLS-CUSUM", border = myborder)
lu.mosum <- mefp(x ~ 1, type = "OLS-MOSUM", h = 0.5, period = 4)
x <- as.vector(stock.returns[, "LU"])
lu.cusum <- monitor(lu.cusum)
lu.mosum <- monitor(lu.mosum)

## pretty plotting
## (needs some work because lm() does not keep "zoo" attributes)
cus.bound <- zoo(c(rep(NA, 27), myborder(28:102)), index(stock.returns))
mos.bound <- as.vector(boundary(dal.mosum))
mos.bound <- zoo(c(rep(NA, 27), mos.bound[1], mos.bound), index(stock.returns))

## Lucent Technologies: CUSUM test
plot(zoo(c(lu.cusum$efpprocess, lu.cusum$process), index(stock.prices)),
  ylim = c(-1, 1) * coredata(cus.bound)[102], xlab = "Time", ylab = "empirical fluctuation process")
abline(0, 0)
abline(v = as.Date("2001-09-10"), lty = 2)
lines(cus.bound, col = 2)
lines(-cus.bound, col = 2)

## Lucent Technologies: MOSUM test
plot(zoo(c(lu.mosum$efpprocess, lu.mosum$process), index(stock.prices)[-(1:14)]),
  ylim = c(-1, 1) * coredata(mos.bound)[102], xlab = "Time", ylab = "empirical fluctuation process")
abline(0, 0)
abline(v = as.Date("2001-09-10"), lty = 2)
lines(mos.bound, col = 2)
lines(-mos.bound, col = 2)

## Delta Air Lines: CUSUM test
plot(zoo(c(dal.cusum$efpprocess, dal.cusum$process), index(stock.prices)),
  ylim = c(-1, 1) * coredata(cus.bound)[102], xlab = "Time", ylab = "empirical fluctuation process")
abline(0, 0)
abline(v = as.Date("2001-09-10"), lty = 2)
lines(cus.bound, col = 2)
lines(-cus.bound, col = 2)

## Delta Air Lines: MOSUM test
plot(zoo(c(dal.mosum$efpprocess, dal.mosum$process), index(stock.prices)[-(1:14)]),
  ylim = range(dal.mosum$process), xlab = "Time", ylab = "empirical fluctuation process")
abline(0, 0)
abline(v = as.Date("2001-09-10"), lty = 2)
lines(mos.bound, col = 2)
lines(-mos.bound, col = 2)
}

\keyword{datasets}
