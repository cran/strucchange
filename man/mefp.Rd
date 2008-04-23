\name{mefp}
\alias{mefp}
\alias{mefp.formula}
\alias{mefp.efp}
\alias{print.mefp}
\alias{monitor}
\title{Monitoring of Empirical Fluctuation Processes}
\description{
  Online monitoring of structural breaks in a linear regression model. A
  sequential fluctuation test based on parameter estimates or OLS residualas
  signals structural breaks.
}
\usage{
mefp(obj, ...)

\method{mefp}{formula}(formula, type = c("OLS-CUSUM", "OLS-MOSUM", "RE", "ME",
    "fluctuation"), data, h = 1, alpha = 0.05,
    functional = c("max", "range"), period = 10,
    tolerance = .Machine$double.eps^0.5, CritvalTable = NULL,
    rescale = NULL, border = NULL, ...)

\method{mefp}{efp}(obj, alpha=0.05, functional = c("max", "range"),
    period = 10, tolerance = .Machine$double.eps^0.5,
    CritvalTable = NULL, rescale = NULL, border = NULL, ...)

monitor(obj, data = NULL, verbose = TRUE)
}
\arguments{
  \item{formula}{a symbolic description for the model to be tested.}
  \item{data}{an optional data frame containing the variables in the model. By
     default the variables are taken from the environment which \code{efp} is
     called from.}
  \item{type}{specifies which type of fluctuation process will be
     computed.}
  \item{h}{(only used for MOSUM/ME processes). A numeric scalar from interval
    (0,1) specifying the size of the data window relative to the sample
    size.}
  \item{obj}{Object of class \code{"efp"} (for \code{mefp}) or
    \code{"mefp"} (for \code{monitor}).}
  \item{alpha}{Significance level of the test, i.e., probability of
    type I error.}
  \item{functional}{Determines if maximum or range of parameter
    differences is used as statistic.}
  \item{period}{(only used for MOSUM/ME processes). Maximum time (relative
    to the history period) that will be monitored. Default is 10 times
    the history period.}
  \item{tolerance}{Tolerance for numeric \code{==} comparisons.}
  \item{CritvalTable}{Table of critical values, this table is
    interpolated to get critical values
    for arbitrary \code{alpha}s. The default depends on the \code{type}
    of fluctuation process (pre-computed tables are available for all
    types). \emph{This argument is under development.}}
  \item{rescale}{If \code{TRUE} the estimates will be standardized by
    the regressor matrix of the corresponding subsample similar to
    Kuan & Chen (1994); if \code{FALSE} the historic regressor matrix will
    be used. The default is to rescale the monitoring processes of type
    \code{"ME"} but not of \code{"RE"}.}
  \item{border}{An optional user-specified border function for the
      empirical process. \emph{This argument is under development.}}
  \item{verbose}{If \code{TRUE}, signal breaks by text output.}
  \item{...}{Currently not used.}
}
\details{
  \code{\link{mefp}} creates an object of class \code{"mefp"} either
  from a model formula or from an object of class \code{"efp"}. In
  addition to the arguments of \code{\link{efp}}, the type of statistic
  and a significance level for the monitoring must be specified. The
  monitoring itself is performed by \code{monitor}, which can be
  called arbitrarily often on objects of class \code{"mefp"}. If new
  data have arrived, then the empirical fluctuation process is computed
  for the new data. If the process crosses the boundaries corresponding
  to the significance level \code{alpha}, a structural break is detected
  (and signaled).

  The typical usage is to initialize the monitoring by creation of an
  object of class \code{"mefp"} either using a formula or an
  \code{"efp"} object. Data available at this stage are considered the
  \emph{history sample}, which is kept fixed during the complete
  monitoring process, and may not contain any structural changes.

  Subsequent calls to \code{monitor} perform a sequential test of the
  null hypothesis of no structural change in new data against the
  general alternative of changes in one or more of the coefficients of
  the regression model.

  The recursive
  estimates test is also called fluctuation test, therefore setting \code{type}
  to \code{"fluctuation"} was used to specify it in earlier versions of
  strucchange. It still can be used now, but will be forced to \code{"RE"}
}
\seealso{\code{\link{plot.mefp}}, \code{\link{boundary.mefp}}}
\references{
  Leisch F., Hornik K., Kuan C.-M. (2000), Monitoring
  Structural Changes with the Generalized Fluctuation Test,
  \emph{Econometric Theory}, \bold{16}, 835--854.

  Zeileis A., Leisch F., Kleiber C., Hornik K. (2005), Monitoring
  Structural Change in Dynamic Econometric Models,
  \emph{Journal of Applied Econometrics}, \bold{20}, 99--121.
  doi:10.1002/jae.776.
}
\examples{
df1 <- data.frame(y=rnorm(300))
df1[150:300,"y"] <- df1[150:300,"y"]+1

## use the first 50 observations as history period
e1 <- efp(y~1, data=df1[1:50,,drop=FALSE], type="ME", h=1)
me1 <- mefp(e1, alpha=0.05)

## the same in one function call
me1 <- mefp(y~1, data=df1[1:50,,drop=FALSE], type="ME", h=1,
              alpha=0.05)

## monitor the 50 next observations
me2 <- monitor(me1, data=df1[1:100,,drop=FALSE])
plot(me2)

# and now monitor on all data
me3 <- monitor(me2, data=df1)
plot(me3)


## Load dataset "USIncExp" with income and expenditure in the US
## and choose a suitable subset for the history period
data("USIncExp")
USIncExp3 <- window(USIncExp, start=c(1969,1), end=c(1971,12))
## initialize the monitoring with the formula interface
me.mefp <- mefp(expenditure~income, type="ME", rescale=TRUE,
                   data=USIncExp3, alpha=0.05)

## monitor the new observations for the year 1972
USIncExp3 <- window(USIncExp, start=c(1969,1), end=c(1972,12))
me.mefp <- monitor(me.mefp)

## monitor the new data for the years 1973-1976
USIncExp3 <- window(USIncExp, start=c(1969,1), end=c(1976,12))
me.mefp <- monitor(me.mefp)
plot(me.mefp, functional = NULL)
}

\concept{CUSUM}
\concept{MOSUM}
\concept{recursive estimates}
\concept{moving estimates}
\concept{fluctuation test}
\concept{monitoring}
\concept{structural change}
\keyword{regression}
