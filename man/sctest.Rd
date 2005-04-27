\name{sctest.formula}
\alias{sctest}
\alias{sctest.formula}
\title{Structural Change Tests}
\description{Performs tests for structural change.}
\usage{
\method{sctest}{formula}(formula, type = , h = 0.15,
    alt.boundary = FALSE, functional = c("max", "range",
    "maxL2", "meanL2"), from = 0.15, to = NULL, point = 0.5,
    asymptotic = FALSE, data, ...)
}

\arguments{
  \item{formula}{a formula describing the model to be tested.}
  \item{type}{a character string specifying the structural change test that ist
     to be performed. Besides the tests types described in \code{\link{efp}}
     and \code{\link{sctest.Fstats}}. The Chow test and the Nyblom-Hansen test
     can be performed by setting type to \code{"Chow"} or \code{"Nyblom-Hansen"},
     respectively.}
  \item{h}{numeric from interval (0,1) specifying the bandwidth. Determins the
     size of the data window relative to sample size (for MOSUM and ME tests
     only).}
  \item{alt.boundary}{logical. If set to \code{TRUE} alternative boundaries
     (instead of the standard linear boundaries) will be used (for CUSUM
     processes only).}
  \item{functional}{indicates which functional should be used to aggregate
     the empirical fluctuation processes to a test statistic.}
 \item{from, to}{numerics. If \code{from} is smaller than 1 they are
     interpreted as percentages of data and by default \code{to} is taken to be
     the 1 - \code{from}. F statistics will be calculated for the observations
     \code{(n*from):(n*to)}, when \code{n} is the number of observations in the
     model. If \code{from} is greater than 1 it is interpreted to be the index
     and \code{to} defaults to \code{n - from}. (for F tests only)}
  \item{point}{parameter of the Chow test for the potential change point.
    Interpreted analogous to the \code{from} parameter. By
     default taken to be \code{floor(n*0.5)} if \code{n} is the  number of
     observations in the model.}
  \item{asymptotic}{logical. If \code{TRUE} the asymptotic (chi-square)
     distribution instead of the exact (F) distribution will be used to compute
     the p value (for Chow test only).}
  \item{data}{an optional data frame containing the variables in the model. By
     default the variables are taken from the environment which
     \code{sctest} is called from.}
  \item{...}{further arguments passed to \code{\link{efp}} or
     \code{\link{Fstats}}.} }

\details{\code{sctest.formula} is mainly a wrapper for \code{\link{sctest.efp}}
and \code{\link{sctest.Fstats}} as it fits an empirical fluctuation process
first or computes the F statistics respectively and subsequently performs the
corresponding test. The Chow test and the Nyblom-Hansen test are available explicitely here.}

\value{an object of class \code{"htest"} containing:
\itemize
\item{statistic}{the test statistic}
\item{p.value}{the corresponding p value}
\item{method}{a character string with the method used}
\item{data.name}{a character string with the data name}}

\seealso{\code{\link{sctest.efp}}, \code{\link{sctest.Fstats}}}

\examples{
## Example 7.4 from Greene (1993), "Econometric Analysis"
## Chow test on Longley data
data(longley)
sctest(Employed ~ Year + GNP.deflator + GNP + Armed.Forces, data = longley,
  type = "Chow", point = 7)

## which is equivalent to segmenting the regression via
fac <- factor(c(rep(1, 7), rep(2, 9)))
fm0 <- lm(Employed ~ Year + GNP.deflator + GNP + Armed.Forces, data = longley)
fm1 <- lm(Employed ~ fac/(Year + GNP.deflator + GNP + Armed.Forces), data = longley)
anova(fm0, fm1)

## estimates from Table 7.5 in Greene (1993)
summary(fm0)
summary(fm1)
}
\keyword{htest}
