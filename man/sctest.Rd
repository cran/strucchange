\name{sctest.formula}
\alias{sctest}
\alias{sctest.formula}
\title{Structural Change Tests}
\description{Performs tests on structural change.}
\usage{
\method{sctest}{formula}(formula, type = <<see below>>, h = 0.15, dynamic = FALSE,
    tol = 1e-7, alt.boundary = FALSE, functional = c("max", "range"),
    from = 0.15, to = NULL, point = floor(0.5*nrow(model.frame(x))),
    asymptotic = FALSE, data, ...)
}

\arguments{
  \item{formula}{a formula describing the model to be tested.}
  \item{type}{a character string specifying the structural change test that ist
     to be performed. Besides the tests types described in \code{\link{efp}}
     and \code{\link{sctest.Fstats}} the Chow test is can be performed by
     setting type to \code{"Chow"}.}
  \item{h}{numeric from interval (0,1) specifying the bandwidth. Determins the
     size of the data window relative to sample size (for MOSUM and ME tests
     only).}
  \item{dynamic}{logical. If \code{TRUE} the lagged observations are included as
     a regressor (for generalized fluctuation tests only).}
  \item{tol}{tolerance when \code{solve} is used}
  \item{alt.boundary}{logical. If set to \code{TRUE} alternative boundaries
     (instead of the standard linear boundaries) will be used (for CUSUM
     processes only).}
  \item{functional}{indicates which functional should be applied to the
     estimates based processes (\code{"fluctuation"} and \code{"ME"}).}
 \item{from, to}{numerics. If \code{from} is smaller than 1 they are
     interpreted as percentages of data and by default \code{to} is taken to be
     the 1 - \code{from}. F statistics will be calculated for the observations
     \code{(n*from):(n*to)}, when \code{n} is the number of observations in the
     model. If \code{from} is greater than 1 it is interpreted to be the index
     and \code{to} defaults to \code{n - from}. (for F tests only)}
  \item{point}{parameter of the Chow test for the potential change point. By
     default taken to be \code{floor(n*0.5)} if \code{n} is the  number of
     observations in the model.}
  \item{asymptotic}{logical. If \code{TRUE} the asymptotic (chi-square)
     distribution instead of the exact (F) distribution will be used to compute
     the p value (for Chow test only).}
  \item{data}{an optional data frame containing the variables in the model. By
     default the variables are taken from the environment which
     \code{sctest} is called from.}
  \item{...}{currently not used.}
}

\details{\code{sctest.formula} is mainly a wrapper for \code{\link{sctest.efp}}
and \code{\link{sctest.Fstats}} as it fits an empirical fluctuation process
first or computes the F statistics respectively and subsequently performs the
corresponding test. The Chow test is available explicitely here.}

\value{an object of class \code{"htest"} containing:
\itemize
\item{statistic}{the test statistic}
\item{p.value}{the corresponding p value}
\item{method}{a character string with the method used}
\item{data.name}{a character string with the data name}}

\author{Achim Zeileis \email{zeileis@ci.tuwien.ac.at}}

\seealso{\code{\link{sctest.efp}}, \code{\link{sctest.Fstats}}}

\examples{
## Load dataset "nhtemp" with average yearly temperatures in New Haven
data(nhtemp)
## plot the data
plot(nhtemp)

## test the model null hypothesis that the average temperature remains constant
## over the years with the Standard CUSUM test
sctest(nhtemp ~ 1)
## with the Chow test (under the alternative that there is a change after the
## 30th observation)
sctest(nhtemp ~ 1, type = "Chow", point = 30)
}
\keyword{htest}
