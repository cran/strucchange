\name{sctest.Fstats}
\alias{sctest.Fstats}
\title{supF-, aveF- and expF-Test}
\description{Performs the supF-, aveF- or expF-test}
\usage{
\method{sctest}{Fstats}(x, type = c("supF", "aveF", "expF"),
    asymptotic = FALSE, ...)
}
\arguments{
 \item{x}{an object of class \code{"Fstats"}.}
 \item{type}{a character string specifying which test will be performed.}
 \item{asymptotic}{logical. Only necessary if \code{x} contains just a single F
     statistic and type is \code{"supF"} or \code{"aveF"}. If then set to
     \code{TRUE} the asymptotic (chi-square) distribution instead of the exact
     (F) distribution will be used to compute the p value.}
 \item{...}{currently not used.}
}

\details{If \code{x} contains just a single F statistic and type is
\code{"supF"} or \code{"aveF"} the Chow test will be performed.

The original GAUSS code for computing the p values of the supF-, aveF- and
expF-test was written by Bruce Hansen and is available from
\url{https://www.ssc.wisc.edu/~bhansen/}. R port by Achim Zeileis.
}

\value{
  An object of class \code{"htest"} containing:
  \item{statistic}{the test statistic,}
  \item{p.value}{the corresponding p value,}
  \item{method}{a character string with the method used,}
  \item{data.name}{a character string with the data name.}
}

\references{
Andrews D.W.K. (1993), Tests for parameter instability and structural
change with unknown change point, \emph{Econometrica}, \bold{61}, 821-856.

Andrews D.W.K., Ploberger W. (1994), Optimal tests when a nuisance parameter
is present only under the alternative, \emph{Econometrica}, \bold{62}, 1383-1414.

Hansen B. (1992), Tests for parameter instability in regressions with I(1)
processes, \emph{Journal of Business & Economic Statistics}, \bold{10}, 321-335.

Hansen B. (1997), Approximate asymptotic p values for structural-change
tests, \emph{Journal of Business & Economic Statistics}, \bold{15}, 60-67. }

\seealso{\code{\link{Fstats}}, \code{\link{plot.Fstats}}}

\examples{
## Load dataset "nhtemp" with average yearly temperatures in New Haven
data(nhtemp)
## plot the data
plot(nhtemp)

## test the model null hypothesis that the average temperature remains
## constant over the years for potential break points between 1941
## (corresponds to from = 0.5) and 1962 (corresponds to to = 0.85)
## compute F statistics
fs <- Fstats(nhtemp ~ 1, from = 0.5, to = 0.85)
## plot the F statistics
plot(fs, alpha = 0.01)
## and the corresponding p values
plot(fs, pval = TRUE, alpha = 0.01)
## perform the aveF test
sctest(fs, type = "aveF")
}
\keyword{htest}
