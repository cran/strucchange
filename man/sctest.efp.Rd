\name{sctest.efp}
\alias{sctest.efp}
\title{Generalized Fluctuation Tests}
\description{Performs a generalized fluctuation test.}
\usage{
\method{sctest}{efp}(x, alt.boundary = FALSE,
    functional = c("max", "range", "maxL2", "meanL2"), ...)
}

\arguments{
  \item{x}{an object of class \code{"efp"}.}
  \item{alt.boundary}{logical. If set to \code{TRUE} alternative boundaries
 (instead of
     the standard linear boundaries) will be used (for CUSUM
     processes only).}
  \item{functional}{indicates which functional should be applied to the
     empirical fluctuation process.}
  \item{...}{currently not used.}
}

\details{The critical values for the MOSUM tests and the ME test are just
tabulated for confidence levels between 0.1 and 0.01, thus the p
value approximations will be poor for other p values. Similarly the
critical values for the maximum and mean squared Euclidian norm (\code{"maxL2"}
and \code{"meanL2"}) are tabulated for confidence levels between 0.2 and 0.005.}

\value{an object of class \code{"htest"} containing:
\itemize
\item{statistic}{the test statistic}
\item{p.value}{the corresponding p value}
\item{method}{a character string with the method used}
\item{data.name}{a character string with the data name}}

\references{Brown R.L., Durbin J., Evans J.M. (1975), Techniques for
testing constancy of regression relationships over time, \emph{Journal of the
Royal Statistal Society}, B, \bold{37}, 149-163.

Chu C.-S., Hornik K., Kuan C.-M. (1995), MOSUM tests for parameter
constancy, \emph{Biometrika}, \bold{82}, 603-617.

Chu C.-S., Hornik K., Kuan C.-M. (1995), The moving-estimates test for
parameter stability, \emph{Econometric Theory}, \bold{11}, 669-720.

Krämer W., Ploberger W., Alt R. (1988), Testing for structural change in
dynamic models, \emph{Econometrica}, \bold{56}, 1355-1369.

Kuan C.-M., Hornik K. (1995), The generalized fluctuation test: A
unifying view, \emph{Econometric Reviews}, \bold{14}, 135 - 161.

Kuan C.-M., Chen (1994), Implementing the fluctuation and moving estimates
tests in dynamic econometric models, \emph{Economics Letters}, \bold{44},
235-239.

Ploberger W., Krämer W. (1992), The CUSUM Test with OLS Residuals,
\emph{Econometrica}, \bold{60}, 271-285.

Zeileis A. (2000), p Values and Alternative Boundaries for CUSUM Tests,
Working Paper 78, SFB "Adaptive Information Systems and Modelling in Economics
and Management Science", Vienna University of Economics,
\url{http://www.wu-wien.ac.at/am/wp00.htm#78}.
}

\seealso{\code{\link{efp}}, \code{\link{plot.efp}}}

\examples{
## Load dataset "nhtemp" with average yearly temperatures in New Haven
data(nhtemp)
## plot the data
plot(nhtemp)

## test the model null hypothesis that the average temperature remains
## constant over the years compute OLS-CUSUM fluctuation process
temp.cus <- efp(nhtemp ~ 1, type = "OLS-CUSUM")
## plot the process with alternative boundaries
plot(temp.cus, alpha = 0.01, alt.boundary = TRUE)
## and calculate the test statistic
sctest(temp.cus)

## compute moving estimates fluctuation process
temp.me <- efp(nhtemp ~ 1, type = "ME", h = 0.2)
## plot the process with functional = "max"
plot(temp.me)
## and perform the corresponding test
sctest(temp.me)
}
\keyword{htest}
