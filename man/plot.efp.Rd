\name{plot.efp}
\alias{plot.efp}
\alias{lines.efp}
\title{Plot Empirical Fluctuation Process}
\description{Plot and lines method for objects of class \code{"efp"}}
\usage{
\method{plot}{efp}(x, alpha = 0.05, alt.boundary = FALSE, boundary = TRUE,
    functional = "max", main = NULL,  ylim = NULL,
    ylab = "empirical fluctuation process", ...)
\method{lines}{efp}(x, functional = "max", ...)
}

\arguments{
  \item{x}{an object of class \code{"efp"}.}
  \item{alpha}{numeric from interval (0,1) indicating the confidence level for
     which the boundary of the corresponding test will be computed.}
  \item{alt.boundary}{logical. If set to \code{TRUE} alternative boundaries
     (instead of the standard linear boundaries) will be plotted (for CUSUM
     processes only).}
  \item{boundary}{logical. If set to \code{FALSE} the boundary will be computed
     but not plotted.}
  \item{functional}{indicates which functional should be applied to the
     process before plotting and which boundaries should be used. If set to \code{NULL}
     a multiple process with boundaries for the \code{"max"} functional is plotted.
     For more details see below.}
  \item{main, ylim, ylab, ...}{high-level \code{\link{plot}} function
 parameters.}
}

\details{Plots are available for the \code{"max"} functional for all process types.
For Brownian bridge type processes the maximum or mean squared Euclidian norm
(\code{"maxL2"} and \code{"meanL2"}) can be used for aggregating before plotting.
No plots are available for the \code{"range"} functional.

Alternative boundaries that are proportional to the standard deviation
of the corresponding limiting process are available for processes with Brownian
motion or Brownian bridge limiting processes.
}

\value{\code{\link{efp}} returns an object of class \code{"efp"} which inherits
from the class \code{"ts"} or \code{"mts"} respectively. The function
\code{\link{plot}} has a method to plot the
empirical fluctuation process; with \code{sctest} the corresponding test for
structural change can be performed.}

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

Ploberger W., Krämer W. (1992), The CUSUM test with OLS residuals,
\emph{Econometrica}, \bold{60}, 271-285.

Zeileis A. (2000), p Values and Alternative Boundaries for CUSUM Tests,
Working Paper 78, SFB "Adaptive Information Systems and Modelling in Economics
and Management Science", Vienna University of Economics,
\url{http://www.wu-wien.ac.at/am/wp00.htm#78}.
}

\seealso{\code{\link{efp}}, \code{\link{boundary.efp}},
\code{\link{sctest.efp}}}

\examples{
## Load dataset "nhtemp" with average yearly temperatures in New Haven
data(nhtemp)
## plot the data
plot(nhtemp)

## test the model null hypothesis that the average temperature remains
## constant over the years
## compute Rec-CUSUM fluctuation process
temp.cus <- efp(nhtemp ~ 1)
## plot the process
plot(temp.cus, alpha = 0.01)
## and calculate the test statistic
sctest(temp.cus)

## compute (recursive estimates) fluctuation process
## with an additional linear trend regressor
lin.trend <- 1:60
temp.me <- efp(nhtemp ~ lin.trend, type = "fluctuation")
## plot the bivariate process
plot(temp.me, functional = NULL)
## and perform the corresponding test
sctest(temp.me)
}
\keyword{hplot}
