\name{efp}
\alias{efp}
\alias{print.efp}
\title{Empirical Fluctuation Process}
\description{Computes an empirical fluctuation process according
to a specified method from the generalized fluctuation test
framework}
\usage{
efp(formula, data, type = <<see below>>, h = 0.15, dynamic = FALSE,
  rescale = TRUE, tol = 1e-7)}
\arguments{
  \item{formula}{a symbolic describtion for the model to be tested.}
  \item{data}{an optional data frame containing the variables in the model. By
     default the variables are taken from the environment which \code{efp} is
     called from.}
  \item{type}{specifies which type of fluctuation process will be
     computed. For details see below.}
  \item{h}{a numeric from interval (0,1) sepcifying the bandwidth. determins the
     size of the data window relative to sample size (for MOSUM and ME processes
     only).}
  \item{dynamic}{logical. If \code{TRUE} the lagged observations are included as
     a regressor.}
  \item{rescale}{logical. If \code{TRUE} the estimates will be standardized by
    the regressor matrix of the corresponding subsample according to Kuan & Chen
    (1994); if \code{FALSE} the whole regressor matrix will be used.
    (only if \code{type} is either \code{"fluctuation"} or \code{"ME"})}
  \item{tol}{tolerance when \code{solve} is used}
}

\details{If \code{type} is one of \code{"Rec-CUSUM"}, \code{"OLS-CUSUM"},
\code{"Rec-MOSUM"} or \code{"OLS-MOSUM"} the function \code{efp} will return a
one-dimensional empiricial process of sums of residuals. Either it will be based
on recursive residuals or on OLS residuals and the process will contain
CUmulative SUMs or MOving SUMs of residuals in a certain data window.
For the MOSUM and ME processes all estimations are done for the
observations in a moving data window, whose size is determined by \code{h} and
which is shifted over the whole sample.

If there is a single structural change point \eqn{t^*}, the standard CUSUM path starts
to depart from its mean 0 at \eqn{t^*}. The OLS-based CUSUM path will have its peak
around \eqn{t^*}. The MOSUM path should have a strong change at \eqn{t^*}.

If \code{type} is either \code{"fluctuation"} or \code{"ME"} a
\emph{k}-dimensional process will be returned, if \emph{k} is the number of
regressors in the model, as it is based on recursive OLS estimates of the
regression coefficients or moving OLS estimates respectively.

Both paths should have a peak around \eqn{t^*} if there is a single structural shift.
}

\value{
  \code{efp} returns an object of class \code{"efp"} which inherits from
  the class \code{"ts"} or \code{"mts"} respectively, to which a string
  with the \code{type} of the process, the number of regressors
  (\code{nreg}), the bandwidth (\code{h}) and the function call
  (\code{call}) are added as attributes. The function \code{\link{plot}}
  has a method to plot the empirical fluctuation process; with
  \code{sctest} the corresponding test on structural change can be
  performed.
}

\references{Brown R.L., Durbin J., Evans J.M. (1975), Techniques for
testing constancy of regression relationships over time, \emph{Journal of the
Royal Statistal Society}, B, \bold{37}, 149-163.

Chu C.-S., Hornik K., Kuan C.-M. (1995), MOSUM tests for parameter
constancy, \emph{Biometrika}, \bold{82}, 603-617.

Chu C.-S., Hornik K., Kuan C.-M. (1995), The moving-estimates test for
parameter stability, \emph{Econometric Theory}, \bold{11}, 669-720.

Kraemer W., Ploberger W., Alt R. (1988), Testing for structural change in
dynamic models, \emph{Econometrica}, \bold{56}, 1355-1369.

Kuan C.-M., Hornik K. (1995), The generalized fluctuation test: A
unifying view, \emph{Econometric Reviews}, \bold{14}, 135 - 161.

Kuan C.-M., Chen (1994), Implementing the fluctuation and moving estimates
tests in dynamic econometric models, \emph{Economics Letters}, \bold{44},
235-239.

Ploberger W., Kraemer W. (1992), The CUSUM test with OLS residuals,
\emph{Econometrica}, \bold{60}, 271-285. }

\author{Achim Zeileis \email{zeileis@ci.tuwien.ac.at}}

\seealso{\code{\link{plot.efp}}, \code{\link{print.efp}},
\code{\link{sctest.efp}}, \code{\link{boundary.efp}}}

\examples{
## Load dataset "nhtemp" with average yearly temperatures in New Haven
data(nhtemp)
## plot the data
plot(nhtemp)

## test the model null hypothesis that the average temperature remains constant
## over the years
## compute OLS-CUSUM fluctuation process
temp.cus <- efp(nhtemp ~ 1, type = "OLS-CUSUM")
## plot the process with alternative boundaries
plot(temp.cus, alpha = 0.01, alt.boundary = TRUE)
## and calculate the test statistic
sctest(temp.cus)


## Load dataset "USIncExp" with income and expenditure in the US
## and choose a suitable subset
data(USIncExp)
USIncExp2 <- window(USIncExp, start=c(1970,1), end=c(1989,12))

## test the null hypothesis that the way the income is spent in expenditure
## does not change over time
## compute moving estimates fluctuation process
me <- efp(expenditure~income, type="ME", data=USIncExp2, h=0.2)
## plot the two dimensional fluctuation process with boundaries
plot(me, functional=NULL)
## and perform the corresponding test
sctest(me)
}
\keyword{regression}
