\name{covHC}
\alias{covHC}
\title{Heteroskedasticity-Consistent Covariance Matrix Estimation}
\description{
 Heteroskedasticity-consistent estimation of the covariance matrix of the
 coefficient estimates in a linear regression model.
}
\usage{
covHC(formula, type = c("HC2", "const", "HC", "HC1", "HC3"), tol = 1e-10,
   data=list())
}
\arguments{
 \item{formula}{a symbolic description for the model to be fitted.}
  \item{type}{a character string specifying the estimation type. For
     details see below.}
 \item{tol}{tolerance when \code{solve} is used}
 \item{data}{an optional data frame containing the variables in the model.
    By default the variables are taken from the environment which \code{covHC}
    is called from.}
}
\details{When \code{type = "const"} constant variances are assumed and
  and \code{covHC} gives the usual estimate of the covariance matrix of
  the coefficient estimates:

  \deqn{\hat \sigma^2 (X^\top X)^{-1}}{sigma^2 (X'X)^{-1}}

  All other methods do not assume constant variances and are suitable in case of
  heteroskedasticity. \code{"HC"} gives White's estimator; for details see the
  references.
}
\value{A matrix containing the covariance matrix estimate.
}
\references{

MacKinnon J. G., White H. (1985),
Some heteroskedasticity-consistent
covariance matrix estimators with improved finite sample properties.
\emph{Journal of Econometrics} \bold{29}, 305-325

}

\seealso{\code{\link{lm}}}

\examples{
## generate linear regression relationship
## with homoskedastic variances
x <- sin(1:100)
y <- 1 + x + rnorm(100)
## compute usual covariance matrix of coefficient estimates
covHC(y~x, type="const")

sigma2 <- sum(residuals(lm(y~x))^2)/98
sigma2 * solve(crossprod(cbind(1,x)))
}
\keyword{htest}
