\name{recresid}
\alias{recresid}
\alias{recresid.default}
\alias{recresid.formula}
\alias{recresid.lm}
\title{Recursive Residuals}
\description{
  A generic function for computing the recursive residuals
  (standardized one step prediction errors) of a linear regression model.
}
\usage{
\method{recresid}{default}(x, y, start = ncol(x) + 1, end = nrow(x),
  tol = sqrt(.Machine$double.eps)/ncol(x), \dots)
\method{recresid}{formula}(formula, data = list(), \dots)
\method{recresid}{lm}(x, data = list(), \dots)
}
\arguments{
  \item{x, y, formula}{specification of the linear regression model:
    either by a regressor matrix \code{x} and a response variable \code{y},
    or by a \code{formula} or by a fitted object \code{x} of class \code{"lm"}.}
  \item{start, end}{integer. Index of the first and last observation, respectively,
    for which recursive residuals should be computed. By default, the maximal
    range is selected.}
  \item{tol}{numeric. A relative tolerance for precision of recursive
    coefficient estimates, see details.}
  \item{data}{an optional data frame containing the variables in the model. By
     default the variables are taken from the environment which \code{recresid} is
     called from. Specifying \code{data} might also be necessary when applying
     \code{recresid} to a fitted model of class \code{"lm"} if this does not
     contain the regressor matrix and the response.}
  \item{\dots}{\emph{currently not used.}}
}
\details{
  Recursive residuals are standardized one-step-ahead prediction errors.
  Under the usual assumptions for the linear regression model they
  are (asymptotically) normal and i.i.d. (see Brown, Durbin, Evans, 1975,
  for details).
  
  The default method computes the initial coefficient estimates via QR
  decomposition, using \code{\link{lm.fit}}. In subsequent steps, the
  updating formula provided by Brown, Durbin, Evans (1975) is employed.
  To avoid numerical instabilities in the first steps (with typically
  small sample sizes), the QR solution is computed for comparison.
  When the relative difference (assessed bey \code{\link{all.equal}})
  between the two solutions falls below \code{tol}, only the updating
  formula is used in subsequent steps.
}
\value{
  A vector containing the recursive residuals.
}

\references{
  Brown R.L., Durbin J., Evans J.M. (1975), Techniques for
  testing constancy of regression relationships over time, \emph{Journal of the
  Royal Statistical Society}, B, \bold{37}, 149-163.
}


\seealso{\code{\link{efp}}}

\examples{
x <- rnorm(100) + rep(c(0, 2), each = 50)
rr <- recresid(x ~ 1)
plot(cumsum(rr), type = "l")

plot(efp(x ~ 1, type = "Rec-CUSUM"))
}
\keyword{regression}
