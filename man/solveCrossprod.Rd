\name{solveCrossprod}
\alias{solveCrossprod}
\title{Inversion of X'X}
\description{Computes the inverse of the cross-product of a matrix X.}
\usage{
solveCrossprod(X, method = c("qr", "chol", "solve"))
}
\arguments{
  \item{X}{a matrix, typically a regressor matrix.}
  \item{method}{a string indicating whether the QR decomposition,
    the Cholesky decomposition or \code{solve} should be used.}
}
\details{Using the Cholesky decomposition of X'X (as computed by \code{crossprod(X)})
is computationally faster and preferred to \code{solve(crossprod(X))}. Using the
QR decomposition of X is slower but should be more accurate.}

\value{a matrix containing the inverse of \code{crossprod(X)}.}

\examples{
X <- cbind(1, rnorm(100))
solveCrossprod(X)
solve(crossprod(X))
}
\keyword{algebra}
