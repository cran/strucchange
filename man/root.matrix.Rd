\name{root.matrix}
\alias{root.matrix}
\title{Root of a Matrix}
\description{Computes the root of a symmetric and positive semidefinite matrix.}
\usage{
root.matrix(X)
}
\arguments{
  \item{X}{a symmetric and positive semidefinite matrix}
}

\value{a symmetric matrix of same dimensions as \code{X}}

\examples{
X <- matrix(c(1,2,2,8), ncol=2)
test <- root.matrix(X)
## control results
X
test \%*\% test
}
\keyword{algebra}
