\name{sctest}
\alias{sctest}

\title{Structural Change Tests}

\description{
Generic function for performing structural change tests.
}

\usage{
sctest(x, \dots)
}

\arguments{
  \item{x}{an object.}
  \item{\dots}{arguments passed to methods.}
}

\details{
\code{sctest} is a generic function for performing/extracting structural
change tests based on various types of objects. The \code{strucchange} package
provides various types of methods.

First, structural change tests based on
F statistics in linear regression models (\code{\link{Fstats}}),
empirical fluctuation processes in linear regression models (\code{\link{efp}}),
and generalized empirical fluctuation processes in parametric models (\code{\link{gefp}})
are available in the corresponding \code{sctest} methods.

Second, convenience interfaces for carrying out structural change tests
in linear regression models and general parametric models are provided in
\code{\link{sctest.formula}} and \code{\link{sctest.default}}, respectively.
}

\value{
  An object of class \code{"htest"} containing:
  \item{statistic}{the test statistic,}
  \item{p.value}{the corresponding p value,}
  \item{method}{a character string with the method used,}
  \item{data.name}{a character string with the data name.}
}

\references{
Zeileis A., Leisch F., Hornik K., Kleiber C. (2002), \code{strucchange}:
An R Package for Testing for Structural Change in Linear Regression Models,
\emph{Journal of Statistical Software}, \bold{7}(2), 1-38.
URL \url{http://www.jstatsoft.org/v07/i02/}.

Zeileis A. (2006), Implementing a Class of Structural Change Tests: An
Econometric Computing Approach. \emph{Computational Statistics & Data Analysis}, 
\bold{50}, 2987--3008. doi:10.1016/j.csda.2005.07.001.
}

\seealso{\code{\link{sctest.formula}}, \code{\link{sctest.default}},
  \code{\link{sctest.Fstats}}, \code{\link{sctest.efp}}, \code{\link{sctest.gefp}}}
