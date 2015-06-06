\name{catL2BB}
\alias{catL2BB}
\alias{ordL2BB}
\alias{ordwmax}

\title{Generators for efpFunctionals along Categorical Variables}

\description{
Generators for \code{efpFunctional} objects suitable for aggregating
empirical fluctuation processes to test statistics along (ordinal)
categorical variables.
}

\usage{
catL2BB(freq)
ordL2BB(freq, nproc = NULL, nrep = 1e5, probs = c(0:84/100, 850:1000/1000), \dots)
ordwmax(freq, algorithm = mvtnorm::GenzBretz(), \dots)
}

\arguments{
  \item{freq}{object specifying the category frequencies for the
    categorical variable to be used for aggregation: either a
    \code{\link{gefp}} object, a \code{\link{factor}}, or a numeric
    vector with either absolute or relative category frequencies.}
  \item{nproc}{numeric. Number of processes used for simulating
    from the asymptotic distribution (passed to \code{\link{efpFunctional}}).
    If \code{feq} is a \code{\link{gefp}} object, then its number of
    processes is used by default.}
  \item{nrep}{numeric. Number of replications used for simulating
    from the asymptotic distribution (passed to \code{\link{efpFunctional}}).}
  \item{probs}{numeric vector specifying for which probabilities 
    critical values should be tabulated.}
  \item{\dots}{further arguments passed to \code{\link{efpFunctional}}.}
  \item{algorithm}{algorithm specification passed to \code{\link[mvtnorm]{pmvnorm}}
    for computing the asymptotic distribution.}
}

\details{
  Merkle, Fan, and Zeileis (2014) discuss three functionals that are 
  suitable for aggregating empirical fluctuation processes along categorical
  variables, especially ordinal variables. The functions \code{catL2BB},
  \code{ordL2BB}, and \code{ordwmax} all require a specification of the
  relative frequencies within each category (which can be computed from
  various specifications, see arguments). All of them employ
  \code{\link{efpFunctional}} (Zeileis 2006) internally to set up an
  object that can be employed with \code{\link{gefp}} fluctuation
  processes.
  
  \code{catL2BB} results in a chi-squared test. This is essentially
  the LM test counterpart to the likelihood ratio test that assesses
  a split into unordered categories.
  
  \code{ordL2BB} is the ordinal counterpart to \code{\link{supLM}}
  where aggregation is done along the ordered categories (rather than
  continuously). The asymptotic distribution is non-standard and needs
  to be simulated via \code{\link[mvtnorm]{rmvnorm}} for every combination
  of frequencies and number of processes. This can be somewhat
  time-consuming, hence it is recommended to store the result of 
  \code{ordL2BB} in case it needs to be applied several \code{\link{gefp}}
  fluctuation processes.

  \code{ordwmax} is a weighted double maximum test based on ideas
  previously suggested by Hothorn and Zeileis (2008) in the context of
  maximally selected statistics. The asymptotic distribution is
  (multivariate) normal and computed by means of \code{\link[mvtnorm]{pmvnorm}}.
}

\value{
  An object of class \code{efpFunctional}.
}

\references{
Hothorn T., Zeileis A. (2008), Generalized Maximally Selected Statistics.
\emph{Biometrics}, \bold{64}, 1263--1269.

Merkle E.C., Fan J., Zeileis A. (2014), Testing for Measurement Invariance with
Respect to an Ordinal Variable. \emph{Psychometrika}, \bold{79}(4), 569--584.
doi:10.1007/S11336-013-9376-7.

Zeileis A. (2006), Implementing a Class of Structural Change Tests: An
Econometric Computing Approach. \emph{Computational Statistics & Data Analysis}, 
\bold{50}, 2987--3008. doi:10.1016/j.csda.2005.07.001.
}

\seealso{\code{\link{efpFunctional}}, \code{\link{gefp}}}

\examples{
## artificial data
set.seed(1)
d <- data.frame(
  x = runif(200, -1, 1),
  z = factor(rep(1:4, each = 50)),
  err = rnorm(200)
)
d$y <- rep(c(0.5, -0.5), c(150, 50)) * d$x + d$err

## empirical fluctuation process
scus <- gefp(y ~ x, data = d, fit = lm, order.by = ~ z)

## chi-squared-type test (unordered LM-type test)
LMuo <- catL2BB(scus)
plot(scus, functional = LMuo)
sctest(scus, functional = LMuo)

## ordinal maxLM test (with few replications only to save time)
maxLMo <- ordL2BB(scus, nrep = 10000)
plot(scus, functional = maxLMo)
sctest(scus, functional = maxLMo)

## ordinal weighted double maximum test
WDM <- ordwmax(scus)
plot(scus, functional = WDM)
sctest(scus, functional = WDM)
}

\keyword{regression}
