\name{efpFunctional}
\alias{efpFunctional}
\alias{simulateBMDist}
\alias{maxBM}
\alias{maxBB}
\alias{maxBMI}
\alias{maxBBI}
\alias{maxL2BB}
\alias{meanL2BB}
\alias{rangeBM}
\alias{rangeBB}
\alias{rangeBMI}
\alias{rangeBBI}
\alias{supLM}
\encoding{latin1}
\title{Functionals for Fluctuation Processes}
\description{Computes an object for aggregating, plotting and testing
empirical fluctuation processes.}
\usage{
efpFunctional(functional = list(comp = function(x) max(abs(x)), time = max),
  boundary = function(x) rep(1, length(x)),
  computePval = NULL, computeCritval = NULL,
  plotProcess = NULL, lim.process = "Brownian bridge",
  nobs = 10000, nrep = 50000, nproc = 1:20, h = 0.5,
  probs = c(0:84/100, 850:1000/1000))
}
\arguments{
  \item{functional}{either a function for aggregating fluctuation processes
    or a list with two functions names \code{"comp"} and \code{"time"}.}
  \item{boundary}{a boundary function.}
  \item{computePval}{a function for computing p values. If neither
    \code{computePval} nor \code{computeCritval} are specified
    critical values are simulated with settings as specified below.}
  \item{computeCritval}{a function for computing critical values. If neither
    \code{computePval} nor \code{computeCritval} are specified
    critical values are simulated with settings as specified below.}
  \item{plotProcess}{a function for plotting the empirical process,
    if set to \code{NULL} a suitable function is set up.}
  \item{lim.process}{a string specifying the limiting process.}
  \item{nobs}{integer specifying the number of observations of each
    Brownian motion simulated.}
  \item{nrep}{integer specifying the number of replications.}
  \item{nproc}{integer specifying for which number of processes
    Brownian motions should be simulated. If set to \code{NULL} only
    \code{nproc = 1} is used and all other values are derived from
    a Bonferroni correction.}
  \item{h}{bandwidth parameter for increment processes.}
  \item{probs}{numeric vector specifying for which probabilities 
    critical values should be tabulated.}
}

\details{
  \code{efpFunctional} computes an object of class \code{"efpFunctional"}
  which should know how to do inference based on empirical fluctuation processes
  (currently only for \code{\link{gefp}} objects and not yet for \code{\link{efp}}
  objects) and how to visualize the corresponding processes.   

  \code{efpFunctional}s for many frequently used test statistics are provided:
  \code{maxBB} for the double maximum function, \code{meanL2BB} for the Cramer-von Mises
  statistic, or \code{rangeBB} for the range statistic. Furthermore, \code{supLM}
  generates an object of class \code{"efpFunctional"} for a certain trimming parameter,
  see the examples. More details can be found in Zeileis (2006).    

  For setting up an \code{efpFunctional}, the functions
  \code{computeStatistic}, \code{computePval}, and \code{plotProcess} need to be
  supplied. These should have the following interfaces:
  \code{computeStatistic} should take a single argument which is the process
    itself, i.e., essentially a n x k matrix where n is the number of
    observations and k the number of processes (regressors).
  \code{computePval} should take two arguments: a scalar test statistic and the
    number of processes k.
  \code{plotProcess} should take two arguments: an object of class \code{"gefp"}
    and \code{alpha} the level of significance for any boundaries or critical
    values to be visualized.
}

\value{
  \code{efpFunctional} returns a list of class \code{"efpFunctional"} with components inlcuding:
  \item{plotProcess}{a function for plotting empirical fluctuation processes,}
  \item{computeStatistic}{a function for computing a test statistic from an empirical fluctuation process,}
  \item{computePval}{a function for computing the corresponding p value,}
  \item{computeCritval}{a function for computing critical values.}
}

\references{
Zeileis A. (2005), A Unified Approach to Structural Change Tests Based on
ML Scores, F Statistics, and OLS Residuals. \emph{Econometric Reviews}, \bold{24},
445--466. doi:10.1080/07474930500406053.

Zeileis A. (2006), Implementing a Class of Structural Change Tests: An
Econometric Computing Approach. \emph{Computational Statistics & Data Analysis}, 
\bold{50}, 2987--3008. doi:10.1016/j.csda.2005.07.001.

Zeileis A., Hornik K. (2007), Generalized M-Fluctuation Tests for Parameter
Instability, \emph{Statistica Neerlandica}, \bold{61}, 488--508.
doi:10.1111/j.1467-9574.2007.00371.x.
}

\seealso{\code{\link{efp}}, \code{\link{efpFunctional}}}

\examples{

data("BostonHomicide")
gcus <- gefp(homicides ~ 1, family = poisson, vcov = kernHAC,
	     data = BostonHomicide)
plot(gcus, functional = meanL2BB)	 
gcus
sctest(gcus, functional = meanL2BB)

y <- rnorm(1000)
x1 <- runif(1000)
x2 <- runif(1000)

## supWald statistic computed by Fstats()
fs <- Fstats(y ~ x1 + x2, from = 0.1)
plot(fs)
sctest(fs)

## compare with supLM statistic
scus <- gefp(y ~ x1 + x2, fit = lm)
plot(scus, functional = supLM(0.1))
sctest(scus, functional = supLM(0.1))

## seatbelt data
data("UKDriverDeaths")
seatbelt <- log10(UKDriverDeaths)
seatbelt <- cbind(seatbelt, lag(seatbelt, k = -1), lag(seatbelt, k = -12))
colnames(seatbelt) <- c("y", "ylag1", "ylag12")
seatbelt <- window(seatbelt, start = c(1970, 1), end = c(1984,12))

scus.seat <- gefp(y ~ ylag1 + ylag12, data = seatbelt)

## double maximum test
plot(scus.seat)
## range test
plot(scus.seat, functional = rangeBB)
## Cramer-von Mises statistic (Nyblom-Hansen test)
plot(scus.seat, functional = meanL2BB)
## supLM test
plot(scus.seat, functional = supLM(0.1))
}
\keyword{regression}
