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

\details{\code{efpFunctional} computes an object of class \code{"efpFunctional"}
  which should know how to do inference based on empirical fluctuation processes
  (currently only for \code{\link{gefp}} objects and not yet for \code{\link{efp}}
  objects) and how to visualize the corresponding processes. In particular, it has
  slots for the functions \code{computeStatistic}, \code{computePval} and \code{plotProcess}.
  These should have the following interfaces:
  {\itemize
    \item{computeStatistic}{should take a single argument which is the process
      itself, i.e., essentially a n x k matrix where n is the number of
      observations and k the number of processes (regressors).}
    \item{computePval}{should take two arguments: a scalar test statistic and the
      number of processes k}
    \item{plotProcess}{should take two arguments: an object of class \code{"gefp"}
      and \code{alpha} the level of significance for any boundaries or critical
      values to be visualized.}}
}

\value{
  \code{efpFunctional} returns a list of class \code{"efpFunctional"} with components inlcuding
  \itemize
  \item{plotProcess}{a function for plotting empirical fluctuation processes,}
  \item{computeStatistic}{a function for computing a test statistic from an empirical fluctuation process,}
  \item{computePval}{a function for computing the corresponding p value,}
  \item{computeCritval}{a function for computing critical values.}
}

\references{
Zeileis A., Hornik K. (2003), Generalized M-Fluctuation Tests for Parameter
Instability, Report 80, SFB "Adaptive Information Systems and Modelling in Economics
and Management Science", Vienna University of Economics,
\url{http://www.wu-wien.ac.at/am/reports.htm#80}.
}

\seealso{\code{\link{efp}}, \code{\link{efpFunctional}}}

\examples{

if("package:sandwich" \%in\% search() || require(sandwich)) {
  data(BostonHomicide)
  gcus <- gefp(homicides ~ 1, family = poisson, vcov = kernHAC,
               data = BostonHomicide)
  plot(gcus, functional = meanL2BB)	   
  gcus
  sctest(gcus, functional = meanL2BB)
}
}
\keyword{regression}