\name{gefp}
\alias{gefp}
\alias{print.gefp}
\alias{sctest.gefp}
\alias{plot.gefp}
\alias{time.gefp}
\alias{print.gefp}
\title{Generalized Empirical M-Fluctuation Processes}
\description{Computes an empirical M-fluctuation process 
from the scores of a fitted model.}
\usage{
gefp(\dots, fit = glm, scores = estfun, vcov = NULL,
  decorrelate = TRUE, sandwich = TRUE, order.by = NULL,
  fitArgs = NULL, parm = NULL, data = list())
}
\arguments{
  \item{\dots}{specification of some model which is passed together
     with \code{data} to the \code{fit} function: \code{fm <- fit(\dots, data = data)}.
     If \code{fit} is set to \code{NULL} the first argument \code{\dots}
     is assumed to be already the fitted model \code{fm}
     (all other arguments in \code{\dots} are ignored and a warning
     is issued in this case).}
  \item{fit}{a model fitting function, typically \code{\link{lm}},
     \code{\link{glm}} or \code{\link[MASS]{rlm}}.}
  \item{scores}{a function which extracts the scores or estimating
     function from the fitted object: \code{scores(fm)}.}
  \item{vcov}{a function to extract the covariance matrix
     for the coefficients of the fitted model:
     \code{vcov(fm, order.by = order.by, data = data)}.}
  \item{decorrelate}{logical. Should the process be decorrelated?}
  \item{sandwich}{logical. Is the function \code{vcov} the sandwich
     estimator or only the middle part?}
  \item{order.by}{Either a vector \code{z} or a formula with a single explanatory
    variable like \code{~ z}. The observations in the model
    are ordered by the size of \code{z}. If set to \code{NULL} (the
    default) the observations are assumed to be ordered (e.g., a
    time series).}
  \item{fitArgs}{List of additional arguments which could be passed to
    the \code{fit} function. Usually, this is not needed and \code{\dots}
    will be sufficient to pass arguments to \code{fit}.}
  \item{parm}{integer or character specifying the component of the estimating
    functions which should be used (by default all components are used).}
 \item{data}{an optional data frame containing the variables in the \code{\dots}
    specification and the \code{order.by} model. By default the variables are
    taken from the environment which \code{gefp} is called from.}
}

\value{
  \code{gefp} returns a list of class \code{"gefp"} with components inlcuding
  \itemize
  \item{process}{the fitted empirical fluctuation process of class
    \code{"zoo"},}
  \item{nreg}{the number of regressors,}
  \item{nobs}{the number of observations,}
  \item{fit}{the fit function used,}
  \item{scores}{the scores function used,}
  \item{fitted.model}{the fitted model.}
}

\references{
Zeileis A., Hornik K. (2003), Generalized M-Fluctuation Tests for Parameter
Instability, Report 80, SFB "Adaptive Information Systems and Modelling in Economics
and Management Science", Vienna University of Economics,
\url{http://www.wu-wien.ac.at/am/reports.htm#80}.

Zeileis A. (2004), Implementing a Class of Structural Change Tests: An
Econometric Computing Approach. Report 7, Department of Statistics and Mathematics,
Wirtschaftsuniversität Wien, Research Report Series,
\url{http://epub.wu-wien.ac.at/dyn/virlib/wp/mediate/epub-wu-01_750.pdf?ID=epub-wu-01_750}.
}

\seealso{\code{\link{efp}}, \code{\link{efpFunctional}}}

\examples{

if("package:sandwich" \%in\% search() || require(sandwich)) {
  data(BostonHomicide)
  gcus <- gefp(homicides ~ 1, family = poisson, vcov = kernHAC,
               data = BostonHomicide)
  plot(gcus, aggregate = FALSE)	   
  gcus
  sctest(gcus)
}

}
\concept{M-fluctuation}
\concept{fluctuation test}
\concept{maximum likelihood scores}
\concept{structural change}
\keyword{regression}
