\name{sctest.default}
\alias{sctest.default}

\title{Structural Change Tests in Parametric Models}

\description{
Performs model-based tests for structural change (or parameter instability)
in parametric models.
}

\usage{
\method{sctest}{default}(x, order.by = NULL, functional = maxBB,
  vcov = NULL, scores = estfun, decorrelate = TRUE, sandwich = TRUE,
  parm = NULL, plot = FALSE, from = 0.1, to = NULL, nobs = NULL,
  nrep = 50000, width = 0.15, xlab = NULL, \dots)
}

\arguments{
  \item{x}{a model object. The model class can in principle be arbitrary
    but needs to provide suitable methods for extracting the \code{scores}
    and associated variance-covariance matrix \code{vcov}.}
  \item{order.by}{either a vector \code{z} or a formula with a single explanatory
    variable like \code{~ z}. The observations in the model
    are ordered by the size of \code{z}. If set to \code{NULL} (the
    default) the observations are assumed to be ordered (e.g., a
    time series).}
  \item{functional}{either a character specification of the functional
    to be used or an \code{\link{efpFunctional}} object. For a list
    of functionals see the details.}
  \item{vcov}{a function to extract the covariance matrix
     for the coefficients of the fitted model:
     \code{vcov(x, order.by = order.by, data = data)}.
     Alternatively, the character string \code{"info"}, for details see
     below.}
  \item{scores}{a function which extracts the scores or estimating
     function from the fitted object: \code{scores(x)}, by default
     this is \code{\link[sandwich]{estfun}}.}
  \item{decorrelate}{logical. Should the process be decorrelated?}
  \item{sandwich}{logical. Is the function \code{vcov} the full sandwich
     estimator or only the meat?}
  \item{parm}{integer or character specifying the component of the estimating
    functions which should be used (by default all components are used).}
  \item{plot}{logical. Should the result of the test also be visualized?}
  \item{from, to}{numeric. In case the \code{functional} is \code{"supLM"}
    (or equivalently \code{"maxLM"}), \code{from} and \code{to} can be
    passed to the \code{\link{supLM}} functional.}
  \item{nobs, nrep}{numeric. In case the \code{functional} is \code{"maxLMo"},
    \code{nobs} and \code{nrep} are passed to the \code{\link{catL2BB}} functional.}
  \item{width}{numeric. In case the \code{functional} is \code{"MOSUM"},
    the bandwidth \code{width} is passed to the \code{\link{maxMOSUM}}
    functional.}
  \item{xlab, \dots}{graphical parameters passed to the plot method (in case
    \code{plot = TRUE}).}
}

\details{
\code{sctest.default} is a convenience interface to \code{\link{gefp}} for
structural change tests (or parameter instability tests) in general
parametric models. It proceeds in the following steps:

\enumerate{
\item The generalized empirical fluctuation process (or score-based CUSUM process)
  is computed via \code{scus <- gefp(x, fit = NULL, \dots)} where \code{\dots}
  comprises the arguments \code{order.by}, \code{vcov}, \code{scores}, \code{decorrelate},
  \code{sandwich}, \code{parm} that are simply passed on to \code{\link{gefp}}.
\item The empirical fluctuation process is visualized (if \code{plot = TRUE}) via
  \code{plot(scus, functional = functional, \dots)}.
\item The empirical fluctuation is assessed by the corresponding significance test
  via \code{sctest(scus, functional = functional)}.
}

The main motivation for prociding the convenience interface is that these three
steps can be easily carried out in one go along with a two convenience options:

\enumerate{
\item By default, the covariance is computed by an outer-product of gradients
  estimator just as in \code{gefp}. This is always available based on the \code{scores}.
  Additionally, by setting \code{vcov = "info"}, the corresponding information
  matrix can be used. Then the average information is assumed to be provided by
  the \code{vcov} method for the model class. (Note that this is only sensible
  for models estimated by maximum likelihood.)
\item Instead of providing the \code{functional} by an \code{\link{efpFunctional}}
  object, the test labels employed by Merkle and Zeileis (2013) and Merkle, Fan,
  and Zeileis (2013) can be used for convenience. Namely, for continuous numeric
  orderings, the following functionals are available:
  \code{functional = "DM"} or \code{"dmax"} provides the double-maximum test (\code{\link{maxBB}}).
  \code{"CvM"} is the Cramer-von Mises functional \code{\link{meanL2BB}}.
  \code{"supLM"} or equivalently \code{"maxLM"} is Andrews' supLM test
  (\code{\link{supLM}}). \code{"MOSUM"} or \code{"maxMOSUM"} is the MOSUM
  functional (\code{\link{maxMOSUM}}), and \code{"range"} is the range
  functional \code{\link{rangeBB}}. Furthermore, several functionals suitable
  for (ordered) categorical \code{order.by} variables are provided:
  \code{"LMuo"} is the unordered LM test (\code{\link{catL2BB}}),
  \code{"WDMo"} is the weighted double-maximum test for ordered variables
  (\code{\link{ordwmax}}), and \code{"maxLMo"} is the maxLM test for
  ordered variables (\code{\link{ordL2BB}}).
}

The theoretical model class is introduced in Zeileis and Hornik (2007) with a
unifying view in Zeileis (2005), especially from an econometric perspective.
Zeileis (2006) introduces the underling computational tools \code{gefp} and
\code{efpFunctional}.

Merkle and Zeileis (2013) discuss the methods in the context of measurement
invariance which is particularly relevant to psychometric models for cross section
data. Merkle, Fan, and Zeileis (2013) extend the results to ordered categorical
variables.

Zeileis, Shah, and Patnaik (2013) provide a unifying discussion in the context
of time series methods, specifically in financial econometrics.
}

\value{
  An object of class \code{"htest"} containing:
  \item{statistic}{the test statistic,}
  \item{p.value}{the corresponding p value,}
  \item{method}{a character string with the method used,}
  \item{data.name}{a character string with the data name.}
}

\references{
Merkle E.C., Zeileis A. (2013), Tests of Measurement Invariance without Subgroups:
A Generalization of Classical Methods. \emph{Psychometrika}, \bold{78}(1), 59--82.
doi:10.1007/S11336-012-9302-4

Merkle E.C., Fan J., Zeileis A. (2013), Testing for Measurement Invariance with
Respect to an Ordinal Variable. \emph{Psychometrika}, Forthcoming.

Zeileis A. (2005), A Unified Approach to Structural Change Tests Based on
ML Scores, F Statistics, and OLS Residuals. \emph{Econometric Reviews}, \bold{24},
445--466. doi:10.1080/07474930500406053.

Zeileis A. (2006), Implementing a Class of Structural Change Tests: An
Econometric Computing Approach. \emph{Computational Statistics & Data Analysis}, 
\bold{50}, 2987--3008. doi:10.1016/j.csda.2005.07.001.

Zeileis A., Hornik K. (2007), Generalized M-Fluctuation Tests for Parameter
Instability, \emph{Statistica Neerlandica}, \bold{61}, 488--508.
doi:10.1111/j.1467-9574.2007.00371.x.

Zeileis A., Shah A., Patnaik I. (2010), Testing, Monitoring, and Dating Structural
Changes in Exchange Rate Regimes, \emph{Computational Statistics and Data Analysis},
Forthcoming. doi:10.1016/j.csda.2009.12.005.
}

\seealso{\code{\link{gefp}}, \code{\link{efpFunctional}}}

\examples{
## Zeileis and Hornik (2007), Section 5.3, Figure 6
data("Grossarl")
m <- glm(cbind(illegitimate, legitimate) ~ 1, family = binomial, data = Grossarl,
  subset = time(fraction) <= 1800)
sctest(m, order.by = 1700:1800, functional = "CvM")
}

\keyword{htest}
