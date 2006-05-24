\name{plot.Fstats}
\alias{plot.Fstats}
\alias{lines.Fstats}
\title{Plot F Statistics}
\description{Plotting method for objects of class \code{"Fstats"}}
\usage{
\method{plot}{Fstats}(x, pval = FALSE, asymptotic = FALSE, alpha = 0.05,
    boundary = TRUE, aveF = FALSE, xlab = "Time", ylab = NULL,
    ylim = NULL, ...)
}

\arguments{
 \item{x}{an object of class \code{"Fstats"}.}
 \item{pval}{logical. If set to \code{TRUE} the corresponding p values instead
     of the original F statistics will be plotted.}
 \item{asymptotic}{logical. If set to \code{TRUE} the asymptotic (chi-square)
     distribution instead of the exact (F) distribution will be used to compute
     the p values (only if \code{pval} is \code{TRUE}).}
 \item{alpha}{numeric from interval (0,1) indicating the confidence level for
     which the boundary of the supF test will be computed.}
 \item{boundary}{logical. If set to \code{FALSE} the boundary will be computed
     but not plotted.}
 \item{aveF}{logical. If set to \code{TRUE} the boundary of the aveF test will
     be plotted. As this is a boundary for the mean of the F statistics rather
     than for the F statistics themselves a dashed line for the mean of the F
     statistics will also be plotted.}
 \item{xlab, ylab, ylim, ...}{high-level \code{\link{plot}} function parameters.}}

\references{
Andrews D.W.K. (1993), Tests for parameter instability and structural
change with unknown change point, \emph{Econometrica}, \bold{61}, 821-856.

Hansen B. (1992), Tests for parameter instability in regressions with I(1)
processes, \emph{Journal of Business & Economic Statistics}, \bold{10}, 321-335.

Hansen B. (1997), Approximate asymptotic p values for structural-change
tests, \emph{Journal of Business & Economic Statistics}, \bold{15}, 60-67. }

\seealso{\code{\link{Fstats}}, \code{\link{boundary.Fstats}},
\code{\link{sctest.Fstats}}}

\examples{
## Load dataset "nhtemp" with average yearly temperatures in New Haven
data("nhtemp")
## plot the data
plot(nhtemp)

## test the model null hypothesis that the average temperature remains
## constant over the years for potential break points between 1941
## (corresponds to from = 0.5) and 1962 (corresponds to to = 0.85)
## compute F statistics
fs <- Fstats(nhtemp ~ 1, from = 0.5, to = 0.85)
## plot the F statistics
plot(fs, alpha = 0.01)
## and the corresponding p values
plot(fs, pval = TRUE, alpha = 0.01)
## perform the aveF test
sctest(fs, type = "aveF")
}
\keyword{hplot}
