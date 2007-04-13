\name{BostonHomicide}
\alias{BostonHomicide}
\title{Youth Homicides in Boston}
\usage{data("BostonHomicide")}
\description{
Data about the number of youth homicides in Boston during the
`Boston Gun Project'---a policing initiative aiming at lowering
homicide victimization among young people in Boston.
}
\format{
A data frame containing 6 monthly time series
and two factors coding seasonality and year, respectively.
 \describe{
    \item{homicides}{time series. Number of youth homicides.}
    \item{population}{time series. Boston population (aged 25-44),
       linearly interpolated from annual data.}
    \item{populationBM}{time series. Population of black males (aged 15-24),
       linearly interpolated from annual data.}
    \item{ahomicides25}{time series. Number of adult homicides (aged 25 and older).}
    \item{ahomicides35}{time series. Number of adult homicides (aged 35-44).}
    \item{unemploy}{time series. Teen unemployment rate (in percent).}
    \item{season}{factor coding the month.}
    \item{year}{factor coding the year.}
}
}

\details{The `Boston Gun Project' is a policing initiative aiming at lowering
youth homicides in Boston. The project began in early 1995 and implemented the
so-called `Operation Ceasefire' intervention which began in the late spring of 1996.

More information is available at:

  \url{http://www.ksg.harvard.edu/criminaljustice/research/bgp.htm}
}

\source{Piehl et al. (2004), Figure 1, Figure 3, and Table 1.

From the table it is not clear how the data should be linearly interpolated.
Here, it was chosen to use the given observations for July of the corresponding
year and then use \code{\link{approx}} with \code{rule = 2}.
}

\references{
Piehl A.M., Cooper S.J., Braga A.A., Kennedy D.M. (2003), Testing for Structural
Breaks in the Evaluation of Programs, \emph{The Review of Economics and Statistics},
\bold{85}(3), 550-558.

Kennedy D.M., Piehl A.M., Braga A.A. (1996), Youth Violence in Boston: Gun Markets,
Serious Youth Offenders, and a Use-Reduction Strategy, \emph{Law and Contemporary Problems},
\bold{59}, 147-183.
}

\examples{
data("BostonHomicide")
attach(BostonHomicide)

## data from Table 1
tapply(homicides, year, mean)
populationBM[0:6*12 + 7]
tapply(ahomicides25, year, mean)
tapply(ahomicides35, year, mean)
population[0:6*12 + 7]
unemploy[0:6*12 + 7]

## model A
## via OLS
fmA <- lm(homicides ~ populationBM + season)
anova(fmA)
## as GLM
fmA1 <- glm(homicides ~ populationBM + season, family = poisson)
anova(fmA1, test = "Chisq")

## model B & C
fmB <- lm(homicides ~ populationBM + season + ahomicides25)
fmC <- lm(homicides ~ populationBM + season + ahomicides25 + unemploy)

detach(BostonHomicide)
}

\keyword{datasets}
