\name{BostonHomicide}
\alias{BostonHomicide}
\title{Youth Homicides in Boston}
\usage{data(BostonHomicide)}
\description{
Data about the number of youth homicides in Boston during the
`Boston Gun Project'---a policing initiative aiming at lowering
homicide victimization among young people in Boston.
}
\format{
A data frame containing 2 monthly time series
and a factor coding seasonality.
 \describe{
    \item{homicides}{time series. Number of youth homicides,}
    \item{population}{time series. Boston population (aged 25-44),
       linearly interpolated from annual data.}
    \item{season}{factor coding the month.}
}
}

\details{The `Boston Gun Project' is a policing initiative aiming at lowering
youth homicides in Boston. The project began in early 1995 and implemented the
so-called `Operation Ceasefire' intervention which began in the late spring of 1996.

More information is available at:

  \url{http://www.ksg.harvard.edu/criminaljustice/research/bgp.htm}
}

\source{Cooper et al. (2001), Figure 1 and Table I.}

\references{
Cooper S.J., Piehl A.M., Braga A.A., Kennedy D.M. (2001), Testing for Structural
Breaks in the Evaluation of Programs, Unpublished Paper, John F. Kennedy School
of Government.

Kennedy D.M., Piehl A.M., Braga A.A. (1996), Youth Violence in Boston: Gun Markets,
Serious Youth Offenders, and a Use-Reduction Strategy, \emph{Law and Contemporary Problems},
\bold{59}, 147-183.
}

\examples{
data(BostonHomicide)

fm <- glm(homicides ~ population + season, data = BostonHomicide,
          family = poisson)
anova(fm, test = "F")
}

\keyword{datasets}
