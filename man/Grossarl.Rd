\name{Grossarl}
\alias{Grossarl}
\title{Marriages and Births in Grossarl}
\usage{data(Grossarl)}
\description{
Data about the number of marriages and illegitimate and legitimate
births in the Austrian Alpine village Grossarl in the 18th and 19th
century.
}
\format{
\code{Grossarl} is a data frame containing 5 annual time series
(1700 - 1899) and 3 factors coding policy interventions.
 \describe{
    \item{marriages}{time series. Number of marriages,}
    \item{illegitimate}{time series. Number of illegitimate births,}
    \item{legitimate}{time series. Number of legitimate births,}
    \item{fraction}{time series. Fraction of illegitimate births,}
    \item{lag.marriages}{time series. Number of marriages in the previous year,}
    \item{politics}{factor coding 4 different political regimes,}
    \item{morals}{factor coding 5 different moral regulations,}
    \item{nuptiality}{factor coding 5 different marriage restrictions.}
}
}

\details{The data frame contains historical demographic data from
Grossarl, a village in the Alpine region of Salzburg, Austria,
during the 18th and 19th century.
During this period, the total population of Grossarl did not vary much on the whole,
with the very exception of the period of the protestant emigrations
in 1731/32.

Especially
during the archbishopric, moral interventions aimed at lowering
the proportion of illegitimate baptisms. For details see the references.}

\source{Parish registers provide the basic
demographic series of baptisms (which is almost equivalent to births
in the study area) and marriages. For
more information see Veichtlbauer et al. (2002).}

\references{
Veichtlbauer O., Hanser E., Zeileis A., Leisch F. (2002),
The Impact Of Policy Interventions on a Pre-Industrial
Population System in the Austrian Alps, forthcoming.

Zeileis A., Veichtlbauer O. (2002), Policy Interventions
Affecting Illegitimacy in Preindustrial Austria:
A Structural Change Analysis, In R. Dutter (ed.),
\emph{Festschrift 50 Jahre Österreichische Statistische Gesellschaft}, 133-146,
Österreichische Statistische Gesellschaft,
\url{http://www.statistik.tuwien.ac.at/oezstat/}.
}

\examples{
data(Grossarl)

## illegitimate births
######################
## lm + MOSUM
plot(Grossarl$fraction)
fm.min <- lm(fraction ~ politics, data = Grossarl)
fm.max <- lm(fraction ~ politics + morals + nuptiality + lag.marriages,
             data = Grossarl)
fm.final <- step(fm.max)
lines(ts(fitted(fm.min), start = 1700), col = 3)
lines(ts(fitted(fm.final), start = 1700), col = 4)
mos.min <- efp(fraction ~ politics, data = Grossarl, type = "OLS-MOSUM")
mos.final <- efp(fraction ~ politics + morals + nuptiality, data = Grossarl,
                 type = "OLS-MOSUM")
plot(mos.min)
lines(mos.final, lty = 2)

## dating
bp <- breakpoints(fraction ~ 1, data = Grossarl, h = 0.1)
summary(bp)
## RSS, BIC, AIC
plot(bp)
plot(0:8, AIC(bp), type = "b")

## probably use 5 (or maybe 6) breakpoints and compare with
## coding of the factors as used by us
##
## politics                   1803      1816 1850
## morals      1736 1753 1771 1803
## nuptiality                 1803 1810 1816      1883
##
## m = 5            1753 1785           1821 1856 1878
## m = 6       1734 1754 1785           1821 1856 1878
##              6    2    5              1    4    3

fm.bp <- lm(fraction ~ breakfactor(breakpoints(bp, breaks = 6)),
            data = Grossarl)

plot(Grossarl$fraction)
lines(fitted(fm.final), col = 3)
lines(fitted(fm.bp), col = 4)



## marriages
############
## lm + MOSUM
plot(Grossarl$marriages)
fm.min <- lm(marriages ~ politics, data = Grossarl)
fm.final <- lm(marriages ~ politics + morals + nuptiality, data = Grossarl)
lines(ts(fitted(fm.min), start = 1700), col = 3)
lines(ts(fitted(fm.final), start = 1700), col = 4)
mos.min <- efp(marriages ~ politics, data = Grossarl, type = "OLS-MOSUM")
mos.final <- efp(marriages ~ politics + morals + nuptiality, data = Grossarl,
                 type = "OLS-MOSUM")
plot(mos.min)
lines(mos.final, lty = 2)

## dating
bp <- breakpoints(marriages ~ 1, data = Grossarl, h = 0.1)
summary(bp)
## RSS, BIC, AIC
plot(bp)
plot(0:8, AIC(bp), type = "b")

## probably use 3 (or maybe 4) breakpoints and compare with
## coding of the factors as used by us
##
## politics                   1803      1816 1850
## morals      1736 1753 1771 1803
## nuptiality                 1803 1810 1816      1883
##
## m = 3       1738                     1813      1875
## m = 4       1738      1794           1814      1875
##              2         4              1         3
fm.bp <- lm(marriages ~ breakfactor(breakpoints(bp, breaks = 4)),
            data = Grossarl)

plot(Grossarl$marriages)
lines(fitted(fm.final), col = 3)
lines(fitted(fm.bp), col = 4)
}

\keyword{datasets}
