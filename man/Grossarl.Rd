\name{Grossarl}
\alias{Grossarl}
\docType{data}
\encoding{latin1}
\title{Marriages, Births and Deaths in Grossarl}
\usage{data("Grossarl")}
\description{
Data about the number of marriages, illegitimate and legitimate
births, and deaths in the Austrian Alpine village Grossarl during the
18th and 19th century.
}
\format{
\code{Grossarl} is a data frame containing 6 annual time series
(1700 - 1899), 3 factors coding policy interventions and 1 vector
with the year (plain numeric).
 \describe{
    \item{marriages}{time series. Number of marriages,}
    \item{illegitimate}{time series. Number of illegitimate births,}
    \item{legitimate}{time series. Number of legitimate births,}
    \item{legitimate}{time series. Number of deaths,}
    \item{fraction}{time series. Fraction of illegitimate births,}
    \item{lag.marriages}{time series. Number of marriages in the previous year,}
    \item{politics}{ordered factor coding 4 different political regimes,}
    \item{morals}{ordered factor coding 5 different moral regulations,}
    \item{nuptiality}{ordered factor coding 5 different marriage restrictions,}
    \item{year}{numeric. Year of observation.}
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

\source{Parish registers provide the basic demographic series of baptisms 
and burials (which is almost equivalent to births and deaths in the study
area) and marriages. For more information see Veichtlbauer et al. (2006).}

\references{
Veichtlbauer O., Zeileis A., Leisch F. (2006),
The Impact Of Policy Interventions on a Pre-Industrial
Population System in the Austrian Alps, forthcoming.

Zeileis A., Veichtlbauer O. (2002), Policy Interventions
Affecting Illegitimacy in Preindustrial Austria:
A Structural Change Analysis, In R. Dutter (ed.),
\emph{Festschrift 50 Jahre Österreichische Statistische Gesellschaft}, 133-146,
Österreichische Statistische Gesellschaft.
}

\examples{
data("Grossarl")

## time series of births, deaths, marriages
###########################################

with(Grossarl, plot(cbind(deaths, illegitimate + legitimate, marriages),
  plot.type = "single", col = grey(c(0.7, 0, 0)), lty = c(1, 1, 3),
  lwd = 1.5, ylab = "annual Grossarl series"))
legend("topright", c("deaths", "births", "marriages"), col = grey(c(0.7, 0, 0)),
  lty = c(1, 1, 3), bty = "n")

## illegitimate births
######################
## lm + MOSUM
plot(Grossarl$fraction)
fm.min <- lm(fraction ~ politics, data = Grossarl)
fm.ext <- lm(fraction ~ politics + morals + nuptiality + marriages,
  data = Grossarl)
lines(ts(fitted(fm.min), start = 1700), col = 2)
lines(ts(fitted(fm.ext), start = 1700), col = 4)
mos.min <- efp(fraction ~ politics, data = Grossarl, type = "OLS-MOSUM")
mos.ext <- efp(fraction ~ politics + morals + nuptiality + marriages,
  data = Grossarl, type = "OLS-MOSUM")
plot(mos.min)
lines(mos.ext, lty = 2)

## dating
bp <- breakpoints(fraction ~ 1, data = Grossarl, h = 0.1)
summary(bp)
## RSS, BIC, AIC
plot(bp)
plot(0:8, AIC(bp), type = "b")

## probably use 5 or 6 breakpoints and compare with
## coding of the factors as used by us
##
## politics                   1803      1816 1850
## morals      1736 1753 1771 1803
## nuptiality                 1803 1810 1816      1883
##
## m = 5            1753 1785           1821 1856 1878
## m = 6       1734 1754 1785           1821 1856 1878
##              6    2    5              1    4    3

## fitted models
coef(bp, breaks = 6)
plot(Grossarl$fraction)
lines(fitted(bp, breaks = 6), col = 2)
lines(ts(fitted(fm.ext), start = 1700), col = 4)


## marriages
############
## lm + MOSUM
plot(Grossarl$marriages)
fm.min <- lm(marriages ~ politics, data = Grossarl)
fm.ext <- lm(marriages ~ politics + morals + nuptiality, data = Grossarl)
lines(ts(fitted(fm.min), start = 1700), col = 2)
lines(ts(fitted(fm.ext), start = 1700), col = 4)
mos.min <- efp(marriages ~ politics, data = Grossarl, type = "OLS-MOSUM")
mos.ext <- efp(marriages ~ politics + morals + nuptiality, data = Grossarl,
  type = "OLS-MOSUM")
plot(mos.min)
lines(mos.ext, lty = 2)

## dating
bp <- breakpoints(marriages ~ 1, data = Grossarl, h = 0.1)
summary(bp)
## RSS, BIC, AIC
plot(bp)
plot(0:8, AIC(bp), type = "b")

## probably use 3 or 4 breakpoints and compare with
## coding of the factors as used by us
##
## politics                   1803      1816 1850
## morals      1736 1753 1771 1803
## nuptiality                 1803 1810 1816      1883
##
## m = 3       1738                     1813      1875
## m = 4       1738      1794           1814      1875
##              2         4              1         3

## fitted models
coef(bp, breaks = 4)
plot(Grossarl$marriages)
lines(fitted(bp, breaks = 4), col = 2)
lines(ts(fitted(fm.ext), start = 1700), col = 4)
}

\keyword{datasets}
