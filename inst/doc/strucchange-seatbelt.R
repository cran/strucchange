## ----include = FALSE----------------------------------------------------------
library("strucchange")
knitr::opts_chunk$set(
  engine = "R",
  collapse = TRUE,
  comment = "##",
  message = FALSE,
  warning = FALSE,
  echo = TRUE,
  fig.height = 5,
  fig.width = 6
)
options(digits = 4)

## ----data---------------------------------------------------------------------
data("UKDriverDeaths", package = "datasets")
plot(log10(UKDriverDeaths))

## ----preprocessing------------------------------------------------------------
library("strucchange")
seatbelt <- UKDriverDeaths |>
  as.zoo() |>
  log10() |>
  lag(c(0, -1, -12)) |>
  setNames(c("y", "ylag1", "ylag12")) |>
  na.trim() |>
  as.ts()

## ----lm-----------------------------------------------------------------------
m <- lm(y ~ ylag1 + ylag12, data = seatbelt)
summary(m)

## ----re-----------------------------------------------------------------------
re <- efp(y ~ ylag1 + ylag12, data = seatbelt, type = "RE")
sctest(re)
plot(re)

## ----fstats-------------------------------------------------------------------
fs <- Fstats(y ~ ylag1 + ylag12, data = seatbelt, from = 0.1)
sctest(fs)
plot(fs)

## ----gefp---------------------------------------------------------------------
scus <- gefp(m, fit = NULL)
sctest(scus, functional = supLM(0.1))
plot(scus, functional = supLM(0.1))

## ----sctest, include=FALSE----------------------------------------------------
sctest(m, functional = supLM(0.1), plot = TRUE)

## ----breakpoints--------------------------------------------------------------
bp <- breakpoints(y ~ ylag1 + ylag12, data = seatbelt, h = 0.1)
summary(bp)
plot(bp)

## ----confint------------------------------------------------------------------
confint(bp, breaks = 2)

## ----coef---------------------------------------------------------------------
coef(bp, breaks = 2)

## ----plot---------------------------------------------------------------------
plot(seatbelt[, "y"], ylab = expression(log[10](casualties)), col = "lightgray", lwd = 2)
lines(fitted(bp, breaks = 2))
lines(confint(bp, breaks = 2))

## ----mefp---------------------------------------------------------------------
sb <- window(seatbelt, start = c(1976, 1), end = c(1982, 12))
me <- mefp(y ~ ylag1 + ylag12, data = sb, type = "ME", h = 0.5)

## ----monitor------------------------------------------------------------------
sb <- window(seatbelt, start = c(1976, 1))
mon <- monitor(me)
plot(mon)

