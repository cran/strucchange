if(!("package:zoo" %in% search() || require(zoo))) warning("could not load package zoo")
if(!("package:sandwich" %in% search() || require(sandwich))) warning("could not load package sandwich")

.First.lib <- function(lib, pkg) {
  if(as.numeric(R.Version()$minor) < 7) {
    autoload("confint", "MASS")
  }
}

maxBB <- efpFunctional(lim.process = "Brownian bridge",
  functional = list(comp = function(x) max(abs(x)), time = max),
  computePval = function(x, nproc = 1)
    pvalue.efp(x, lim.process = "Brownian bridge", functional = "max", alt.boundary = FALSE, k = nproc))

maxBBI <- efpFunctional(lim.process = "Brownian bridge increments",
  functional = list(comp = function(x) max(abs(x)), time = max),
  computePval = function(x, nproc = 1)
    pvalue.efp(x, lim.process = "Brownian bridge increments", functional = "max", alt.boundary = FALSE, k = nproc))

maxBM <- efpFunctional(lim.process = "Brownian motion",
  functional = list(comp = function(x) max(abs(x)), time = max),
  boundary = function(x) 1 + 2*x,
  computePval = function(x, nproc = 1)
    pvalue.efp(x, lim.process = "Brownian motion", functional = "max", alt.boundary = FALSE, k = nproc))

maxBMI <- efpFunctional(lim.process = "Brownian motion increments",
  functional = list(comp = function(x) max(abs(x)), time = max),
  computePval = function(x, nproc = 1)
    pvalue.efp(x, lim.process = "Brownian motion increments", functional = "max", alt.boundary = FALSE, k = nproc))


rangeBB <- efpFunctional(lim.process = "Brownian bridge",
  functional = list(time = function(x) diff(range(x)), comp = max),
  computePval = function(x, nproc = 1)
    pvalue.efp(x, lim.process = "Brownian bridge", functional = "range", alt.boundary = FALSE, k = nproc))

rangeBBI <- efpFunctional(lim.process = "Brownian bridge increments",
  functional = list(time = function(x) diff(range(x)), comp = max),
  computePval = function(x, nproc = 1)
    pvalue.efp(x, lim.process = "Brownian bridge increments", functional = "range", alt.boundary = FALSE, k = nproc))

rangeBM <- efpFunctional(lim.process = "Brownian motion",
  functional = list(time = function(x) diff(range(x)), comp = max),
  computePval = function(x, nproc = 1)
    pvalue.efp(x, lim.process = "Brownian motion", functional = "range", alt.boundary = FALSE, k = nproc))

rangeBMI <- efpFunctional(lim.process = "Brownian motion increments",
  functional = list(time = function(x) diff(range(x)), comp = max),
  computePval = function(x, nproc = 1)
    pvalue.efp(x, lim.process = "Brownian motion increments", functional = "range", alt.boundary = FALSE, k = nproc))


maxL2BB <- efpFunctional(lim.process = "Brownian bridge",
  functional = list(comp = function(x) sum(x^2), time = max),
  computePval = function(x, nproc = 1)
    pvalue.efp(x, lim.process = "Brownian bridge", functional = "maxL2", alt.boundary = FALSE, k = nproc))

meanL2BB <- efpFunctional(lim.process = "Brownian bridge",
  functional = list(comp = function(x) sum(x^2), time = mean),
  computePval = function(x, nproc = 1)
    pvalue.efp(x, lim.process = "Brownian bridge", functional = "meanL2", alt.boundary = FALSE, k = nproc))

supLM <- function(from = 0.15, to = NULL) {
  if(is.null(to)) to <- 1 - from
  stopifnot(from < 1, to < 1, from < to)
  lambda <- ((1-from)*to)/(from*(1-to))    
    
  wmax <- function(x) {
    n <- length(x)
    n1 <- floor(from * n)
    n2 <- floor(to * n)
    tt <- seq(along = x)/n
    x <- x[n1:n2]
    tt <- tt[n1:n2]
    x <- x/(tt * (1-tt))
    return(max(x))
  }

  computePval <- function(x, nproc = 1)
    pvalue.Fstats(x, type = "supF", k = nproc, lambda = lambda)
  
  computeCritval <- function(alpha, nproc = 1)
    uniroot(function(y) {computePval(y, nproc = nproc) - alpha}, c(0, 1000))$root

  boundary <- function(x) rep(1, length(x))

  plotProcess <- function(x, alpha = 0.05, aggregate = TRUE,
    xlab = NULL, ylab = NULL, main = x$type.name, ylim = NULL, ...)
  {
    n <- x$nobs
    n1 <- floor(from * n)
    n2 <- floor(to * n)
    tt <- (1:n)/n

    bound <- computeCritval(alpha = alpha, nproc = NCOL(x$process)) * boundary(0:n/n)
    if(is.null(xlab)) {
      if(!is.null(x$order.name)) xlab <- x$order.name
        else xlab <- "Time"
    }

    if(aggregate) {
      proc <- zoo(rowSums(as.matrix(x$process)^2), time(x))
      proc <- proc[-1]
      bound <- zoo(bound, time(x)[-1])
      tt <- zoo(tt, time(x)[-1])
      
      proc <- proc[n1:n2,]
      bound <- bound[n1:n2,]
      tt <- tt[n1:n2,]
      proc <- proc/(tt * (1-tt))
      
      if(is.null(ylab)) ylab <- "empirical fluctuation process"
      if(is.null(ylim)) ylim <- range(c(range(proc), range(bound)))
    
      plot(proc, xlab = xlab, ylab = ylab, main = main, ylim = ylim, ...)
      abline(0, 0)
      lines(bound, col = 2)	    
    } else {
      if(is.null(ylim) & NCOL(x$process) < 2) ylim <- range(c(range(x$process), range(bound), range(-bound)))
      if(is.null(ylab) & NCOL(x$process) < 2) ylab <- "empirical fluctuation process"

      panel <- function(x, ...)
      {
        lines(x, ...)
        abline(0, 0)
      }
      plot(x$process, xlab = xlab, ylab = ylab, main = main, panel = panel, ylim = ylim, ...)
    }
  }
  
  efpFunctional(lim.process = "Brownian bridge",
    functional = list(comp = function(x) sum(x^2), time = wmax),
    boundary = boundary,
    computePval = computePval,
    computeCritval = computeCritval,
    plotProcess = plotProcess)
}
