if(!("package:zoo" %in% search() || require(zoo))) warning("Could not load package zoo")

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

