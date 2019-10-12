recresid <- function(x, ...)
{
  UseMethod("recresid")
}

recresid.formula <- function(formula, data = list(), ...)
{
    mf <- model.frame(formula, data = data)
    y <- model.response(mf)
    modelterms <- terms(formula, data = data)
    X <- model.matrix(modelterms, data = data)
    rr <- recresid(X, y, ...)
    return(rr)
}

recresid.lm <- function(x, data = list(), ...)
{
    X <- if(is.matrix(x$x)) x$x else model.matrix(terms(x), model.frame(x))
    y <- if(is.vector(x$y)) x$y else model.response(model.frame(x))
    rr <- recresid(X, y, ...)
    return(rr)
}


## R version of recresid.
recresid_r <- function(x, y, start = ncol(x) + 1, end = nrow(x),
  tol = sqrt(.Machine$double.eps)/ncol(x), qr.tol = 1e-7, ...)
{
  n <- end
  q <- start - 1
  k <- ncol(x)
  rval <- rep(0, n - q)

  ## convenience function to replace NAs with 0s in coefs
  coef0 <- function(obj) {
    cf <- obj$coefficients
    ifelse(is.na(cf), 0, cf)
  }
  Xinv0 <- function(obj) {
    qr <- obj$qr
    rval <- matrix(0, ncol = k, nrow = k)
    wi <- qr$pivot[1:qr$rank]
    rval[wi,wi] <- chol2inv(qr$qr[1:qr$rank, 1:qr$rank, drop = FALSE])
    rval
  }

  ## initialize recursion
  y1 <- y[1:q]
  fm <- lm.fit(x[1:q, , drop = FALSE], y1, tol = qr.tol, ...)
  X1 <- Xinv0(fm)	 
  betar <- coef0(fm)
  xr <- as.vector(x[q+1,])
  fr <- as.vector((1 + (t(xr) %*% X1 %*% xr)))
  rval[1] <- (y[q+1] - t(xr) %*% betar)/sqrt(fr)

  ## check recursion agains full QR decomposition?
  check <- TRUE

  if((q+1) < n)
  {
    for(r in ((q+2):n))
    {
        ## check for NAs in coefficients
	nona <- all(!is.na(fm$coefficients))
    
	## recursion formula
        X1 <- X1 - (X1 %*% outer(xr, xr) %*% X1)/fr
        betar <- betar + X1 %*% xr * rval[r-q-1] * sqrt(fr)

	## full QR decomposition
	if(check) {
	  y1 <- y[1:(r-1)]
	  fm <- lm.fit(x[1:(r-1), , drop = FALSE], y1, tol = qr.tol, ...)
	  nona <- nona & all(!is.na(betar)) & all(!is.na(fm$coefficients))
	  ## keep checking?
	  if(nona && isTRUE(all.equal(as.vector(fm$coefficients), as.vector(betar), tol = tol))) check <- FALSE 
	  X1 <- Xinv0(fm)
          betar <- coef0(fm)
        }
        
        ## residual
        xr <- as.vector(x[r,])
	fr <- as.vector((1 + (t(xr) %*% X1 %*% xr)))
	rval[r-q] <- (y[r] - sum(xr * betar, na.rm = TRUE))/sqrt(fr)
    }
  }
  return(rval)
}


## C version of recresid.
recresid_c <- function(x, y, start = ncol(x) + 1, end = nrow(x),
  tol = sqrt(.Machine$double.eps)/ncol(x), ...) 
{
  n <- end
  q <- start - 1
  k <- ncol(x)
  rval <- rep(0, n - q)

  ## convenience function to replace NAs with 0s in coefs
  coef0 <- function(obj) {
    cf <- obj$coefficients
    ifelse(is.na(cf), 0, cf)
  }
  Xinv0 <- function(obj) {
    qr <- obj$qr
    rval <- matrix(0, ncol = k, nrow = k)
    wi <- qr$pivot[1:qr$rank]
    rval[wi, wi] <- chol2inv(qr$qr[1:qr$rank, 1:qr$rank, drop = FALSE])
    rval
  }

  ## initialize recursion
  y1 <- y[1:q]
  fm <- lm.fit(x[1:q, , drop = FALSE], y1)
  X1 <- Xinv0(fm)
  betar <- coef0(fm)
  xr <- as.vector(x[q + 1, ])
  fr <- as.vector((1 + (t(xr) %*% X1 %*% xr)))
  rval[1] <- (y[q + 1] - t(xr) %*% betar)/sqrt(fr)
  check <- TRUE

  ## fallback function to be called from C
  fallback <- function(r, fm, betar, check) {
    nona <- all(!is.na(fm$coefficients))
    y1 <- y[1:(r - 1)]
    fm <- lm.fit(x[1:(r - 1), , drop = FALSE], y1)
    nona <- nona & all(!is.na(betar)) & all(!is.na(fm$coefficients))
    if(nona && isTRUE(all.equal(as.vector(fm$coefficients), as.vector(betar), tol = tol))) 
      check <- FALSE
    return(list("fm" = fm, "X1" = Xinv0(fm), "betar" = coef0(fm), "check" = check))
  }

  rho <- new.env()

  ## run the C version of recresid
  if((q + 1) < n) {
    rval <- .Call("recresid", as.integer(q + 2), as.integer(n), X1, xr, fr,
      betar, rval, x, y, check, fallback, fm, rho, PACKAGE = "strucchange")
  }

  return(rval)
}


## Default wrapper for recresid.
recresid.default <- function(x, y, start = ncol(x) + 1, end = nrow(x),
  tol = sqrt(.Machine$double.eps)/ncol(x), qr.tol = 1e-7, engine = c("R", "C"), ...)
{
  ## checks and data dimensions
  stopifnot(start > ncol(x) & start <= nrow(x))
  stopifnot(end >= start & end <= nrow(x))

  engine <- match.arg(engine, c("R", "C"))
  ## FIXME: C++ available in strucchangeArmadillo on R-Forge

  recresid_fun <- switch(engine,
    "R" = recresid_r,
    "C" = recresid_c
    ## "C++" = strucchangeArmadillo::recresid_cpp
  )

  return(recresid_fun(x = x, y = y, start = start, end = end, tol = tol, qr.tol = qr.tol, ...))
}

