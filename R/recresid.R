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

recresid.default <- function(x, y, ...)
{
    n <- nrow(x)
    q <- ncol(x)
    w <- rep(0,(n-q))
    Xr1 <- x[1:q,,drop = FALSE]
    xr <- as.vector(x[q+1,])
    X1 <- solveCrossprod(Xr1)
    fr <- as.vector((1 + (t(xr) %*% X1 %*% xr)))
    betar <- X1 %*%t(Xr1)%*% y[1:q]
    w[1] <- (y[q+1] - t(xr) %*% betar)/sqrt(fr)

    if((q+1) < n)
    {
      for(r in ((q+2):n))
      {
          X1 <- X1 - (X1 %*% outer(xr, xr) %*% X1)/fr
	  betar <- betar + X1 %*% xr * w[r-q-1]*sqrt(fr)
	  xr <- as.vector(x[r,])
          fr <- as.vector((1 + (t(xr) %*% X1 %*% xr)))
          w[r-q] <- (y[r] - t(xr) %*% betar)/sqrt(fr)
      }
    }
    return(w)
}

