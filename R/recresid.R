recresid <- function(x, ..., tol = 1e-7)
{
  UseMethod("recresid")
}

recresid.formula <- function(formula, data=list(), ..., tol = 1e-7)
{
    mf <- model.frame(formula, data = data)
    y <- model.response(mf)
    modelterms <- terms(formula, data = data)
    X <- model.matrix(modelterms, data = data)
    rr <- recresid(X, y, tol = tol)
    return(rr)
}

recresid.lm <- function(x, data=list(), ..., tol = 1e-7)
{
    mf <- model.frame(x, data = data)
    y <- model.response(mf)
    modelterms <- terms(x, data = data)
    X <- model.matrix(modelterms, data = data)
    rr <- recresid(X, y, tol = tol)
    return(rr)
}

recresid.default <- function(x, y, ..., tol = 1e-7)
{
    n <- nrow(x)
    q <- ncol(x)
    w <- rep(0,(n-q))
    Xr1 <- x[1:q,,drop = FALSE]
    xr <- as.vector(x[q+1,])
    X1 <- solve(t(Xr1)%*%Xr1, tol=tol)
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

