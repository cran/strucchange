root.matrix <- function(X)
{
    if((ncol(X)==1)&&(nrow(X)==1)) return(sqrt(X))
    else
    {
        X.eigen <- eigen(X, symmetric=TRUE)
        if(any(X.eigen$values < 0)) stop("matrix is not positive semidefinite")
        sqomega <- sqrt(diag(X.eigen$values))
        V <- X.eigen$vectors
        return(V%*%sqomega%*%t(V))
    }
}

solveCrossprod <- function(X, method = c("qr", "chol")) {
  if(match.arg(method) == "qr") chol2inv(qr.R(qr(X)))
  else chol2inv(chol(crossprod(X)))
}

covHC <- function(formula, type = c("HC2", "const", "HC", "HC1", "HC3"), data = list())
{
  mf <- model.frame(formula, data = data)
  y <- model.response(mf)
  X <- model.matrix(formula, data = data)
  n <- nrow(X)
  k <- ncol(X)
  Q1 <- solveCrossprod(X)
  res <- lm.fit(X,y)$residuals
  sigma2 <- var(res)*(n-1)/(n-k)
  type <- match.arg(type)

  if( type == "const") {
    V <- sigma2 * Q1 }
  else
  {
    if(type == "HC2")
    {
      diaghat <- 1 - diag(X %*% Q1 %*% t(X))
      res <- res/sqrt(diaghat)
    }
    if(type == "HC3")
    {
      diaghat <- 1 - diag(X %*% Q1 %*% t(X))
      res <- res/diaghat
      Xu <- as.vector(t(X) %*% res)
    }
    VX <- apply(X, 2, function(x) x * res)
    if(type %in% c("HC", "HC1", "HC2")) {V <- Q1 %*% t(VX) %*% VX %*% Q1}
    if(type == "HC1") {V <- V * (n/(n-k))}
    if(type == "HC3") {V <- Q1 %*% (t(VX) %*% VX - (outer(Xu,Xu) /n)) %*% Q1 * (n-1)/n}
  }
  return(V)
}
