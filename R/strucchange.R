efp <- function(formula, data=list(),
                type = c("Rec-CUSUM", "OLS-CUSUM", "Rec-MOSUM",
                "OLS-MOSUM", "fluctuation", "ME"), h = 0.15,
                dynamic = FALSE, rescale = TRUE, tol = 1e-7)
{
    mf <- model.frame(formula, data = data)
    y <- model.response(mf)
    modelterms <- terms(formula, data = data)
    X <- model.matrix(modelterms, data = data)
    n <- nrow(X)
    if(dynamic)
    {
      Xnames <- colnames(X)
      X <- cbind(y[1:(n-1)],X[2:n,])
      colnames(X) <- c("lag", Xnames)
      y <- y[-1]
      n <- n-1
    }
    k <- ncol(X)
    type <- match.arg(type)

    ## recursive residuals

    rec.res <- function(X, y, tol = 1e-7)
    {
        n <- nrow(X)
        q <- ncol(X)
        w <- rep(0,(n-q))
        for(r in ((q+1):n))
        {
            Xr1 <- X[1:(r-1),]
            xr <- as.vector(X[r,])
            X1 <- solve(t(Xr1)%*%Xr1, tol=tol)
            fr <- sqrt(1 + (t(xr) %*% X1 %*% xr))
            wr <- t(xr)%*% X1 %*%t(Xr1)%*% y[1:(r-1)]
            w[r-q] <- (y[r] - wr)/fr
        }
        return(w)
    }


    retval <- list(process = NULL,
                   type = type,
                   nreg = k,
                   nobs = n,
                   call = match.call(),
                   formula = formula,
                   par = NULL,
                   type.name = NULL,
                   coef = NULL,
                   Q12 = NULL,
                   datatsp = NULL,
		   rescale = rescale)

    switch(type,

           ## empirical process of Standard CUSUM model

           "Rec-CUSUM" = {
               w <- rec.res(X, y, tol = tol)
               sigma <- sqrt(var(w))
               process <- cumsum(c(0,w))/(sigma*sqrt(n-k))
               if(is.ts(data))
                   process <- ts(process, end = end(data),
                                 frequency = frequency(data))
               else
               {
               if(is.ts(y))
                   process <- ts(process, end = end(y),
                                 frequency = frequency(y))
               }
               retval$type.name <- "Standard CUSUM test"
           },

           ## empirical process of OLS-based CUSUM model

           "OLS-CUSUM" = {
               fm <- lm.fit(X,y)
               e <- fm$residuals
               sigma <- sqrt(sum(e^2)/fm$df.residual)
               process <- cumsum(c(0,e))/(sigma*sqrt(n))
               if(is.ts(data))
                   process <- ts(process, end = end(data),
                                 frequency = frequency(data))
               else
               {
               if(is.ts(y))
                   process <- ts(process, end = end(y),
                                 frequency = frequency(y))
               }
               retval$type.name <- "OLS-based CUSUM test"
           },

           ## empirical process of Recursive MOSUM model

           "Rec-MOSUM" = {
               w <- rec.res(X, y, tol = tol)
               nw <- n - k
               nh <- floor(nw*h)
               process <- rep(0, (nw-nh))
               for(i in 0:(nw-nh))
               {
                   process[i+1] <- sum(w[(i+1):(i+nh)])
               }
               sigma <- sqrt(var(w)*(nw-1)/(nw-k))
               process <- process/(sigma*sqrt(nw))
               if(is.ts(data))
                   process <- ts(process, end = time(data)[(n-floor(0.5 + nh/2))],
                                 frequency = frequency(data))
               else
               {
                 if(is.ts(y))
                   process <- ts(process,
                                 end = time(y)[(n-floor(0.5 + nh/2))],
                                 frequency = frequency(y))
                 else
                   process <- ts(process,
                                 end = (n-floor(0.5 + nh/2))/n,
                                 frequency=n)
               }
               retval$par <- h
               retval$type.name <- "Recursive MOSUM test"
           },

           ## empirical process of OLS-based MOSUM model

           "OLS-MOSUM" = {
               fm <- lm.fit(X,y)
               e <- fm$residuals
               sigma <- sqrt(sum(e^2)/fm$df.residual)
               nh <- floor(n*h)
               process <- rep(0, n-nh+1)
               for(i in 0:(n-nh))
               {
                   process[i+1] <- sum(e[(i+1):(i+nh)])
               }
               process <- process/(sigma*sqrt(n))
               if(is.ts(data))
                   process <- ts(process, end = time(data)[(n-floor(0.5 + nh/2))],
                                 frequency = frequency(data))
               else
               {
                 if(is.ts(y))
                   process <- ts(process,
                                 end = time(y)[(n-floor(0.5 + nh/2))],
                                 frequency = frequency(y))
                 else
                   process <- ts(process,
                                 end = (n-floor(0.5 + nh/2))/n,
                                 frequency=n)
               }
               retval$par <- h
               retval$type.name <- "OLS-based MOSUM test"
           },

           ## empirical process of recursive estimates fluctuation

           "fluctuation" = {
               m.fit <- lm.fit(X,y)
               beta.hat <- m.fit$coefficients
               sigma <- sqrt(sum(m.fit$residual^2)/m.fit$df.residual)
               process <- matrix(rep(0,k*(n-k+1)), nrow=k)
               Q12 <- root.matrix(crossprod(X))/sqrt(n)
               if(rescale)
               {
                   for(i in k:(n-1))
                   {
                       Qi12 <- root.matrix(crossprod(X[1:i,]))/sqrt(i)
                       process[,(i-k+1)] <- Qi12 %*%
                               (lm.fit(as.matrix(X[1:i,]), y[1:i])$coefficients - beta.hat)
                   }
               }
               else
               {
                   for(i in k:(n-1))
                   {
                       process[,(i-k+1)] <- Q12 %*% (lm.fit(as.matrix(X[1:i,]),
                                                           y[1:i])$coefficients - beta.hat)
                   }
               }
               process <- t(cbind(0, process))*matrix(rep((k-1):n,k),
                                                      ncol=k)/(sigma*sqrt(n))
               colnames(process) <- colnames(X)
               if(is.ts(data))
                   process <- ts(process, end = end(data),
                                 frequency = frequency(data))
               else
               {
                 if(is.ts(y))
                   process <- ts(process, end = end(y), frequency = frequency(y))
                 else
                   process <- ts(process, start = 0, frequency = nrow(process) - 1)
               }
               retval$Q12 <- Q12
               retval$type.name <- "Fluctuation test (recursive estimates test)"
           },

           ## empirical process of moving estimates fluctuation

           "ME" = {
               m.fit <- lm.fit(X,y)
               beta.hat <- m.fit$coefficients
               sigma <- sqrt(sum(m.fit$residual^2)/m.fit$df.residual)
               nh <- floor(n*h)
               process <- matrix(rep(0,k*(n-nh+1)), nrow=k)
               Q12 <- root.matrix(crossprod(X))/sqrt(n)
               if(rescale)
               {
                   for(i in 0:(n-nh))
                   {
                       Qnh12 <- root.matrix(crossprod(X[(i+1):(i+nh),]))/sqrt(nh)
                       process[, i+1] <-  Qnh12 %*% (lm.fit(
                                             as.matrix(X[(i+1):(i+nh),]), y[(i+1):(i+nh)])$coefficients - beta.hat)
                   }
               }
               else
               {
                   for(i in 0:(n-nh))
                   {
                       process[, i+1] <- Q12 %*% (lm.fit(as.matrix(X[(i+1):(i+nh),]),
		         y[(i+1):(i+nh)])$coefficients - beta.hat)
                   }
               }
               process <- nh*t(process)/(sqrt(n)*sigma)
               colnames(process) <- colnames(X)
               if(is.ts(data))
                   process <- ts(process, end = time(data)[(n-floor(0.5 + nh/2))],
                                 frequency = frequency(data))
               else
               {
                 if(is.ts(y))
                   process <- ts(process, end = time(y)[(n-floor(0.5 + nh/2))],
                                 frequency = frequency(y))
                 else
                   process <- ts(process, end = (n-floor(0.5 + nh/2))/n, frequency = n)
               }
               retval$par <- h
               retval$Q12 <- Q12
               retval$type.name <- "ME test (moving estimates test)"
           })

    if(!is.ts(process))
        process <- ts(process, start = 0, frequency = (length(process)-1))

    retval$process <- process

    if(is.ts(data))
        retval$datatsp <- tsp(data)
    else if(is.ts(y))
        retval$datatsp <- tsp(y)
    else
        retval$datatsp <- c(0, 1, n)

    m.fit <- lm.fit(X,y)
    retval$coefficients <- coefficients(m.fit)
    retval$sigma <-  sqrt(sum(m.fit$residual^2)/m.fit$df.residual)
    class(retval) <- c("efp")
    return(retval)
}


plot.efp <- function(x, alpha = 0.05, alt.boundary = FALSE, boundary = TRUE,
                     functional = "max", main = NULL,  ylim = NULL,
                     ylab = "empirical fluctuation process", ...)
{
    bound <- boundary(x, alpha = alpha, alt.boundary = alt.boundary)
    pos <- FALSE

    if(is.null(main)){
        if(alt.boundary & (x$type %in% c("Rec-CUSUM", "OLS-CUSUM"))){
            main <- paste(x$type.name, "with alternative boundaries")
        }
        else{
            main <- x$type.name
        }
    }

    z <- x$process
    switch(x$type,
           "fluctuation" = {
               if(!is.null(functional) && (functional == "max"))
               {
                   z <- ts(apply(abs(x$process), 1, 'max'),
                           start = start(x$process),
                           frequency = frequency(x$process))
                   pos <- TRUE
               }
           },

           "ME" = {
               if(!is.null(functional) && (functional == "max"))
               {
                   z <- ts(apply(abs(x$process), 1, 'max'),
                           start = start(x$process),
                           frequency = frequency(x$process))
                   pos <- TRUE
               }
           })

    ymax <- max(c(z, bound))
    if(pos)
        ymin <- 0
    else
        ymin <- min(c(z, -bound))
    if(is.null(ylim)) ylim <- c(ymin, ymax)
    if(boundary)
        panel <- function(y, ...)
        {
            lines(y, ...)
            lines(bound, col=2)
            lines(-bound, col=2)
            abline(0,0)
        }
    else
        panel <- function(y, ...)
        {
            lines(y, ...)
            abline(0,0)
        }
    if(any(attr(z, "class") == "mts"))
        plot(z, main = main, panel = panel, ...)
    else
    {
        plot(z, main = main, ylab = ylab, ylim = ylim, ...)
        if(boundary)
        {
            lines(bound, col=2)
            if(!pos) lines(-bound, col=2)
        }
        abline(0,0)
    }
}

pvalue.efp <- function(x, type, alt.boundary, functional = "max", h = NULL, k =
 NULL)
{
  type <- match.arg(type,
    c("Rec-CUSUM", "OLS-CUSUM", "Rec-MOSUM", "OLS-MOSUM", "fluctuation", "ME"))
  functional <- match.arg(functional, c("max", "range"))
  switch(type,

  "Rec-CUSUM" = {
  if(alt.boundary)
  {
      pval <- c(1, 0.997, 0.99, 0.975, 0.949, 0.912, 0.864, 0.806, 0.739, 0.666,
           0.589, 0.512, 0.437, 0.368, 0.307, 0.253, 0.205, 0.163, 0.129, 0.100,
           0.077, 0.058, 0.043, 0.032, 0.024, 0.018, 0.012, 0.009, 0.006, 0.004,
           0.003, 0.002, 0.001, 0.001, 0.001)
      critval <- (10:44)/10
      p <- approx(critval, pval, x, rule=2)$y
    }
    else
    {
      p <- ifelse(x < 0.3, 1 - 0.1465*x,
       2*(1-pnorm(3*x) +
 exp(-4*(x^2))*(pnorm(x)+pnorm(5*x)-1)-exp(-16*(x^2))*(1-pnorm(x))))
    }
  },

  "OLS-CUSUM" = {
    if(alt.boundary)
    {
      pval <- c(1, 1, 0.997, 0.99, 0.977, 0.954, 0.919, 0.871, 0.812, 0.743,
       0.666, 0.585, 0.504, 0.426, 0.353, 0.288, 0.230, 0.182, 0.142, 0.109,
       0.082, 0.062, 0.046, 0.034, 0.025, 0.017, 0.011, 0.008, 0.005, 0.004,
       0.003, 0.002, 0.001, 0.001, 0.0001)
      critval <- (12:46)/10
      p <- approx(critval, pval, x, rule=2)$y
    }
    else
    {
      p <- ifelse(x < 0.48, 1 - 0.1147*x,
             2*(exp(-2*x^2)-exp(-8*x^2)))
    }
  },

  "Rec-MOSUM" = {
    crit.table <- cbind(c(3.2165, 2.9795, 2.8289, 2.7099, 2.6061, 2.5111, 2.4283, 2.3464, 2.2686, 2.2255),
                        c(3.3185, 3.0894, 2.9479, 2.8303, 2.7325, 2.6418, 2.5609, 2.4840, 2.4083, 2.3668),
                        c(3.4554, 3.2368, 3.1028, 2.9874, 2.8985, 2.8134, 2.7327, 2.6605, 2.5899, 2.5505),
                        c(3.6622, 3.4681, 3.3382, 3.2351, 3.1531, 3.0728, 3.0043, 2.9333, 2.8743, 2.8334),
                        c(3.8632, 3.6707, 3.5598, 3.4604, 3.3845, 3.3102, 3.2461, 3.1823, 3.1229, 3.0737),
                        c(4.1009, 3.9397, 3.8143, 3.7337, 3.6626, 3.5907, 3.5333, 3.4895, 3.4123, 3.3912))
    tablen <- dim(crit.table)[2]
    tableh <- (1:10)*0.05
    tablep <- c(0.2, 0.15, 0.1, 0.05, 0.025, 0.01)
    tableipl <- numeric(tablen)
    for(i in (1:tablen)) tableipl[i] <- approx(tableh, crit.table[,i], h, rule = 2)$y
    p <- approx(c(0,tableipl), c(1,tablep), x, rule = 2)$y
  },

  "OLS-MOSUM" = {
    if(k>6) k <- 6
    crit.table <- get("sc.me")[(((k-1)*10+1):(k*10)), ]
    tablen <- dim(crit.table)[2]
    tableh <- (1:10)*0.05
    tablep <- c(0.1, 0.05, 0.025, 0.01)
    tableipl <- numeric(tablen)
    for(i in (1:tablen)) tableipl[i] <- approx(tableh, crit.table[,i], h, rule = 2)$y
    p <- approx(c(0,tableipl), c(1,tablep), x, rule = 2)$y
  },

  "fluctuation" = {
    switch(functional,
    "max" = {
      p <- ifelse(x<0.1, 1,
      {
        summand <- function(a,b)
        {
          exp(-2*(a^2)*(b^2))*(-1)^a
        }
        p <- 0
        for(i in 1:100) p <- p + summand(i,x)
        1-(1+2*p)^k
      })
    },
    "range" = {
      p <- ifelse(x<0.4,1,
      {
        p <- 0
        for(i in 1:10) p <- p + (4*i^2*x^2 - 1) * exp(-2*i^2*x^2)
        1-(1-2*p)^k
      })
    })
  },

  "ME" = {
    switch(functional,
    "max" = {
      if(k>6) k <- 6
      crit.table <- get("sc.me")[(((k-1)*10+1):(k*10)), ]
      tablen <- dim(crit.table)[2]
      tableh <- (1:10)*0.05
      tablep <- c(0.1, 0.05, 0.025, 0.01)
      tableipl <- numeric(tablen)
      for(i in (1:tablen)) tableipl[i] <- approx(tableh, crit.table[,i], h, rule = 2)$y
      p <- approx(c(0,tableipl), c(1,tablep), x, rule = 2)$y
    },
    "range" = {
      p <- ifelse(x<0.53,1,
      {
        p <- 0
        for(i in 1:10) p <- p + (-1)^(i-1) * i * pnorm(-i*x)
        1-(1-8*p)^k
      })
    })
  })
  return(p)
}


sctest.efp <- function(x, alt.boundary = FALSE, functional = c("max", "range"),
 ...)
{

    k <- x$nreg
    h <- x$par
    type <- x$type
    functional <- match.arg(functional)
    dname <- paste(deparse(substitute(x)))
    METHOD <- x$type.name
    x <- x$process

    switch(type,

           "Rec-CUSUM" = {
               j <- (1:(length(x)-1))/(length(x)-1)
               x <- x[-1]
               if(alt.boundary)
               {
                   STAT <- max(abs(x/sqrt(j)))
                   METHOD <- paste(METHOD, "with alternative boundaries")
               }
               else
               {
                   STAT <- max(abs(x/(1 + 2*j)))
               }
               names(STAT) <- "S"
           },

           "OLS-CUSUM" = {
               if(alt.boundary)
               {
                   j <- (1:(length(x)-2))/(length(x)-1)
                   x <- x[-1*c(1,length(x))]
                   STAT <- max(abs(x/sqrt(j*(1-j))))
                   METHOD <- paste(METHOD, "with alternative boundaries")
               }
               else
               {
                   STAT <- max(abs(x))
               }
               names(STAT) <- "S0"
           },

           "Rec-MOSUM" = {
               STAT <- max(abs(x))
               names(STAT) <- "M"
           },

           "OLS-MOSUM" = {
               STAT <- max(abs(x))
               names(STAT) <- "M0"
           },

           "fluctuation" = {
               switch(functional,
                      "max" = {
                          STAT <- max(abs(x))
                      },
                      "range" = {
                          myrange <- function(y) diff(range(y))
                          if(any(class(x)=="mts")) x <- x[-1,]
                          else x <- as.matrix(x[-1])
                          STAT <- max(apply(x,2,myrange))
                          METHOD <- paste(METHOD, "with range norm")
                      })
               names(STAT) <- "FL"
           },

           "ME" = {
               switch(functional,
                      "max" = {
                          STAT <- max(abs(x))
                      },
                      "range" = {
                          myrange <- function(y) diff(range(y))
                          STAT <- max(apply(x,2,myrange))
                          METHOD <- paste(METHOD, "with range norm")
                      })
               names(STAT) <- "ME"
           })
    PVAL <- pvalue.efp(STAT, type, alt.boundary,
                       functional = functional, h = h, k = k)
    RVAL <- list(statistic = STAT, p.value = PVAL,
                 method = METHOD, data.name = dname)
    class(RVAL) <- "htest"
    return(RVAL)
}

Fstats <- function(formula, from = 0.15, to = NULL, data = list(),
   cov.type = c("const","HC", "HC1"), tol=1e-7)
{
  mf <- model.frame(formula, data = data)
  y <- model.response(mf)
  modelterms <- terms(formula, data = data)
  X <- model.matrix(modelterms, data = data)
  k <- ncol(X)
  n <- length(y)
  e <- lm.fit(X,y)$residuals
  cov.type <- match.arg(cov.type)

  if(is.ts(data)){
      ytime <- time(data)
      yfreq <- frequency(data)
  }
  else if(is.ts(y)){
      ytime <- time(y)
      yfreq <- frequency(y)
  }

  ts.eps <- getOption("ts.eps")

  if(length(from)==2)
  {
    from <- which(abs(ytime-(from[1]+(from[2]-1)/yfreq)) < ts.eps)
    if(!is.null(to)) to <- which(abs(ytime-(to[1]+(to[2]-1)/yfreq)) < ts.eps)
  }
  else
  if(from < 1)
  {
    from <- floor(from*n)
    if(!is.null(to)) to <- floor(to*n)
  }
  if(is.null(to)) to <- n - from
  if(from < (k+1))
  {
    from <- k+1
    warning("'from' changed (was too small)")
  }
  if(to > (n-k-1))
  {
    to <- n-k-1
    warning("'to' changed (was too large)")
  }
  if(from <= to)
    point <- (from:to)
  else
    stop("inadmissable change points: 'from' is larger than 'to'")

  sume2 <- sum(e^2)
  lambda <- ((n-from)*to)/(from*(n-to))
  np <- length(point)
  stats <- rep(0,np)
  for(i in 1:np)
  {
    X1 <- as.matrix(X[(1:point[i]),])
    X2 <- as.matrix(X[((point[i]+1):n),])
    fm1 <- lm.fit(X1,y[1:point[i]])
    fm2 <- lm.fit(X2,y[((point[i]+1):n)])
    u <- c(fm1$residuals, fm2$residuals)
    sigma2 <- (sum(u^2))/(n-2*k)
    if(cov.type == "const") {
      stats[i] <- (sume2-sum(u^2))/sigma2}
    else {
      allX <- cbind(X1, matrix(rep(0, point[i]*k), ncol=k))
      allX <- rbind(allX, cbind(X2, X2))
      beta2 <- fm2$coefficients - fm1$coefficients
      Q1 <- solve(crossprod(allX), tol=tol)

      VX <- apply(allX, 2, function(x) x * e)
      V <- Q1 %*% t(VX) %*% VX %*% Q1
      if(cov.type == "HC1") {V <- V * (n/(n-k))}

      stats[i] <- as.vector(t(beta2) %*% solve(V[-(1:k),-(1:k)],
                    tol=tol) %*% beta2)
     }
  }
  if(is.ts(data)){
      stats <- ts(stats, start = time(data)[from], frequency = frequency(data))
      datatsp <- tsp(data)
  }
  else if(is.ts(y)){
      stats <- ts(stats, start = time(y)[from], frequency = frequency(y))
      datatsp <- tsp(y)
  }
  else{
      stats <- ts(stats, start = from/n, frequency = n)
      datatsp <- c(0, 1, n)
  }

  retval <- list(Fstats = stats,
                 nreg = k,
                 nobs = n,
                 par = lambda,
                 call = match.call(),
                 formula = formula,
                 datatsp = datatsp)

  class(retval) <- "Fstats"
  return(retval)
}

print.efp <- function(x, ...)
{
    cat("\nEmpirical Fluctuation Process:", x$type.name, "\n\n")
    cat("Call: ")
    print(x$call)
    cat("\n")
}

print.Fstats <- function(x, ...)
{
    cat("\nF statistics \n\n")
    cat("Call: ")
    print(x$call)
    cat("\n")
}

plot.Fstats <- function(x, pval = FALSE, asymptotic = FALSE,
                        alpha = 0.05, boundary = TRUE, aveF = FALSE,
                        xlab = "Time", ylab = NULL,
                        ylim = NULL, ...)
{
    k <- x$nreg
    n <- x$nobs
    bound <- boundary(x, alpha = alpha, pval = pval, aveF = aveF, asymptotic =
                      asymptotic)
    x <- x$Fstats

    if(pval)
    {
        if(asymptotic)
            x <- 1 - pchisq(x, k)
        else
            x <- 1 - pf(x, k, (n-2*k))
        if(is.null(ylab)) ylab <- "p values"
    }
    else
        if(is.null(ylab)) ylab <- "F statistics"
    if(is.null(ylim)) ylim <- c(0, max(c(x,bound)))
    plot(x, xlab = xlab, ylab = ylab, ylim = ylim, ...)
    abline(0,0)
    if(boundary)
    {
        lines(bound, col=2)
        if(aveF) lines(ts(rep(mean(x),length(x)),start=start(x),
                       frequency = frequency(x)),lty=2)
    }
}

sctest.Fstats <- function(x, type = c("supF", "aveF", "expF"), asymptotic = FALSE, ...)
{
    dname <- paste(deparse(substitute(x)))
    type <- match.arg(type)
    switch(type,
           supF = {
               STATISTIC <- max(x$Fstats)
               names(STATISTIC) <- "sup.F"
               METHOD <- "supF test"
           },
           aveF = {
               STATISTIC <- mean(x$Fstats)
               names(STATISTIC) <- "ave.F"
               METHOD <- "aveF test"
           },
           expF = {
               STATISTIC <- log(mean(exp(0.5*x$Fstats)))
               names(STATISTIC) <- "exp.F"
               METHOD <- "expF test"
           })
    if((x$par == 1) & !(type == "expF") & !asymptotic)
    {
        METHOD <- "Chow test"
        PVAL <- 1 - pf(STATISTIC, k, (n-2*k))
    }
    else
        PVAL <- pvalue.Fstats(STATISTIC, type = type,
                              k=x$nreg, lambda=x$par)

    RVAL <- list(statistic = STATISTIC, p.value = PVAL,
                 method = METHOD, data.name = dname)
    class(RVAL) <- "htest"
    return(RVAL)
}

sctest <- function(x, ...)
{
    UseMethod("sctest")
}

sctest.formula <- function(formula, type = c("Rec-CUSUM", "OLS-CUSUM",
  "Rec-MOSUM", "OLS-MOSUM", "fluctuation", "ME", "Chow", "supF", "aveF",
  "expF"), h = 0.15, alt.boundary = FALSE, functional = c("max", "range"),
  from = 0.15, to = NULL, point = 0.5, asymptotic = FALSE, data = list(), ...)
{
  type <- match.arg(type)
  dname <- paste(deparse(substitute(formula)))
  if(type == "Chow")
  {
    mf <- model.frame(formula, data = data)
    y <- model.response(mf)
    modelterms <- terms(formula, data = data)
    X <- model.matrix(modelterms, data = data)
    METHOD <- "Chow test"
    k <- ncol(X)
    n <- length(y)

    if(is.ts(data)){
        ytime <- time(data)
        yfreq <- frequency(data)
    }
    else if(is.ts(y)){
        ytime <- time(y)
        yfreq <- frequency(y)
    }

    ts.eps <- getOption("ts.eps")

    if(length(point)==2) {
      point <- which(abs(ytime-(point[1]+(point[2]-1)/yfreq)) < ts.eps)
    }
    else {
      if(point < 1)
      {
        point <- floor(point*n)
      }
    }

    if(!((point>k) & (point<(n-k)))) stop("inadmissable change point")
    e <- lm.fit(X,y)$residuals
    u <- c(lm.fit(as.matrix(X[(1:point),]),y[1:point])$residuals,
      lm.fit(as.matrix(X[((point+1):n),]),y[((point+1):n)])$residuals)
    STATISTIC <- ((sum(e^2)-sum(u^2))/k)/((sum(u^2))/(n-2*k))
    names(STATISTIC) <- "F"
    if(asymptotic)
    {
      STATISTIC <- STATISTIC * k
      PVAL <- 1 - pchisq(STATISTIC, k)
    }
    else
      PVAL <- 1 - pf(STATISTIC, k, (n-2*k))
    RVAL <- list(statistic = STATISTIC, p.value = PVAL, method = METHOD, data.name = "formula")
    class(RVAL) <- "htest"
  }
  else if(any(type == c("Rec-CUSUM", "OLS-CUSUM",
  "Rec-MOSUM", "OLS-MOSUM", "fluctuation", "ME")))
  {
    process <- efp(formula, type = type, h = h, data = data, ...)
    RVAL <- sctest(process, alt.boundary = alt.boundary, functional = functional)
  }
  else if(any(type == c("supF", "aveF", "expF")))
  {
    fs <- Fstats(formula, from = from, to = to, data=data, ...)
    RVAL <- sctest(fs, type = type)
  }
  RVAL$data.name <- dname
  return(RVAL)
}

boundary <- function(x, ...)
{
    UseMethod("boundary")
}

boundary.efp <- function(x, alpha = 0.05, alt.boundary = FALSE, ...)
{
    pos <- FALSE
    k <- x$nreg
    h <- x$par

    bound <- uniroot(function(y) {pvalue.efp(y, type=x$type,
        alt.boundary=alt.boundary, h=h, k=k) - alpha}, c(0,20))$root
    switch(x$type,
           "Rec-CUSUM" = {
               if(alt.boundary)
                   bound <- sqrt((0:(length(x$process)-1))/(length(x$process)-1))*bound
               else
                   bound <- bound + (2*bound*(0:(length(x$process)-1))/(length(x$process)-1))
           },
           "OLS-CUSUM" = {
               if(alt.boundary)
               {
                   j <- (0:(length(x$process)-1))/(length(x$process)-1)
                   bound <- sqrt(j*(1-j))*bound
               }
               else
                   bound <- rep(bound,length(x$process))
           },
           "Rec-MOSUM" = { bound <- rep(bound, length(x$process))},
           "OLS-MOSUM" = { bound <- rep(bound, length(x$process))},
           "fluctuation" = { bound <- rep(bound, length(x$process))},
           "ME" = { bound <- rep(bound, length(x$process))})

    bound <- ts(bound, end = end(x$process), frequency = frequency(x$process))
    return(bound)
}

boundary.Fstats <- function(x, alpha = 0.05, pval = FALSE, aveF =
                            FALSE, asymptotic = FALSE, ...)
{
    if(aveF)
    {
      myfun <-  function(y) {pvalue.Fstats(y, type="ave", x$nreg, x$par) - alpha}
      upper <- 20
    }
    else
    {
      myfun <-  function(y) {pvalue.Fstats(y, type="sup", x$nreg, x$par) - alpha}
      upper <- 40
    }
    bound <- uniroot(myfun, c(0,upper))$root
    if(pval)
    {
        if(asymptotic)
            bound <- 1 - pchisq(bound, x$nreg)
        else
            bound <- 1 - pf(bound, x$nreg, (x$nobs-2*x$nreg))
    }
    bound <- ts(bound,
                start = start(x$Fstats),
                end = end(x$Fstats),
                frequency = frequency(x$Fstats))
    return(bound)
}

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

lines.efp <- function(x, ...)
{
    proc <- x$process
    if(any(class(proc)=="mts"))
    {
        proc <- ts(apply(abs(x$process), 1, 'max'),
                   start = start(x$process), frequency = frequency(x$process))
    }
    lines(proc, ...)
}

lines.Fstats <- function(x, ...)
{
    lines(x$Fstats, ...)
}

covHC <- function(formula, type=c("HC2", "const", "HC", "HC1", "HC3"),
 data=list())
{
  mf <- model.frame(formula, data = data)
  y <- model.response(mf)
  X <- model.matrix(formula, data = data)
  n <- nrow(X)
  k <- ncol(X)
  Q1 <- solve(crossprod(X))
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
