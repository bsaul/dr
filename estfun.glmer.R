estfun.glmer <- function(x, grad.method, ...)
{
  xmat   <- getME(x, 'X')
  resp   <- getME(x, 'y')
  parms  <- unlist(getME(x, c('beta', 'theta')))
  nparms <- length(parms)
  clust  <- getME(x, 'flist')
  family <- gfit@resp$family$family
  
  if(length(getME(x, 'theta')) > 1){
    stop('Not yet sure how to handle > 1 random effect')
  }
  
  if(family == 'binomial'){
    # Logistic-Normal Model
    integrand <- function(b, Y, X, parms){
      lc <- outer(X %*% parms[-nparms], b, '+')
      h  <- apply(lc, 3, function(x) dbinom(Y, 1, plogis(x) ) )
      hh <- apply(h, 2, prod)
      hh * dnorm(b, mean = 0, sd = parms[nparms])
    }
    
    objective.fun <- function(parms, Y, X){
      log(integrate(integrand, lower = -Inf, upper = Inf, parms = parms, 
                    Y = Y, X = X)$value)
    }
  } else {
    stop("binomial family is alls we got")
  }
  

  
  out <- by(cbind(resp, xmat), clust, simplify = F, FUN = function(x) {
    x <- as.matrix(x)
    numDeriv::grad(objective.fun, x = parms, 
                   Y = x[, 1], X = x[, -1], method = grad.method)
  })
  
  do.call('rbind', out)
}

