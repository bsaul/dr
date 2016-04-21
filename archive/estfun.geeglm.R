library(geepack)
library(magrittr)
data("ohio")

fit <- geeglm(resp ~ age + smoke + age:smoke, id=id, data=ohio,
      family=binomial, corstr="independence", scale.fix= F)


ee.gaussian <- function(x, w, r){
  x %*% diag(w) %*% r
}

ee.binomial <- function(x, lc, w, f, r, psi){
  D <- apply(x, 2, function(i) i * exp(lc)/((1+exp(lc))^2) )
  # This only works for independence working correlation matrices
  V <- psi * diag(f * (1 - f), ncol = length(f) )/length(f)
  D %*% V %*% diag(w) %*% r
}  

ee <- function(x, w, lc, f, r, psi, family){
  switch(family,
         gaussian = ee.gaussian(x, w, r),
         binomial = ee.binomial(x, w, lc, f, r, psi))
}

estfun.geeglm <- function(x, ...)
{
  xmat   <- model.matrix(x)
  r      <- x$residuals
  f      <- x$fitted.values
  w      <- x$weights
  lc     <- x$linear.predictor
  clust  <- x$id
  family <- x$family$family
  link   <- x$family$link
  psi    <- summary(x)$dispersion$Estimate
  
  if(x$corstr != 'independence'){
    stop("only independence working correlation is supported at this time")
  }
  
  m <- data.frame(f = f, r = r, w = w, lc = lc, xmat)
  out <- lapply(split(m, clust), FUN = function(x){
    x <- as.matrix(x)
    ee(x = x[, 5:ncol(x)], w = x[ ,'w'], f = x[ , 'f'],
       r = x[ , 'r'], lc = x[, 'lc'], psi = psi,
       family = family) %>% t() 
    }) %>% 
    unlist() %>%
    matrix(ncol = ncol(xmat), byrow = T, 
           dimnames = list(id = unique(clust), 
                           variable = dimnames(xmat)[[2]]))

  return(out)
}

estfun.geeglm(fit)
