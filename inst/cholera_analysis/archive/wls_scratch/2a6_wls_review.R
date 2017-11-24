grab_psiFUN2 <- function(object, data, weights = 1, ...){
  
  X  <- stats::model.matrix(object$formula, data = data)
  
  Y  <- as.numeric(stats::model.frame(grab_response_formula(object), data = data)[[1]])
  n  <- length(Y)
  p  <- length(stats::coef(object))
  phi    <- as.numeric(summary(object)$dispersion[1])
  W      <- weights
  family <- object$family$family
  link   <- object$family$link
  invlnk <- object$family$linkinv
  family_link <- paste(family, link, sep = '_')
  
  stopifnot(length(W) == 1 | length(W) == n)
  if(length(W) == 1){
    W <- rep(W, n)
  }
  
  function(theta){
    lp <- X %*% theta # linear predictor
    f  <- as.numeric(invlnk(lp))  # fitted values
    r  <- Y - f       # residuals
    
    ### TODO: this is cludgy and needs to be reworked to be more general
    if(family_link == 'gaussian_identity'){
      D <- X
      V <- phi * diag(1, nrow = n, ncol = n)
    } else if(family_link %in% c('binomial_logit', 'quasibinomial_logit')){
      D <- apply(X, 2, function(x) x * exp(lp)/((1+exp(lp))^2) )
      V <- phi * diag(f * (1 - f), ncol = length(f) )/length(f)
    }
    
    t(D) %*% solve(V) %*% diag(W, nrow = n, ncol = n) %*% (r)
  }
}
