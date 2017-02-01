
#### IPW estimator ####

#------------------------------------------------------------------------------#
#' 
#' @export
#------------------------------------------------------------------------------#


weight_estimator <- function(A, X, lower = -Inf, upper = Inf)
{
  X <- as.matrix(X)
  f <- function(theta, alpha){
    vapply(alpha, function(x){
      w <- try(integrate(integrand, lower = lower, upper = upper,
                         theta = theta, alpha = x, response = A, xmatrix = X),
               silent = TRUE)
      if(is(w, 'try-error')){
        NA
      } else {
        1/w$value
      } }, numeric(1))
  }
  f
  # memoise::memoise(f)
}

#------------------------------------------------------------------------------#
#' 
#' @export
#------------------------------------------------------------------------------#


make_ipw_estimator <- function(Y, A, X_treatment, ...){
  w <- weight_estimator(A = A, X = X_treatment, ...)
  pi_t <- pi_term(A = A)
  
  function(theta, alpha, a = NULL){
    Ia <- if(is.null(a)) 1 else (A == a) * 1
    mean(Y * Ia) * w(theta, alpha)  / {if(!is.null(a)) dbinom(a, 1, alpha) else 1}
  }
}

#------------------------------------------------------------------------------#
#### Outcome Estimator ####
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#' 
#' @export
#------------------------------------------------------------------------------#

make_otc_estimator <- function(X_outcome, rhs_formula_outcome, inv_link, ...){
  
  X_expanded <- expand_outcome_frame(X_outcome, rhs_formula_outcome)
  
  function(theta, alpha, a = NULL){
    n1 <- nrow(X_outcome)
    n2 <- n1 - {if(is.null(a)) 0 else 1}
    a <- if(is.null(a)) 0:1 else a
    
    x <- X_expanded %>%
      filter_(~ A %in% a) %>%
      mutate_(pi_value = ~ dbinom(sum_a, n2, prob = alpha),
              fitted   = ~ as.numeric(inv_link(model.matrix(rhs_formula_outcome, data = .) %*% theta )),
              estimate = ~ fitted * pi_value) 
    
    sum(tapply(tapply(x$estimate, paste(x$A, x$ID), sum), 
               rep(a, each = n1), sum)) / n1
  }
}

#------------------------------------------------------------------------------#
#### Doubly Robust Estimator ####
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#' 
#' @export
#------------------------------------------------------------------------------#


make_dr_term1 <- function(X, inv_link){
  X <- as.matrix(X)
  f <- function(theta){
    inv_link(X %*% theta)
  }
  memoise::memoise(f)
}

#------------------------------------------------------------------------------#
#' 
#' @export
#------------------------------------------------------------------------------#

make_dbr_estimator <- function(Y, A, X_outcome, X_treatment, rhs_formula_outcome, inv_link){
  w <- weight_estimator2(A = A, X = X_treatment)
  dr_term1 <- make_dr_term1(X_outcome, inv_link = inv_link)
  dr_term2 <- make_otc_estimator(X_outcome, 
                                 rhs_formula_outcome = rhs_formula_outcome, 
                                 inv_link = inv_link)
  
  q_treatment <- ncol(X_treatment) + 1
  
  function(theta, alpha, a = NULL){
    Ia <- if(is.null(a)) 1 else (A == a) * 1
    Ybar <- mean(Ia * (Y - dr_term1(theta[(q_treatment + 1):length(theta)] ) ) )
    term1 <- Ybar * w(theta[1:q_treatment], alpha) / (dbinom(a, 1, alpha) * !is.null(a) * 1)
    term2 <- dr_term2(theta[(q_treatment + 1):length(theta)], alpha, a)
    term1 + term2
  }
}
