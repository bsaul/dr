#------------------------------------------------------------------------------#
#' Makes IPW, OTC, and DBR estimators for group-level data
#' @export
#------------------------------------------------------------------------------#

dr_estimators <- function(data, t_model, o_model){
  
  Y      <- geex::get_response(formula(o_model), data = data)
  A      <- geex::get_response(A ~ 1, data = data)
  X_t    <- geex::get_design_matrix(geex::get_fixed_formula(t_model), data = data)
  X_o    <- as.data.frame(geex::get_design_matrix(geex::get_fixed_formula(o_model), data = data))
  X_o_ex <- expand_outcome_frame(X_o, geex::get_fixed_formula(o_model))
  inv_link_o <- family(o_model)$linkinv
  rhs_o <- geex::get_fixed_formula(o_model)
  n_ <- nrow(X_o)  
  
  ## components for IPW estimator
  ip_fun    <- weight_estimator(A = A, X = X_t, randomization = 2/3)
  
  ## components for DBR estimator
  dr_term1_fun <- make_dr_term1(X_o, inv_link = inv_link_o)
  
  ## indices
  p_t <- ncol(X_t) + 1
  p_o <- ncol(X_o)
  p   <- p_t + p_o
  index_t <- 1:p_t
  index_o <- (p_t + 1):(p_t + p_o)
  
  function(theta, alpha){

    ### IPW estimator ###
    ipw <- ip_fun(theta = theta[index_t], alpha = alpha)
    ipw_ce0 <- mean(Y * (A == 0)) * ipw / dbinom(0, 1, alpha)
    ipw_ce1 <- mean(Y * (A == 1)) * ipw / dbinom(1, 1, alpha)
    ipw_ce  <- mean(Y) * ipw 
    
    ### OTC estimator ###
    fitted <- as.numeric(
      inv_link_o(model.matrix(rhs_o, data = X_o_ex) %*% theta[index_o] ))
    
    pi_term_a <- vapply(alpha, function(x) {dbinom(X_o_ex$sum_a, n_ - 1, x)}, 
                        numeric(nrow(X_o_ex)))
    pi_term   <- vapply(alpha, function(x) {dbinom(X_o_ex$sum_a, n_    , x)}, 
                        numeric(nrow(X_o_ex)))
    
    estimates_a <- apply(pi_term_a, 2, function(col) col * fitted)
    estimates   <- apply(pi_term,   2, function(col) col * fitted)
    
    part1_a <- apply(estimates_a, 2, function(col) { 
      tapply(col, paste(X_o_ex$A, X_o_ex$ID), sum) 
    }) 
    
    part1   <- apply(estimates, 2, function(col) { 
      x <- tapply(col, paste(X_o_ex$ID, X_o_ex$A), sum) 
      tapply(x, rep(unique(X_o_ex$ID), each = 2), mean)
    }) 
    
    otc_ce_a <- apply(part1_a , 2, function(col) {
      tapply(col, rep(0:1, each = n_), sum)} ) / n_
    
    otc_ce0 <- otc_ce_a[1, ]
    otc_ce1 <- otc_ce_a[2, ]
    otc_ce  <- apply(part1, 2, mean)
    
    ## DBR estimator ###
    fY <- dr_term1_fun(theta[index_o])
    Ybar0 <- mean((A == 0) * (Y - fY ) )
    Ybar1 <- mean((A == 1) * (Y - fY ) )
    Ybar  <- mean(Y - fY )
    
    term1_0 <- Ybar0 * ipw / dbinom(0, 1, alpha)
    term1_1 <- Ybar1 * ipw / dbinom(1, 1, alpha)
    term1   <- Ybar * ipw 
    
    dbr_ce0 <- term1_0 + otc_ce0
    dbr_ce1 <- term1_1 + otc_ce1
    dbr_ce  <- term1 + otc_ce
    
    c(ipw_ce0 = ipw_ce0,
      ipw_ce1 = ipw_ce1,
      ipw_ce  = ipw_ce,
      otc_ce0 = otc_ce0 ,
      otc_ce1 = otc_ce1,
      otc_ce  = otc_ce,
      dbr_ce0 = dbr_ce0,
      dbr_ce1 = dbr_ce1,
      dbr_ce  = dbr_ce)
  }
}
#### IPW estimator ####

#------------------------------------------------------------------------------#
#' 
#' @export
#------------------------------------------------------------------------------#


weight_estimator <- function(A, X, lower = -Inf, upper = Inf, randomization = 1)
{
  X <- as.matrix(X)
  f <- function(theta, alpha){
    vapply(alpha, function(x){
      w <- try(integrate(integrand, lower = lower, upper = upper,
                         theta = theta, alpha = x, response = A, xmatrix = X, 
                         randomization = randomization),
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
