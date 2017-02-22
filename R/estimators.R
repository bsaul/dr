#------------------------------------------------------------------------------#
#' Makes IPW, OTC, and DBR estimators for group-level data
#' @export
#------------------------------------------------------------------------------#

dr_estimators <- function(data, t_model, o_model, randomization){
  
  ## component data
  Y      <- geex::get_response(formula(o_model), data = data)
  A      <- geex::get_response(A ~ 1, data = data)
  X_t    <- geex::get_design_matrix(geex::get_fixed_formula(t_model), data = data)
  X_o    <- as.data.frame(geex::get_design_matrix(geex::get_fixed_formula(o_model), data = data))
  X_o_ex <- expand_outcome_frame(X_o, geex::get_fixed_formula(o_model))
  N      <- nrow(X_o)  
  rhs_o  <- geex::get_fixed_formula(o_model)
  inv_link_o <- family(o_model)$linkinv
  
  ## components for IPW estimator
  ip_fun    <- weight_estimator(A = A, X = X_t, randomization = randomization)

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
    ipw_ce0 <- mean(Y * (A == 0)) * ipw / (1 - alpha)
    ipw_ce1 <- mean(Y * (A == 1)) * ipw / alpha
    ipw_ce  <- mean(Y) * ipw 
    
    ### OTC estimator ###
    # compute fitted value for expanded data.frame
    mu <- as.numeric(
      inv_link_o(model.matrix(rhs_o, data = X_o_ex) %*% theta[index_o] ))
    
    pi_term_a <- vapply(alpha, function(x) {dbinom(X_o_ex$sum_a, N - 1, x)}, 
                        numeric(nrow(X_o_ex)))
    pi_term   <- vapply(alpha, function(x) {dbinom(X_o_ex$sum_a, N    , x)}, 
                        numeric(nrow(X_o_ex)))
    
    estimates_a <- apply(pi_term_a, 2, function(col) col * mu)
    estimates   <- apply(pi_term,   2, function(col) col * mu)
    
    part1_a <- apply(estimates_a, 2, function(col) { 
      tapply(col, paste(X_o_ex$A, X_o_ex$ID), sum) 
    }) 
    
    part1   <- apply(estimates, 2, function(col) { 
      x <- tapply(col, paste(X_o_ex$ID, X_o_ex$A), sum) 
      tapply(x, rep(unique(X_o_ex$ID), each = 2), mean)
    }) 
    
    otc_ce_a <- apply(part1_a , 2, function(col) {
      tapply(col, rep(0:1, each = N), sum)} )
    
    otc_ce0 <- otc_ce_a[1, ]/N
    otc_ce1 <- otc_ce_a[2, ]/N
    otc_ce  <- apply(part1, 2, sum)/N
    
    ## DBR estimator ###
    fY <- dr_term1_fun(theta[index_o])
    Ybar0 <- sum((A == 0) * (Y - fY) )
    Ybar1 <- sum((A == 1) * (Y - fY) )
    Ybar  <- sum(Y - fY)
    
    term1_0 <- Ybar0 * ipw / (1 - alpha)
    term1_1 <- Ybar1 * ipw / alpha
    term1   <- Ybar * ipw 
    
    dbr_ce0 <- (term1_0 + otc_ce0*N)/N
    dbr_ce1 <- (term1_1 + otc_ce1*N)/N
    dbr_ce  <- (term1 + otc_ce*N)/N
    
    ## Hajek corrections ###
    
    ## IPW
    Nhat0 <- sum(A == 0) * ipw / (1 - alpha)
    Nhat1 <- sum(A == 1) * ipw / alpha
    
    hjk_ipw_ce0 <- (sum(Y * (A == 0)) * ipw / (1 - alpha )) / Nhat0
    hjk_ipw_ce1 <- (sum(Y * (A == 1)) * ipw / alpha) / Nhat1
    hjk_ipw_ce0 <- ifelse(is.nan(hjk_ipw_ce0) | is.infinite(hjk_ipw_ce0), 0, hjk_ipw_ce0)
    hjk_ipw_ce1 <- ifelse(is.nan(hjk_ipw_ce1) | is.infinite(hjk_ipw_ce1), 0, hjk_ipw_ce1)
    
    ## DBR
    hjk_dbr_ce0 <- (term1_0 + otc_ce0*N)/Nhat0
    hjk_dbr_ce1 <- (term1_1 + otc_ce1*N)/Nhat1
    hjk_dbr_ce0 <- ifelse(is.nan(hjk_dbr_ce0) | is.infinite(hjk_dbr_ce0), 0, hjk_dbr_ce0)
    hjk_dbr_ce1 <- ifelse(is.nan(hjk_dbr_ce1) | is.infinite(hjk_dbr_ce1), 0, hjk_dbr_ce1)
    
    c(ipw_ce0 = ipw_ce0,
      ipw_ce1 = ipw_ce1,
      ipw_ce  = ipw_ce,
      otc_ce0 = otc_ce0 ,
      otc_ce1 = otc_ce1,
      otc_ce  = otc_ce,
      dbr_ce0 = dbr_ce0,
      dbr_ce1 = dbr_ce1,
      dbr_ce  = dbr_ce,
      hjk_ipw_ce0 = hjk_ipw_ce0,
      hjk_ipw_ce1 = hjk_ipw_ce1,
      hjk_ipw_ce  = NA,
      hjk_dbr_ce0 = hjk_dbr_ce0, 
      hjk_dbr_ce1 = hjk_dbr_ce1,
      hjk_dbr_ce  = NA)
  }
}


#------------------------------------------------------------------------------#
#' IP weight estimator
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
               silent = FALSE)
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
#' Makes the first term in the DR estimator
#' @export
#------------------------------------------------------------------------------#

make_dr_term1 <- function(X, inv_link){
  X <- as.matrix(X)
  f <- function(theta){
    inv_link(X %*% theta)
  }
  memoise::memoise(f)
}