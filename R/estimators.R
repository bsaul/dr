#------------------------------------------------------------------------------#
#' Makes IPW estimator for group-level data
#' @export
#------------------------------------------------------------------------------#

ipw_estimator <- function(data, models, randomization, hajek = FALSE, ...){
  
  ## component data
  comp <- extract_model_info(model = models, data = data, 'ipw')
  Y <- comp$Y
  A <- comp$A
  N <- comp$N
  ## components for IPW estimator
  ip_fun <- weight_estimator(
    A = comp$A, 
    X = comp$X_t, 
    randomization = randomization)

  index_t <- 1:comp$p

  function(theta, alpha){
    
    ### IPW estimator ###
    ipw <- ip_fun(theta = theta[index_t], alpha = alpha)
    ipw_0 <- sum(Y * (A == 0)) * ipw / (1 - alpha)
    ipw_1 <- sum(Y * (A == 1)) * ipw / alpha
    # ipw_ce  <- mean(Y) * ipw 
    
    if(hajek){
      Nhat0 <- sum(A == 0) * ipw / (1 - alpha)
      Nhat1 <- sum(A == 1) * ipw / alpha
      out <- c(ipw_0/Nhat0, ipw_1/Nhat1)
      names(out) <- paste0(c('ipw_hjk_Y0_', 'ipw_hjk_Y1_' ), alpha)
    } else {
      out <- c(ipw_0/N, ipw_1/N)
      names(out) <- paste0(c('ipw_Y0_', 'ipw_Y1_' ), alpha)
    }
    out
  }
}



#------------------------------------------------------------------------------#
#' Makes OTC (outcome based) estimator for group-level data
#' @export
#------------------------------------------------------------------------------#

otc_estimator <- function(data, models, randomization, ...){
  
  ## component data
  comp <- extract_model_info(model = models, data = data, 'otc')
  MM   <- model.matrix(comp$rhs_o, data = comp$X_o_ex) 
  N    <- comp$N
  ## indices
  index_o <- 1:comp$p
  function(theta, alpha){

    ### OTC estimator ###
    # compute fitted value for expanded data.frame
    mu <- as.numeric(comp$inv_link_o(MM %*% theta[index_o]))
    
    # compute pi term per number treated in group per subject
    pi_term_a <- vapply(alpha, function(x) { # vapply so it can work over a vector of alphas
      dbinom(comp$X_o_ex$sum_a, comp$N - 1, x)}, 
      numeric(nrow(comp$X_o_ex)))
    
    # mulitply mu_ij by the pi term rowwise
    # apply over columns so that it can work over a vector of alphas
    piXmu_a <- apply(pi_term_a, 2, function(col) col * mu)
    
    # sum within levels of A (0:1) WITHIN subjects
    piXmu_a <- apply(piXmu_a, 2, function(col) { 
      tapply(col, paste(comp$X_o_ex$A, comp$X_o_ex$ID), sum) 
    }) 
    
    # sum within levels of A (0:1) ACROSS subjects
    otc_ce_a <- apply(piXmu_a , 2, function(col) {
      tapply(col, rep(0:1, each = N), sum)} )
    
    otc_ce0 <- otc_ce_a[1, ]/N
    otc_ce1 <- otc_ce_a[2, ]/N
    
    ## 'Overall' marginal mean...
    # pi_term   <- vapply(alpha, function(x) {dbinom(comp$X_o_ex$sum_a, N    , x)}, 
    #                     numeric(nrow(X_o_ex)))
    # piXmu   <- apply(pi_term,   2, function(col) col * mu)
    # part1   <- apply(estimates, 2, function(col) { 
    #   x <- tapply(col, paste(X_o_ex$ID, X_o_ex$A), sum) 
    #   tapply(x, rep(unique(X_o_ex$ID), each = 2), mean)
    # }) 
    # 
    # otc_ce  <- apply(part1, 2, sum)/N
    
    x <- c(otc_ce0, otc_ce1)
    names(x) <- paste0(c('otc_Y0_', 'otc_Y1_'), alpha)
    x
  }
}
#------------------------------------------------------------------------------#
#' Makes DBR estimators for group-level data
#' @export
#------------------------------------------------------------------------------#

dbr_estimator <- function(data, models, randomization, ...){
  
  ## component data
  comp <- extract_model_info(model = models, data = data, 'dbr')
  Y <- comp$Y
  A <- comp$A
  N <- comp$N
  ## components for IPW part
  ip_fun <- weight_estimator(
    A = comp$A, 
    X = comp$X_t, 
    randomization = randomization)
  ## components for OTC part
  dr_term1_fun <- make_dr_term1(
    comp$X_o, 
    inv_link = comp$inv_link_o)

  dr_term2_fun <- otc_estimator(data, models, randomization)

  ## indices
  p_t <- comp$p_t
  p_o <- comp$p_o
  p   <- p_t + p_o
  index_t <- 1:p_t
  index_o <- (p_t + 1):(p_t + p_o)
  
  function(theta, alpha){

    fY      <- dr_term1_fun(theta[index_o])
    ipw     <- ip_fun(theta[index_t], alpha)
    Ybar0   <- sum((A == 0) * (Y - fY) )
    Ybar1   <- sum((A == 1) * (Y - fY) )
    term1_0 <- Ybar0 * ipw / (1 - alpha)
    term1_1 <- Ybar1 * ipw / alpha
    term2   <- dr_term2_fun(theta[index_o], alpha)
    
    dbr_ce0 <- (term1_0 + term2[1]*N)/N
    dbr_ce1 <- (term1_1 + term2[2]*N)/N
    
    ## Overall marginal means... ###
    # term1   <- Ybar * ipw 
    # Ybar  <- sum(Y - fY)

    x <- c(dbr_ce0, dbr_ce1)
    names(x) <- paste0(c('dbr_Y0_', 'dbr_Y1_'), alpha)
    x
  }
}

#------------------------------------------------------------------------------#
#' IP weight estimator
#' @export
#------------------------------------------------------------------------------#
holding <- function(){
  
  ## Hajek corrections
  Nhat0 <- sum(A == 0) * ipw / (1 - alpha)
  Nhat1 <- sum(A == 1) * ipw / alpha
  dbr_hjk_ce0 <- (term1_0 + term2[1]*N)/Nhat0
  dbr_hjk_ce1 <- (term1_1 + term2[2]*N)/Nhat1
  ## Hajek corrections
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