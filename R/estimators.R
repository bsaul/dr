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
      out <- c(
        ifelse(sum(A == 0) > 0, ipw_0/Nhat0, 0), 
        ifelse(sum(A == 1) > 0, ipw_1/Nhat1, 0))
      names(out) <- paste0(rep(c('ipw_hjk_Y0_', 'ipw_hjk_Y1_'), each = length(alpha)), alpha)
    } else {
      out <- c(ipw_0/N, ipw_1/N)
      names(out) <- paste0(rep(c('ipw_Y0_', 'ipw_Y1_'), each = length(alpha)), alpha)
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
  index_o <- 1:comp$p
 
  N    <- comp$N
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
    names(x) <- paste0(rep(c('otc_Y0_', 'otc_Y1_'), each = length(alpha)), alpha)
    x
  }
}
#------------------------------------------------------------------------------#
#' Makes DBR estimators for group-level data
#' @export
#------------------------------------------------------------------------------#

dbr_estimator <- function(data, models, randomization, hajek = FALSE, ...){
  
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
    
    dbr_0 <- (term1_0 + term2[1:length(alpha)]*N)
    dbr_1 <- (term1_1 + term2[(length(alpha) + 1):length(term2)]*N)
    
    if(hajek){
      Nhat0 <- sum(A == 0) * ipw / (1 - alpha)
      Nhat1 <- sum(A == 1) * ipw / alpha
      out <- c(
        ifelse(sum(A == 0) > 0, term1_0/Nhat0 + term2[1], 0),
        ifelse(sum(A == 1) > 0, term1_1/Nhat1 + term2[2], 0))
      names(out) <- paste0(rep(c('dbr_hjk_Y0_', 'dbr_hjk_Y1_'), each = length(alpha)), alpha)
    } else {
      out <- c(dbr_0/N, dbr_1/N)
      names(out) <- paste0(rep(c('dbr_Y0_', 'dbr_Y1_'), each = length(alpha)), alpha)
    }
    out
  }
}

#------------------------------------------------------------------------------#
#' Makes WLS (DBR) estimators for group-level data
#' (i.e.) regression based
#' @export
#------------------------------------------------------------------------------#

wls_dbr_estimator <- function(data, models, randomization, regression_type = 'wls', ...){
  
  ## component data
  comp <- extract_model_info(
    model = models, 
    data = data, 
    estimator_type = 'wls_dbr',
    regression_type = regression_type)
  
  X_ex_0 <- comp$X_o_ex %>% filter(A == 0)
  X_ex_1 <- comp$X_o_ex %>% filter(A == 1)
  MM_0   <- model.matrix(comp$rhs_o_reg_0, data = X_ex_0)
  MM_1   <- model.matrix(comp$rhs_o_reg_1, data = X_ex_1)

  index_t   <- 1:comp$p_t
  N    <- comp$N
  
  function(theta, alpha){
    # stopifnot(length(alpha) == 1)
    # theta should be ordere: theta_0_alpha1, theta_0_alpha2, ..., theta_1_alpha1, theta_1_alpha2, ... 
    ce0 <- ce1 <- numeric(length(alpha))
    index0 <- (comp$p_t + 1):(comp$p_t + comp$p_o_0)
    index1 <- (comp$p_t + (length(alpha)*comp$p_o_0) + 1):(comp$p_t + (length(alpha)*comp$p_o_0) + comp$p_o_1)
    ### Regression-based DRR estimator ###
    for(k in 1:length(alpha)){
      if(k > 1){
        index0 <- index0 + comp$p_o_0
        index1 <- index1 + comp$p_o_1
      }

      # compute fitted value for expanded data.frame
      mu_0 <- as.numeric(comp$inv_link_o(MM_0 %*% theta[index0]))
      mu_1 <- as.numeric(comp$inv_link_o(MM_1 %*% theta[index1]))
      # compute pi term per number treated in group per subject
      pi_term_a <- dbinom(X_ex_0$sum_a, comp$N - 1, alpha[k])
      
      # mulitply mu_ij by the pi term rowwise
      piXmu_a_0 <- mu_0 * pi_term_a
      piXmu_a_1 <- mu_1 * pi_term_a
      
      # sum within levels of A (0:1) WITHIN subjects
      piXmu_a <- tapply(
        X     = c(piXmu_a_0, piXmu_a_1), 
        INDEX = paste(rep(0:1, each = nrow(X_ex_0)), c(X_ex_0$ID, X_ex_1$ID)), 
        FUN   = sum)
      
      # sum within levels of A (0:1) ACROSS subjects
      dbr2_ce_a <- tapply(
        X     = piXmu_a, 
        INDEX = rep(0:1, each = N), 
        FUN   = sum)
      
      ce0[k] <- dbr2_ce_a[1]/N
      ce1[k] <- dbr2_ce_a[2]/N
      
    }

    x <- c(ce0, ce1) 
    label0 <- paste0(regression_type, '_dbr_Y0_')
    label1 <- paste0(regression_type, '_dbr_Y1_')
    names(x) <- paste0(rep(c(label0, label1), each = length(alpha)), alpha)
    x
  }
}

#------------------------------------------------------------------------------#
#' Makes regression with propensity (DBR) estimators for group-level data
#' (i.e.) regression based
#' @export
#------------------------------------------------------------------------------#

pcov_dbr_estimator <- function(data, models, randomization, ...){
  
  ## component data
  comp <- extract_model_info(
    model = models, 
    data = data, 
    estimator_type = 'pcov_dbr',
    regression_type = 'pcov')
  
  A  <- comp$A
  N  <- comp$N
  lnkinv <- comp$inv_link_o
  
  ip_fun <- weight_estimator(
    A = comp$A, 
    X = comp$X_t, 
    randomization = randomization)
  
  dr_term1_fun_0 <- make_dr_term1(
    comp$X_o_reg_0, 
    inv_link = comp$inv_link_o)
  
  dr_term1_fun_1 <- make_dr_term1(
    comp$X_o_reg_1, 
    inv_link = comp$inv_link_o)
  
  MM_0   <- model.matrix(comp$rhs_o_reg_0, data = comp$X_o_reg_0)
  MM_1   <- model.matrix(comp$rhs_o_reg_1, data = comp$X_o_reg_1)
  
  index_t   <- 1:comp$p_t
  index_o_0 <- (comp$p_t + 1):(comp$p_t + comp$p_o_0)
  index_o_1 <- (comp$p_t + comp$p_o_0 + 1):(comp$p_t + comp$p_o_0 + comp$p_o_1)
  
  function(theta, alpha){
    stopifnot(length(alpha) == 1)
    
    ## Method 1 ##
    K <- 10 # number of resamples for A_tilde
    out_0 <- out_1 <- numeric(K)
    for(k in 1:K){
      # Sample an A_tilde for each k
      A_tilde     <- rbinom(N, 1, prob = alpha)
      
      # Replace A with A_tilde
      new_data    <- data
      new_data$A  <- A_tilde
      
      # With A_tilde, compute pi/f PER subject after setting 
      # a given subject's a_ij to 0.
      
      hold_0 <- hold_1 <- numeric(nrow(data))
      # compute ip weight and m_ij PER subject
      for(j in 1:nrow(data)){
        dt_0 <- dt_1 <- data
        A_new_0 <- A_new_1 <- data$A
        A_new_0[j] <- 0 # set A_ij to a_ij = 0
        A_new_1[j] <- 1 # set A_ij to a_ij = 1
        
        ip_fun_0 <- weight_estimator(
          A = A_new_0,
          X = comp$X_t)
        
        ip_fun_1 <- weight_estimator(
          A = A_new_1,
          X = comp$X_t)
        
        # IPW has pi = prod_i^n, so to remove contribution of jth subject,
        # divide by (1 - alpha)^(1 - 0)
        dt_0$ipw[j] <- ip_fun_0(theta[index_t], alpha)/( 1 - alpha)
        dt_1$ipw[j] <- ip_fun_1(theta[index_t], alpha)/(alpha)
        
        
        j_data_0    <- dt_0[j, , drop = FALSE]
        j_data_0$fA <- mean(A_new_0)
        j_data_1    <- dt_1[j, , drop = FALSE]
        j_data_1$fA <- mean(A_new_1)
        j_row_0     <- model.matrix(comp$rhs_o_reg_0, data = j_data_0)
        j_row_1     <- model.matrix(comp$rhs_o_reg_1, data = j_data_1)
        
        m_ij_0    <- lnkinv(j_row_0 %*% theta[index_o_0])
        m_ij_1    <- lnkinv(j_row_1 %*% theta[index_o_1])
        hold_0[j] <- m_ij_0
        hold_1[j] <- m_ij_1
      }
      
      out_0[k] <- mean(hold_0)
      out_1[k] <- mean(hold_1)
    }
    
    x <- c(mean(out_0), mean(out_1))
    x
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
}

#------------------------------------------------------------------------------#
#' Create vector of IP weights
#' @export
#------------------------------------------------------------------------------#


make_ipw_vector <- function(fulldata, models, group, randomization = 1, alpha){
  splitdt <- split(fulldata, fulldata[[group]])
  
  lapply(splitdt, function(x){
    comp <- extract_model_info(model = models, data = x, 'dbr')
    Y   <- comp$Y
    A   <- comp$A
    N   <- comp$N
    X_o <- comp$X_o
    f_o <- terms.formula(comp$rhs_o)
    L   <- model.matrix(drop.terms(f_o, attr(f_o, 'term.labels') == 'A'),
                        data = x)
    
    ## components for IPW part
    ip_fun <- weight_estimator(
      A = A, 
      X = comp$X_t, 
      randomization = randomization)
    
    IPW <- ip_fun(unlist(lme4::getME(models$t_model, c('beta', 'theta'))), alpha = alpha)
    IPW / dbinom(A, 1, prob = alpha)
    
  }) %>%
    unlist()
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