#------------------------------------------------------------------------------#
#' Estimating equations for IPW, OTC, and DBR estimators
#' @export
#------------------------------------------------------------------------------#
lan_eefun <- function(data, t_model, o_model){
  
  Y <- geex::get_response(formula(o_model), data = data)
  A <- geex::get_response(formula(t_model), data = data)
  X_t <- geex::get_design_matrix(geex::get_fixed_formula(t_model), data = data)
  X_o <- as.data.frame(geex::get_design_matrix(geex::get_fixed_formula(o_model), data = data))
  X_o_ex <- expand_outcome_frame(X_o, geex::get_fixed_formula(o_model))
  inv_link_o <- family(o_model)$linkinv
  rhs_o <- geex::get_fixed_formula(o_model)
  
  ## components for IPW estimator
  ip_fun    <- weight_estimator(A = A, X = X_t)
  score_fun_t <- geex::make_eefun(t_model, data = data, 
                                  numderiv_opts = list(method = 'simple'))
  
  ## components for OTC estimator
  n_ <- nrow(X_o)  
  score_fun_o <- geex::make_eefun(o_model, data = data)
  
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
    scores_t <- score_fun_t(theta[index_t])
    
    ipw <- ip_fun(theta = theta[index_t], alpha = alpha)
    ipw_ce0 <- mean(Y * (A == 0)) * ipw / dbinom(0, 1, alpha)
    ipw_ce1 <- mean(Y * (A == 1)) * ipw / dbinom(1, 1, alpha)
    ipw_ce  <- mean(Y) * ipw 
    
    ### OTC estimator ###
    scores_o <- score_fun_o(theta[index_o])
    
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
    
    nalpha <- length(alpha)
    
    ### Estimating Equations ###
    c(scores_t, 
      scores_o,
      ipw_ce0 - theta[p + 1 + (0 * nalpha)],
      ipw_ce1 - theta[p + 1 + (1 * nalpha)],
      ipw_ce  - theta[p + 1 + (2 * nalpha)],
      otc_ce0 - theta[p + 1 + (3 * nalpha)],
      otc_ce1 - theta[p + 1 + (4 * nalpha)],
      otc_ce  - theta[p + 1 + (5 * nalpha)],
      dbr_ce0 - theta[p + 1 + (6 * nalpha)],
      dbr_ce1 - theta[p + 1 + (7 * nalpha)],
      dbr_ce  - theta[p + 1 + (8 * nalpha)])
  }
}