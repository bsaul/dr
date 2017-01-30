#### IPW eefun ####

#------------------------------------------------------------------------------#
#' IPW estimator
#' @export
#------------------------------------------------------------------------------#
ipw_eefun <- function(data, t_model, o_model){
  
  Y <- geex::get_response(formula(o_model), data = data)
  A <- geex::get_response(formula(t_model), data = data)
  X_t <- geex::get_design_matrix(geex::get_fixed_formula(t_model), data = data)
  X_o <- as.data.frame(geex::get_design_matrix(geex::get_fixed_formula(o_model), data = data))
  print(X_o)
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
  
  ## indices
  p_t <- ncol(X_t) + 1
  p_o <- ncol(X_o)
  index_t <- 1:p_t
  index_o <- (p_t + 1):(p_t + p_o)
  p  <- p_t + p_o
  
  function(theta, alpha1, alpha2, a1, a2){
    
    ### IPW estimator ###
    scores_t <- score_fun_t(theta[index_t])
    
    ipw1 <- ip_fun(theta = theta[index_t], alpha = alpha1)
    ipw2 <- ip_fun(theta = theta[index_t], alpha = alpha2)
    
    Ia1 <- if(is.null(a1)) 1 else (A == a1) * 1
    Ia2 <- if(is.null(a2)) 1 else (A == a2) * 1
   
    ipw_ce1 <- mean(Y * Ia1) * ipw1 / {if(!is.null(a1)) dbinom(a1, 1, alpha1) else 1}
    ipw_ce2 <- mean(Y * Ia2) * ipw2 / {if(!is.null(a2)) dbinom(a2, 1, alpha2) else 1}
    
    ### OTC estimator ###
    scores_o <- score_fun_o(theta[index_o])
    
    nn1_o <- n_ - {if(is.null(a1)) 0 else 1}
    nn2_o <- n_ - {if(is.null(a2)) 0 else 1}
    a1_o <- if(is.null(a1)) 0:1 else a1
    a2_o <- if(is.null(a2)) 0:1 else a2
    
    x1 <- X_o_ex %>%
      filter_(~ A %in% a1_o) %>%
      mutate_(pi_value1 = ~ dbinom(sum_a, nn1_o, prob = alpha1),
              fitted   = ~ as.numeric(inv_link_o(model.matrix(rhs_o, data = .) %*% theta[index_o] )),
              estimate = ~ fitted * pi_value1) 
    
    x2 <- X_o_ex %>%
      filter_(~ A %in% a2_o) %>%
      mutate_(pi_value2 = ~ dbinom(sum_a, nn1_o, prob = alpha2),
              fitted   = ~ as.numeric(inv_link_o(model.matrix(rhs_o, data = .) %*% theta[index_o] )),
              estimate = ~ fitted * pi_value2) 
    
    otc_ce1 <- 
      sum(tapply(tapply(x1$estimate, paste(x1$A, x1$ID), sum), 
               rep(a1_o, each = n_), sum)) / n_
    
    otc_ce2 <- 
      sum(tapply(tapply(x2$estimate, paste(x2$A, x2$ID), sum), 
                 rep(a2_o, each = n_), sum)) / n_
    
    ### Estimating Equations ###
    c(scores_t, 
      scores_o,
      ipw_ce1 - theta[p + 1],
      ipw_ce2 - theta[p + 2],
      otc_ce1 - theta[p + 3],
      otc_ce2 - theta[p + 4])
  }
}