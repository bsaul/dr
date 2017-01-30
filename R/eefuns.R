#### IPW eefun ####

#------------------------------------------------------------------------------#
#' IPW estimator
#' @export
#------------------------------------------------------------------------------#
ipw_eefun <- function(data, t_model, o_model){
  
  Y <- geex::get_response(formula(o_model), data = data)
  A <- geex::get_response(formula(t_model), data = data)
  X <- geex::get_design_matrix(get_fixed_formula(t_model), data = data)
  
  ip_fun    <- weight_estimator(A = A, X = X)
  score_fun <- geex::make_eefun(t_model, data = data, numderiv_opts = list(method = 'simple'))
  
  function(theta, alpha1, alpha2, a1, a2){
    p <- length(theta) - 2
    
    scores <- score_fun(theta[1:p])
    
    ipw1 <- ip_fun(theta = theta[1:p], alpha = alpha1)
    ipw2 <- ip_fun(theta = theta[1:p], alpha = alpha2)
    
    Ia1 <- if(is.null(a1)) 1 else (A == a1) * 1
    Ia2 <- if(is.null(a2)) 1 else (A == a2) * 1
   
    ce1 <- mean(Y * Ia1) * ipw1 / {if(!is.null(a1)) dbinom(a1, 1, alpha1) else 1}
    ce2 <- mean(Y * Ia2) * ipw2 / {if(!is.null(a2)) dbinom(a2, 1, alpha2) else 1}
    
    c(scores, 
      ce1 - theta[p + 1],
      ce2 - theta[p + 2])
  }
}