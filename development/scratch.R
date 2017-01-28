## 
library(inferference)

vaccinesim

theta_t <- unlist(lme4::getME(models$model_treatment, c('beta', 'theta')))
Y   <- get_response(form_o, group_data)
A   <- get_response(form_t, group_data)
X_t <- get_design_frame(form_rhs_t, group_data)
X_o <- get_design_frame(form_rhs_o, group_data)

make_ipw_estimator(Y = Y, A = A, X_treatment = X_t)

ipw_eefun <- function(data, t_model, o_model){
  
  Y <- get_response(o_model)
  A <- get_response(t_model)
  X <- get_design_frame(get_fixed_formula(t_model), data = data)
  
  ip_fun    <- weight_estimator(A = A, X = X)
  score_fun <- make_ee(t_model)
  
  function(theta, alpha1, alpha2, a1, a2){
    p <- length(theta) - 2
    
    scores <- score_fun(theta[1:p])
    
    ipw1 <- ip_fun(theta = theta[1:p], alpha = alpha1)
    ipw2 <- ip_fun(theta = theta[1:p], alpha = alpha2)
    
    Ia1 <- if(is.null(a1)) 1 else (A == a) * 1
    Ia2 <- if(is.null(a2)) 1 else (A == a) * 1
    
    ce1 <- mean(Y * Ia1) * ipw1 / {if(!is.null(a1)) dbinom(a1, 1, alpha1) else 1}
    ce2 <- mean(Y * Ia2) * ipw2 / {if(!is.null(a2)) dbinom(a2, 1, alpha2) else 1}

    c(scores, 
      ce1 - theta[p + 1],
      ce2 - theta[p + 2])
  }
}