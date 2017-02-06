#------------------------------------------------------------------------------#
#' Estimating equations for IPW, OTC, and DBR estimators
#' @export
#------------------------------------------------------------------------------#

dr_eefun <- function(data, t_model, o_model){
  
  score_fun_t <- geex::make_eefun(t_model, data = data, 
                                  numderiv_opts = list(method = 'simple'))
  score_fun_o <- geex::make_eefun(o_model, data = data)
  estimators <- dr_estimators(data = data, t_model = t_model, o_model = o_model)
  
  ## indices
  X_t <- geex::get_design_matrix(geex::get_fixed_formula(t_model), data = data)
  X_o <- geex::get_design_matrix(geex::get_fixed_formula(o_model), data = data)
  p_t <- ncol(X_t) + 1
  p_o <- ncol(X_o)
  p   <- p_t + p_o
  index_t <- 1:p_t
  index_o <- (p_t + 1):(p_t + p_o)
  
  function(theta, alpha){
    
    index_target <- (p + 1):(length(theta))
    ### Non-target parameters ###
    scores_t <- score_fun_t(theta[index_t])
    scores_o <- score_fun_o(theta[index_o])
    
    ### Target parameters ###
    target <- estimators(theta, alpha = alpha)

    ### Estimating Equations ###
    c(scores_t, 
      scores_o,
      target - theta[index_target])
  }
}