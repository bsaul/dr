#------------------------------------------------------------------------------#
#' Estimating equations for IPW, OTC, and DBR estimators
#' @export
#------------------------------------------------------------------------------#

generic_eefun <- function(data, models, randomization, estimator_type, hajek = FALSE){
  
  comp <- extract_model_info(models = models, data = data, estimator_type)
  A <- comp$A
  p <- comp$p
  
  ## Create estimating equation functions for nontarget parameters
  if(estimator_type %in% c('ipw', 'dbr')){
    score_fun_t <- geex::make_eefun(
      models$t_model, 
      data = data, 
      numderiv_opts = list(method = 'simple'))
    p_t <- comp$p_t
    index_t <- 1:p_t # index for nontarget propensity model parameters
  }
  
  if(estimator_type %in% c('otc', 'dbr')){
    score_fun_o <- geex::make_eefun(
      models$o_model, 
      data = data)
    
    p_o <- comp$p_o
    
    if(estimator_type == 'dbr'){
      index_o <- (p_t + 1):(p_t + p_o) # index for nontarget propensity model parameters
    } else {
      index_o <- 1:comp$p_o
    }
  }
  
  ## Create estimating equation functions for target parameters
  estimatorFUN <- match.fun(paste0(estimator_type, '_estimator'))
  estimator <- estimatorFUN(
    data          = data, 
    models        = models,
    randomization = randomization,
    hajek         = hajek
  )
  
  ### PSI function ###
  function(theta, alpha){
    
    index_target <- (p + 1):(length(theta))
    ### Non-target parameters ###
    if(estimator_type == 'dbr'){
      scores_t <- score_fun_t(theta[index_t])
      scores_o <- score_fun_o(theta[index_o])
      scores   <- c(scores_t, scores_o)
    } else if(estimator_type == 'ipw'){
      scores <- score_fun_t(theta[index_t])
    } else if(estimator_type == 'otc'){
      scores <- score_fun_o(theta[index_o])
    }
    
    ### Target parameters ###
    grp_est <- estimator(theta[1:p], alpha = alpha)
   
    ### Estimating Equations ###
    if(hajek){
      c(scores,
        # cancel out contribution of group estimate AND target parameter
        (sum(A == 0) > 0 ) * (grp_est[1] - theta[p + 1]), 
        (sum(A == 1) > 0 ) * (grp_est[2] - theta[p + 2]))
    } else {
      c(scores,
        grp_est - theta[index_target])
    }
  }
}