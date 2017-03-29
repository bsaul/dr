#------------------------------------------------------------------------------#
#' Estimating equations for IPW, OTC, and DBR estimators
#' @export
#------------------------------------------------------------------------------#

generic_eefun <- function(
  data, 
  models, 
  randomization, 
  estimator_type, 
  regression_type = 'none',
  hajek = FALSE){
  
  comp <- extract_model_info(models = models, data = data, 
                             estimator_type = estimator_type,
                             regression_type = regression_type)
  A <- comp$A
  p <- comp$p
  
  ## Create estimating equation functions for nontarget parameters
  if(estimator_type %in% c('ipw', 'dbr', 'reg_dbr')){
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
 
  if(estimator_type == 'reg_dbr'){
    if(regression_type == 'wls'){
      score_fun_reg_0 <- geex::make_eefun(
        models$wls_model_0, 
        data = data)
      score_fun_reg_1 <- geex::make_eefun(
        models$wls_model_1, 
        data = data)
    } else if(regression_type == 'pcov'){
      score_fun_reg_0 <- geex::make_eefun(
        models$pcov_model_0, 
        data = data)
      score_fun_reg_1 <- geex::make_eefun(
        models$pcov_model_1, 
        data = data)
    }

    
    index_t <- 1:comp$p_t
    index_o_0 <- (comp$p_t + 1):(comp$p_t + comp$p_o_0)
    index_o_1 <- (comp$p_t + comp$p_o_0 + 1):(comp$p_t + comp$p_o_0 + comp$p_o_1)
  }
   
  ## Create estimating equation functions for target parameters
  estimatorFUN <- match.fun(paste0(estimator_type, '_estimator'))
  estimator <- estimatorFUN(
    data          = data, 
    models        = models,
    randomization = randomization,
    hajek         = hajek,
    regression_type = regression_type
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
    } else if(estimator_type == 'reg_dbr'){
      scores_t <- score_fun_t(theta[index_t])
      scores_reg_0 <- score_fun_reg_0(theta[index_o_0])
      scores_reg_1 <- score_fun_reg_1(theta[index_o_1])
      scores   <- c(scores_t, scores_reg_0, scores_reg_1)
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