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
  if(estimator_type %in% c('ipw', 'dbr', 'wls_dbr')){
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
 
  if(estimator_type == 'wls_dbr'){
    score_funs_0 <- lapply(models$wls_model_0, function(x) {
        geex::make_eefun(x, data = data)})
    score_funs_1 <- lapply(models$wls_model_1, function(x) {
        geex::make_eefun(x, data = data)})

    
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
    
    index_target <- ((length(theta) - (length(alpha) * 2)) + 1):(length(theta))
    
    ### Non-target parameters ###
    if(estimator_type == 'dbr'){
      scores_t <- score_fun_t(theta[index_t])
      scores_o <- score_fun_o(theta[index_o])
      scores   <- c(scores_t, scores_o)
    } else if(estimator_type == 'ipw'){
      scores <- score_fun_t(theta[index_t])
    } else if(estimator_type == 'otc'){
      scores <- score_fun_o(theta[index_o])
    } else if(estimator_type == 'wls_dbr'){
      scores_t <- score_fun_t(theta[index_t])
      index0 <- (comp$p_t + 1):(comp$p_t + comp$p_o_0)
      index1 <- (comp$p_t + (length(alpha)*comp$p_o_0) + 1):(comp$p_t + (length(alpha)*comp$p_o_0) + comp$p_o_1)
      
      lapply(seq_along(score_funs_0), function(k){
        if(k > 1){
          index0 <<- index0 + comp$p_o_0
        }
        score_funs_0[[k]](theta[index0])
      }) %>% unlist() -> scores_0
      
      lapply(seq_along(score_funs_1), function(k){
        if(k > 1){
          index1 <<- index1 + comp$p_o_1
        }
        score_funs_1[[k]](theta[index1])
      }) %>% unlist() -> scores_1

      scores <- c(scores_t, scores_0, scores_1)
    }
    
    ### Target parameters ###
    grp_est <- estimator(theta[-index_target], alpha = alpha)
   
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