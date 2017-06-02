#------------------------------------------------------------------------------#
#' Estimating equations for WLS scores
#' @export
#------------------------------------------------------------------------------#

make_eefun_wls <- function(t_model, wls_model, data, a, randomization = 1){

  X  <- model.matrix(wls_model$formula, data = data)
  Y  <- as.numeric(model.frame(geex::get_response_formula(wls_model), data = data)[[1]])
  A  <- geex::get_response(A ~ 1, data = data)
  X_t<- geex::get_design_matrix(geex::get_fixed_formula(t_model), data)
  n  <- length(Y)
  p_o <- length(coef(wls_model))
  p_t <- ncol(X_t) + 1
  phi    <- as.numeric(summary(wls_model)$dispersion[1])
  
  ## components for IPW estimator
  ip_fun <- weight_estimator(
    A = A, 
    X = X_t, 
    randomization = randomization)

  family <- wls_model$family$family
  link   <- wls_model$family$link
  invlnk <- wls_model$family$linkinv
  family_link <- paste(family, link, sep = '_')
  
  # stopifnot(length(W) == 1 | length(W) == n)
  # if(length(W) == 1){
  #   W <- rep(W, n)
  # }
  index_tt <- 1:(p_t)
  index_oo <- (p_t + 1):(p_t + p_o)
  
  function(theta, alpha){
    ipw <- ip_fun(theta = theta[index_tt], alpha = alpha)
    W   <- I(A == a) * ipw
    
    lp <- X %*% theta[index_oo] # linear predictor
    f  <- as.numeric(invlnk(lp))  # fitted values
    r  <- Y - f       # residuals
    
    ### TODO: this is cludgy and needs to be reworked to be more general
    if(family_link == 'gaussian_identity'){
      D <- X
      V <- phi * diag(1, nrow = n, ncol = n)
    } else if(family_link == 'binomial_logit'){
      D <- apply(X, 2, function(x) x * exp(lp)/((1+exp(lp))^2) )
      V <- phi * diag(f * (1 - f), ncol = length(f) )/length(f)
    }
    
    t(D) %*% solve(V) %*% diag(W, nrow = n, ncol = n) %*% (r)
  }
}

#------------------------------------------------------------------------------#
#' Estimating equations for IPW, OTC, DBR, and WLS_DBR estimators
#' @export
#------------------------------------------------------------------------------#

generic_eefun <- function(
  data, 
  models, 
  randomization = 1, 
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
      make_eefun_wls(t_model = models$t_model, wls_model = x, 
                     data = data, a = 0, randomization = randomization)})
    score_funs_1 <- lapply(models$wls_model_1, function(x) {
      make_eefun_wls(t_model = models$t_model, wls_model = x, 
                     data = data, a = 1, randomization = randomization)})

    
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
        score_funs_0[[k]](theta[c(index_t, index0)], alpha[k])
      }) %>% unlist() -> scores_0
      
      lapply(seq_along(score_funs_1), function(k){
        if(k > 1){
          index1 <<- index1 + comp$p_o_1
        }
        score_funs_1[[k]](theta[c(index_t, index1)], alpha[k])
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