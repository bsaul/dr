#------------------------------------------------------------------------------#
#' Estimating equations for WLS scores
#' @export
#------------------------------------------------------------------------------#

wls_score_maker <- function(theta, X, Y, phi, ipw, A, a, family_link, invlnk){
  W   <- I(A == a) * ipw
  
  lp <- X %*% theta # linear predictor
  f  <- as.numeric(invlnk(lp))  # fitted values
  r  <- Y - f       # residuals
  n  <- length(f)
  ### TODO: this is cludgy and needs to be reworked to be more general
  if(family_link == 'gaussian_identity'){
    D <- X
    V <- phi * diag(1, nrow = n, ncol = n)
  } else if(family_link %in% c('binomial_logit', 'quasibinomial_logit')){
    D <- apply(X, 2, function(x) x * exp(lp)/((1+exp(lp))^2) )
    V <- phi * diag(f * (1 - f), ncol = length(f) )/length(f)

    
  }
  
  t(D) %*% solve(V) %*% diag(W, nrow = n, ncol = n) %*% (r)
  # t(D) %*% diag(W, nrow = n, ncol = n) %*% (r)
}

#------------------------------------------------------------------------------#
#' Estimating equations for WLS scores
#' @export
#------------------------------------------------------------------------------#

make_eefun_wls <- function(data, models, randomization = 1){

  # Treatment (IPW part)
  
  A  <- geex::grab_response(A ~ 1, data = data)
  X_t<- geex::grab_design_matrix(
    data = data, 
    rhs_formula = geex::grab_fixed_formula(models$t_model))
  
  n  <- length(A)
  p_t <- ncol(X_t) + 1
  index_tt <- 1:(p_t)
  
  score_fun_t <- geex::grab_psiFUN(
    models$t_model, 
    data = data, 
    numderiv_opts = list(method = 'simple'))
  ## components for IPW estimator
  ip_fun <- weight_estimator(
    A = A, 
    X = X_t, 
    randomization = randomization)

  # WLS models
  wls0 <- models$wls_model_0
  wls1 <- models$wls_model_1
  wls  <- append(wls0, wls1)
  length0 <- length(wls0)
  
  coef_lengths <- unlist(lapply(wls, function(x){ length(coef(x))}))
  index_starts <- (cumsum(coef_lengths) - coef_lengths[1]) + (p_t + 1)
  alphas <- as.numeric(names(wls))

  
  score_parts <- lapply(seq_along(wls), function(j) {
    x     <- wls[[j]]
    X     <- model.matrix(formula(x), data = data)
    Y     <- as.numeric(model.frame(geex::grab_response_formula(x), data = data)[[1]])
    p_o   <- length(coef(x))
    phi   <- suppressWarnings(as.numeric(summary(x)$dispersion[1]))
    a     <- if(j <= length0) 0 else 1
    alpha <- alphas[j]
    index <- index_starts[j]:(index_starts[j] + (p_o - 1))
    
    list(X = X, Y = Y, p_o = p_o, phi = phi, alpha = alpha, a = a, 
         index = index)
  })
  
  family <- wls0[[1]]$family$family
  link   <- wls0[[1]]$family$link
  invlnk <- wls0[[1]]$family$linkinv
  family_link <- paste(family, link, sep = '_')
  
  function(theta){
    scores_t <- score_fun_t(theta[index_tt])
    
    
    scores_wls <- lapply(score_parts, function(x){
      denom <- (x$alpha^x$a) * (1- x$alpha)^(1-x$a)
      ipw <- ip_fun(theta = theta[index_tt], alpha = x$alpha)/denom
      
      wls_score_maker(
        theta = theta[x$index],
        X     = x$X,
        Y     = x$Y,
        phi   = x$phi,
        ipw   = ipw,
        A     = A,
        a     = x$a,
        invlnk = invlnk,
        family_link = family_link)
    }) %>% unlist()
    
    #Estimating functions
    c(scores_t, scores_wls)
  }
}

#------------------------------------------------------------------------------#
#' Estimator used in wls estFUN
#' @export
#------------------------------------------------------------------------------#


wls_dbr_estimator_estFUN <- function(data, models, randomization, regression_type = 'wls', ...){
  
  ## component data
  comp <- extract_model_info(
    model = models, 
    data = data, 
    estimator_type = 'wls_dbr',
    regression_type = regression_type)
  
  Y <- comp$Y
  A <- comp$A
  ## components for IPW estimator
  ip_fun <- weight_estimator(
    A = comp$A, 
    X = comp$X_t, 
    randomization = randomization)
  
  ## components for OTC part
  dr_term1_fun_0 <- make_dr_term1(
    comp$X_o_reg_0, 
    inv_link = comp$inv_link_o)
  
  dr_term1_fun_1 <- make_dr_term1(
    comp$X_o_reg_1, 
    inv_link = comp$inv_link_o)
  
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
      
      
      fY_0    <- dr_term1_fun_0(theta[index0])
      fY_1    <- dr_term1_fun_1(theta[index1])
      ipw     <- ip_fun(theta[index_t], alpha[k])
      Ybar0   <- sum((A == 0) * (Y - fY_0) )
      Ybar1   <- sum((A == 1) * (Y - fY_1) )
      term1_0 <- Ybar0 * ipw / (1 - alpha[k])
      term1_1 <- Ybar1 * ipw / alpha[k]
      
      
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
      
      ce0[k] <- (dbr2_ce_a[1] + term1_0)/N 
      ce1[k] <- (dbr2_ce_a[2] + term1_1)/N
      
    }
    
    x <- c(ce0, ce1) 
    label0 <- paste0(regression_type, '_dbr_Y0_')
    label1 <- paste0(regression_type, '_dbr_Y1_')
    names(x) <- paste0(rep(c(label0, label1), each = length(alpha)), alpha)
    x
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
    score_fun_t <- geex::grab_psiFUN(
      models$t_model, 
      data = data, 
      numderiv_opts = list(method = 'simple'))
    p_t <- comp$p_t
    index_t <- 1:p_t # index for nontarget propensity model parameters
  }
  
  if(estimator_type %in% c('otc', 'dbr')){
    score_fun_o <- geex::grab_psiFUN(
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
    score_funs <- make_eefun_wls(models = models, 
                   data   = data,
                   randomization = randomization)
   
  }
   
  ## Create estimating equation functions for target parameters
  estimatorFUN <- match.fun(paste0(estimator_type, '_estimator'))
  if(estimator_type == "wls_dbr"){
    estimatorFUN <- wls_dbr_estimator_estFUN
  }
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
      scores <- score_funs(theta[-index_target])
      # scores <- c(scores_t, scores_0, scores_1)
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
