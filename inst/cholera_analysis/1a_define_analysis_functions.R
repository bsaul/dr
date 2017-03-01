#------------------------------------------------------------------------------#
#   Title: Functions for analysis of cholera data with IPW, OTC,
#          and DR estimators
#  Author: B. Saul
#    Date: 2017-02-27
# Purpose: 
#------------------------------------------------------------------------------#

estimate_cholera_parms <- function(data, allocations, model_args, ...){
  
  ## Component models
  m <- make_models(model_args = model_args, data = data)
  
  ## Grab component model parameter estimates
  theta_t <- unlist(lme4::getME(m$t_model, c('beta', 'theta')))
  theta_o <- coef(m$o_model)
  
  temp <- list(eeFUN = generic_eefun, splitdt = split(data, data$group))
  
  ## TEMPORARY!! (hopefully) ##
  estimator_args <- list(
    ipw = list(type      = 'ipw',
               hajek     = FALSE,
               theta     = theta_t),
    otc = list(type      = 'otc',
               hajek     = FALSE,
               theta     = theta_o),
    dbr = list(type      = 'dbr',
               hajek     = FALSE,
               theta     = c(theta_t, theta_o))
    # ,
    # ipw_hjk = list(type  = 'ipw',
    #            hajek     = TRUE,
    #            theta     = theta_t),
    # dbr_hjk = list(type  = 'dbr',
    #            hajek     = TRUE,
    #            theta     = c(theta_t, theta_o))
  )
  
  all <- lapply(estimator_args, function(eargs){
    ## Estimate parameters for each allocation
    hold <- lapply(allocations, function(allocation){
      
      temp <- append(temp, list(ee_args = list(alpha = allocation)))
      p <- length(eargs$theta)
      make_estimator_fun <- match.fun(paste0(eargs$type, '_estimator'))
      
      ## BEGIN  Point estimates ##
      target <- lapply(temp$splitdt, function(grp_dt){
        # Create estimator function
        estimator <- make_estimator_fun(
          data          = grp_dt, 
          models        = m,
          randomization = 2/3,
          hajek         = eargs$hajek)
        # Evaluate estimator function
        estimator(eargs$theta, alpha = allocation)
      })
      
      target <- target %>% list_matrix() %>% apply(., 2, mean)
      print(target)
      ## END Point estimates ##
      
      ## BEGIN VCOV estimates ##
      mats <- geex::compute_matrices(
        geex_list        = temp,
        theta            = c(eargs$theta, target),
        numDeriv_options = list(method = 'simple'),
        models           = m,
        randomization    = 2/3,
        estimator_type   = eargs$type,
        hajek            = eargs$hajek)
      
      Sigma <- try(geex::compute_sigma(mats$A, mats$B), silent = TRUE)
      if(is(Sigma, 'try-error')){
        std_error <- NA
      } else {
        std_error <- sqrt(diag(Sigma)[(p + 1):(p + 2)])
      }
      ## END VCOV estimates ##
      
      # Convert to data_frame
      data_frame(
        method    = eargs$type,
        hajek     = eargs$hajek,
        a         = 0:1,
        alpha     = allocation,
        estimate  = target,
        std_error = std_error)
    }) # END lapply per allocation
    bind_rows(hold) # combine estimates for allocations
  }) # END lapply per estimator
  
  bind_rows(all) # combine estimates for estimators
}