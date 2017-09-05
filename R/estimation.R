#------------------------------------------------------------------------------#
#' Step 0 of estimation
#' 
#' Creates treatment and outcome models
#' 
#' @param data dataset
#' @param model_args list of model arguments
#' @param ... additional arguments
#' @export
#------------------------------------------------------------------------------#

est_step0 <- function(data, model_args, ...){
  ## Component models
  m <- make_models(model_args = model_args[c('t_model', 'o_model')], 
                   data = data)
}

#------------------------------------------------------------------------------#
#' Step 1 of estimation
#' 
#' Creates  estimator arguments and models
#' 
#' @export
#------------------------------------------------------------------------------#

est_step1 <- function(data, step0, model_args, allocations, randomization, ...){
  
  m <- step0 
  ## Grab component model parameter estimates
  theta_t <- unlist(lme4::getME(m$t_model, c('beta', 'theta')))
  theta_o <- coef(m$o_model)
  
  m$wls_model_0 <- m$wls_model_1 <- vector('list', length(allocations))
  names(m$wls_model_0) <- names(m$wls_model_1) <- allocations
  
  hold <- lapply(seq_along(allocations), function(k){
    
    ipwv <- make_ipw_vector(fulldata = data, 
                            models  = m, 
                            group   = 'group', 
                            alpha   = allocations[k],
                            randomization = 1)
    # Add IP weights to dataset
    data <- data %>%
      mutate_(
        ipw  =~ ipwv,
        ipw0 =~ ipw * (A == 0),
        ipw1 =~ ipw * (A == 1)
      )
    ipw0 <- data[['ipw0']]
    ipw1 <- data[['ipw1']]
    A    <- data[['A']]
    
    # Fit regression-based models
    wls_model_0 <- try(glm(
      formula = model_args[['wls_model_0']][['formula']], 
      family  = model_args[['wls_model_0']][['options']][['family']],
      weights = ipw0,
      data    = data))
    
    wls_model_1 <- try(glm(
      formula = model_args[['wls_model_1']][['formula']], 
      family  = model_args[['wls_model_1']][['options']][['family']],
      weights = ipw1,
      data    = data))
    
    # in some data generating situations the glm model fails to converge,
    # this makes a way for the lapply over different estimators to skip
    # this situation, but to keep track that it failed.
    skipwls <<- FALSE
    if(is(wls_model_0, 'try-error') | is(wls_model_1, 'try-error')){
      skipwls <<- TRUE
    } else {
      m$wls_model_0[[k]] <<- wls_model_0
      m$wls_model_1[[k]] <<- wls_model_1
    }
  }) 
  
  theta_wls <- c(theta_t,
                 unlist(lapply(m$wls_model_0, coef)),
                 unlist(lapply(m$wls_model_1, coef))) 
  ## TEMPORARY!! (hopefully) ##
  estimator_args <- list(
    ipw = list(type      = 'ipw',
               hajek     = FALSE,
               theta     = theta_t,
               regtyp    = 'none',
               skipit    = FALSE),
    otc = list(type      = 'otc',
               hajek     = FALSE,
               theta     = theta_o,
               regtyp    = 'none',
               skipit    = FALSE),
    dbr = list(type      = 'dbr',
               hajek     = FALSE,
               theta     = c(theta_t, theta_o),
               regtyp    = 'none',
               skipit    = FALSE),
    wls_dbr = list(type  = 'wls_dbr',
                   theta = if(skipwls == TRUE) NA else theta_wls,
                   hajek = FALSE,
                   regtyp    = 'wls',
                   skipit    = skipwls))
  list(
    estimator_args = estimator_args,
    models         = m)
}
