#------------------------------------------------------------------------------#
#   Title: Functions for analysis of cholera data with IPW, OTC,
#          and DR estimators
#  Author: B. Saul
#    Date: 2017-02-27
# Purpose: 
#------------------------------------------------------------------------------#

estimate_cholera_parms_step0 <- function(data, model_args, allocations, ...){
  
  ## Component models
  m <- make_models(model_args = model_args[c('t_model', 'o_model')], data = data)
}

estimate_cholera_parms_step1 <- function(data, models, model_args, allocations, randomization, ...){
  
  ## Grab component model parameter estimates
  theta_t <- unlist(lme4::getME(models$t_model, c('beta', 'theta')))
  theta_o <- coef(models$o_model)
  # Make estimator args
  m <- models
  
  m$wls_model_0 <- m$wls_model_1 <- m$pcov_model_0 <- m$pcov_model_1 <-  vector('list', length(allocations))
  names(m$wls_model_0) <-names(m$wls_model_1) <- names(m$pcov_model_0) <- names(m$pcov_model_1) <- allocations
  
  hold <- lapply(seq_along(allocations), function(k){
    ipwv <- make_ipw_vector(fulldata = data, 
                            models  = models, 
                            group   = 'group', 
                            alpha   = allocations[k],
                            randomization = randomization)
    
    # Add IP weights to dataset
    data <- data %>%
      mutate_(
        ipw  =~ ipwv,
        ipw0 =~ ipw * (A == 0),
        ipw1 =~ ipw * (A == 1)
      )
    ipw0 <- data[['ipw0']]
    ipw1 <- data[['ipw1']]
    A <- data[['A']]
    
    # Fit regression-based models
    m$wls_model_0[[k]] <<- try(glm(
      formula = model_args[['wls_model_0']][['formula']],
      family  = model_args[['wls_model_0']][['options']][['family']],
      weights = ipw0,
      data    = data))
    
    m$wls_model_1[[k]] <<- try(glm(
      formula = model_args[['wls_model_1']][['formula']],
      family  = model_args[['wls_model_1']][['options']][['family']],
      weights = ipw1,
      data    = data))
    # 
    # m$pcov_model_0[[k]] <<- try(glm(
    #   formula = model_args[['pcov_model_0']][['formula']],
    #   family  = model_args[['pcov_model_0']][['options']][['family']],
    #   weights = (A == 0) * 1,
    #   data    = data))
    # 
    # m$pcov_model_1[[k]] <<- try(glm(
    #   formula = model_args[['pcov_model_1']][['formula']],
    #   family  = model_args[['pcov_model_1']][['options']][['family']],
    #   weights = (A == 1) * 1,
    #   data    = data))
    
  })
 
  theta_wls <- c(theta_t,
                 unlist(lapply(m$wls_model_0, coef)),
                 unlist(lapply(m$wls_model_1, coef))) 

  estimator_args <- list(
    ipw = list(type      = 'ipw',
               hajek     = FALSE,
               theta     = theta_t,
               regtyp    = 'none',
               skipit    = TRUE),
    otc = list(type      = 'otc',
               hajek     = FALSE,
               theta     = theta_o,
               regtyp    = 'none',
               skipit    = TRUE),
    dbr = list(type      = 'dbr',
               hajek     = FALSE,
               theta     = c(theta_t, theta_o),
               regtyp    = 'none',
               skipit    = TRUE),
    wls_dbr = list(type    = 'wls_dbr',
                   theta   =  theta_wls,
                   hajek   = FALSE,
                   regtyp  = 'wls',
                   skipit  = FALSE)
    # pcov_dbr = list(type   = 'pcov_dbr',
    #                 theta  = c(theta_t, theta_pcov_0, theta_pcov_1),
    #                 hajek  = FALSE,
    #                 regtyp = 'pcov',
    #                 skipit = TRUE)
  )
  
  list(
     estimator_args = estimator_args,
     models         = m)
}


estimate_cholera_parms_step2 <- function(data, allocations, models, model_args, randomization, compute_se = TRUE, ...){  

  geexList <- list(eeFUN = generic_eefun, splitdt = split(data, data$group))
  
  ## Estimate parameters for each allocation
  hold <- lapply(allocations, function(allocation){
    
    temp <- estimate_cholera_parms_step1(data          = data, 
                                         models        = models, 
                                         model_args    = model_args, 
                                         allocations   = allocation, 
                                         randomization = randomization)
    estimator_args <- temp$estimator_args
    m <- temp$models
    
    all <- lapply(estimator_args, function(eargs){
      ## Estimate parameters for each allocation

      if(eargs$skip == FALSE){
        geexList <- append(geexList, list(ee_args = list(alpha = allocation)))
        p        <- length(eargs$theta)
        n_allocation <- length(allocation)
        make_estimator_fun <- match.fun(paste0(eargs$type, '_estimator'))

        ## BEGIN  Point estimates ##
        target <- lapply(geexList$splitdt, function(grp_dt){
          # Create estimator function
          estimator <- make_estimator_fun(
            data          = grp_dt,
            models        = m,
            randomization = randomization,
            hajek         = eargs$hajek,
            regression_type = eargs$regtyp)
          # Evaluate estimator function
          estimator(eargs$theta, alpha = allocation)
        })
        

        target <- target %>% list_matrix() %>% apply(., 2, mean)
        print(eargs$type)
        print(target)
        ## END Point estimates ##
        ## BEGIN VCOV estimates ##
        if(compute_se){
          mats <- geex::compute_matrices(
            geex_list        = geexList,
            theta            = c(eargs$theta, target),
            derivFUN_control = list(method = 'simple'),
            models           = m,
            randomization    = randomization,
            estimator_type   = eargs$type,
            hajek            = eargs$hajek,
            regression_type  = eargs$regtyp)

          Sigma <- try(geex::compute_sigma(mats$A, mats$B), silent = TRUE)
          print(Sigma)
          if(is(Sigma, 'try-error')){
            Sigma <- NA
          }
        } else {
          Sigma <- NA
        }

        ## END VCOV estimates ##
      } else if(eargs$skipit == TRUE){
        n_allocation <- length(allocation)
        target <- NA
        Sigma <- NA
      }

      ## END VCOV estimates ##

      list(
        alpha    = allocation,
        estimate = target,
        vcov     = Sigma
      )

    }) # END lapply per estimator
    all
  }) # END lapply per allocation
  hold
}
