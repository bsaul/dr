
est_sim <- function(simdt, allocations, model_args, compute_se, ...){
  
  ## Component models
  m <- make_models(model_args = model_args[c('t_model', 'o_model')], 
                   data = simdt)

  ## Grab component model parameter estimates
  theta_t <- unlist(lme4::getME(m$t_model, c('beta', 'theta')))
  theta_o <- coef(m$o_model)

  hold <- lapply(allocations, function(allocation){
    
    ipw <- make_ipw_vector(fulldata = simdt, models = m, group = 'group', 
                           alpha = allocation)
    A <- simdt$A
    # Add IP weights to dataset
    simdt <- simdt %>%
      mutate_(
        ipw0 =~ ipw * (A == 0),
        ipw1 =~ ipw * (A == 1)
      )
    ipw0 <- simdt$ipw0
    ipw1 <- simdt$ipw1
    #print(model_args)
    # Fit regression-based models
    wls_model_0 <- try(glm(
      formula = model_args[['wls_model_0']][['formula']], 
      family  = model_args[['wls_model_0']][['options']][['family']],
      weights = ipw0,
      data    = simdt))

    wls_model_1 <- try(glm(
      formula = model_args[['wls_model_1']][['formula']], 
      family  = model_args[['wls_model_1']][['options']][['family']],
      weights = ipw1,
      data    = simdt))
   
    # in some data generating situations the glm model fails to converge,
    # this makes a way for the lapply over different estimators to skip
    # this situation, but to keep track that it failed.
    if(is(wls_model_0, 'try-error') | is(wls_model_1, 'try-error')){
      skipwls <- TRUE
      theta_wls_0 <- NA
      theta_wls_1 <- NA
    } else {
      skipwls <- FALSE
      m$wls_model_0 <- wls_model_0
      m$wls_model_1 <- wls_model_1
      theta_wls_0 <- coef(m$wls_model_0)
      theta_wls_1 <- coef(m$wls_model_1)
    }
    
    pcov_model_0 <- try(glm(
      formula = model_args[['pcov_model_0']][['formula']], 
      family  = model_args[['pcov_model_0']][['options']][['family']],
      weights = (A == 0) * 1,
      data    = simdt))
    
    pcov_model_1 <- try(glm(
      formula = model_args[['pcov_model_1']][['formula']], 
      family  = model_args[['pcov_model_1']][['options']][['family']],
      weights = (A == 1) * 1,
      data    = simdt))
    
    if(is(pcov_model_0, 'try-error') | is(pcov_model_1, 'try-error')){
      skippcov <- TRUE
      theta_pcov_0 <- NA
      theta_pcov_1 <- NA
    } else {
      skippcov <- FALSE
      m$pcov_model_0  <- pcov_model_0
      m$pcov_model_1  <- pcov_model_1
      theta_pcov_0 <- coef(m$pcov_model_0)
      theta_pcov_1 <- coef(m$pcov_model_1)
    }

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
      wls_dbr = list(type  = 'reg_dbr',
                     theta = if(skipwls == TRUE) {NA} else {
                                    c(theta_t, theta_wls_0, theta_wls_1)},
                     hajek = FALSE,
                     regtyp    = 'wls',
                     skipit    = skipwls),
      pcov_dbr = list(type  = 'reg_dbr',
                     theta = if(skippcov == TRUE) {NA} else {
                                    c(theta_t, theta_pcov_0, theta_pcov_1)},
                     hajek = FALSE,
                     regtyp    = 'pcov',
                     skipit    = skippcov)
    )
    
    
    temp <- list(eeFUN = generic_eefun, splitdt = split(simdt, simdt$group))
    
    
    all <- lapply(estimator_args, function(eargs){
    ## Estimate parameters for each allocation

       if(eargs$skip == FALSE){
         temp <- append(temp, list(ee_args = list(alpha = allocation)))
         p    <- length(eargs$theta)
         n_allocation <- length(allocation)
         make_estimator_fun <- match.fun(paste0(eargs$type, '_estimator'))
         
         ## BEGIN  Point estimates ##
         target <- lapply(temp$splitdt, function(grp_dt){
           # Create estimator function
           estimator <- make_estimator_fun(
             data          = grp_dt, 
             models        = m,
             randomization = 1,
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
             geex_list        = temp,
             theta            = c(eargs$theta, target),
             numDeriv_options = list(method = 'simple'),
             models           = m,
             randomization    = 1,
             estimator_type   = eargs$type,
             hajek            = eargs$hajek,
             regression_type  = eargs$regtyp)

           Sigma <- try(geex::compute_sigma(mats$A, mats$B), silent = TRUE)
           if(is(Sigma, 'try-error')){
             std_error <- NA
           } else {
             std_error <- sqrt(diag(Sigma)[(p + 1):(p + (n_allocation*2))])
           }
         } else {
           std_error <- NA
         }

         ## END VCOV estimates ##
       } else if(eargs$skipit == TRUE){
         n_allocation <- length(allocation)
         target <- NA
         std_error <- NA
       }

        
        # Convert to data_frame
        data_frame(
          method    = eargs$type,
          regression_type = eargs$regtyp,
          hajek     = eargs$hajek,
          a         = rep(0:1, each = n_allocation),
          alpha     = rep(allocation, times = n_allocation),
          estimate  = target,
          std_error = std_error
          )
    }) # END lapply per allocation
    bind_rows(all) # combine estimates for estimators
  }) # END lapply per estimator
  
  bind_rows(hold) # combine estimates for allocations
}

arg_maker <- function(scenario, nsims){
  temp <- scenarios[scenario, ]
  list(
    m     = temp$m,
    ni    = temp$ni,
    beta  = temp$beta[[1]],
    theta = temp$theta[[1]],
    gamma = temp$gamma[[1]],
    nsims = nsims
  )
}

do_sim <- function(scenario_no, allocations, all_model_args, compute_se, ...){
  simdt <- do.call(gen_sim, args = arg_maker(scenario_no, nsims = 1))
  on.exit(rm(simdt))
  out <- lapply(seq_along(all_model_args), function(j){
    est_sim(
      simdt       = simdt[[1]],
      allocations = allocations, 
      model_args  = all_model_args[[j]],
      compute_se  = compute_se) %>%
      mutate_(
        sid = ~scenario_no,
        model_spec = ~names(all_model_args)[j])
  }) 

  bind_rows(out)
}

do_scenarios <- function(nsims, scenario_nos, allocations, all_model_args,
                         compute_se, 
                         .parallel = TRUE, ...){
  
  lapply(scenario_nos, function(s){
    plyr::ldply(1:nsims, .parallel = .parallel, function(b){
      do_sim(s, allocations = allocations, all_model_args, compute_se = compute_se) %>%
       mutate_(simid =~ b,
               sid   =~ s)
    })
  })
}
