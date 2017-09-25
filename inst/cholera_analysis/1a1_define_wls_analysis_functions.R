#------------------------------------------------------------------------------#
#   Title: Functions for analysis of cholera data with IPW, OTC,
#          and DR estimators
#  Author: B. Saul
#    Date: 2017-02-27
# Purpose: 
#------------------------------------------------------------------------------#
wls_estFUN <- function(data, wls_formula, ipw_formula, invlnk, randomization){
  
  L <- grab_design_matrix(data, wls_formula)
  L_ex <- expand_outcome_frame(data, geex::grab_fixed_formula(wls_formula))
  L_ex_0 <- L_ex %>% filter(A == 0)
  L_ex_1 <- L_ex %>% filter(A == 1)
  MM_0 <- grab_design_matrix(L_ex_0, wls_formula)
  MM_1 <- grab_design_matrix(L_ex_1, wls_formula)
  X <- grab_design_matrix(data, ipw_formula)
  Y <- grab_response(data, wls_formula)
  A <- grab_response(data, ipw_formula)
  n <- length(Y)
  p <- ncol(L)

  ipwFUN <- weight_estimator(
    A = A, 
    X = X, 
    randomization = randomization)
  
  function(theta, alpha, ipw_theta){
    nalpha <- length(alpha)
    index0 <- 1:p
    index1 <- (p*nalpha + 1):(p*nalpha + p)
    
    wls0_ee <- lapply(alpha, function(x){
      ipw <- ipwFUN(ipw_theta, alpha = x)/(1 - x)
      W  <- diag(ipw * (A == 0), nrow = n, ncol = n)
      ee <- t(L) %*% W %*% (Y - invlnk(L %*% theta[index0]))
      
      index0 <<- index0 + p
      return(ee)
    }) %>% unlist()
    
    wls1_ee <- lapply(alpha, function(x){
      ipw <- ipwFUN(ipw_theta, alpha = x)/(x)
      W  <- diag(ipw * (A == 1), nrow = n, ncol = n)
      ee <- t(L) %*% W %*% (Y - invlnk(L %*% theta[index1]))
      index1 <<- index1 + p
      return(ee)
    }) %>% unlist()
    
    ce0 <- ce1 <- numeric(nalpha)
    ### Regression-based DRR estimator ###
    for(k in 1:length(alpha)){
      if(k == 1){
        index0 <- 1:p
        index1 <- (p*nalpha + 1):(p*nalpha + p)
      }
      if(k > 1){
        index0 <- index0 + p
        index1 <- index1 + p
      }
      
      # compute fitted value for expanded data.frame
      mu_0 <- as.numeric(invlnk(MM_0 %*% theta[index0]))
      mu_1 <- as.numeric(invlnk(MM_1 %*% theta[index1]))
      # compute pi term per number treated in group per subject
      pi_term_a <- dbinom(L_ex_0$sum_a, n - 1, alpha[k])
      
      # mulitply mu_ij by the pi term rowwise
      piXmu_a_0 <- mu_0 * pi_term_a
      piXmu_a_1 <- mu_1 * pi_term_a
      
      # sum within levels of A (0:1) WITHIN subjects
      piXmu_a <- tapply(
        X     = c(piXmu_a_0, piXmu_a_1), 
        INDEX = paste(rep(0:1, each = nrow(L_ex_0)), c(L_ex_0$ID, L_ex_0$ID)), 
        FUN   = sum)
      
      # sum within levels of A (0:1) ACROSS subjects
      wls_ce_a <- tapply(
        X     = piXmu_a, 
        INDEX = rep(0:1, each = n), 
        FUN   = sum)
      
      ce0[k] <- wls_ce_a[1]/n
      ce1[k] <- wls_ce_a[2]/n
    }
    
    ntheta <- length(theta)
    
    c(wls0_ee, wls1_ee, 
      ce0 - theta[(ntheta - 2*nalpha + 1):(ntheta - nalpha)], 
      ce1 - theta[(ntheta - nalpha + 1):ntheta])
    
  }
}

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
        ipw0 =~ ipw * (A == 0)/(1 - allocations[k]),
        ipw1 =~ ipw * (A == 1)/allocations[k]
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
  })
 
  theta_wls <- c(theta_t,
                 unlist(lapply(m$wls_model_0, coef)),
                 unlist(lapply(m$wls_model_1, coef))) 

  # theta_wls <- c(unlist(lapply(m$wls_model_0, function(x) c(theta_t, coef(x)))),
  #                unlist(lapply(m$wls_model_1, function(x) c(theta_t, coef(x))))) 
                        
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
               regtyp    = 'TRUE',
               skipit    = FALSE),
    wls_dbr = list(type    = 'wls_dbr',
                   theta   =  theta_wls,
                   hajek   = FALSE,
                   regtyp  = 'wls',
                   skipit  = FALSE)
  )
  
  list(
     estimator_args = estimator_args,
     models         = m)
}


estimate_cholera_parms_step2_wls <- function(data, allocations, models, model_args, randomization, compute_se = TRUE, ...){  

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
        # geexList <- append(geexList, list(ee_args = list(alpha = allocation)))
        p            <- length(eargs$theta)
        n_allocation <- length(allocation)


        ## BEGIN  Point estimates ##
        # start wls targets at DBR estimates
        wls_target_start <- lapply(geexList$splitdt, function(grp_dt){
          # Create estimator function
          estimator <- dbr_estimator(
            data          = grp_dt,
            models        = m,
            randomization = randomization,
            hajek         = eargs$hajek,
            regression_type = eargs$regtyp)
          # Evaluate estimator function
          estimator(eargs$theta, alpha = allocation)
        })
        

        theta_t <- unlist(lme4::getME(m$t_model, c('beta', 'theta')))
        wls_target_start <- wls_target_start %>% list_matrix() %>% apply(., 2, mean)
        
        wls_beta_start <- glm(y_obs ~ fA + age + rivkm, data = choleradt, 
                         family = binomial) %>% coef()
        
        wls_start <- c(rep(wls_beta_start, times = 4), wls_target_start)
        
        print("starting geex")
        ## END Point estimates ##
        ## BEGIN VCOV estimates ##
        if(compute_se){
          geexOut <- geex::m_estimate(
            estFUN = wls_estFUN,
            data   = data,
            units  = 'group',
            root_control = setup_root_control(start = wls_start),
            deriv_control = setup_deriv_control(method = "simple"),
            compute_roots = TRUE,
            outer_args = list(randomization    = randomization,
                              wls_formula = y_obs ~ fA + age + rivkm,
                              ipw_formula = B ~ age + rivkm,
                              invlnk = plogis),
            inner_args = list(alpha = allocation, ipw_theta = theta_t)
          )
          
          target <- roots(geexOut)
          Sigma <- vcov(geexOut)
          Ai    <- geexOut@sandwich_components@.A_i
          Bi    <- geexOut@sandwich_components@.B_i
          
          print(target)
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
        Ai <- NA
        Bi <- NA
      }

      ## END VCOV estimates ##

      list(
        alpha    = allocation,
        estimate = target,
        vcov     = Sigma,
        Ai       = Ai,
        Bi       = Bi
      )

    }) # END lapply per estimator
    all
  }) # END lapply per allocation
  hold
}
