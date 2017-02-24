gen_data <- function(m, ni, gamma, theta, beta){
  n <- m * ni
  data_frame(
    m   = m,
    ni  = ni,
    group = rep(1:m, each = ni),
    Z1  = rnorm(n, sd = 1),
    Z1_abs = abs(Z1),
    Z2  = rbinom(n, size = 1, prob = .5),
    b   = rep(rnorm(m, sd = theta), each = ni),
    p   = as.numeric(plogis(cbind(1, Z1_abs, Z2) %*% gamma + b)),
    A   = rbinom(n, size = 1, prob = p),
    fA  = rep(tapply(A, group, mean), each = ni),
    Y   = as.numeric(cbind(1, A, fA, Z1_abs, Z2) %*% beta)
  )
}

gen_sims <- function(nsims, ...){
  dots <- list(...)
  args <- dots[pmatch(names(dots), names(formals((gen_data))))]
  x    <- as.call(append(list(gen_data), args))
  replicate(nsims, eval(x), simplify = FALSE)
}

est_sims <- function(sims, allocations, model_args, .parallel = TRUE){
  plyr::llply(seq_along(sims), .parallel = .parallel, .progress = 'text', 
              function(i){
    
    x <- sims[[i]]
    ## Component models
    m <- make_models(model_args = model_args, data = x)
  
    ## Grab component model parameter estimates
    theta_t <- unlist(lme4::getME(m$tmodel, c('beta', 'theta')))
    theta_o <- coef(m$omodel)

    temp <- list(eeFUN = dr_eefun, splitdt = split(x, x$group))
    
    ## Estimate parameters for each allocation
    hold <- lapply(allocations, function(allocation){
      temp <- append(temp, list(ee_args = list(alpha = allocation)))
      
      # Point estimates
      target <- lapply(temp$splitdt, function(grp_dt){
        # Create estimator function
        estimator <- dr_estimators(
          grp_dt, 
          t_model = m$tmodel, 
          o_model = m$omodel,
          randomization = 1)
        
        # Evaluate estimator function
        estimator(c(theta_t, theta_o), alpha = allocation)
      }) %>% 
        list_matrix() %>% 
        apply(., 2, mean)
      
      # vcov estimates
      mats <- geex::compute_matrices(
        temp,
        theta   = c(theta_t, theta_o, target),
        numDeriv_options = list(method = 'simple'),
        t_model = m$tmodel,
        o_model = m$omodel,
        randomization = 1)
      
      Ainv <- try(solve(Matrix(mats$A)), silent = TRUE)
      if(!is(Ainv, 'try-error')){
        Sigma <-Ainv %*% mats$B %*% t(Ainv)
        std_error <- sqrt(diag(Sigma)[-(1:(length(theta_t) + length(theta_o)))])
      } else {
        std_error <- NA
      }

      # Convert to data_frame
      data_frame(
        # method    = rep(c('ipw', 'otc', 'dbr', 'hjk_ipw', 'hjk_dbr'), each = 2),
        method    = rep(c('ipw', 'otc', 'dbr'), each = 2),
        a         = rep(c(0, 1), 3),
        alpha     = allocation,
        estimate  = target,
        std_error = std_error,
        simID     = i
      )
    })
    
    bind_rows(hold) # combine estimates for allocations from single simulation
  }) %>% 
    bind_rows() # combine results for all simulations
}
