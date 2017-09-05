
est_step2 <- function(data, allocations, model_args, randomization, compute_se = TRUE, verbose = TRUE,...){  
   models <- dr::est_step0(data, model_args)
   geexList <- list(eeFUN = generic_eefun, splitdt = split(data, data$group))
  
   hold <- lapply(allocations, function(allocation){
     temp <- dr::est_step1(
        data          = data, 
        step0         =  models, 
        model_args    = model_args, 
        allocations   = allocation, 
        randomization = randomization)
     estimator_args <- temp$estimator_args
     m <- temp$models
    
    all <- lapply(estimator_args, function(eargs){
    ## Estimate parameters for each allocation

       if(eargs$skip == FALSE){
         geexList <- append(geexList, list(ee_args = list(alpha = allocation)))
         p    <- length(eargs$theta)

         n_allocation <- length(allocation)
         make_estimator_fun <- match.fun(paste0(eargs$type, '_estimator'))
         
         ## BEGIN  Point estimates ##
         target <- lapply(geexList$splitdt, function(grp_dt){
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
      
         ## END Point estimates ##
         ## BEGIN VCOV estimates ##
         if(compute_se){
           mats <- geex::compute_matrices(
             geex_list        = geexList,
             theta            = c(eargs$theta, target),
             derivFUN_control = list(method = 'simple'),
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
         n_allocation <- length(allocations)
         target <- NA
         std_error <- NA
       }

      if(verbose == TRUE){
        print(eargs$type)
        print(paste('target estimates: ', paste(round(target, 4), collapse = ' ')))
        print(paste('std_err estimates: ', paste(round(std_error, 4), collapse = ' ')))
      }
        
        # Convert to data_frame
        data_frame(
          method    = eargs$type,
          regression_type = eargs$regtyp,
          hajek     = eargs$hajek,
          a         = rep(0:1, each = n_allocation),
          alpha     = rep(allocation, times = 2),
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

do_sim <- function(scenario_no, allocations, all_model_args, compute_se, verbose = verbose, ...){
  simdt <- do.call(gen_sim, args = arg_maker(scenario_no, nsims = 1))
  on.exit(rm(simdt))
  out <- lapply(seq_along(all_model_args), function(j){
    est_step2(
      data        = simdt[[1]],
      allocations = allocations, 
      model_args  = all_model_args[[j]],
      compute_se  = compute_se,
      verbose     = verbose ) %>%
      mutate_(
        sid = ~scenario_no,
        model_spec = ~names(all_model_args)[j])
  }) 

  bind_rows(out)
}

do_scenarios <- function(nsims, scenario_nos, allocations, all_model_args,
                         compute_se, verbose,
                         .parallel = TRUE, ...){
  
  lapply(scenario_nos, function(s){
    plyr::l_ply(1:nsims, .parallel = .parallel, function(b){
      if(verbose){print(b)}
      x <- do_sim(s, allocations = allocations, all_model_args, compute_se = compute_se,
             verbose = verbose) %>%
       mutate_(simid =~ b,
               sid   =~ s)
      save(x, file = paste0('inst/experiments/ex30/ex30_results/results_', b, '.rda'))
    })
  })
}
