#------------------------------------------------------------------------------#
#' 
#' @export
#------------------------------------------------------------------------------#

estimation <- function(treatment_formula, 
                       outcome_formula,
                       outcome_method = geepack::geeglm,
                       outcome_method_opts = list(family = gaussian, id = quote(group)),
                       this_data,
                       allocations,
                       target_a)
{
  
  estimator_types <-  c('ipw', 'otc', 'dbr')
  counterfactual_grid <- expand.grid(estimator_type = estimator_types, 
                                     alpha = allocations, a = target_a,
                                     stringsAsFactors = F)
  #### Fit nuisance models ####
  model_args = list(
    model_treatment = list(
      formula = treatment_formula,
      method  = lme4::glmer,
      options = list(family = binomial)
    ) ,
    model_outcome = list(
      formula = outcome_formula,
      method  = outcome_method,
      options = outcome_method_opts
    ) )
  
  print('making models')
  models <- sandwichShop::make_models(model_args = model_args, data = this_data)
  
  treatment_warnings <- model_warnings(models$model_treatment)
  outcome_warnings   <- model_warnings(models$model_outcome)
  
  print('models made')
  
  #### Grab nuisance parameter estimates ####
  theta_t <- unlist(lme4::getME(models$model_treatment, c('beta', 'theta')))
  theta_o <- coef(models$model_outcome)
  theta_a <- c(theta_t, theta_o)
  # theta   <- list(ipw = theta_t, out = theta_o, dr = theta_a)
  theta <- data_frame(estimator_type = estimator_types,
                      theta = list(theta_t, theta_o, theta_a))
  
  #### Grab full formulas ####
  form_t <- formula(models$model_treatment)
  form_o <- formula(models$model_outcome)
  
  #### Grab rhs formulas ####
  form_rhs_t <- get_fixed_formula(models$model_treatment)
  form_rhs_o <- get_fixed_formula(models$model_outcome)
  
  #### Outcome Inv Link ####
  otc_inv_link <- models$model_outcome$family$linkinv
  
  #### Estimating Equations for each nuisance model ####  
  print('making estimating eqns')
  ee_t <- estfun(models$model_treatment, grad_method = 'simple')
  ee_o <- estfun(models$model_outcome, grad_method = 'simple')
  print('estimating equations made')
  
  #### Derivatives of Estimating Eqns for each nuisance model (U matrix) ####  
  #   bread_t <- bread(models$model_treatment)
  #   bread_o <- bread(models$model_outcome)
  
  #### Create list of group level data  ####
  split_data <- split(this_data, this_data$group, drop = FALSE)
  m <- length(split_data)
  
  print('forming frame')
  #### Create data_frame of estimator functions ####
  frame <- lapply(split_data, function(group_data) {
    Y   <- geex::get_response(form_o, group_data)
    A   <- geex::get_response(form_t, group_data)
    X_t <- geex::get_design_frame(form_rhs_t, group_data)
    X_o <- geex::get_design_frame(form_rhs_o, group_data)
    
    estimators <- list(make_ipw_estimator(Y = Y, A = A, X_treatment = X_t),
                       make_otc_estimator(X_outcome = X_o, 
                                          rhs_formula_outcome = form_rhs_o,
                                          inv_link = otc_inv_link),
                       make_dbr_estimator(Y = Y, A = A, 
                                          X_outcome = X_o, 
                                          X_treatment = X_t, 
                                          rhs_formula_outcome = form_rhs_o,
                                          inv_link = otc_inv_link))
    U21_funcs <- list(make_ipw_U21_estimator(Y = Y, A = A, X_treatment = X_t, theta_t = theta_t),
                      make_otc_U21_estimator(A = A,
                                             X_outcome = X_o,
                                             rhs_formula_outcome = form_rhs_o,
                                             theta_o = theta_o,
                                             inv_link = otc_inv_link),
                      make_dbr_U21_estimator(Y = Y, A = A,
                                             X_outcome = X_o,
                                             X_treatment = X_t,
                                             rhs_formula_outcome = form_rhs_o,
                                             theta_t = theta_t,
                                             theta_o = theta_o,
                                             inv_link = otc_inv_link))
    data_frame(n_i = length(Y),
               estimator_type = c('ipw', 'otc', 'dbr'),
               estimator = estimators,
               U21_func  = U21_funcs)
  }) %>% bind_rows() 
  
  frame %<>% 
    left_join(theta, by = 'estimator_type') %>%
    left_join(counterfactual_grid, by = 'estimator_type')
  
  #### Evaluate group-level estimators ####
  frame %<>%
    rowwise() %>%
    mutate(estimate = evaluate_df_function(estimator, theta = theta, alpha = alpha, a = a),
           U21      = list(evaluate_df_function(U21_func, alpha = alpha, a = a) ),
           psi_func = list(evaluate_df_function(psi, estimator)) ) %>%
    ungroup() %>%
    group_by(estimator_type, alpha, a) %>%
    mutate_(target_estimate = ~ sum(estimate, na.rm = T)/m,
            psi_val  = ~ estimate - target_estimate) %>%
    ungroup() 
  
  print('frame created')
  print('evaluating results')
  results <- plyr::ddply(frame, plyr::.(estimator_type, alpha, a), .progress = 'text', function(x){
    if(x$estimator_type[1] == 'ipw'){
      ee <- cbind(ee_t, x$psi_val)
      # bb <- bread_t
    } else if(x$estimator_type[1] == 'otc') {
      ee <- cbind(ee_o, x$psi_val)
      # bb <- bread_o
    } else if(x$estimator_type[1] == 'dbr') {
      ee <- cbind(ee_t, ee_o, x$psi_val)
      # bb <- Matrix::bdiag(bread_t, bread_o)
    }
    
    V <- crossprod(ee)/nrow(ee)
    #     U2 <- c(apply(list_matrix(x$U21), 2, mean), 1)
    #     U <- rbind_fill_zero(list(bb, matrix(-U2, ncol = length(U2))))
    #     sigma <- (solve(U) %*% V %*% t(solve(U)))/nrow(ee)
    #     std.error <- sqrt(sigma[nrow(sigma), ncol(sigma)])
    
    U21 <- apply(list_matrix(x$U21), 2, mean)
    U21 <- matrix(-U21, ncol = length(U21))
    V11 <- V[1:(nrow(V) - 1), 1:(ncol(V) - 1)]
    V21 <- V[nrow(V), 1:(ncol(V) - 1)]
    V22 <- V[nrow(V), ncol(V)]
    var.est <- ( ( (U21 - 2 * V21) %*% solve(V11) %*% t(U21) ) + V22) / nrow(ee)
    std.error <- as.numeric(sqrt(var.est))
    # std.error <- NA
    df_out <- data_frame(estimate = x$target[1], 
                         std.error = std.error,
                         treatment_warnings = treatment_warnings,
                         outcome_warnings   = outcome_warnings)
    return(df_out)
  })
  
  #   return(results)
  results
  # frame
}

#------------------------------------------------------------------------------#
#' 
#' @export
#------------------------------------------------------------------------------#


estimate_sims <- function(sims,
                          formula_treatment,
                          formula_outcome,
                          progress = 'text',
                          parallel = FALSE)
{
  plyr::ddply(sims, plyr::.(simID), .progress = progress, .parallel = parallel,
              function(x){
                m <- try(estimation(treatment_formula = formula_treatment, 
                                    outcome_formula   = formula_outcome, 
                                    this_data =  x,
                                    allocations = alphas, 
                                    target_a = 1))
                # print(m)
                if(is(m,  'try-error')){
                  data_frame(error = 'error')
                } else{
                  out <- m
                  out$simID <- x$simID[1]
                  out$error <- ''
                  out
                }
              })
}

