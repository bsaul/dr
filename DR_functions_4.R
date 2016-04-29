#------------------------------------------------------------------------------#
#   Title: Doubly Robust Functions
#  Author: B. Saul
#    Date: 2016-04-28
# Purpose: functions for IPW and DR estimators
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
#### Utilities ####
#------------------------------------------------------------------------------#
get_fixed_formula <- function(model_object){
  formula(model_object, fixed.only = TRUE)[-2]
}

get_design_frame <- function(rhs_formula, data){
  as.data.frame(model.matrix(rhs_formula, data))
}

get_response <- function(formula, data){
  model.response(model.frame(formula, data = data))
}

list_matrix <- function(this_list)
{
  ulist <- unlist(this_list)
  m <- length(this_list)
  p <- length(ulist)/m
  matrix(ulist, nrow = m, ncol = p, byrow = T)
}

rbind_fill_zero <- function(this_list)
{
  if(!is.list(this_list)){
    this_list <- list(this_list)
  }
  p <- max(unlist(lapply(this_list, ncol)))
  out <- lapply(this_list, function(x){
    xp <- ncol(x)
    if(xp < p ){
      cbind(x,matrix(0, nrow = nrow(x), ncol = p - xp))
    } else {
      x
    }
  })
  do.call('rbind', out)
}

evaluate_df_function <- function(f, ...){
  f(...)
}

#------------------------------------------------------------------------------#
#### Internals ####
# functions used in estimators
#------------------------------------------------------------------------------#

integrand <- function(b, response, xmatrix, theta){
  lp <- outer(xmatrix %*% theta[-length(theta)], b, '+')
  h  <- apply(lp, 3, function(x) dbinom(response, 1, plogis(x) ) )
  hh <- apply(h, 2, prod)
  hh * dnorm(b, mean = 0, sd = theta[length(theta)])
}

pi_term <- function(A){
  f <- function(alpha){
    prod(dbinom(x = A, 1, prob = alpha))
  }
  memoise::memoise(f)
}

grids <- new.env()

expand_grid_n <- function(n)
{
  lookup <-as.character(n)
  if(!exists(lookup, envir = grids)){
    grids[[lookup]] <- expand.grid(ID = 1:n, sum_a = 0:n, A = 0:1)
  }
  return(grids[[lookup]])
}

expand_grid_n2 <- function(n1, n2, a)
{
  lookup <- paste(n1, n2, a, sep = '_')
  if(!exists(lookup, envir = grids)){
    grids[[lookup]] <- expand.grid(ID = 1:n1, sum_a = 0:n2, 
                                   A = {if(is.null(a)) 0:1 else a })
  }
  return(grids[[lookup]])
}


# microbenchmark::microbenchmark(
#   expand_grid_n(100),
#   expand.grid(ID = 1:100, sum_a = 0:100, A = 0:1))

expand_outcome_frame <- function(X_outcome, formula_outcome){
  new_grid <- expand_grid_n(nrow(X_outcome))
  f <- function(theta){
    X_expanded <- X_outcome %>%
      mutate(ID = row_number()) %>%
      select(-A) %>%
      full_join(new_grid, by = 'ID')  
    
    new_frame <- get_design_frame(formula_outcome, X_expanded) 
    fit_func <- make_dr_term1(new_frame)

    X_expanded <- X_expanded %>%
      ungroup() %>%
      mutate_(fitted_value = ~fit_func(theta))
  }
 memoise::memoise(f)
}

#------------------------------------------------------------------------------#
#### IPW Estimator ####
#------------------------------------------------------------------------------#

weight_estimator <- function(A, X)
{
  X <- as.matrix(X)
  f <- function(theta){
    1/integrate(integrand, lower = -Inf, upper = Inf,
                theta = theta, response = A, xmatrix = X)$value
  }
  memoise::memoise(f)
}

make_ipw_estimator <- function(Y, A, X_treatment, ...){
  w <- weight_estimator(A = A, X = X_treatment)
  pi_t <- pi_term(A = A)
  
  function(theta, alpha, a = NULL){
    Ia <- if(is.null(a)) 1 else (A == a) * 1
    mean(Y * Ia) * w(theta) * pi_t(alpha)  / (dbinom(a, 1, alpha) * !is.null(a) * 1)  
  }
}

#------------------------------------------------------------------------------#
#### Outcome Estimator ####
#------------------------------------------------------------------------------#

make_otc_estimator <- function(X_outcome, rhs_formula_outcome, ...){
  function(theta, alpha, a = NULL){
    n1 <- nrow(X_outcome) 
    n2 <- n1 - {if(is.null(a)) 0 else 1}
    
    X_outcome %>%
      mutate(ID = row_number()) %>%
      # Remove all treatment variables so that model.matrix correctly generates
      # interactions
      select(-contains('A')) %>%
      # Generate the relevant values of A and 
      # all possible sum(a_i) for each subject
      full_join(expand_grid_n2(n1, n1 - 1, a), by = "ID") %>%
      # group_by_(~ID, ~A) %>%
      # Compute fA and pi for each sum(a_i) a
      mutate_(fA = ~ sum_a / n1,
              pi_value = ~ dbinom(sum_a, n2, prob = alpha) ) %>%
      ungroup() %>%
      # Compute fitted values
      mutate_(fits = ~ as.numeric(model.matrix(rhs_formula_outcome, data = .) %*% theta ) ) %>%
      # Sum by individual to compute term2
      group_by_(~A, ~ID) %>%
      summarize_(estimate = ~sum(fits * pi_value)) %>%
      group_by_(~A) %>%
      summarize_(estimate = ~mean(term2)) %>%
      summarize_(estimate = ~mean(term2)) %>%
      as.numeric()
  }
}

#------------------------------------------------------------------------------#
#### Doubly Robust Estimator ####
#------------------------------------------------------------------------------#

make_dr_term1 <- function(X){
  X <- as.matrix(X)
  f <- function(theta){
    X %*% theta
  }
  memoise::memoise(f)
}

make_dbr_estimator <- function(Y, A, X_outcome, X_treatment, formula_outcome){
  w <- weight_estimator(A = A, X = X_treatment)
  pi_t <- pi_term(A = A)
  dr_term1 <- make_dr_term1(X_outcome)
  dr_term2 <- make_outcome_estimator(X_outcome, formula_outcome = formula_outcome)
  
  q_treatment <- ncol(X_treatment) + 1

  function(theta, alpha, a = NULL){
    Ia <- if(is.null(a)) 1 else (A == a) * 1
    Ybar <- mean(Ia * (Y - dr_term1(theta[(q_treatment + 1):length(theta)] ) ) )
    term1 <- Ybar * w(theta[1:q_treatment]) * pi_t(alpha) / (dbinom(a, 1, alpha) * !is.null(a) * 1)
    term2 <- dr_term2(theta[(q_treatment + 1):length(theta)], alpha, a)
    term1 + term2
  }
}

#------------------------------------------------------------------------------#
#### Estimating Equations ####
#------------------------------------------------------------------------------#

psi <- function(estimator){
  function(theta, alpha, a){
    q <- length(theta)
    estimator(theta[1:(q-1)], alpha, a) - theta[q]
  }
}

#------------------------------------------------------------------------------#
#### Estimation ####
#------------------------------------------------------------------------------#

estimation <- function(treatment_formula, 
                       outcome_formula,
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
      method  = geepack::geeglm,
      options = list(family = gaussian, id = quote(group))
    ) )
  
  models <- sandwichShop::make_models(model_args = model_args, data = this_data)
  
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
  
#   #### Estimating Equations for each nuisance model ####  
#   ee_t <- estfun(models$model_treatment)
#   ee_o <- estfun(models$model_outcome)
#   
#   #### Derivatives of Estimating Eqns for each nuisance model (U matrix) ####  
#   bread_t <- bread(models$model_treatment)
#   bread_o <- bread(models$model_outcome)
  
  #### Create list of group level data  ####
  split_data <- split(this_data, this_data$group, drop = FALSE)
  
  #### Create data_frame of estimator functions ####
  frame <- lapply(split_data, function(group_data) {
    Y   <- get_response(form_o, group_data)
    A   <- get_response(form_t, group_data)
    X_t <- get_design_frame(form_rhs_t, group_data)
    X_o <- get_design_frame(form_rhs_o, group_data)
    
    estimators <- list(make_ipw_estimator(Y = Y, A = A, X_treatment = X_t),
                       make_otc_estimator(X_outcome = X_o, 
                                         formula_outcome = form_rhs_o),
                       make_dbr_estimator(Y = Y, A = A, 
                                                  X_outcome = X_o, 
                                                  X_treatment = X_t, 
                                                  formula_outcome = form_rhs_o))
    data_frame(estimator_type = c('ipw', 'otc', 'dbr'),
               estimator = estimators)
  }) %>% bind_rows() 
  
  frame %<>% 
    left_join(theta, by = 'estimator_type') %>%
    left_join(counterfactual_grid, by = 'estimator_type')
  #### Evaluate group-level estimators ####
  
  frame %<>%
    group_by(estimator_type) %>%
    rowwise() %>%
    mutate(estimate = evaluate_df_function(estimator, theta = theta, alpha = alpha, a = a))

  #### Evaluate population targets ####
  
  #### Evaluate psi functions ####
  # psis <- 
  #### Evaluate psi' functions ####
  
  #### Evaluate variance estimates ####
  
  return(frame)
}
