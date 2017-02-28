#------------------------------------------------------------------------------#
#   Title: Experiment 01: speeding up the outcome estimator functions
#  Author: B. Saul
#    Date: 2016-04-30
# Purpose: develop faster outcome estimators
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
#### Utilities ####
#------------------------------------------------------------------------------#
get_fixed_formula_00 <- function(model_object){
  formula(model_object, fixed.only = TRUE)[-2]
}

get_design_frame_00 <- function(rhs_formula, data){
  as.data.frame(model.matrix(rhs_formula, data))
}

get_response_00 <- function(formula, data){
  model.response(model.frame(formula, data = data))
}

#------------------------------------------------------------------------------#
#### Internals - copied from DR_functions_4.R ####
# functions used in estimators
#------------------------------------------------------------------------------#

grids <- new.env()

expand_grid_n_00 <- function(n1, n2)
{
  lookup <- paste(n1, n2, sep = '_')
  if(!exists(lookup, envir = grids)){
    grids[[lookup]] <- expand.grid(ID = 1:n1, sum_a = 0:n2, 
                                   A = 0:1) %>%
      mutate_(fA = ~ sum_a / n1) 
  }
  return(grids[[lookup]])
}


expand_outcome_frame_00 <- function(X_outcome, rhs_formula_outcome){
  f <- function(theta){
    n <- nrow(X_outcome)
    X_outcome %>%
      mutate_(ID = ~ row_number()) %>%
      # Remove treatment variables so that they are updated 
      select(-A, -fA) %>%
      # Generate the relevant values of A and 
      # all possible sum(a_i) for each subject
      full_join(expand_grid_n_00(n, n - 1), by = "ID") %>%
      # Compute fitted values
      mutate_(fits = ~ as.numeric(model.matrix(rhs_formula_outcome, data = .) %*% theta ) ) 
  }
  memoise::memoise(f)
}

expand_outcome_frame_01 <- function(X_outcome, rhs_formula_outcome){
  n <- nrow(X_outcome)
  X_outcome %>%
    mutate_(ID = ~ row_number()) %>%
    # Remove treatment variables so that they are updated 
    select(-A, -fA) %>%
    # Generate the relevant values of A and 
    # all possible sum(a_i) for each subject
    full_join(expand_grid_n_00(n, n - 1), by = "ID") 
}
# expand_outcome_frame_01 <- memoise::memoise(expand_outcome_frame_01)

#------------------------------------------------------------------------------#
#### Outcome Estimator ####
#------------------------------------------------------------------------------#

make_otc_estimator_00 <- function(X_outcome, rhs_formula_outcome, ...){
  
  X_expanded <- expand_outcome_frame_00(X_outcome, rhs_formula_outcome)
  
  function(theta, alpha, a = NULL){
    n <- nrow(X_outcome) - {if(is.null(a)) 0 else 1}
    a <- if(is.null(a)) 0:1 else a
    
    X_expanded(theta) %>%
      filter_(~ A %in% a) %>%
      mutate_(pi_value = ~ dbinom(sum_a, n, prob = alpha)) %>%
      # Sum by individual to compute term2
      group_by_(~A, ~ID) %>%
      summarize_(estimate = ~sum(fits * pi_value)) %>%
      group_by_(~A) %>%
      summarize_(estimate = ~mean(estimate)) %>%
      summarize_(estimate = ~mean(estimate)) %>%
      as.numeric()
  }
}

make_otc_estimator_01 <- function(X_outcome, rhs_formula_outcome, ...){
  
  X_expanded <- expand_outcome_frame_01(X_outcome, rhs_formula_outcome)
  
  function(theta, alpha, a = NULL){
    n <- nrow(X_outcome) - {if(is.null(a)) 0 else 1}
    a <- if(is.null(a)) 0:1 else a
    
    X_expanded %>%
      filter_(~ A %in% a) %>%
      mutate_(pi_value = ~ dbinom(sum_a, n, prob = alpha),
              fits = ~ as.numeric(model.matrix(rhs_formula_outcome, data = .) %*% theta )) %>%
      # Sum by individual to compute term2
      group_by_(~A, ~ID) %>%
      summarize_(estimate = ~sum(fits * pi_value)) %>%
      group_by_(~A) %>%
      summarize_(estimate = ~mean(estimate)) %>%
      summarize_(estimate = ~mean(estimate)) %>%
      as.numeric()
  }
}

make_otc_estimator_02 <- function(X_outcome, rhs_formula_outcome, ...){
  
  X_expanded <- expand_outcome_frame_01(X_outcome, rhs_formula_outcome)
  
  function(theta, alpha, a = NULL){
    n <- nrow(X_outcome) - {if(is.null(a)) 0 else 1}
    a <- if(is.null(a)) 0:1 else a
    
    x <- X_expanded %>%
      filter_(~ A %in% a) %>%
      mutate_(pi_value = ~ dbinom(sum_a, n, prob = alpha),
              fitted   = ~ as.numeric(model.matrix(rhs_formula_outcome, data = .) %*% theta ),
              estimate = ~ fitted * pi_value) 
    
    # This is returns an incorrect value
    sum(tapply(x$estimate, paste(x$A, x$ID), sum)) / n  
  }
}

make_otc_estimator_03 <- function(X_outcome, rhs_formula_outcome, ...){
  
  X_expanded <- expand_outcome_frame_01(X_outcome, rhs_formula_outcome)
  
  function(theta, alpha, a = NULL){
    n1 <- nrow(X_outcome)
    n2 <- n1 - {if(is.null(a)) 0 else 1}
    a <- if(is.null(a)) 0:1 else a
    
    x <- X_expanded %>%
      filter_(~ A %in% a) %>%
      mutate_(pi_value = ~ dbinom(sum_a, n2, prob = alpha),
              fitted   = ~ as.numeric(model.matrix(rhs_formula_outcome, data = .) %*% theta ),
              estimate = ~ fitted * pi_value) 
    
    sum(tapply(tapply(x$estimate, paste(x$A, x$ID), sum), 
           rep(a, each = n1), sum)) / n1
  }
}

#------------------------------------------------------------------------------#
                          #### Begin Experiments ####

#------------------------------------------------------------------------------#
#### Experiment 01_00 ####
# profile functions to see where bottle necks are
#------------------------------------------------------------------------------#

ex01_00_1 <- function(){
  X <- get_design_frame_00(rhs_outcome, g1_dt)
  f <- make_otc_estimator_00(X, rhs_outcome)
  f(theta_o, alpha = .5, a = 1)
}

ex01_00_2 <- function(){
  lapply(split_dt, function(group_data){
    X <- get_design_frame_00(rhs_outcome, group_data)
    f <- make_otc_estimator_00(X, rhs_outcome)
    f(theta_o, alpha = .5, a = 1)
  })
}

ex01_00_3 <- function(){
  X <- get_design_frame_00(rhs_outcome, g1_dt)
  f <- make_otc_estimator_00(X, rhs_outcome)
  numDeriv::grad(f, x = theta_o, alpha = .5, a = 1)
}

#------------------------------------------------------------------------------#
#### Experiment 01_01 ####
# profile functions to see where bottle necks are
#------------------------------------------------------------------------------#

ex01_01_1 <- function(){
  X <- get_design_frame_00(rhs_outcome, g1_dt)
  f <- make_otc_estimator_01(X, rhs_outcome)
  f(theta_o, alpha = .5, a = 1)
}

ex01_01_2 <- function(){
  lapply(split_dt, function(group_data){
    X <- get_design_frame_00(rhs_outcome, group_data)
    f <- make_otc_estimator_01(X, rhs_outcome)
    f(theta_o, alpha = .5, a = 1)
  })
}

ex01_01_3 <- function(){
  X <- get_design_frame_00(rhs_outcome, g1_dt)
  f <- make_otc_estimator_01(X, rhs_outcome)
  numDeriv::grad(f, x = theta_o, alpha = .5, a = 1)
}

#------------------------------------------------------------------------------#
#### Experiment 01_02 ####
# use base R functions to take means/sums
#------------------------------------------------------------------------------#

ex01_02_1 <- function(){
  X <- get_design_frame_00(rhs_outcome, g1_dt)
  f <- make_otc_estimator_02(X, rhs_outcome)
  f(theta_o, alpha = .5, a = 1)
}

ex01_02_2 <- function(){
  lapply(split_dt, function(group_data){
    X <- get_design_frame_00(rhs_outcome, group_data)
    f <- make_otc_estimator_02(X, rhs_outcome)
    f(theta_o, alpha = .5, a = 1)
  })
}

ex01_02_3 <- function(){
  X <- get_design_frame_00(rhs_outcome, g1_dt)
  f <- make_otc_estimator_02(X, rhs_outcome)
  numDeriv::grad(f, x = theta_o, alpha = .5, a = 1)
}


#------------------------------------------------------------------------------#
#### Experiment 01_03 ####
# Follow up on 02 to get accurate values
#------------------------------------------------------------------------------#
ex01_03_1 <- function(){
  X <- get_design_frame_00(rhs_outcome, g1_dt)
  f <- make_otc_estimator_03(X, rhs_outcome)
  f(theta_o, alpha = .5, a = 1)
}

ex01_03_2 <- function(){
  lapply(split_dt, function(group_data){
    X <- get_design_frame_00(rhs_outcome, group_data)
    f <- make_otc_estimator_03(X, rhs_outcome)
    f(theta_o, alpha = .5, a = 1)
  })
}

ex01_03_3 <- function(){
  X <- get_design_frame_00(rhs_outcome, g1_dt)
  f <- make_otc_estimator_03(X, rhs_outcome)
  numDeriv::grad(f, x = theta_o, alpha = .5, a = 1)
}
