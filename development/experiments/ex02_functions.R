#------------------------------------------------------------------------------#
#   Title: Experiment 02: speeding up the dr estimator functions
#  Author: B. Saul
#    Date: 2016-04-30
# Purpose: develop faster dr estimators
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

weight_estimator <- function(A, X)
{
  X <- as.matrix(X)
  f <- function(theta){
    1/integrate(integrand, lower = -Inf, upper = Inf,
                theta = theta, response = A, xmatrix = X)$value
  }
  memoise::memoise(f)
}

make_dr_term1 <- function(X){
  X <- as.matrix(X)
  f <- function(theta){
    X %*% theta
  }
  memoise::memoise(f)
}

expand_outcome_frame <- function(X_outcome, rhs_formula_outcome){
  f <- function(theta){
    n <- nrow(X_outcome)
    X_outcome %>%
      mutate_(ID = ~ row_number()) %>%
      # Remove treatment variables so that they are updated 
      select(-A, -fA) %>%
      # Generate the relevant values of A and 
      # all possible sum(a_i) for each subject
      full_join(expand_grid_n(n, n - 1), by = "ID") %>%
      # Compute fitted values
      mutate_(fits = ~ as.numeric(model.matrix(rhs_formula_outcome, data = .) %*% theta ) ) 
  }
  memoise::memoise(f)
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
#### Doubly Robust Estimator ####
#------------------------------------------------------------------------------#




make_dbr_estimator_00 <- function(Y, A, X_outcome, X_treatment, rhs_formula_outcome){
  w <- weight_estimator(A = A, X = X_treatment)
  pi_t <- pi_term(A = A)
  dr_term1 <- make_dr_term1(X_outcome)
  dr_term2 <- make_otc_estimator_03(X_outcome, rhs_formula_outcome = rhs_formula_outcome)
  
  q_treatment <- ncol(X_treatment) + 1
  
  function(theta, alpha, a = NULL){
    Ia <- if(is.null(a)) 1 else (A == a) * 1
    Ybar <- mean(Ia * (Y - dr_term1(theta[(q_treatment + 1):length(theta)] ) ) )
    term1 <- Ybar * w(theta[1:q_treatment]) * pi_t(alpha) / (dbinom(a, 1, alpha) * !is.null(a) * 1)
    term2 <- dr_term2(theta[(q_treatment + 1):length(theta)], alpha, a)
    term1 + term2
  }
}


make_dbr_U21_estimator_01 <- function(Y, A, X_outcome, 
                                      X_treatment, 
                                      rhs_formula_outcome){
  w <- weight_estimator(A = A, X = X_treatment)
  ws <- w(theta_t)
  wd <- numDeriv::grad(w, x = theta_t)
  pi_t <- pi_term(A = A)
  dr_term1 <- make_dr_term1(X_outcome)(theta_o)
  n <- nrow(X_outcome)
  function(alpha, a = NULL){
    pi_t1 <- pi_t(alpha) / {if(!is.null(a)) dbinom(a, 1, alpha) else 1}
    # U21 corresponding to treatment parameters
    Ia <- if(is.null(a)) 1 else (A == a) * 1
    Ybar <- mean(Ia * (Y - dr_term1) ) 
    U21_t <- Ybar * wd * pi_t1
    
    # U21 corresponding to  outcome parameters
    xmat <- model.matrix(rhs_formula_outcome, data = X_outcome)
    tt <- 1 - (Ia * pi_t1 * ws)
    
    U21_o <- apply(xmat, 2, function(col) {
     sum(col * tt)/n
    })
    
    c(U21_t, U21_o)
  }
}


make_dbr_U21_estimator_02 <- function(Y, A, X_outcome, W, Wd, rhs_formula_outcome){
  pi_t <- pi_term(A = A)
  dr_term1 <- make_dr_term1(X_outcome)(theta_o)
  n <- nrow(X_outcome)
  function(alpha, a = NULL){
    pi_t1 <- pi_t(alpha) / {if(!is.null(a)) dbinom(a, 1, alpha) else 1}
    # U21 corresponding to treatment parameters
    Ia <- if(is.null(a)) 1 else (A == a) * 1
    Ybar <- mean(Ia * (Y - dr_term1) ) 
    U21_t <- Ybar * Wd * pi_t1
    
    # U21 corresponding to  outcome parameters
    xmat <- model.matrix(rhs_formula_outcome, data = X_outcome)
    tt <- 1 - (Ia * pi_t1 * W)
    
    U21_o <- apply(xmat, 2, function(col) {
      sum(col * tt)/n
    })
    
    c(U21_t, U21_o)
  }
}


make_dbr_U21_estimator_03 <- function(Y, A, X_outcome, X_treatment, rhs_formula_outcome){
  w <- weight_estimator(A = A, X = X_treatment)
  ws <- w(theta_t)
  wd <- numDeriv::grad(w, x = theta_t, method = 'simple')
  pi_t <- pi_term(A = A)
  dr_term1 <- make_dr_term1(X_outcome)(theta_o)
  n <- nrow(X_outcome)
  function(alpha, a = NULL){
    pi_t1 <- pi_t(alpha) / {if(!is.null(a)) dbinom(a, 1, alpha) else 1}
    # U21 corresponding to treatment parameters
    Ia <- if(is.null(a)) 1 else (A == a) * 1
    Ybar <- mean(Ia * (Y - dr_term1) ) 
    U21_t <- Ybar * wd * pi_t1
    
    # U21 corresponding to  outcome parameters
    xmat <- model.matrix(rhs_formula_outcome, data = X_outcome)
    tt <- 1 - (Ia * pi_t1 * ws)
    
    U21_o <- apply(xmat, 2, function(col) {
      sum(col * tt)/n
    })
    
    c(U21_t, U21_o)
  }
}
#------------------------------------------------------------------------------#
                          #### Begin Experiments ####

#------------------------------------------------------------------------------#
#### Experiment 02_00 ####
# profile functions to see where bottle necks are
#------------------------------------------------------------------------------#

ex02_00_1 <- function(){
  X_o <- get_design_frame_00(rhs_outcome, g1_dt)
  X_t <- get_design_frame_00(get_fixed_formula(m_A), g1_dt)
  A   <- get_response(formula(m_A), g1_dt)
  Y   <- get_response(formula(m_Y), g1_dt)
#   w   <- weight_estimator(A, X_t)
#   w   <- w(theta_t)
  f <- make_dbr_estimator_00(Y, A, X_o, X_t, rhs_o)
  f(theta_a, alpha = .5, a = 1)
}

ex02_00_2 <- function(){
  X_o <- get_design_frame_00(rhs_outcome, g1_dt)
  X_t <- get_design_frame_00(get_fixed_formula(m_A), g1_dt)
  A   <- get_response(formula(m_A), g1_dt)
  Y   <- get_response(formula(m_Y), g1_dt)
  #   w   <- weight_estimator(A, X_t)
  #   w   <- w(theta_t)
  f <- make_dbr_estimator_00(Y, A, X_o, X_t, rhs_o)
  numDeriv::grad(f, x = theta_a, alpha = .5, a = 1)
}


ex02_00_3 <- function(){
  lapply(split_dt, function(g_dt){
    X_o <- get_design_frame_00(rhs_outcome, g1_dt)
    X_t <- get_design_frame_00(get_fixed_formula(m_A), g1_dt)
    A   <- get_response(formula(m_A), g1_dt)
    Y   <- get_response(formula(m_Y), g1_dt)
    f <- make_dbr_estimator_00(Y, A, X_o, X_t, rhs_o)
    numDeriv::grad(f, x = theta_a, alpha = .5, a = 1)
  })
}

ex02_00_4 <- function(){
  X_o <- get_design_frame_00(rhs_outcome, g1_dt)
  X_t <- get_design_frame_00(get_fixed_formula(m_A), g1_dt)
  A   <- get_response(formula(m_A), g1_dt)
  Y   <- get_response(formula(m_Y), g1_dt)
  #   w   <- weight_estimator(A, X_t)
  #   w   <- w(theta_t)
  f <- make_dbr_estimator_00(Y, A, X_o, X_t, rhs_o)
  numDeriv::grad(f, x = theta_a, alpha = .5, a = 1, method = 'simple')
}


#------------------------------------------------------------------------------#
#### Experiment 02_01 ####
# create a function that directly computes the closed form of the U21 matrix
#------------------------------------------------------------------------------#

ex02_01_1 <- function(){
  X_o <- get_design_frame_00(rhs_outcome, g1_dt)
  X_t <- get_design_frame_00(get_fixed_formula(m_A), g1_dt)
  A   <- get_response(formula(m_A), g1_dt)
  Y   <- get_response(formula(m_Y), g1_dt)
  f <- make_dbr_U21_estimator_01(Y, A, X_o, X_t, rhs_o)
  f(alpha = .5, a = 1)
}

ex02_01_2 <- function(){
  lapply(split_dt, function(g_dt){
  X_o <- get_design_frame_00(rhs_outcome, g1_dt)
  X_t <- get_design_frame_00(get_fixed_formula(m_A), g1_dt)
  A   <- get_response(formula(m_A), g1_dt)
  Y   <- get_response(formula(m_Y), g1_dt)
  f <- make_dbr_U21_estimator_01(Y, A, X_o, X_t, rhs_o)
  f(alpha = .5, a = 1)
  })
}

ex02_01_3 <- function(){
  X_o <- get_design_frame_00(rhs_outcome, g1_dt)
  X_t <- get_design_frame_00(get_fixed_formula(m_A), g1_dt)
  A   <- get_response(formula(m_A), g1_dt)
  Y   <- get_response(formula(m_Y), g1_dt)
  f <- make_dbr_U21_estimator_01(Y, A, X_o, X_t, rhs_o)
  f(alpha = .25, a = NULL)
}


ex02_01_4 <- function(){
  lapply(split_dt, function(g_dt){
    X_o <- get_design_frame_00(rhs_outcome, g1_dt)
    X_t <- get_design_frame_00(get_fixed_formula(m_A), g1_dt)
    A   <- get_response(formula(m_A), g1_dt)
    Y   <- get_response(formula(m_Y), g1_dt)
    f <- make_dbr_U21_estimator_01(Y, A, X_o, X_t, rhs_o)
    list(f(alpha = .55, a = 1), f(alpha = .25, a = 1))
  })
}

#------------------------------------------------------------------------------#
#### Experiment 02_02 ####
# create a function that directly computes the closed form of the U21 matrix;
# but passing in the Weights and derivatives precomputed
#------------------------------------------------------------------------------#


ex02_02_1 <- function(){
  X_t <- get_design_frame_00(get_fixed_formula(m_A), g1_dt)
  A   <- get_response(formula(m_A), g1_dt)
  wf  <- weight_estimator(A, X_t)
  W   <- wf(theta_t)
  Wd  <- numDeriv::grad(wf, x = theta_t)
  X_o <- get_design_frame_00(rhs_outcome, g1_dt)
  A   <- get_response(formula(m_A), g1_dt)
  Y   <- get_response(formula(m_Y), g1_dt)
  f <- make_dbr_U21_estimator_02(Y, A, X_o, W, Wd, rhs_o)
  f(alpha = .5, a = 1)
}


ex02_02_2 <- function(){
  weights <- lapply(split_dt, function(g_dt){
    X_t <- get_design_frame_00(get_fixed_formula(m_A), g1_dt)
    A   <- get_response(formula(m_A), g1_dt)
    wf  <- weight_estimator(A, X_t)
    list(W = wf(theta_t), Wd = numDeriv::grad(wf, x = theta_t))
  })
  
  lapply(seq_along(split_dt), function(i){
    X_o <- get_design_frame_00(rhs_outcome, split_dt[[i]])
    A   <- get_response(formula(m_A), split_dt[[i]])
    Y   <- get_response(formula(m_Y), split_dt[[i]])
    W   <- weights[[i]][[1]]
    Wd  <-  weights[[i]][[2]]
    
    f <- make_dbr_U21_estimator_02(Y, A, X_o, W, Wd, rhs_o)
    f(alpha = .5, a = 1)
  })
}

#------------------------------------------------------------------------------#
#### Experiment 02_03 ####
# modify 01 by using method = 'simple' to observe speed gains
#------------------------------------------------------------------------------#

ex02_03_1 <- function(){
  X_o <- get_design_frame_00(rhs_outcome, g1_dt)
  X_t <- get_design_frame_00(get_fixed_formula(m_A), g1_dt)
  A   <- get_response(formula(m_A), g1_dt)
  Y   <- get_response(formula(m_Y), g1_dt)
  f <- make_dbr_U21_estimator_03(Y, A, X_o, X_t, rhs_o)
  f(alpha = .5, a = 1)
}

ex02_03_2 <- function(){
  lapply(split_dt, function(g_dt){
    X_o <- get_design_frame_00(rhs_outcome, g1_dt)
    X_t <- get_design_frame_00(get_fixed_formula(m_A), g1_dt)
    A   <- get_response(formula(m_A), g1_dt)
    Y   <- get_response(formula(m_Y), g1_dt)
    f <- make_dbr_U21_estimator_03(Y, A, X_o, X_t, rhs_o)
    list(f(alpha = .55, a = 1), f(alpha = .25, a = 1))
  })
}

