#------------------------------------------------------------------------------#
#   Title: Doubly Robust Functions
#  Author: B. Saul
#    Date: 2016-06-13
# Purpose: functions for IPW and DR estimators
#          - updates outcome functions to handle different link functions
#          - phutz with how the weights are estimator for the ipw
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

# integrand <- function(b, response, xmatrix, theta){
#   lp <- outer(xmatrix %*% theta[-length(theta)], b, '+')
#   h  <- apply(lp, 3, function(x) dbinom(response, 1, plogis(x) ) )
#   print(h)
#   hh <- apply(h, 2, prod)
#   print(hh)
#   hh * dnorm(b, mean = 0, sd = theta[length(theta)])
# }

integrand2 <- function(b, response, xmatrix, theta, alpha){
  # print(b)
  lp <- outer(xmatrix %*% theta[-length(theta)], b, '+')
  plp <- plogis(lp)
  h  <- apply(plp, 3, function(x) response * log(x/alpha) + (1 - response) * log((1 - x)/(1 - alpha)))
  # h  <- apply(plp, 3, function(x) response * log(x) + (1 - response) * log(1 - x))
  # print(h)
  hh <- apply(h, 2, function(x) exp(sum(x)))
  # print(hh)
  hh * dnorm(b, mean = 0, sd = theta[length(theta)])
}

pi_term <- function(A){
  f <- function(alpha){
    prod(dbinom(x = A, 1, prob = alpha))
  }
  memoise::memoise(f)
}

grids <- new.env()

expand_grid_n <- function(n1, n2)
{
  lookup <- paste(n1, n2, sep = '_')
  if(!exists(lookup, envir = grids)){
    grids[[lookup]] <- expand.grid(ID = 1:n1, sum_a = 0:n2, 
                                   A = 0:1) %>%
      mutate_(fA = ~ sum_a / n1) 
  }
  return(grids[[lookup]])
}

expand_outcome_frame <- function(X_outcome, rhs_formula_outcome){
  n <- nrow(X_outcome)
  X_outcome %>%
    mutate_(ID = ~ row_number()) %>%
    # Remove treatment variables so that they are updated 
    select(-A, -fA) %>%
    # Generate the relevant values of A and 
    # all possible sum(a_i) for each subject
    full_join(expand_grid_n(n, n - 1), by = "ID") 
}


#------------------------------------------------------------------------------#
#### IPW Estimator ####
#------------------------------------------------------------------------------#

# weight_estimator <- function(A, X)
# {
#   X <- as.matrix(X)
#   f <- function(theta){
#     1/integrate(integrand, lower = -Inf, upper = Inf,
#                 theta = theta, response = A, xmatrix = X)$value
#   }
#   memoise::memoise(f)
# }

weight_estimator2 <- function(A, X, lower = -5, upper = 5)
{
  X <- as.matrix(X)
  f <- function(theta, alpha){
    1/integrate(integrand2, lower = lower, upper = upper,
                theta = theta, alpha = alpha, response = A, xmatrix = X)$value
  }
  f
  # memoise::memoise(f)
}


# make_ipw_estimator <- function(Y, A, X_treatment, ...){
#   w <- weight_estimator(A = A, X = X_treatment)
#   pi_t <- pi_term(A = A)
#   
#   function(theta, alpha, a = NULL){
#     Ia <- if(is.null(a)) 1 else (A == a) * 1
#     mean(Y * Ia) * w(theta) * pi_t(alpha)  / {if(!is.null(a)) dbinom(a, 1, alpha) else 1}
#   }
# }

make_ipw_estimator <- function(Y, A, X_treatment, ...){
  w <- weight_estimator2(A = A, X = X_treatment, ...)
  pi_t <- pi_term(A = A)
  
  function(theta, alpha, a = NULL){
    Ia <- if(is.null(a)) 1 else (A == a) * 1
    mean(Y * Ia) * w(theta, alpha)  / {if(!is.null(a)) dbinom(a, 1, alpha) else 1}
  }
}

# make_ipw_U21_estimator <- function(Y, A, 
#                                    X_treatment, 
#                                    theta_t)
# {
#   w <- weight_estimator(A = A, X = X_treatment)
#   wd <- numDeriv::grad(w, x = theta_t, method = 'simple')
#   pi_t <- pi_term(A = A)
#   function(alpha, a = NULL){
#     pi_t1 <- pi_t(alpha) / {if(!is.null(a)) dbinom(a, 1, alpha) else 1}
#     # U21 corresponding to treatment parameters
#     Ia <- if(is.null(a)) 1 else (A == a) * 1
#     Ybar <- mean(Ia * Y) 
#     U21_t <- Ybar * wd * pi_t1
#     return(U21_t)
#   }
# }

make_ipw_U21_estimator <- function(Y, A, 
                                   X_treatment, 
                                   theta_t)
{
  w <- weight_estimator(A = A, X = X_treatment)
  function(alpha, a = NULL){
    wd <- numDeriv::grad(w, x = theta_t, alpha = alpha, method = 'simple')
    wd <- wd / {if(!is.null(a)) dbinom(a, 1, alpha) else 1}
    # U21 corresponding to treatment parameters
    Ia <- if(is.null(a)) 1 else (A == a) * 1
    Ybar <- mean(Ia * Y) 
    U21_t <- Ybar * wd
    return(U21_t)
  }
}

#------------------------------------------------------------------------------#
#### Outcome Estimator ####
#------------------------------------------------------------------------------#

make_otc_estimator <- function(X_outcome, rhs_formula_outcome, inv_link, ...){
  
  X_expanded <- expand_outcome_frame(X_outcome, rhs_formula_outcome)
  
  function(theta, alpha, a = NULL){
    n1 <- nrow(X_outcome)
    n2 <- n1 - {if(is.null(a)) 0 else 1}
    a <- if(is.null(a)) 0:1 else a
    
    x <- X_expanded %>%
      filter_(~ A %in% a) %>%
      mutate_(pi_value = ~ dbinom(sum_a, n2, prob = alpha),
              fitted   = ~ as.numeric(inv_link(model.matrix(rhs_formula_outcome, data = .) %*% theta )),
              estimate = ~ fitted * pi_value) 
    
    sum(tapply(tapply(x$estimate, paste(x$A, x$ID), sum), 
               rep(a, each = n1), sum)) / n1
  }
}

make_otc_U21_estimator <- function(X_outcome, 
                                   rhs_formula_outcome,
                                   theta_o,
                                   inv_link,
                                   ...)
{
  function(alpha, a = NULL){
      f <- make_otc_estimator(X_outcome, rhs_formula_outcome, inv_link = inv_link)
      numDeriv::grad(f, x = theta_o, method = 'simple', alpha = alpha, a = a)
  }
}


# make_otc_U21_estimator <- function(A, 
#                                    X_outcome, 
#                                    rhs_formula_outcome,
#                                    theta_o,
#                                    inv_link)
# {
#   pi_t <- pi_term(A = A)
#   dr_term1 <- make_dr_term1(X_outcome, inv_link = inv_link)(theta_o)
#   n <- nrow(X_outcome)
#   function(alpha, a = NULL){
#     pi_t1 <- pi_t(alpha) / {if(!is.null(a)) dbinom(a, 1, alpha) else 1}
#     Ia <- if(is.null(a)) 1 else (A == a) * 1
#     # U21 corresponding to  outcome parameters
#     xmat <- model.matrix(rhs_formula_outcome, data = X_outcome)
#     
#     U21_o <- apply(xmat, 2, function(col) {
#       sum(col)/n
#     })
#     
#     return(U21_o)
#   }
# }

#------------------------------------------------------------------------------#
#### Doubly Robust Estimator ####
#------------------------------------------------------------------------------#

make_dr_term1 <- function(X, inv_link){
  X <- as.matrix(X)
  f <- function(theta){
    inv_link(X %*% theta)
  }
  memoise::memoise(f)
}

make_dbr_estimator <- function(Y, A, X_outcome, X_treatment, rhs_formula_outcome, inv_link){
  w <- weight_estimator2(A = A, X = X_treatment)
  dr_term1 <- make_dr_term1(X_outcome, inv_link = inv_link)
  dr_term2 <- make_otc_estimator(X_outcome, 
                                 rhs_formula_outcome = rhs_formula_outcome, 
                                 inv_link = inv_link)
  
  q_treatment <- ncol(X_treatment) + 1

  function(theta, alpha, a = NULL){
    Ia <- if(is.null(a)) 1 else (A == a) * 1
    Ybar <- mean(Ia * (Y - dr_term1(theta[(q_treatment + 1):length(theta)] ) ) )
    term1 <- Ybar * w(theta[1:q_treatment], alpha) / (dbinom(a, 1, alpha) * !is.null(a) * 1)
    term2 <- dr_term2(theta[(q_treatment + 1):length(theta)], alpha, a)
    term1 + term2
  }
}


make_dbr_U21_estimator <- function(Y, A, 
                                   X_outcome, 
                                   X_treatment, 
                                   rhs_formula_outcome,
                                   theta_t,
                                   theta_o,
                                   inv_link,
                                   ...)
{
  f <- make_dbr_estimator(Y = Y, A = A, 
                          X_outcome = X_outcome, 
                          X_treatment = X_treatment,
                          rhs_formula_outcome = rhs_formula_outcome, 
                          inv_link = inv_link)
  theta <- c(theta_t, theta_o)
  function(alpha, a = NULL){
    numDeriv::grad(f, x = theta, method = 'simple', alpha = alpha, a = a)
  }
}

# make_dbr_U21_estimator <- function(Y, A, 
#                                    X_outcome, 
#                                    X_treatment, 
#                                    rhs_formula_outcome,
#                                    theta_t,
#                                    theta_o,
#                                    inv_link)
# {
#   w <- weight_estimator(A = A, X = X_treatment)
#   ws <- w(theta_t)
#   wd <- numDeriv::grad(w, x = theta_t, method = 'simple')
#   pi_t <- pi_term(A = A)
#   dr_term1 <- make_dr_term1(X_outcome, inv_link = inv_link)(theta_o)
#   n <- nrow(X_outcome)
#   
#   function(alpha, a = NULL){
#     pi_t1 <- pi_t(alpha) / {if(!is.null(a)) dbinom(a, 1, alpha) else 1}
#     # U21 corresponding to treatment parameters
#     Ia <- if(is.null(a)) 1 else (A == a) * 1
#     Ybar <- mean(Ia * (Y - dr_term1) ) 
#     U21_t <- Ybar * wd * pi_t1
#     
#     # U21 corresponding to  outcome parameters
#     xmat <- model.matrix(rhs_formula_outcome, data = X_outcome)
#     tt <- 1 - (Ia * pi_t1 * ws)
#     
#     U21_o <- apply(xmat, 2, function(col) {
#       sum(col * tt)/n
#     })
#     
#     return(c(U21_t, U21_o))
#   }
# }
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
    Y   <- get_response(form_o, group_data)
    A   <- get_response(form_t, group_data)
    X_t <- get_design_frame(form_rhs_t, group_data)
    X_o <- get_design_frame(form_rhs_o, group_data)
    
    estimators <- list(make_ipw_estimator(Y = Y, A = A, X_treatment = X_t),
                       make_otc_estimator(X_outcome = X_o, 
                                         rhs_formula_outcome = form_rhs_o,
                                         inv_link = otc_inv_link),
                       make_dbr_estimator(Y = Y, A = A, 
                                          X_outcome = X_o, 
                                          X_treatment = X_t, 
                                          rhs_formula_outcome = form_rhs_o,
                                          inv_link = otc_inv_link))
#     U21_funcs <- list(make_ipw_U21_estimator(Y = Y, A = A, X_treatment = X_t, theta_t = theta_t),
#                       make_otc_U21_estimator(A = A,
#                                              X_outcome = X_o, 
#                                              rhs_formula_outcome = form_rhs_o,
#                                              theta_o = theta_o,
#                                              inv_link = otc_inv_link),
#                       make_dbr_U21_estimator(Y = Y, A = A, 
#                                              X_outcome = X_o, 
#                                              X_treatment = X_t, 
#                                              rhs_formula_outcome = form_rhs_o,
#                                              theta_t = theta_t,
#                                              theta_o = theta_o,
#                                              inv_link = otc_inv_link))
    data_frame(n_i = length(Y),
               estimator_type = c('ipw', 'otc', 'dbr'),
               estimator = estimators
               # U21_func  = U21_funcs
               )
  }) %>% bind_rows() 
  
  frame %<>% 
    left_join(theta, by = 'estimator_type') %>%
    left_join(counterfactual_grid, by = 'estimator_type')
  
  #### Evaluate group-level estimators ####
  frame %<>%
    rowwise() %>%
    mutate(estimate = evaluate_df_function(estimator, theta = theta, alpha = alpha, a = a)
           # U21      = list(evaluate_df_function(U21_func, alpha = alpha, a = a) )
           # psi_func = list(evaluate_df_function(psi, estimator))
           ) %>%
    ungroup() %>%
    group_by(estimator_type, alpha, a) %>%
    mutate_(target_estimate = ~ sum(estimate, na.rm = T)/m,
            psi_val  = ~ estimate - target_estimate) %>%
    ungroup() 
  
  print('frame created')
  print('evaluating results')
  results <- plyr::ddply(frame, plyr::.(estimator_type, alpha, a), .progress = 'text', function(x){
#     if(x$estimator_type[1] == 'ipw'){
#       ee <- cbind(ee_t, x$psi_val)
#       # bb <- bread_t
#     } else if(x$estimator_type[1] == 'otc') {
#       ee <- cbind(ee_o, x$psi_val)
#       # bb <- bread_o
#     } else if(x$estimator_type[1] == 'dbr') {
#       ee <- cbind(ee_t, ee_o, x$psi_val)
#       # bb <- Matrix::bdiag(bread_t, bread_o)
#     }
#   
#     V <- crossprod(ee)/nrow(ee)
# #     U2 <- c(apply(list_matrix(x$U21), 2, mean), 1)
# #     U <- rbind_fill_zero(list(bb, matrix(-U2, ncol = length(U2))))
# #     sigma <- (solve(U) %*% V %*% t(solve(U)))/nrow(ee)
# #     std.error <- sqrt(sigma[nrow(sigma), ncol(sigma)])
#     
#     U21 <- apply(list_matrix(x$U21), 2, mean)
#     U21 <- matrix(-U21, ncol = length(U21))
#     V11 <- V[1:(nrow(V) - 1), 1:(ncol(V) - 1)]
#     V21 <- V[nrow(V), 1:(ncol(V) - 1)]
#     V22 <- V[nrow(V), ncol(V)]
#     var.est <- ( ( (U21 - 2 * V21) %*% solve(V11) %*% t(U21) ) + V22) / nrow(ee)
#     std.error <- as.numeric(sqrt(var.est))
    std.error <- NA
    df_out <- data_frame(estimate = x$target[1], 
                         std.error = NA)
    return(df_out)
  })
          
#   return(results)
  # results
  frame
}

#------------------------------------------------------------------------------#
#### Run on simulations ####
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
      out
    }
  })
}

