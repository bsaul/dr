#------------------------------------------------------------------------------#
#   Title: Doubly Robust Functions
#  Author: B. Saul
#    Date: 2016-04-20
# Purpose: functions for IPW and DR estimators
#------------------------------------------------------------------------------#

#### Functions ####
get_fixed_formula <- function(model_object){
  formula(model_object, fixed.only = TRUE)[-2]
}

get_design_frame <- function(rhs_formula, data){
  as.data.frame(model.matrix(rhs_formula, data))
}

get_response <- function(formula, data){
  model.response(model.frame(formula, data = data))
}

integrand <- function(b, response, xmatrix, theta){
  lp <- outer(xmatrix %*% theta[-length(theta)], b, '+')
  h  <- apply(lp, 3, function(x) dbinom(response, 1, plogis(x) ) )
  hh <- apply(h, 2, prod)
  hh * dnorm(b, mean = 0, sd = theta[length(theta)])
}

weight_estimator <- function(A, X)
{
  X <- as.matrix(X)
  function(theta){
    1/integrate(integrand, lower = -Inf, upper = Inf,
                theta = theta, response = A, xmatrix = X)$value
  }
}

pi_term <- function(A){
  function(alpha){
    prod(dbinom(x = A, 1, prob = alpha))
  }
}

make_dr_term1 <- function(X, theta){
  X <- as.matrix(X)
  X %*% theta
}

make_outcome_estimator <- function(X_outcome, formula_outcome, theta, ...){
  function(alpha, a = NULL){
    n <- nrow(X_outcome) - {if(is.null(a)) 0 else 1}
    
    X_outcome %>%
      mutate(ID = row_number()) %>%
      select(-A) %>%
      # Generate all possible sum(a_i) for each subject
      merge(expand.grid(sum_a = 0:n,
                        A = {if(is.null(a)) 0:1 else a }),
            all = T) %>%
      group_by_(~ID) %>%
      # Compute pi and p_ij for each sum(a_i)
      mutate_(fA = ~ sum_a/n(),
              pi = ~ dbinom(sum_a, n, prob = alpha)  )  %>%
      ungroup()  %>%
      # Compute mu_ij for each a_i per subject
      mutate_(mu =~ as.numeric(model.matrix(formula_outcome[-2], data = .) %*% theta ) ) %>%
      # Sum by individual to compute term2
      group_by_(~ID) %>%
      summarize_(term2 = ~sum(mu * pi)) %>%
      summarize_(mean = ~mean(term2)) %>%
      as.numeric()
  }
}

make_ipw_estimator <- function(Y, A, W, ...){
  pi_t <- pi_term(A = A)

  f <- function(alpha, a = NULL){
      Ia <- if(is.null(a)) 1 else (A == a) * 1
      da <- if(is.null(a)) 1 else dbinom(a, 1, alpha)
      mean(Y * Ia) * W * pi_t(alpha)  / da  
  }
  Vectorize(f)
}


check_a <- function(A, a){
  if(is.null(a)) 1 else (A == a) * 1
}

make_ipw_estimator2 <- function(Y, A, W, ...){
  pi_t <- pi_term(A = A)
  
  f <- function(alpha1, a1 = NULL, alpha2 = NULL, a2 = NULL, gFun = `-`){
    Ia1 <- check_a(A, a1)
    Ia2 <- check_a(A, a2)
    da1 <- if(is.null(a1)) 1 else dbinom(a1, 1, alpha1)
    da2 <- if(is.null(a2)) 1 else dbinom(a2, 1, alpha2)
   
    g1 <- mean(Y * Ia1) * W * pi_t(alpha1)  / da1
    
    g2 <- if(!is.null(alpha2)) mean(Y * Ia2) * W * pi_t(alpha2) / da2 else 0
    return(gFun(g1, g2))
  }
}

make_dr_estimator <- function(Y, A, X_outcome, W, formula_outcome, theta){
  pi_t <- pi_term(A = A)
  dr_term1 <- make_dr_term1(X_outcome, theta)
  dr_term2 <- make_outcome_estimator(X_outcome, formula_outcome = formula_outcome, theta)

  function(alpha, a = NULL){
    Ia <- if(is.null(a)) 1 else (A == a) * 1
    Ybar <- mean(Ia * (Y - dr_term1) )
    term1 <- Ybar * W * pi_t(alpha) / (dbinom(a, 1, alpha) * !is.null(a) * 1)
    term2 <- dr_term2(alpha, a)
    term1 + term2
  }
}

make_dr_estimator_deriv <- function(X_outcome, W, A, ...){
  pi_t <- pi_term(A)
  function(alpha, a = NULL){
    tt <- 1 - pi_t(alpha) * W / {if(is.null(a)) 1 else dbinom(a, 1, alpha)}
    
    apply(X_outcome, 2, function(col) mean(col * tt) )
  }
}

psi <- function(estimate, target){
  estimate - target
}

make_group_estimators <- function(split_data,
                                  estimator,
                                  formula_outcome  = NULL)
{
  out <- lapply(split_data, function(group_data)
  {
    X_treatment <- group_data[['X_treatment']]
    X_outcome <- group_data[['X_outcome']]
    A <- group_data[['treatment']]
    Y <- group_data[['outcome']]
    eval(call(estimator, Y = Y, A = A,
              X_treatment = X_treatment,
              X_outcome = X_outcome,
              formula_outcome = formula_outcome))
  } )
  return(out)
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

#### Estimation ####
estimate_sims <- function(formula_treatment, formula_outcome)
{
  plyr::dlply(DRsims, plyr::.(simID), .progress = 'text', .parallel = TRUE, function(this_data)
  {  
    model_args = list(
      model_treatment = list(
        formula = formula_treatment,
        method  = lme4::glmer,
        options = list(family = binomial)
      ) ,
      model_outcome = list(
        formula = formula_outcome,
        method  = geepack::geeglm,
        options = list(family = gaussian, id = quote(group))
      ) )
    
    models <- make_models(model_args = model_args, data = this_data)
    
    theta_treatment <- unlist(lme4::getME(models$model_treatment, c('beta', 'theta')))
    theta_outcome   <- coef(models$model_outcome)
    theta_both      <- c(theta_treatment, theta_outcome)
    
    rhs_treatment <- get_fixed_formula(models$model_treatment)
    rhs_outcome   <- get_fixed_formula(models$model_outcome)
    #### Create list of group level data  ####
    
    split_data <- split(this_data, this_data$group)
    
    weight_funcs <- lapply(split_data, function(group_data){
      XX <- get_design_frame(rhs_treatment, group_data)
      A  <- get_response(formula(models$model_treatment), group_data)
      weight_estimator(A, XX)
    })
    
    weights <- lapply(weight_funcs, function(f) f(theta_treatment))
    weightd <- lapply(weight_funcs, function(f) numDeriv::grad(f, x = theta_treatment, method = 'simple'))
    
    ee_treatment <- estfun(models$model_treatment)
    bread_treatment <- bread(models$model_treatment)
    ee_outcome <- estfun(models$model_outcome)
    bread_outcome <- bread(models$model_outcome)
    
    lapply(alphas, function(target_alpha){
      #### IPW ####
      
      ipw_estimates <- lapply(seq_along(split_data), function(i){
        W  <- weights[[i]]
        A  <- get_response(formula(models$model_treatment), split_data[[i]])
        Y  <- get_response(formula(models$model_outcome), split_data[[i]])
        f <- make_ipw_estimator(Y, A, W)
        f(alpha = target_alpha, a = target_a)
      })
      
      ipw_estimates_U <- lapply(seq_along(split_data), function(i){
        W  <- weightd[[i]]
        A  <- get_response(formula(models$model_treatment), split_data[[i]])
        Y  <- get_response(formula(models$model_outcome), split_data[[i]])
        f  <- make_ipw_estimator(Y, A, W)
        f(alpha = target_alpha, a = target_a)
      })
      
      target <- ipw_estimates %>%
        list_matrix() %>%
        mean()
      
      target_U <- ipw_estimates_U %>% list_matrix() %>% apply(2, mean)
      
      psi_val <- lapply(ipw_estimates, function(estimate) {
        psi(estimate, target)
      }) %>% list_matrix()
      
      ee   <- cbind(ee_treatment, psi_val)
      V    <- crossprod(ee)/nrow(ee)
      U <- rbind_fill_zero(list(bread_treatment, 
                                matrix(-c(target_U, 1), ncol = length(target_U) + 1)))
      
      sigma <- (solve(U) %*% V %*% t(solve(U)))/nrow(ee)
      
      target.sigma <- sigma[nrow(sigma), ncol(sigma)]
      std.error <- sqrt(target.sigma)
      
      #output dataframe
      ipw <- data.frame(estimator = 'ipw',
                        alpha     = target_alpha,
                        a         = target_a,
                        estimate  = target, 
                        std.error = std.error,
                        conf.high = target + 1.96 * std.error,
                        conf.low  = target - 1.96 * std.error,
                        stringsAsFactors = F)
      
      #### Outcome ####
      
      outcome_estimates <- lapply(split_data, function(group_data){
        ff <- formula(models$model_outcome)
        XX <- get_design_frame(rhs_outcome, group_data)
        f <- make_outcome_estimator(XX, ff, theta = theta_outcome)
        f(alpha = target_alpha, a = target_a)
      })
      
      outcome_estimates_U <- lapply(split_data, function(group_data){
        XX <- get_design_frame(rhs_outcome, group_data)
        z1 <- apply(XX, 2, mean)
        c(z1, 1)
      }) %>% list_matrix() %>%
        apply(2, mean)
      
      target <- outcome_estimates %>%
        list_matrix() %>%
        mean()
      
      psi_val <- lapply(outcome_estimates, function(estimate) {
        psi(estimate, target)
      }) %>% list_matrix()

      ee <- cbind(ee_outcome, psi_val)
      V  <- crossprod(ee)/nrow(ee)
      U <- rbind_fill_zero(list(bread_outcome, 
                                matrix(-outcome_estimates_U, 
                                       ncol = length(outcome_estimates_U))))
      
      sigma <- (solve(U) %*% V %*% t(solve(U)))/nrow(ee)
      
      target.sigma <- sigma[nrow(sigma), ncol(sigma)]
      std.error <- sqrt(target.sigma)
      
      outcome <- data.frame(estimator = 'outcome',
                            alpha     = target_alpha,
                            a         = target_a,
                            estimate  = target, 
                            std.error = std.error,
                            conf.high = target + 1.96 * std.error,
                            conf.low  = target - 1.96 * std.error,
                            stringsAsFactors = F)
      
      #### DR #### 
      
      dr_estimates <- lapply(seq_along(split_data), function(i){
        ff <- formula(models$model_outcome)
        XX <- get_design_frame(rhs_outcome, split_data[[i]])
        W  <- weights[[i]]
        A  <- get_response(formula(models$model_treatment), split_data[[i]])
        Y  <- get_response(formula(models$model_outcome), split_data[[i]])
        f  <- make_dr_estimator(Y, A, XX, W = W, formula= ff, theta = theta_outcome)
        f(alpha = target_alpha, a = target_a)
      })
      
      dr_estimates_ipwU <- lapply(seq_along(split_data), function(i){
        W  <- weightd[[i]]
        A  <- get_response(formula(models$model_treatment), split_data[[i]])
        Y  <- get_response(formula(models$model_outcome), split_data[[i]])
        XX_outcome <- get_design_frame(rhs_outcome, split_data[[i]])
        mu <- make_dr_term1(XX_outcome, theta_outcome)
        f  <- make_ipw_estimator(Y - mu, A, W)
        f(alpha = target_alpha, a = target_a)
      })
      
      dr_estimates_U <- lapply(seq_along(split_data), function(i){
        W  <- weights[[i]]
        ipw_U <- dr_estimates_ipwU[[i]]
        XX <- get_design_frame(rhs_outcome, split_data[[i]])
        A  <- get_response(formula(models$model_treatment), split_data[[i]])
        f  <- make_dr_estimator_deriv(XX, W, A)
        c(ipw_U, f(alpha = target_alpha, a = target_a), 1)
      }) %>% list_matrix() %>%
        apply(2, mean)
      
      target <- dr_estimates %>%
        list_matrix() %>%
        mean()
      
      psi_val <- lapply(dr_estimates, function(estimate) {
        psi(estimate, target)
      }) %>% list_matrix()
      
      ee <- cbind(ee_treatment, ee_outcome, psi_val)
      V  <- crossprod(ee)/nrow(ee)
      bread_m <- Matrix::bdiag(bread_treatment, bread_outcome)
      
      U <- rbind_fill_zero(list(bread_m, 
                                matrix(-dr_estimates_U, 
                                       ncol = length(dr_estimates_U))))
      
      sigma <- (solve(U) %*% V %*% t(solve(U)))/nrow(ee)
      
      target.sigma <- sigma[nrow(sigma), ncol(sigma)]
      std.error <- sqrt(target.sigma)
      
      dr <- data.frame(estimator = 'dr',
                       alpha     = target_alpha,
                       a         = target_a,
                       estimate  = target, 
                       std.error = std.error,
                       conf.high = target + 1.96 * std.error,
                       conf.low  = target - 1.96 * std.error,
                       stringsAsFactors = F)
      
      bind_rows(ipw, outcome, dr)
    }) %>% bind_rows() %>%
      mutate_(simID = ~this_data$simID[1])
  })
} %>%
  bind_rows()
