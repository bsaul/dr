#------------------------------------------------------------------------------#
#   Title: Doubly Robust Functions
#  Author: B. Saul
#    Date: 2016-04-20
# Purpose: functions for IPW and DR estimators
#------------------------------------------------------------------------------#

#### Functions ####

integrand <- function(b, response, xmatrix, theta){
  lp <- outer(xmatrix %*% theta[-length(theta)], b, '+')
  h  <- apply(lp, 3, function(x) dbinom(response, 1, plogis(x) ) )
  hh <- apply(h, 2, prod)
  hh * dnorm(b, mean = 0, sd = theta[length(theta)])
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

pi_term <- function(A){
  f <- function(alpha){
    prod(dbinom(x = A, 1, prob = alpha))
  }
}

make_dr_term1 <- function(X){
  X <- as.matrix(X)
  function(theta){
    X %*% theta
  }
}

make_outcome_estimator <- function(X_outcome, formula_outcome, ...){
  function(theta, alpha, a = NULL){
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
      mutate_(mu =~ as.numeric(model.matrix(formula_outcome[-2], data = .) %*% theta )  ) %>%
      # Sum by individual to compute term2
      group_by_(~ID) %>%
      summarize_(term2 = ~sum(mu * pi)) %>%
      summarize_(mean = ~mean(term2)) %>%
      as.numeric()
  }
}

make_ipw_estimator <- function(Y, A, X_treatment, ...){
  w <- weight_estimator(A = A, X = X_treatment)
  pi_t <- pi_term(A = A)

  function(theta, alpha, a = NULL){
    Ia <- if(is.null(a)) 1 else (A == a) * 1
    mean(Y * Ia) * w(theta) * pi_t(alpha)  / (dbinom(a, 1, alpha) * !is.null(a) * 1)  
  }
}


make_dr_estimator <- function(Y, A, X_outcome, X_treatment, formula_outcome){
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

psi <- function(group_estimator){
  function(theta, alpha, a){
    q <- length(theta)
    group_estimator(theta[1:(q-1)], alpha, a) - theta[q]
  }
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

simestimation <- function(treatment_formula, 
                          outcome_formula,
                          simdata)
{
  
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
  
  models <- make_models(model_args = model_args, data = simdata)
  
  theta_treatment <- unlist(lme4::getME(models$model_treatment, c('beta', 'theta')))
  theta_outcome   <- coef(models$model_outcome)
  theta_observed  <- c(theta_treatment, theta_outcome)
  
  #### Create list of group level data  ####
  groups         <- simdata$group
  outcome        <- simdata$Y
  xmat_treatment <- as.data.frame(model.matrix(models$model_treatment))
  xmat_outcome   <- as.data.frame(model.matrix(models$model_outcome))
  treatment      <- model.response(model.frame(models$model_treatment))
  
  
  frame <- list(outcome = outcome, 
                treatment = treatment,
                X_treatment = xmat_treatment,
                X_outcome   = xmat_outcome) %>%
    lapply(., function(x) split(x, groups, drop = FALSE))
  
  split_frame <- Map(list, outcome      = frame[['outcome']],
                     treatment    = frame[['treatment']],
                     X_treatment  = frame[['X_treatment']],
                     X_outcome    = frame[['X_outcome']])
  
  #### Create list of group level estimators ####
  dr_estimators  <- make_group_estimators(split_frame,
                                          estimator = 'make_dr_estimator',
                                          formula_outcome = formula(models$model_outcome))
  
  
  ipw_estimators  <- make_group_estimators(split_frame,
                                           estimator = 'make_ipw_estimator')
  
  outcome_estimators  <- make_group_estimators(split_frame,
                                               estimator = 'make_outcome_estimator',
                                               formula_outcome = formula(models$model_outcome))
  
  #   dr_estimators  <- lapply(dr_estimators, function(f) memoise::memoise(f))
  #   ipw_estimators <- lapply(ipw_estimators, function(f) memoise::memoise(f))
  #   ipw_estimators <- lapply(outcome_estimators, function(f) memoise::memoise(f))
  
  models2 <- list(models[[1]], models[[2]], models)
  estimators <- list(ipw = ipw_estimators, outcome = outcome_estimators, dr = dr_estimators)
  thetas     <- list(theta_treatment, theta_outcome, theta_observed)
  temp       <- lapply(models, bread)
  breads     <- list(temp[[1]], temp[[2]], Matrix::bdiag(temp))
  temp       <- lapply(models, estfun)
  efuns      <- list(temp[[1]], temp[[2]], cbind(temp[[1]], temp[[2]]))
  
  alphas   <- c(0.1, 0.5, 0.9)
  target_a <- 1
  
  #### Estimate target parameters ####
  results <- lapply(seq_along(estimators),
                    function(i){
                      
                      funcs <- estimators[[i]]
                      theta <- thetas[[i]]
                      model <- models2[[i]]
                      type  <- names(estimators[i])
                      ee_0  <- efuns[[i]]
                      bb_0  <- breads[[i]]
                      
                      lapply(alphas, function(talpha){
                        
                        ## Collect values 
                        target <- lapply(funcs, function(f) {
                          f(theta = theta,
                            alpha = talpha,
                            a     = target_a)
                        }) %>%
                          list_matrix() %>%
                          mean()
                        
                        psi_observed <- lapply(funcs, function(f) {
                          f2 <- psi(f)
                          f2(theta = c(theta, target), 
                             alpha = talpha, 
                             a     = target_a)
                        }) %>%
                          list_matrix()
                        
                        psi_prime_target <- lapply(funcs, function(f) {
                          f2 <- psi(f)
                          numDeriv::grad(f2, 
                                         x = c(theta, target), 
                                         alpha = talpha, 
                                         a = target_a)
                        }) %>%
                          list_matrix() %>% 
                          apply(2, mean) %>% 
                          matrix(nrow = 1)
                        
                        ## Compile values
                        ee   <- cbind(ee_0, psi_observed)
                        V    <- crossprod(ee)/nrow(ee)
                        U <- as.matrix(rbind_fill_zero(list(bb_0, -psi_prime_target)))
                        sigma <- (solve(U) %*% V %*% t(solve(U)))/nrow(ee)
                        
                        target.sigma <- sigma[nrow(sigma), ncol(sigma)]
                        std.error <- sqrt(target.sigma)
                        
                        #output dataframe
                        data.frame(estimator = type,
                                   alpha     = talpha,
                                   a         = target_a,
                                   estimate  = target, 
                                   std.error = std.error,
                                   conf.high = target + 1.96 * std.error,
                                   conf.low  = target - 1.96 * std.error )
                      }) %>%
                        bind_rows()
                    }) %>%
    bind_rows() %>%
    mutate_(simID = simdata$simID[1])
  return(results)
}

simestimation_all <- function(sims,
                              treatment_formula, 
                              outcome_formula,
                              progress = 'none', 
                              parallel = FALSE)
{
  sims %>% 
  {plyr::dlply(., plyr::.(simID), .progress = progress, .parallel = parallel,
               function(sim) { 
                 simestimation(treatment_formula, outcome_formula, simdata =  sim) 
               } ) }
}