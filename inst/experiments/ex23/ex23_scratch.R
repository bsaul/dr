dbr_estimator_test <- function(data, models, randomization, ...){
  
  ## component data
  comp <- extract_model_info(
    model = models, 
    data = data, 
    estimator_type = 'reg_dbr',
    regression_type = 'pcov')
  A  <- comp$A
  N  <- comp$N
  lnkinv <- comp$inv_link_o
  
  ip_fun <- weight_estimator(
    A = comp$A, 
    X = comp$X_t, 
    randomization = randomization)
  
  dr_term1_fun_0 <- make_dr_term1(
    comp$X_o_reg_0, 
    inv_link = comp$inv_link_o)
  
  dr_term1_fun_1 <- make_dr_term1(
    comp$X_o_reg_1, 
    inv_link = comp$inv_link_o)

  MM_0   <- model.matrix(comp$rhs_o_reg_0, data = comp$X_o_reg_0)
  MM_1   <- model.matrix(comp$rhs_o_reg_1, data = comp$X_o_reg_1)
  
  index_t   <- 1:comp$p_t
  index_o_0 <- (comp$p_t + 1):(comp$p_t + comp$p_o_0)
  index_o_1 <- (comp$p_t + comp$p_o_0 + 1):(comp$p_t + comp$p_o_0 + comp$p_o_1)

  function(theta, alpha){
    stopifnot(length(alpha) == 1)
    
    ### Term 1 ###
    fY_0    <- dr_term1_fun_0(theta[index_o_0])
    fY_1    <- dr_term1_fun_1(theta[index_o_1])
    ipw     <- ip_fun(theta[index_t], alpha)
    term1_0 <- (A == 0) * fY_0 * ipw / (1 - alpha)
    term1_1 <- (A == 1) * fY_1 * ipw / alpha
    # term1_0 <- (A == 0) *  ipw / (1 - alpha)
    # term1_1 <- (A == 1) *  ipw / alpha
    # 
    ### Term 2 ###
    nsamples <- 25
    hold_0  <- hold_1 <- matrix(NA, nrow = N, ncol = nsamples)
    for(i in 1:nsamples){
      A_tilde <- rbinom(N, 1, prob = alpha)
      fA <- sum(A_tilde)/N
      fA_0 <- (sum(A_tilde) - (A_tilde == 1))/N # remove 1 when A_ij == 1
      fA_1 <- (sum(A_tilde) + (A_tilde == 0))/N # add 1 when A_ij == 0
      
      ipw0 <- ipw1 <- numeric(N)
      for(j in 1:N){
        A_tilde_0 <-  A_tilde_1 <- A_tilde
        A_tilde_0[j] <- 0
        A_tilde_1[j] <- 1
        
        # print(A_tilde_0)
        ip_fun_0 <- weight_estimator(
          A = A_tilde_0, 
          X = comp$X_t, 
          randomization = randomization)
        
        ip_fun_1 <- weight_estimator(
          A = A_tilde_1, 
          X = comp$X_t, 
          randomization = randomization)
        
        ipw0[j] <- ip_fun_0(theta_t, alpha)/(1 - alpha)
        ipw1[j] <- ip_fun_1(theta_t, alpha)/alpha
      }
      
      # print(cbind(ipw0, ipw1))

      MM_0[, 'fA']  <- fA
      MM_1[, 'fA']  <- fA
      MM_0[, 'ipw'] <- ipw0
      MM_1[, 'ipw'] <- ipw1
      # print(MM_0)
      # print(MM_1)
      hold_0[ , i] <- lnkinv(MM_0 %*% theta[index_o_0])
      hold_1[ , i] <- lnkinv(MM_1 %*% theta[index_o_1])
    }

    term2_0 <- apply(hold_0, 1, mean)
    term2_1 <- apply(hold_1, 1, mean)
    
    # print(cbind(term1_0, term2_0))
    
    ### Sum ###
    dbr2_ce0 <- sum(term1_0) / N
    dbr2_ce1 <- sum(term2_0) / N
    # dbr2_ce0 <- sum(-term1_0 + term2_0) / N
    # dbr2_ce1 <- sum(-term1_1 + term2_1) / N

    x <- c(dbr2_ce0, dbr2_ce1) 
    label0 <- paste0('_dbr_Y0_')
    label1 <- paste0('_dbr_Y1_')
    names(x) <- paste0(rep(c(label0, label1), each = length(alpha)), alpha)
    x
  }
}

dbr_estimator_test2 <- function(data, models, randomization, ...){
  
  ## component data
  comp <- extract_model_info(
    model = models, 
    data = data, 
    estimator_type = 'reg_dbr',
    regression_type = 'pcov')
  A <- comp$A
  Y <- comp$Y
   
  ip_fun <- weight_estimator(
    A = comp$A, 
    X = comp$X_t, 
    randomization = randomization)

  dr_term1_fun_0 <- make_dr_term1(
    comp$X_o_reg_0, 
    inv_link = comp$inv_link_o)
  
  dr_term1_fun_1 <- make_dr_term1(
    comp$X_o_reg_1, 
    inv_link = comp$inv_link_o)
  
  index_t <- 1:comp$p_t
  index_o_0 <- (comp$p_t + 1):(comp$p_t + comp$p_o_0)
  index_o_1 <- (comp$p_t + comp$p_o_0 + 1):(comp$p_t + comp$p_o_0 + comp$p_o_1)
  N    <- comp$N
  
  function(theta, alpha){
    stopifnot(length(alpha) == 1)
    
    fY_0    <- dr_term1_fun_0(theta[index_o_0])
    fY_1    <- dr_term1_fun_1(theta[index_o_1])
    ipw     <- ip_fun(theta[index_t], alpha)
    Ybar0   <- sum((A == 0) * (Y - fY_0) )
    Ybar1   <- sum((A == 1) * (Y - fY_1) )
    term1_0 <- Ybar0 * ipw / (1 - alpha)
    term1_1 <- Ybar1 * ipw / alpha
  
    dbr2_ce0 <- (term1_0) / N
    dbr2_ce1 <- (term1_1) / N
    
    x <- c(dbr2_ce0, dbr2_ce1) 
    label0 <- paste0('_dbr_Y0_')
    label1 <- paste0('_dbr_Y1_')
    names(x) <- paste0(rep(c(label0, label1), each = length(alpha)), alpha)
    x
  }
}



###### 

which_m <- 3
test <- do.call(gen_sim, args = arg_maker(9, nsims = 1))[[1]]

m <- make_models(model_args = margs[[which_m]][c('t_model', 'o_model')], 
                 data = test)

## Grab component model parameter estimates
theta_t <- unlist(lme4::getME(m$t_model, c('beta', 'theta')))
theta_o <- coef(m$o_model)

ipwv <- make_ipw_vector(fulldata = test, models = m, group = 'group', 
                        alpha = 0.5)

A <- test$A
# Add IP weights to dataset
test <- test %>%
  mutate_(
    ipw  =~ ipwv,
    ipw0 =~ (A == 0) * ipw,
    ipw1 =~ (A == 1) * ipw
  )
# Fit regression-based models
wls_model_0 <- try(glm(
  formula = margs[[which_m]][['wls_model_0']][['formula']], 
  family  = margs[[which_m]][['wls_model_0']][['options']][['family']],
  weights = (A == 0) * ipw0,
  data    = test))

wls_model_1 <- try(glm(
  formula = margs[[which_m]][['wls_model_1']][['formula']], 
  family  = margs[[which_m]][['wls_model_1']][['options']][['family']],
  weights = (A == 1) * ipw1,
  data    = test))


if(is(wls_model_0, 'try-error') | is(wls_model_1, 'try-error')){
  skipwls <- TRUE
  theta_wls_0 <- NA
  theta_wls_1 <- NA
} else {
  skipwls <- FALSE
  m$wls_model_0 <- wls_model_0
  m$wls_model_1 <- wls_model_1
  theta_wls_0 <- coef(m$wls_model_0)
  theta_wls_1 <- coef(m$wls_model_1)
}

pcov_model_0 <- try(glm(
  formula = margs[[which_m]][['pcov_model_0']][['formula']], 
  family  = margs[[which_m]][['pcov_model_0']][['options']][['family']],
  weights = (A == 0) * 1 ,
  data    = test))

pcov_model_1 <- try(glm(
  formula = margs[[which_m]][['pcov_model_1']][['formula']], 
  family  = margs[[which_m]][['pcov_model_1']][['options']][['family']],
  weights = (A == 1) * 1,
  data    = test))


if(is(pcov_model_0, 'try-error') | is(pcov_model_1, 'try-error')){
  skipwls <- TRUE
  theta_pcov_0 <- NA
  theta_pcov_0 <- NA
} else {
  skipwls <- FALSE
  m$pcov_model_0 <- pcov_model_0
  m$pcov_model_1 <- pcov_model_1
  theta_pcov_0 <- coef(m$pcov_model_0)
  theta_pcov_1 <- coef(m$pcov_model_1)
}

# theta_pcov_0[1] <- 1
# theta_pcov_1[1] <- 1
# theta_pcov_0[2:5] <- 0
# # theta_pcov_1[2:5] <- 0
# theta_pcov_0[5] <- 0
# theta_pcov_1[5] <- 0
test_split <- split(test, f = test$group)
theta <- c(theta_t, theta_pcov_0, theta_pcov_1)
rm(A)
ff <- dbr_estimator_test(data = test_split[[1]], models = m, randomization = 1)
ff(theta, .5)
# 
# ff2 <- dbr_estimator_test2(data = test_split[[1]], models = m, randomization = 1)
# ff2(theta, .55)

# mean(model.matrix(pcov_model_0) %*% theta_pcov_0)
# mean(model.matrix(pcov_model_1) %*% theta_pcov_1)

target <- lapply(test_split, function(grp_dt){
  # Create estimator function
  estimator <- dbr_estimator_test(
    data          = grp_dt, 
    models        = m,
    randomization = 1,
    regression_type = 'pcov')
  # Evaluate estimator function
  x <- estimator(theta, alpha = 0.55)
  # print(x)
  x
})

target %>% list_matrix() %>% apply(., 2, mean)


# ipw_target <- lapply(test_split, function(grp_dt){
#   # Create estimator function
#   estimator <- ipw_estimator(
#     data          = grp_dt, 
#     models        = m,
#     randomization = 1)
#   # Evaluate estimator function
#   x <- estimator(theta, alpha = 0.55)
#   # print(x)
#   x
# })
# 
# ipw_target %>% list_matrix() %>% apply(., 2, mean)


target2 <- lapply(test_split, function(grp_dt){
  # Create estimator function
  estimator <- dbr_estimator_test2(
    data          = grp_dt, 
    models        = m,
    randomization = 1,
    regression_type = 'pcov')
  # Evaluate estimator function
  x <- estimator(theta, alpha = 0.5)
  x
})


target2 %>% list_matrix() %>% apply(., 2, mean)







####
test_MM <- lapply(test_split, function(grp){
  MM <- model.matrix(margs[[which_m]][['pcov_model_0']][['formula']], data = grp)
  Y  <- grp$Y
  A  <- grp$A
  N  <- length(Y)
  function(theta){
    t(MM) %*% (Y - MM %*% theta)
  }
})

testG <- function(theta){
  hold <- lapply(test_MM, function(f) f(theta))
  apply(simplify2array(hold), 1, sum)
}

rootSolve::multiroot(testG, start = theta_pcov_0)
theta_pcov_0


