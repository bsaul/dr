dbr_estimator_test <- function(data, models, randomization, ...){
  
  ## component data
  comp <- extract_model_info(
    model = models, 
    data = data, 
    estimator_type = 'reg_dbr',
    regression_type = 'pcov')
  A <- comp$A
  
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
  
  X_ex_0 <- comp$X_o_ex %>% filter(A == 0)
  X_ex_1 <- comp$X_o_ex %>% filter(A == 1)

  MM_0   <- model.matrix(comp$rhs_o_reg_0, data = comp$X_o_reg_0)
  MM_1   <- model.matrix(comp$rhs_o_reg_1, data = comp$X_o_reg_1)
  
  index_t <- 1:comp$p_t
  index_o_0 <- (comp$p_t + 1):(comp$p_t + comp$p_o_0)
  index_o_1 <- (comp$p_t + comp$p_o_0 + 1):(comp$p_t + comp$p_o_0 + comp$p_o_1)
  N    <- comp$N
  
  function(theta, alpha){
    stopifnot(length(alpha) == 1)
    
    ### Regression-based DRR estimator ###
    fY_0      <- dr_term1_fun_0(theta[index_o_0])
    fY_1      <- dr_term1_fun_1(theta[index_o_1])
    ipw     <- ip_fun(theta[index_t], alpha)
    Ybar0   <- sum((A == 0) * (fY_0) )
    Ybar1   <- sum((A == 1) * (fY_1) )
    term1_0 <- Ybar0 * ipw / (1 - alpha)
    term1_1 <- Ybar1 * ipw / alpha

    A_tilde <- rbinom(N, 1, prob = alpha)
    fA <- sum(A_tilde)/N
    MM_0[, 'fA'] <- fA
    MM_1[, 'fA'] <- fA
    mu_0 <- as.numeric(comp$inv_link_o(MM_0 %*% theta[index_o_0]))
    mu_1 <- as.numeric(comp$inv_link_o(MM_1 %*% theta[index_o_1]))
    
    term2_0 <- sum(mu_0)
    term2_1 <- sum(mu_1)
    
    dbr2_ce0 <- (-term1_0 + term2_0) / N
    dbr2_ce1 <- (-term1_1 + term2_1) / N
    
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

ipw0 <- make_ipw_vector(fulldata = test, models = m, group = 'group', 
                        a = 0, alpha = 0.5)
ipw1 <- make_ipw_vector(fulldata = test, models = m, group = 'group', 
                        a = 1, alpha = 0.5)
A <- test$A
# Add IP weights to dataset
test <- test %>%
  mutate_(
    ipw0 =~ ipw0,
    ipw1 =~ ipw1
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

rm(ipw0, ipw1)

pcov_model_0 <- try(glm(
  formula = margs[[which_m]][['pcov_model_0']][['formula']], 
  family  = margs[[which_m]][['pcov_model_0']][['options']][['family']],
  weights = (A == 0) * 1,
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


test_split <- split(test, f = test$group)
theta <- c(theta_t, theta_pcov_0, theta_pcov_1)
rm(A)
ff <- dbr_estimator_test(data = test_split[[1]], models = m, randomization = 1)
ff(theta, .5)


target <- lapply(test_split, function(grp_dt){
  # Create estimator function
  estimator <- dbr_estimator_test(
    data          = grp_dt, 
    models        = m,
    randomization = 1,
    regression_type = 'pcov')
  # Evaluate estimator function
  x <- estimator(theta, alpha = 0.5)
  print(x)
  x
})


target %>% list_matrix() %>% apply(., 2, mean)


