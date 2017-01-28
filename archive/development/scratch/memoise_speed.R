

m_A <- lme4::glmer(tru_propen, family = binomial, data = DRsims)
theta_A <- unlist(getME(m_A, c('beta', 'theta')))
xmat_treatment <- as.data.frame(model.matrix(m_A))
treatment      <- model.response(model.frame(m_A))


frame <- list(outcome   = DRsims$Y,
              treatment = treatment,
              X_treatment = xmat_treatment) %>%
  lapply(., function(x) split(x, DRsims$group, drop = FALSE))

split_frame <- Map(list,
                   outcome      = frame[['outcome']],
                   treatment    = frame[['treatment']],
                   X_treatment  = frame[['X_treatment']])

#### 1


funcs <- lapply(split_frame, function(x){
  f <- weight_estimator(A = x$treatment, X = x$X_treatment)
})

ipw_funcs <- lapply(split_frame, function(x){
  f <- make_ipw_estimator(Y = x$outcome, A = x$treatment, X = x$X_treatment)
})


system.time(lapply(ipw_funcs, function(f) f(theta, .25, 0)))
system.time(lapply(ipw_funcs, function(f) f(theta, .5, 0)))
system.time(lapply(ipw_funcs, function(f) f(theta, .5, 1)))

system.time(lapply(ipw_funcs, function(f) numDeriv::grad(f, x = theta, alpha = .25, a = 0, method = 'simple')))
system.time(lapply(ipw_funcs, function(f) numDeriv::grad(f, x = theta, alpha = .5, a = 0, method = 'simple')))
system.time(lapply(ipw_funcs, function(f) numDeriv::grad(f, x = theta, alpha = .5, a = 1, method = 'simple')))
#### 3


funcs2 <- lapply(split_frame, function(x){
  f <- weight_estimator2(A = x$treatment, X = x$X_treatment)
})

ipw_funcs2 <- lapply(split_frame, function(x){
  f <- make_ipw_estimator2(Y = x$outcome, A = x$treatment, X = x$X_treatment)
})



system.time(lapply(ipw_funcs2, function(f) numderiv::grad(f, x = theta, alpha = .25, a = 0)))
system.time(lapply(ipw_funcs2, function(f) f(theta, .5, 0)))
system.time(lapply(ipw_funcs2, function(f) f(theta, .5, 1)))

system.time(lapply(ipw_funcs2, function(f) numDeriv::grad(f, x = theta, alpha = .25, a = 0, method = 'simple')))
system.time(lapply(ipw_funcs2, function(f) numDeriv::grad(f, x = theta, alpha = .5, a = 0, method = 'simple')))
system.time(lapply(ipw_funcs2, function(f) numDeriv::grad(f, x = theta, alpha = .5, a = 1, method = 'simple')))

###### Outcome  