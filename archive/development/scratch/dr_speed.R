

m_A <- lme4::glmer(tru_propen, family = binomial, data = DRsims)
m_Y <- geepack::geeglm(tru_outcome, family = gaussian, data = DRsims, id = ID)
theta1 <- coef(m_Y)
theta2 <- unlist(getME(m_A, c('beta', 'theta')))



xmat_outcome <- as.data.frame(model.matrix(m_Y))
treatment      <- model.response(model.frame(m_Y))
xmat_treatment <- as.data.frame(model.matrix(m_A))
treatment      <- model.response(model.frame(m_A))
rhs_outcome <- get_fixed_formula(m_Y)

frame <- list(outcome   = DRsims$Y,
              treatment = treatment,
              X_treatment = xmat_treatment,
              X_outcome = xmat_outcome) %>%
  lapply(., function(x) split(x, DRsims$group, drop = FALSE))

split_frame <- Map(list,
                   treatment    = frame[['treatment']],
                   X_treatment  = frame[['X_treatment']],
                   outcome      = frame[['outcome']],
                   X_outcome  = frame[['X_outcome']])

####Y, A, X_outcome, X_treatment, formula_outcome
dr_funcs <- lapply(split_frame, function(x){
  make_dr_estimator(Y = x$outcome, A = x$treatment, X_outcome = x$X_outcome, 
                    X_treatment = x$X_treatment, rhs_outcome)
})

system.time({test <- lapply(dr_funcs, function(f) f(c(theta2, theta1), alpha =.5, a = 1))})
system.time({test <- lapply(dr_funcs, function(f) f(c(theta2, theta1), alpha =.5, a = 1))})
system.time({test <- lapply(dr_funcs, function(f) f(c(theta2, theta1), alpha =.5, a = 0))})
system.time({test <- lapply(dr_funcs, function(f) f(c(theta2, theta1), alpha =.25, a = 0))})
