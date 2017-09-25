ipws <- make_ipw_vector(choleradt, 
                models = models0,
                group  = 'group',
                randomization = 2/3,
                alpha = .3)
ipw0 <- ipws * (choleradt$A == 0)

models1 <- estimate_cholera_parms_step1(
  data = choleradt,
  models = models0,
  model_args = analysis_model_args,
  allocations = c(.4, .3),
  randomization = 2/3
)

model_args <- analysis_model_args

# wls_model_test <- glm(
#   as.integer(y_obs) ~ fA + age + rivkm,
#   weights = ipws * (choleradt$A == 0),
#   data    = choleradt,
#   family  = quasibinomial(link = "logit")
# )

test_ee <- make_eefun_wls(
  models = models1$models,
  data = choleradt %>% filter(group == "100"),
  randomization = 2/3)

test_est <- wls_dbr_estimator_estFUN(
  models = models1$models,
  data = choleradt %>% filter(group == "100"),
  randomization = 2/3
)

test_est(c(theta_t, theta_wls), c(.4, .3))

theta_t <- unlist(lme4::getME(models1$models$t_model, c('beta', 'theta')))
theta_wls <- lapply(append(models1$models$wls_model_0, models1$models$wls_model_1),
                    function(x) coef(x)) %>% unlist()

length(theta_wls)
test_ee(c(rep(1, 4), theta_wls))


test_ee2 <- generic_eefun(
  data = choleradt %>% filter(group == "100"),
  models = models1$models,
  randomization = 2/3,
  estimator_type = 'wls_dbr',
  regression_type = 'wls'
)

x <- c(rep(1, 4), theta_wls, 0, 0, 0, 0)

A1 <- numDeriv::jacobian(test_ee2, x = x, alpha = c(.4, .3))
A1

solve(A1) %*% t(solve(A1))

results_hold <- results

results2[[1]][[1]]$wls_dbr$vcov %>% diag()
results2[[1]][[1]]$wls_dbr$Ai[[1]]
results_hold[[1]][[1]]$wls_dbr$vcov %>% diag()
results_hold[[1]][[1]]$wls_dbr$Ai[[1]]
