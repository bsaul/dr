
library(dr)
library(dplyr)
exdt <- gen_data(20, 10, c(0.1, 0.2, 0, .2), 1, c(2, 2, 1, -1.5, 2, -3))

mods <- make_models(model_args = list(
    t_model = 
      list(method = lme4::glmer,
           formula = A ~ Z1_abs + Z2 + Z1_abs*Z2 + (1|group),
           options = list(family = binomial(link = 'logit'))),
    o_model =
      list(method  = geepack::geeglm,
           formula = Y ~ A + fA + Z1_abs + Z2 + Z1_abs*Z2,
           options = list(
             family  = gaussian(link = 'identity'),
             id      = quote(group)))),
    data = exdt)

theta_tt <- unlist(lme4::getME(mods$t_model, c('beta', 'theta')))
theta_oo <- coef(mods$o_model)

fun0 <- dbr_estimator(
  data = exdt %>% filter(group == 1),
  models = mods,
  randomization = 1)

fun0(c(theta_tt, theta_oo, 0, 0), alpha = .5)

fun_ee0 <- generic_eefun(
  data = exdt %>% filter(group == 1),
  models = mods,
  estimator_type = 'dbr',
  randomization = 1)
fun_ee0(c(theta_tt, theta_oo, 0, 0), alpha = .5)

mods$wls_model_0 <- glm(Y ~ fA + Z1_abs + Z2 + Z1_abs*Z2, 
                        family  = gaussian(link = 'identity'),
                        weights = (exdt$A == 0) * make_ipw_vector(exdt, mods, 'group', a = 0, alpha= .5),
                        data    = exdt)
mods$wls_model_1 <- glm(Y ~ fA + Z1_abs + Z2 + Z1_abs*Z2, 
                      family  = gaussian(link = 'identity'),
                      weights = (exdt$A == 1) * make_ipw_vector(exdt, mods, 'group', a = 1, alpha= .5),
                      data    = exdt)

fun1 <- reg_dbr_estimator(
  data = exdt %>% filter(group == 1),
  models = mods,
  regression_type = 'wls',
  randomization = 1)

fun1(c(theta_tt, coef(mods$wls_model_0), coef(mods$wls_model_1), 0, 0), alpha = .5)


fun_ee1 <- generic_eefun(
  data = exdt %>% filter(group == 1),
  models = mods,
  estimator_type = 'reg_dbr',
  randomization = 1,
  regression_type = 'wls')
fun_ee1(c(theta_tt, coef(mods$wls_model_0), coef(mods$wls_model_1), 0, 0), alpha = .5)

####  regression-on-propensity

exdt <- exdt %>%
  mutate_(
    ipw0 =~ make_ipw_vector(exdt, mods, 'group', a = 0, alpha= .5),
    ipw1 =~ make_ipw_vector(exdt, mods, 'group', a = 1, alpha= .5)
  )

mods$pcov_model_0 <- glm(Y ~ fA + Z1_abs + Z2 + Z1_abs*Z2 + ipw0, 
                        family  = gaussian(link = 'identity'),
                        weights = (exdt$A == 0) * 1,
                        data    = exdt)
mods$pcov_model_1 <- glm(Y ~ fA + Z1_abs + Z2 + Z1_abs*Z2 + ipw1, 
                        family  = gaussian(link = 'identity'),
                        weights = (exdt$A == 1) * 1,
                        data    = exdt)
fun2 <-  reg_dbr_estimator(
  data = exdt %>% filter(group == 1),
  models = mods,
  regression_type = 'pcov',
  randomization = 1)

fun2(c(theta_tt, coef(mods$pcov_model_0), coef(mods$pcov_model_1), 0, 0), alpha = .5)

fun_ee2 <- generic_eefun(
  data = exdt %>% filter(group == 1),
  models = mods,
  estimator_type = 'reg_dbr',
  randomization = 1,
  regression_type = 'pcov')
fun_ee2(c(theta_tt, coef(mods$pcov_model_0), coef(mods$pcov_model_1), 0, 0), alpha = .5)

