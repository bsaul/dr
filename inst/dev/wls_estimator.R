
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


wls_preprocess <- function(data, models, randomization = 1, a, alpha){
  splitdt <- split(data, data$group)
  lapply(splitdt, function(x){
    comp <- extract_model_info(model = models, data = x, 'dbr')
    Y   <- comp$Y
    A   <- comp$A
    N   <- comp$N
    X_o <- comp$X_o
    f_o <- terms.formula(comp$rhs_o)
    L   <- model.matrix(drop.terms(f_o, attr(f_o, 'term.labels') == 'A'),
                        data = x)

    ## components for IPW part
    ip_fun <- weight_estimator(
      A = comp$A, 
      X = comp$X_t, 
      randomization = randomization)
    
    IPW <- ip_fun(unlist(lme4::getME(models$t_model, c('beta', 'theta'))), alpha = alpha)
    (A == a) * IPW / dbinom(A, 1, prob = alpha)

  }) %>%
    unlist()

}

mods$wls_model_0 <- glm(Y ~ fA + Z1_abs + Z2 + Z1_abs*Z2, 
                        family  = gaussian(link = 'identity'),
                        weights = wls_preprocess(exdt, mods, a = 0, alpha= .6),
                        data    = exdt)
mods$wls_model_1 <- glm(Y ~ fA + Z1_abs + Z2 + Z1_abs*Z2, 
                      family  = gaussian(link = 'identity'),
                      weights = wls_preprocess(exdt, mods, a = 1, alpha= .6),
                      data    = exdt)

theta_oo <- unlist(lme4::getME(mods$t_model, c('beta', 'theta')))

fun1 <- wls_dbr_estimator(data = exdt %>% filter(group == 1),
                  models = mods,
                  randomization = 1)

fun1 <- otc_estimator(data = exdt %>% filter(group == 1),
                          models = mods,
                          randomization = 1)

fun1(c(theta_oo, coef(mods$wls_model_0), coef(mods$wls_model_1), 0, 0), alpha = .5)
fun1(c(coef(mods$o_model), 0, 0), alpha = .5)

fun_ee <- generic_eefun(
  data = exdt %>% filter(group == 1),
  models = mods,
  estimator_type = 'wls_dbr',
  randomization = 1)
fun_ee(c(theta_oo, coef(mods$wls_model_0), coef(mods$wls_model_1), 0, 0), alpha = .5)
