library(dplyr)
library(dr)
library(geex)
alphas <- list(c(.3, .4))

load( pipe( 'ssh saulb@diamond.bios.unc.edu "cat /home/groups/projects/mhudgens/emch/data/R_data/emch_analysis_data.Rdata"' ))
source('inst/cholera_analysis/1a_define_analysis_functions.R')

choleradt <- analysis_c %>%
  group_by(group) %>%
  mutate(fA = mean(A)) %>%
  filter(n() < 100) %>%
  ungroup()



analysis_model_args1 <- list(
  t_model = 
    list(method = lme4::glmer,
         formula = B ~ age + rivkm + (1|group),
         options = list(family = binomial(link = 'logit'), nAGQ = 25L)),
  o_model =
    list(method  = geepack::geeglm,
         formula = y_obs ~ A + fA + age + rivkm,
         options = list(
           family  = binomial(link = 'logit'),
           id      = quote(group))),
  wls_model_0 =
    list(method = stats::glm,
         formula = as.integer(y_obs) ~ fA + age + rivkm,
         options = list(family = quasibinomial(link = 'logit'))),
  wls_model_1 =
    list(method = stats::glm,
         formula = as.integer(y_obs) ~ fA + age + rivkm,
         options = list(family = quasibinomial(link = 'logit')))
)

models0 <- estimate_cholera_parms_step0(
  data        = choleradt,
  model_args  = analysis_model_args1
)

models1 <- estimate_cholera_parms_step1(
  data = choleradt,
  models = models0,
  model_args = analysis_model_args1,
  allocations = c(.3, .4),
  randomization = 2/3
)


#####
# Compare point estimates using V1 and V2 of WLS estimator

splitdt <- split(choleradt, f = choleradt$group)
target_1 <- lapply(splitdt, function(grp_dt){
  # Create estimator function
  estimator <- wls_dbr_estimator(
    data          = grp_dt,
    models        = models1$models,
    randomization = 2/3,
    hajek         = FALSE,
    regression_type = "wls")
  # Evaluate estimator function
  estimator(models1$estimator_args$wls_dbr$theta, alpha = c(.3, .4))
})

target_2 <- lapply(splitdt, function(grp_dt){
  # Create estimator function
  estimator <- wls_dbr_estimator_estFUN(
    data          = grp_dt,
    models        = models1$models,
    randomization = 2/3,
    hajek         = FALSE,
    regression_type = "wls")
  # Evaluate estimator function
  estimator(models1$estimator_args$wls_dbr$theta, alpha = c(.3, .4))
})


target_1 <- target_1 %>% list_matrix() %>% apply(., 2, mean)
target_2 <- target_2 %>% list_matrix() %>% apply(., 2, mean)

target_1
target_2

####
# 

basis_wls <- create_basis(
  estFUN = make_eefun_wls,
  data   = choleradt,
  units  = 'group',
  outer_args = list(models = models1$models, randomization = 2/3)
)

basis_wls@.GFUN(models1$estimator_args$wls_dbr$theta)

basis_g <- create_basis(
  estFUN = generic_eefun,
  data   = choleradt,
  units  = 'group',
  outer_args = list(models = models1$models, randomization = 2/3, estimator_type = 'wls_dbr',
                    regression_type = 'wls')
)

models1$models$wls_model_0[[1]]$weights

####
# Check estimating equation for treatment model

test_eefun <- function(data, model){
  f <- geex::grab_psiFUN(
    model, 
    data = data, 
    numderiv_opts = list(method = 'Richardson'))
  function(theta){
    f(theta)
  }
}

basis_t <- create_basis(
  estFUN = test_eefun,
  data   = choleradt,
  units  = 'group',
  outer_args = list(model = models1$models$t_model)
)

basis_t@.GFUN(theta = models1$estimator_args$ipw$theta)

#####
#

ipwv <- make_ipw_vector(fulldata = choleradt, 
                        models  = models1$models, 
                        group   = 'group', 
                        alpha   = .3,
                        randomization = 2/3)

####
# Check estimating equation for GLM model w/o weights using cholera data ####

tester1 <- glm(formula = y_obs ~ fA + age + rivkm,
              family = binomial(link = 'logit'), data = choleradt)
    
tester2 <- glm(formula = y_obs ~ fA + age + rivkm,
               family = binomial(link = 'logit'), data = choleradt,
               control = glm.control(epsilon = 1e-12))



test_eefun_glm <- function(data, model){
  f <- grab_psiFUN(object = model, data = data)
  function(theta){
    f(theta)
  }
}

basis_tester1 <- create_basis(
  estFUN = test_eefun_glm,
  data   = choleradt,
  units  = 'group',
  outer_args = list(model = tester1)
)

basis_tester2 <- create_basis(
  estFUN = test_eefun_glm,
  data   = choleradt,
  units  = 'group',
  outer_args = list(model = tester2)
)


# These should be zero
basis_tester1@.GFUN(theta = coef(tester1))
basis_tester2@.GFUN(theta = coef(tester2))

# What if I use rootSolve to find roots rather than glm()?
mtester1 <- m_estimate(
  estFUN = test_eefun_glm,
  data   = choleradt,
  units  = 'group',
  root_control = setup_root_control(start = coef(tester1)),
  outer_args = list(model = tester1)
)

# The coefficient estimates are different...
coef(tester1)
roots(mtester1)

mtester1@basis@.GFUN(coef(tester1)) # Bad
mtester1@basis@.GFUN(roots(mtester1)) # Good

# Well, does test_eefun_glm even work? ####
# library(geepack)
# data("vaccinesim", package = "inferference")
# glm_test <- geeglm(A ~ X1 + X2, id = group, data = vaccinesim, family = binomial())
# 
# basis_tester_glm <- create_basis(
#   estFUN = test_eefun_glm,
#   data   = vaccinesim,
#   units  = 'group',
#   outer_args = list(model = glm_test)
# )
# 
# basis_tester_glm@.GFUN(coef(glm_test))

######
#

wls_estFUN <- function(data, wls_formula, ipw_formula, randomization,
                       ipw_theta){
  L <- grab_design_matrix(data, wls_formula)
  X <- grab_design_matrix(data, ipw_formula)
  Y <- grab_response(data, wls_formula)
  A <- grab_response(data, ipw_formula)
  n <- length(Y)
  ipwFUN <- weight_estimator(
    A = A, 
    X = X, 
    randomization = randomization)
  ipw <- ipwFUN(ipw_theta, alpha = .3)/(1 - .3)
  
  function(theta){
    W <- diag(ipw * (A == 0), nrow = n, ncol = n)
    t(L) %*% W %*% (Y - plogis(L %*% theta))
  }
}


mtest_wls <- m_estimate(
  estFUN = wls_estFUN,
  data   = choleradt,
  units  = 'group',
  compute_vcov = FALSE,
  root_control = setup_root_control(start =models1$estimator_args$wls_dbr$theta[c(5:8)]),
  outer_args = list(wls_formula = y_obs ~ fA + age + rivkm, 
                    ipw_formula = B ~ age + rivkm,
                    randomization = 2/3,
                    ipw_theta = models1$estimator_args$wls_dbr$theta[c(1:4)])
)

roots(mtest_wls)
mtest_wls@basis@.GFUN(roots(mtest_wls))

models1$estimator_args$wls_dbr$theta[c(5:8)]

ipwv <- make_ipw_vector(fulldata = choleradt, 
                        models  = models1$models, 
                        group   = 'group', 
                        alpha   = .3,
                        randomization = 2/3)
ipw0 <- ipwv/(1 - .3) * (choleradt$A == 0)

test_glm <- glm(y_obs ~ fA + age + rivkm, data = choleradt, family = quasibinomial(),
    weights = ipw0)

roots(mtest_wls)
coef(test_glm)
######## 

mtest_wls <- m_estimate(
  estFUN = make_eefun_wls,
  data   = choleradt,
  units  = 'group',
  compute_vcov = FALSE,
  root_control = setup_root_control(start =models1$estimator_args$wls_dbr$theta),
  outer_args = list(models = models1$models, randomization = 2/3)
)

mtest_wls@basis@.GFUN(roots(mtest_wls))
roots(mtest_wls)

