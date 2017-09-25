#------------------------------------------------------------------------------#
#   Title: Analysis of cholera data with IPW, OTC, and DR estimators
#  Author: B. Saul
#    Date: 2016-05-20
# Purpose:
#     Log:
# 20170623 - running all but WLS estimator; will do seperately as I'm getting
#            an error from the WLS estimator when estimating the variance
#------------------------------------------------------------------------------#

library(dplyr)
library(dr)
library(geex)
alphas <- lapply(seq(.3, .6, by = .02), function(x) sort(c(.4, x)))
load( pipe( 'ssh saulb@diamond.bios.unc.edu "cat /home/groups/projects/mhudgens/emch/data/R_data/emch_analysis_data.Rdata"' ))


# Prepare data
choleradt <- analysis_c %>%
  group_by(group) %>%
  mutate(fA = mean(A)) %>%
  filter(n() < 100) %>% # for computational speed
  ungroup()

# Component model arguments
analysis_model_args <- list(
  t_model = 
    list(method = lme4::glmer,
         formula = B ~ age + rivkm + (1|group),
         options = list(family = binomial(link = 'logit'))),
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

# Treatment and outcome models
models0 <- estimate_cholera_parms_step0(
  data        = choleradt,
  model_args  = analysis_model_args
)

# Results for all but WLS ####
# (setting skipit = TRUE for wls in estimator_args)
source("inst/cholera_analysis/1a_define_analysis_functions.R")
results_nonwls <- lapply(seq_along(alphas), function(i){
  attempt <- try(estimate_cholera_parms_step2(
    data        = choleradt,
    allocations = alphas[i],
    models      = models0,
    model_args  = analysis_model_args,
    randomization  = 2/3,
    compute_se  = TRUE
  ))
  if(is(attempt, 'try-error')){
    NULL
  } else{
    attempt
  }
})

save(results_nonwls, file = paste0('inst/cholera_analysis/cholera_results_', Sys.Date(),'.rda'))

# Results for WLS ####
# (setting skipit = FALSE for wls in estimator_args)
source("inst/cholera_analysis/1a1_define_wls_analysis_functions.R")
results_wls <- lapply(seq_along(alphas), function(i){
  attempt <- try(estimate_cholera_parms_step2_wls(
    data        = choleradt,
    allocations = alphas[i],
    models      = models0,
    model_args  = analysis_model_args,
    randomization  = 2/3,
    compute_se  = TRUE
  ))
  if(is(attempt, 'try-error')){
    NULL
  } else{
    attempt
  }
})

save(results_wls, file = paste0('inst/cholera_analysis/cholera_results_wls_', Sys.Date(),'.rda'))
  
