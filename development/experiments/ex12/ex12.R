#------------------------------------------------------------------------------#
#   Title: Development experiement 12 
#  Author: B. Saul
#    Date: 2017-02-15
# Purpose: Understand and pinpoint bias observed in previous experiments
#------------------------------------------------------------------------------#
library(dplyr)
library(dr)
library(doMC)
registerDoMC(4)
source('development/experiments/ex12/ex12_funs.R')

nsims <- 300

## Generate scenarios ##
scenarios <- data_frame(
  sid = 1:12,
  gamma = rep(list(c(0.1, 0.2, 0),
                   c(2.1, 0.4, -4.2),
                   c( -1, 0.5, -1), 
                   c(0  , 0.75, 0)),
              times = 3),
  theta = rep(list(.3, .3, .3, 1.5), 3),
  beta  = rep(list(c(2, 0, 0, -1.5, 2), 
                   c(2, 2, 0, -1.5, 2),
                   c(2, 2, 1, -1.5, 2)), 
              each = 4)
)

## Oracle ## 
oracle <- expand.grid(
  sid    = 1:12,
  alpha = c(0.1, 0.5, 0.9), 
  a     = c(0, 1, NA)) %>%
  ## sqrt(2/pi) is mean of half normal distr when sigma = 1
  mutate(truth = 3 + (-1.5*sqrt(2/pi)) + I(sid > 4) * a * 2 + I(sid > 8) * alpha * 1) %>%
  mutate(truth = ifelse(!is.na(a), truth,
                        3 + (-1.5*sqrt(2/pi)) + I(sid > 4) * 1 + I(sid > 8) * alpha * 1))
  
# to do: add truth for marginal means

## Component model arguments ##
tcor_ocor_margs <- list(
  tmodel = 
    list(method = lme4::glmer,
         formula = A ~ Z1_abs + Z2 + (1|group),
         options = list(family = binomial(link = 'logit'))),
  omodel =
    list(method  = geepack::geeglm,
         formula = Y ~ A + fA + Z1_abs + Z2,
         options = list(
           family  = gaussian(link = 'identity'),
           id      = quote(group)))
  )

tmis_ocor_margs <- list(
  tmodel = 
    list(method = lme4::glmer,
         formula = A ~ Z1 + Z2 + (1|group),
         options = list(family = binomial(link = 'logit'))),
  omodel =
    list(method  = geepack::geeglm,
         formula = Y ~ A + fA + Z1_abs + Z2,
         options = list(
           family  = gaussian(link = 'identity'),
           id      = quote(group)))
)

tcor_omis_margs <- list(
  tmodel = 
    list(method = lme4::glmer,
         formula = A ~ Z1_abs + Z2 + (1|group),
         options = list(family = binomial(link = 'logit'))),
  omodel =
    list(method  = geepack::geeglm,
         formula = Y ~ A + fA + Z1 + Z2,
         options = list(
           family  = gaussian(link = 'identity'),
           id      = quote(group)))
)

tmis_omis_margs <- list(
  tmodel = 
    list(method = lme4::glmer,
         formula = A ~ Z1 + Z2 + (1|group),
         options = list(family = binomial(link = 'logit'))),
  omodel =
    list(method  = geepack::geeglm,
         formula = Y ~ A + fA + Z1 + Z2,
         options = list(
           family  = gaussian(link = 'identity'),
           id      = quote(group)))
)

margs <- list(
  tcor_ocor = tcor_ocor_margs,
  tmis_ocor = tmis_ocor_margs,
  tcor_omis = tcor_omis_margs,
  tmis_omis = tmis_omis_margs 
)

## Generate simulations ##
simulations <- Map(gen_sims, 
    gamma = scenarios$gamma, 
    theta = scenarios$theta,
    beta  = scenarios$beta,
    m     = 300,
    ni    = 20,
    nsims = nsims)

## Compute estimates ##
estimates <- lapply(seq_along(simulations), function(k) {
  lapply(seq_along(margs), function(j){
    est_sims(
      simulations[[k]], 
      allocations = c(0.1, 0.5, 0.9), 
      model_args = margs[[j]]) %>%
      mutate_(
        sid = ~k,
        model_spec = ~names(margs)[j])
  }) %>% bind_rows()
}) %>% bind_rows() 
  
estimates <- estimates %>% #select(-bias, -failed, -truth) %>%
  group_by(model_spec, sid, method, a, alpha) %>%
  left_join(oracle, by = c('sid', 'a', 'alpha')) %>%
  mutate(bias      = estimate - truth,
         failed    = is.na(estimate) | is.infinite(estimate))

## Compute results ## 
results <- estimates  %>%
  summarise(mean_est  = mean(estimate * !failed, na.rm = TRUE),
            mean_bias = mean(bias * !failed, na.rm = TRUE),
            failures  = sum(failed),
            ese       = sd(estimate))

save(estimates, scenarios, oracle, results, 
     file = 'development/experiments/ex12/ex12_results.rda')

