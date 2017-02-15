#------------------------------------------------------------------------------#
#   Title: Development experiement 11 
#  Author: B. Saul
#    Date: 2017-02-15
# Purpose: Understand and pinpoint bias observed in previous experiments
#------------------------------------------------------------------------------#
library(dplyr)
library(dr)
library(doMC)
registerDoMC(4)

# Generate scenarios
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

# Oracle 
oracle <- expand.grid(
  sid    = 1:12,
  alpha = c(0.1, 0.5, 0.9), 
  a     = c(0, 1, NA)) %>%
  mutate(truth = 3 + I(sid > 4) * a * 2 + I(sid > 8) * alpha * 1)

# Component model arguments
margs <- list(
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

# Generate simulations
simulations <- Map(gen_sims, 
    gamma = scenarios$gamma, 
    theta = scenarios$theta,
    beta  = scenarios$beta,
    m     = 300,
    ni    = 20,
    nsims = 250)

# 
estimates <- lapply(seq_along(simulations), function(k) {
  est_sims(
    simulations[[k]], 
    allocations = c(0.1, 0.5, 0.9), 
    model_args = margs) %>%
    mutate(sid = k)
}) %>% bind_rows()

results <- estimates %>% group_by(sid, method, a, alpha) %>%
  left_join(oracle, by = c('sid', 'a', 'alpha')) %>%
  mutate(bias      = estimate - truth) %>%
  summarise(mean_est  = mean(estimate),
            mean_bias = mean(bias),
            ese       = sd(estimate))

