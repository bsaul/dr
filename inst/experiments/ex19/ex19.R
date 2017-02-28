#------------------------------------------------------------------------------#
#   Title: Development experiment 19
#  Author: B. Saul
#    Date: 2017-02-25
# Purpose: Understand and pinpoint bias observed in previous experiments
#------------------------------------------------------------------------------#
library(dplyr)
library(dr)
library(doMC)
registerDoMC(4)
source('development/experiments/ex19/ex19_funs.R')
source('development/experiments/ex19/ex19_settings.R')
nsims <- 500
which_scenario <- 5
## Generate simulations ##
simulations <- Map(gen_sim, 
    gamma = scenarios$gamma, 
    theta = scenarios$theta,
    beta  = scenarios$beta,
    m     = scenarios$m,
    ni    = scenarios$ni,
    nsims = nsims
)



## Compute estimates ##

# ptm <- proc.time()
estimates <- plyr::ldply(1:nsims, .parallel = TRUE, function(b){
  do_sim(which_scenario, c(0.1, 0.5, 0.9), margs) %>%
    mutate_(simid = ~ b)
}) %>%
  bind_rows()


# estimates <- plyr::ldply(
#   simulations[[which_scenarios]], .parallel = TRUE, .progress = 'text', function(sim) {
#   lapply(seq_along(margs), function(j){
#     est_sim(
#       simdt = sim,
#       allocations = c(0.1, 0.5, 0.9), 
#       model_args = margs[[j]]) %>%
#       mutate_(
#         sid = ~which_scenarios,
#         model_spec = ~names(margs)[j])
#   }) %>% bind_rows()
# })
# #}) %>% bind_rows() 
# # proc.time() - ptm
  
estimates <- estimates %>% #select(-bias, -failed, -truth) %>%
  group_by(model_spec, sid, method, hajek, a, alpha) %>%
  left_join(oracle, by = c('sid', 'a', 'alpha')) %>%
  mutate(bias      = estimate - truth,
         failed    = is.na(estimate) | is.infinite(estimate) | is.nan(estimate),
         conf.low  = estimate - 1.96 * std_error,
         conf.high = estimate + 1.96 * std_error,
         covered   = conf.low < truth & truth < conf.high)

## Compute results ## 
results <- estimates  %>%
  summarise(mean_est  = mean(estimate * !failed, na.rm = TRUE),
            mean_bias = mean(bias * !failed, na.rm = TRUE),
            coverage  = mean(covered * !failed, na.rm = TRUE),
            failures  = sum(failed),
            ese       = sd(estimate * !failed, na.rm = TRUE))

save(estimates, scenarios, oracle, results, nsims, margs,
     file = 'development/experiments/ex19/ex19_results.rda')

