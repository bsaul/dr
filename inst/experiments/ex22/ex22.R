#------------------------------------------------------------------------------#
#   Title: Development ex22
#  Author: B. Saul
#    Date: 2017-02-25
# Purpose: Understand and pinpoint bias observed in previous experiments
#------------------------------------------------------------------------------#
library(dplyr)
library(dr)
library(doMC)
registerDoMC(4)
source('inst/experiments/ex22/ex22_funs.R')
source('inst/experiments/ex22/ex22_settings.R')
nsims <- 350
which_scenarios <- 9

## Compute estimates ##
allocations <- c(0.1, 0.5, 0.9)
estimates <- do_scenarios(nsims, which_scenarios, allocations, 
                          all_model_args = margs,
                          compute_se = FALSE,
                          .parallel = TRUE) %>%
   bind_rows()

estimates <- estimates %>% #select(-bias, -failed, -truth) %>%
  group_by(model_spec, sid, method, regression_type, hajek, a, alpha) %>%
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
     file = 'inst/experiments/ex22/ex22_results.rda')

