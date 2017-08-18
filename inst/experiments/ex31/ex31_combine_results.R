library(dplyr)

load('inst/experiments/ex28/ex28_results.rda')
estimates <- estimates %>% mutate(simid = simid + 700) %>%
  select(-truth)

# estimates <- data_frame()
for(i in 1:700){
  a <- try(load(file = paste0('inst/experiments/ex30/ex30_results/results_', i, '.rda')),
           silent = TRUE)
  if(!is(a, 'try-error')){
    estimates <- bind_rows(estimates, x)
  }
}

estimates <- estimates %>% #select(-bias, -failed, -truth) %>%
  group_by(model_spec, sid, method, regression_type, hajek, a, alpha) %>%
  left_join(oracle, by = c('sid', 'a', 'alpha')) %>%
  mutate(bias      = estimate - truth,
         failed    = is.na(estimate) | is.infinite(estimate) | is.nan(estimate),
         conf.low  = estimate - 1.96 * std_error,
         conf.high = estimate + 1.96 * std_error,
         covered   = conf.low < truth & truth < conf.high)

results <- estimates  %>%
  summarise(mean_est  = mean(estimate * !failed, na.rm = TRUE),
            mean_bias = mean(bias * !failed, na.rm = TRUE),
            coverage  = mean(covered * !failed, na.rm = TRUE),
            failures  = sum(failed),
            ase       = mean(std_error),
            ese       = sd(estimate * !failed, na.rm = TRUE),
            nsims     = n())
