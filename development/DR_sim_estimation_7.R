#------------------------------------------------------------------------------#
#   Title: Doubly Robust simulation estimation
#  Author: B. Saul
#    Date: 2016-06-13
# Purpose: estimation for Lan's DR simulations using functions_V3
#          - this tries to calculate the weights differently
#------------------------------------------------------------------------------#

library(inferference)
library(lme4)
library(dplyr)
library(geepack)
library(magrittr)
library(sandwich)
library(sandwichShop)
library(doMC)
registerDoMC(4)

alphas <- c(0.1, 0.5, 0.9)

tru_outcome <- Y ~ Z1 + Z2 + Z3 + Z4 + A + fA + A*Z1 + fA*Z2
tru_treatment  <- A ~ Z1 + Z2 + Z3 + Z4 + (1|group) 
mis_outcome <- Y ~ X1 + X2 + X3 + X4 + A + fA + A*X1 + fA*X2
mis_treatment  <- A ~ X1 + X2 + X3 + X4 + (1|group)

# DRsims <- DRsims %>% filter(simID == 1) # for testing

target_a <- 1

truth <- data.frame(alpha = alphas, a = target_a, truth = 2+0.5*target_a+6*(alphas*3/4))

### Update to write results to results folder.
# add a comment to each results file
# temp <- estimate_sims(DRsims, tru_treatment, tru_outcome)

tru_tru <- estimate_sims(DRsims %>% filter(simID < 100), 
                         formula_treatment = tru_treatment, 
                         formula_outcome = tru_outcome, 
                         parallel = TRUE) %>%
  left_join(truth, by = c('alpha', 'a')) %>%
  mutate(conf.low  = estimate - 1.96 * std.error,
         conf.high = estimate + 1.96 * std.error, 
         bias      = estimate - truth,
         covered   = truth > conf.low & truth < conf.high * 1)

save(tru_tru, file = paste0('results/', experimentID, '_', Sys.Date(), '.rda'))


#### Results ####

# r1 <- tru_tru %>%
#   group_by(estimator_type, alpha) %>%
#   summarise(mean   = mean(estimate),
#             bias   = mean(bias),
#             coverage = mean(covered),
#             ase    = mean(std.error),
#             ese    = sd (estimate)) %>%
#   mutate(mu = 'True',
#          pi = 'True') 
