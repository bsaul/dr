#------------------------------------------------------------------------------#
#   Title: Doubly Robust simulation estimation
#  Author: B. Saul
#    Date: 2016-04-26
# Purpose: estimation for Lan's DR simulations
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

source(file = 'DR_functions_3.R')

alphas <- c(0.1, 0.5, 0.9)

tru_outcome <- Y ~ Z1 + Z2 + Z3 + Z4 + A + fA + A*Z1 + fA*Z2
tru_treatment  <- A ~ -1 + Z1 + Z2 + Z3 + Z4 + (1|group) 
mis_outcome <- Y ~ X1 + X2 + X3 + X4 + A + fA + A*X1 + fA*X2
mis_treatment  <- A ~ -1 + X1 + X2 + X3 + X4 + (1|group)

target_a = 1
## True values
truth <- data.frame(alpha = alphas, a = target_a, truth = 2+0.5*target_a+6*(alphas*3/4))

DRsims <- sims_500x_m500_n4 %>%  filter(simID < 2)

# Treatment_Outcome
tru_tru <- estimate_sims(tru_treatment, tru_outcome) %>%
  left_join(truth, by = c('alpha', 'a')) %>%
  mutate(bias = estimate - truth,
         covered = truth > conf.low & truth < conf.high * 1)
 
save(tru_tru, file = paste0('results/tru_tru', Sys.Date(), '.rda'))

mis_tru <- estimate_sims(mis_treatment, tru_outcome)%>%
  left_join(truth, by = c('alpha', 'a')) %>%
  mutate(bias = estimate - truth,
         covered = truth > conf.low & truth < conf.high * 1)

save(mis_tru, file = paste0('results/mis_tru', Sys.Date(), '.rda'))

tru_mis <- estimate_sims(tru_treatment, mis_outcome)%>%
  left_join(truth, by = c('alpha', 'a')) %>%
  mutate(bias = estimate - truth,
         covered = truth > conf.low & truth < conf.high * 1)

save(tru_mis, file = paste0('results/tru_mis', Sys.Date(), '.rda'))

mis_mis <- estimate_sims(mis_treatment, mis_outcome)%>%
  left_join(truth, by = c('alpha', 'a')) %>%
  mutate(bias = estimate - truth,
         covered = truth > conf.low & truth < conf.high * 1)

save(mis_mis, file = paste0('results/mis_mis', Sys.Date(), '.rda'))

r1 <- tru_tru %>%
  group_by(estimator, alpha) %>%
  summarise(mean   = mean(estimate),
            bias   = mean(bias),
            coverage = mean(covered),
            ase    = mean(std.error)) %>%
  mutate(mu = 'True',
         pi = 'True') 

r2 <- tru_mis %>%
  group_by(estimator, alpha) %>%
  summarise(mean   = mean(estimate),
            bias   = mean(bias),
            coverage = mean(covered)) %>%
  mutate(pi = 'True', mu = 'Mis')

r3 <- mis_tru %>%
  group_by(estimator, alpha) %>%
  summarise(mean   = mean(estimate),
            bias   = mean(bias),
            coverage = mean(covered)) %>%
  mutate(pi = 'Mis', mu = 'True')

r4 <- mis_mis%>%
  group_by(estimator, alpha) %>%
  summarise(mean   = mean(estimate),
            bias   = mean(bias),
            coverage = mean(covered)) %>%
  mutate(pi = 'True', mu = 'True')


results <- rbind(r1, r2, r3, r4) %>%
  arrange(type)

save.image()


