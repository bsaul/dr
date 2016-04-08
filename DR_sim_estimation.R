#------------------------------------------------------------------------------#
#   Title: Doubly Robust simulation estimation
#  Author: B. Saul
#    Date: 2016-03-07
# Purpose: estimation for Lan's DR simulations
#------------------------------------------------------------------------------#

library(inferference)
library(dplyr)
library(geepack)
library(magrittr)
library(sandwich)
library(sandwichShop)
library(doMC)
registerDoMC(4)
source(file = 'DR_functions.R')

DRsims <- sims_500x_m500_n4
alphas <- c(0.1, 0.5, 0.9)

tru_outcome <- Y ~ Z1 + Z2 + Z3 + Z4 + A + fA + A*Z1 + fA*Z2
tru_propen  <- Y | A ~ -1 + Z1 + Z2 + Z3 + Z4 + (1|group) | group
mis_outcome <- Y ~ X1 + X2 + X3 + X4 + A + fA + A*X1 + fA*X2
mis_propen  <- Y | A ~ -1 + X1 + X2 + X3 + X4 + (1|group) | group  

# Outcome_Propensity
tru_tru <- dr_results(DRsims, tru_outcome, tru_propen)
mis_tru <- dr_results(DRsims, mis_outcome, tru_propen)
tru_mis <- dr_results(DRsims, tru_outcome, mis_propen)
mis_mis <- dr_results(DRsims, mis_outcome, mis_propen)

## True values
truth <- data.frame(alpha1 = alphas, truth = 2+0.5*1+6*(alphas*3/4))
#save.image()
# inter_size=(n_i-1)
# bar_y_alpha_1=2+0.5*alpha_1+6*(alpha_1*inter_size/n_i)*inter
# bar_y1_alpha_1=2+0.5*1+6*(alpha_1*inter_size/n_i)*inter
# bar_y0_alpha_1=2+0.5*0+6*(alpha_1*inter_size/n_i)*inter

r1 <- tru_tru %>%
  group_by(type, alpha1) %>%
  summarise(mean   = mean(estimate)) %>%
  mutate(mu = 'True',
         pi = 'True') 

r2 <- tru_mis %>%
  group_by(type, alpha1) %>%
  summarise(mean   = mean(estimate)) %>%
  mutate(mu = 'True',
         pi = 'Mis')

r3 <- mis_tru %>%
  group_by(type, alpha1) %>%
  summarise(mean   = mean(estimate)) %>%
  mutate(mu = 'Mis',
         pi = 'True')

r4 <- mis_mis%>%
  group_by(type, alpha1) %>%
  summarise(mean   = mean(estimate)) %>%
  mutate(mu = 'Mis',
         pi = 'Mis') 

results <- rbind(r1, r2, r3, r4) %>%
  left_join(truth, by = 'alpha1') %>%
  mutate(bias = mean - truth) %>%
  filter(!(type == 'ipw' & mu == pi)) %>%
  arrange(type)

save.image()
  