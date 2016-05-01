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

source(file = 'DR_functions_5.R')

alphas <- c(0.1, 0.5, 0.9)

tru_outcome <- Y ~ Z1 + Z2 + Z3 + Z4 + A + fA + A*Z1 + fA*Z2
tru_treatment  <- A ~ -1 + Z1 + Z2 + Z3 + Z4 + (1|group) 
mis_outcome <- Y ~ X1 + X2 + X3 + X4 + A + fA + A*X1 + fA*X2
mis_treatment  <- A ~ -1 + X1 + X2 + X3 + X4 + (1|group)

DRsims <- sims_500x_m500_n4 %>% filter(simID < 3)

truth <- data.frame(alpha = alphas, a = target_a, truth = 2+0.5*target_a+6*(alphas*3/4))

### Update to write results to results folder.
# add a comment to each results file
temp <- estimate_sims(DRsims, tru_treatment, tru_outcome)


