#------------------------------------------------------------------------------#
#   Title: Experiment 03: comparing the speed of inferference to DR_functions_5
#  Author: B. Saul
#    Date: 2016-05-1
# Purpose: comparing the speed of inferference to DR_functions_5
#------------------------------------------------------------------------------#
library(inferference)
library(lme4)
library(dplyr)
library(geepack)
library(magrittr)
library(sandwich)
library(sandwichShop)

DRsims <- sims_500x_m500_n4 %>% filter(simID < 2)
tru_outcome <- Y ~ Z1 + Z2 + Z3 + Z4 + A + fA + A*Z1 + fA*Z2
tru_treatment  <- A ~ Z1 + Z2 + Z3 + Z4 + (1|group) 

# This isn't exactly a fair comparison since interference computes all causal 
# contrasts but it should give some sense in comparing the key functions

system.time(t1 <- estimation(tru_treatment, tru_outcome, filter(DRsims, simID == 1),
                               allocations = c(.1, .5, .9), target_a = 1))

# user  system elapsed 
# 81.953   0.497  82.723 



system.time(t2 <- interference(Y | A ~ Z1 + Z2 + Z3 + Z4 + (1|group) | group,
                         data = filter(DRsims, simID == 1),
                         allocations = c(.1, .5, .9),
                         method = 'simple'))


# user  system elapsed 
# 87.672   0.683  88.718 

# compare esimates for ipw
t1 %>% filter(estimator_type == 'ipw')
t2$estimates %>% filter(effect == 'outcome', trt1 == 1) %>% select(alpha1, trt1, estimate, std.error)
