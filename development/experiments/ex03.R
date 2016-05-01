#------------------------------------------------------------------------------#
#   Title: Experiment 03: comparing the speed of inferference to DR_functions_5
#  Author: B. Saul
#    Date: 2016-05-1
# Purpose: comparing the speed of inferference to DR_functions_5
#------------------------------------------------------------------------------#

library(lme4)
library(dplyr)
library(geepack)
library(magrittr)
library(sandwich)
library(sandwichShop)
library(lineprof)

# This isn't exactly a fair comparison since interference computes all causal 
# contrasts but it should give some sense in comparing the key functions

system.time(estimation(tru_treatment, tru_outcome, filter(DRsims, simID == 1),
                               allocations = c(.1, .5, .9), target_a = 1))

# user  system elapsed 
# 81.953   0.497  82.723 

system.time(interference(Y | A ~ -1 + Z1 + Z2 + Z3 + Z4 + (1|group) | group,
                         data = filter(DRsims, simID == 1),
                         allocations = c(.1, .5, .9),
                         method = 'simple'))
# user  system elapsed 
# 87.672   0.683  88.718 