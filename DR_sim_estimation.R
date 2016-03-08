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

alphas <- c(0.1, 0.5, 0.9)

DRresults <- DRsims %>% 
   filter(simID <= 2) %>%
  {plyr::dlply(., plyr::.(simID), function(sim){
    this_sim <- dr_ipw(formula_outcome = Y ~ Z1 + Z2 + Z3 + Z4 + A + fA + A*Z1 + fA*Z2,
           formula_interference = Y | A ~ -1 + Z1 + Z2 + Z3 + Z4 +(1|group) | group,
           method_outcome  = geepack::geeglm,
           method_opts_outcome = list(id = quote(group), family = gaussian),
           method_opts_interference = list(allocations = alphas,
                                           method = 'simple', 
                                           runSilent = T),
           dr_term2_function = dr_term2,
           data = sim)
    dr_ipw_estimate(this_sim) %>%
      mutate_(simID = ~ sim$simID[1])
  })} %>%
  bind_rows()

DRresults
# 
# temp <- dr_ipw(formula_outcome = Y ~ A + X1,
#                formula_interference = Y | A ~ X1 + (1|group) | group,
#                method_outcome  = geepack::geeglm,
#                method_opts_outcome = list(id = quote(group), family = gaussian),
#                method_opts_interference = list(allocations = alphas,
#                                                method = 'simple', 
#                                                runSilent = T),
#                dr_term2_function = dr_term2,
#                data = DRsims %>% filter(simID == 2))
# rm(temp)
