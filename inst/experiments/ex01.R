#------------------------------------------------------------------------------#
#   Title: Experiment 01: speeding up the outcome estimator
#  Author: B. Saul
#    Date: 2016-04-30
# Purpose: develop faster outcome estimators
#------------------------------------------------------------------------------#

library(lme4)
library(dplyr)
library(geepack)
library(magrittr)
library(sandwich)
library(sandwichShop)
library(lineprof)

DRsims <- sims_500x_m500_n4 %>%  filter(simID < 2)

tru_outcome <- Y ~ Z1 + Z2 + Z3 + Z4 + A + fA + A*Z1 + fA*Z2

m_Y <- geeglm(tru_outcome, data = DRsims, id = ID)

theta_o <- coef(m_Y)
rhs_outcome <- formula(m_Y)[-2]

split_dt <- split(DRsims, DRsims$group)
g1_dt <- split_dt[[1]]

#------------------------------------------------------------------------------#
#### Experiment 01-00 ####
# profile functions to see where bottle necks are
# Currently evaluation of otc_estimator takes ~ 0.02s on a single dataset;
# this means ~ 9s for 500 groups. For grad(method = 'Richardson'), 
# this means ~ 1.2 sec per group or 10 minutes.
#------------------------------------------------------------------------------#

source('development/experiments/ex01_functions.R')

system.time(ex01_00_1())
system.time(ex01_00_2())
system.time(ex01_00_3())

ex01_00_1_profile <- lineprof(ex01_00_1())
shine(ex01_00_1_profile)

#------------------------------------------------------------------------------#
#### Experiment 01-01 ####
# move mutate(fitted) from expand_frame to make_otc.
#------------------------------------------------------------------------------#

source('development/experiments/ex01_functions.R')

system.time(ex01_01_1())
system.time(ex01_01_2())
system.time(ex01_01_3())

# Slightly faster in terms of evaluating the estimator
microbenchmark::microbenchmark(ex01_01_1(), ex01_00_1())

# Shaves off ~ 1/3 of the time in terms of taking derivator of estimator
microbenchmark::microbenchmark(ex01_01_3(), ex01_00_3())

#------------------------------------------------------------------------------#
#### Experiment 01-02 ####
# use base R functions to take means/sums
#------------------------------------------------------------------------------#

source('development/experiments/ex01_functions.R')

# Get recursion error when expand_outcome_frame_01 is memoised
system.time(ex01_02_1())
# Faster ~ 30% in terms of evaluating the estimator compared to 01
# Hoever, is it as accurate?
microbenchmark::microbenchmark(ex01_02_1(), ex01_01_1(), ex01_00_1())


#  Roughly 50% faster in terms of taking derivatives
system.time(ex01_02_3())
microbenchmark::microbenchmark(ex01_02_3(), ex01_01_3(), ex01_00_3())

#------------------------------------------------------------------------------#
#### Experiment 01-03 ####
# Follow up on 02 to get accurate values
#------------------------------------------------------------------------------#

source('development/experiments/ex01_functions.R')

# Compare value from 01 and 03
identical(ex01_01_1(), ex01_03_1())

system.time(ex01_03_1())

# Shaves ~2 seconds off all the groups
system.time(ex01_01_2())
system.time(ex01_03_2())

# Faster ~ 40% in terms of evaluating the estimator compared to 01
microbenchmark::microbenchmark(ex01_03_1(), ex01_01_1(), ex01_00_1())

#  Roughly 50% faster in terms of taking derivatives
system.time(ex01_03_3())
microbenchmark::microbenchmark(ex01_03_3(), ex01_01_3(), ex01_00_3(), times = 50L)
