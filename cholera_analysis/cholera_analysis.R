#------------------------------------------------------------------------------#
#   Title: Analysis of cholera data with IPW, OTC, and DR estimators
#  Author: B. Saul
#    Date: 2016-05-20
# Purpose: 
#------------------------------------------------------------------------------#

library(inferference)
library(dplyr)
library(sandwich)
library(sandwichShop)
# Load necessary functions
source('R/functions_v2.R', echo = T, max.deparse.length = 5000)

# alphas <- c(.3, .45)
alphas <- seq(.3, .6, by = .01)


# choleradt <- vaccinesim %>%
#   group_by(group) %>%
#   mutate(fA = mean(A))
# 
# results <- estimation(treatment_formula = A ~ X1 + X2 + (1|group), 
#                       outcome_formula = y ~ A + fA + A*fA + X1 + X2,
#                       outcome_method = glm,
#                       outcome_method_opts = list(family = binomial),
#                       this_data = choleradt,
#                       allocations = alphas,
#                       target_a = 0)

load( pipe( 'ssh saulb@diamond.bios.unc.edu "cat /home/groups/projects/mhudgens/emch/data/R_data/emch_analysis_data.Rdata"' ))

choleradt <- analysis_c %>%
  group_by(group) %>%
  mutate(fA = mean(A)) %>%
  filter(n() < 1074)

results <- estimation(treatment_formula = A ~ age + rivkm + (1|group), 
                      outcome_formula = y_obs ~ A + fA + A*fA + age + rivkm,
                      outcome_method = glm,
                      outcome_method_opts = list(family = binomial),
                      this_data = choleradt,
                      allocations = alphas,
                      target_a = 0) 

results %>%
  filter(estimator_type == 'dbr') %>%
  arrange(estimate)

results %>%
  filter(alpha == 0.3)

results %>%
  filter(n_i > 1074, alpha == 0.3) %>%
  View()
  


target <- results %>%
  group_by(estimator_type, alpha) %>%
  summarize(estimate = first(target_estimate)) %>%
  mutate(estimate = estimate * 1000) 

with(target,{
  dbr <- target[estimator_type == 'dbr', ] 
  otc <- target[estimator_type == 'otc', ] 
  ipw <- target[estimator_type == 'ipw', ] 
  
  plot(alpha, estimate, type = 'n')
  lines(dbr$alpha, dbr$estimate)
  lines(otc$alpha, otc$estimate, col = 'blue')
  lines(ipw$alpha, ipw$estimate, col = 'red')
  legend(.55, 5.5, c('dbr', 'otc', 'ipw'), col = c('black', 'blue', 'red'),
         lty = 1)
  }
)

  


# test_glmer <- lme4::glmer(A ~ age + rivkm + (1|group), data = choleradt,
#                           family = binomial)
# test_geeglm <- geepack::geeglm(y_obs ~ A + fA + A*fA + age + rivkm, data = choleradt, 
#                                family = binomial, id = group)
# 
# test_glm <- glm(y_obs ~ A + fA + A*fA + age + rivkm, data = choleradt, 
#                                family = binomial)
