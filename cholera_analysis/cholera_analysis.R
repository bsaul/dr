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

alphas <- c(.3, .45, .6)



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
  mutate(fA = mean(A))

results <- estimation(treatment_formula = A ~ age + rivkm + (1|group), 
                      outcome_formula = y_obs ~ A + fA + A*fA + age + rivkm,
                      outcome_method = glm,
                      outcome_method_opts = list(family = binomial),
                      this_data = choleradt,
                      allocations = alphas,
                      target_a = 0) 

results %>%
  mutate(estimate = estimate * 1000)
results$estimate * 1000

otc_1 <- results %>% filter(estimator_type == 'otc', alpha == .3)


system.time(
otc_1 %>% 
  filter(row_number() < 10) %>%
  rowwise() %>%
  mutate(estimate = evaluate_df_function(estimator, theta = theta, a = a, alpha = alpha))
)
6300/60

with(results, plot(x = alpha, y = estimate))


test_glmer <- lme4::glmer(A ~ age + rivkm + (1|group), data = choleradt,
                          family = binomial)
test_geeglm <- geepack::geeglm(y_obs ~ A + fA + A*fA + age + rivkm, data = choleradt, 
                               family = binomial, id = group)

test_glm <- glm(y_obs ~ A + fA + A*fA + age + rivkm, data = choleradt, 
                               family = binomial)
