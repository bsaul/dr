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
# source('R/functions_v2.R', echo = T, max.deparse.length = 5000)
source('R/functions_v4.R', echo = T, max.deparse.length = 5000)
# alphas <- c(.3, .45, .6)
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
  mutate(fA = mean(A)) 
# %>%
  # filter(n() < 1074)

results_a0 <- estimation(treatment_formula = A ~ age + rivkm + (1|group), 
                      outcome_formula = y_obs ~ A + fA + A*fA + age + rivkm,
                      outcome_method = glm,
                      outcome_method_opts = list(family = binomial),
                      this_data = choleradt,
                      allocations = alphas,
                      target_a = 0)

results_a1 <- estimation(treatment_formula = A ~ age + rivkm + (1|group), 
                      outcome_formula = y_obs ~ A + fA + A*fA + age + rivkm,
                      outcome_method = glm,
                      outcome_method_opts = list(family = binomial),
                      this_data = choleradt,
                      allocations = alphas,
                      target_a = 1) 


save(results_a0, results_a1, file = 'cholera_analysis/cholera_results.rda')

results %>%
  filter(estimator_type == 'ipw') %>%
  arrange(-n_i)

results %>%
  filter(alpha == 0.3)

results %>%
  filter(n_i > 1074, alpha == 0.3) %>%
  View()
  


target_a0 <- results_a0 %>%
  group_by(estimator_type, alpha) %>%
  summarize(estimate = first(target_estimate)) %>%
  mutate(estimate = estimate * 1000) 

target_a1 <- results_a1 %>%
  group_by(estimator_type, alpha) %>%
  summarize(estimate = first(target_estimate)) %>%
  mutate(estimate = estimate * 1000) 

target <- left_join(target_a0, target_a1, by = c('estimator_type', 'alpha'))

target <- target %>%
  group_by(estimator_type) %>%
  mutate(y0      = estimate.x,
         y1      = estimate.y,
         de_diff = estimate.x - estimate.y,
         de_rat  = 1 - (estimate.x - estimate.y),
         y0.4    = estimate.x[alpha == .4],
         ie_diff = y0.4 - estimate.x,
         ie_rat  = 1 - (y0.4/estimate.x),
         te_diff = y0.4 - estimate.y,
         te_rat  = 1 - (y0.4/estimate.y))


plot_effect <- function(var, ylim, xlab, ylab){
  with(target, {
    dbr <- target[estimator_type == 'dbr', ] 
    otc <- target[estimator_type == 'otc', ] 
    ipw <- target[estimator_type == 'ipw', ] 
    
    plot(dbr$alpha, dbr[[var]], 
         type = 'n', bty = 'l', 
         ylim = ylim,
         xlab = xlab,
         ylab = ylab,
         mgp  = c(2, .5, 0))
    lines(dbr$alpha, dbr[[var]])
    lines(otc$alpha, otc[[var]], lty = 2)
    lines(ipw$alpha, ipw[[var]], lty = 4)
  })
}

pdf(file = 'cholera_analysis/cholera_ce_diff_plots_V001.pdf', width = 7, height = 3)
par(mfrow=c(1,3))
plot_effect('de_diff', c(-2, 8), expression(alpha), expression(widehat("DE")*"("*alpha*")"))
legend(.45, 6.5, c('dbr', 'otc', 'ipw'),  lty = c(1, 2, 4), bty = 'n')
plot_effect('ie_diff', c(-2, 8), expression(alpha), expression(widehat("IE")*"("*0.4*","*alpha*"'"*")"))
plot_effect('te_diff', c(-2, 8), expression(alpha), expression(widehat("TE")*'('*0.4*","*alpha*"'"*")"))
dev.off()


pdf(file = 'cholera_analysis/cholera_ce_diff_plots_V001.pdf', width = 7, height = 3)
par(mfrow=c(1,3))
plot_effect('de_diff', c(-2, 8), expression(alpha), expression(widehat("DE")*"("*alpha*")"))
legend(.45, 6.5, c('dbr', 'otc', 'ipw'),  lty = c(1, 2, 4), bty = 'n')
plot_effect('ie_diff', c(-2, 8), expression(alpha), expression(widehat("IE")*"("*0.4*","*alpha*"'"*")"))
plot_effect('te_diff', c(-2, 8), expression(alpha), expression(widehat("TE")*'('*0.4*","*alpha*"'"*")"))
dev.off()

pdf(file = 'cholera_analysis/cholera_y_plots_V001.pdf', width = 7, height = 3)
par(mfrow=c(1,3))
plot_effect('y0', c(-2, 8), expression(alpha), expression(widehat("Y")*"(0,"*alpha*")"))
legend(.45, 8, c('dbr', 'otc', 'ipw'),  lty = c(1, 2, 4), bty = 'n')
plot_effect('y1', c(-2, 8), expression(alpha), expression(widehat("Y")*"(1,"*alpha*")"))
dev.off()
  





# test_glmer <- lme4::glmer(A ~ age + rivkm + (1|group), data = choleradt,
#                           family = binomial)
# test_geeglm <- geepack::geeglm(y_obs ~ A + fA + A*fA + age + rivkm, data = choleradt, 
#                                family = binomial, id = group)
# 
# test_glm <- glm(y_obs ~ A + fA + A*fA + age + rivkm, data = choleradt, 
#                                family = binomial)
