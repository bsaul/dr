


x <- tru_tru[[1]]$ipw$weightd[ , , 1] 

y <- tapply(fitted(tru_tru[[1]]$outcome), tru_tru[[1]]$outcome$geese$id, sum)

z1 <- tru_tru[[1]]$ipw$Upart$marginal_outcomes$groups[ , , 1]

z2 <- apply(x, 2, function(col) col * y)

z1 - z2
#------------------------------------------------------------------------------#
#### OLD ####
#------------------------------------------------------------------------------#
library(inferference)
library(doMC)
registerDoMC(4)

#------------------------------------------------------------------------------#
#### Get sims working for IPW first ####
#------------------------------------------------------------------------------#

alphas <- c(0.1, 0.5, 0.9)

sims_400x_m500_n4 <- sims_400x_m500_n4 %>%
  mutate(
#     pA1 = plogis(-Z1 + 2*Z2 - 1.25*Z3 - 0.1*Z4 + b),
#     pA2 = plogis(0 - Z1 + .05 * Z2 - .05 * Z3 - .05 * Z4 + b),
    Y1 =  2 - 1.5 * Z1 - 2.7*Z2 + 3*Z3 - Z4  + 0.5*A1 + 6*fA1 + A1*Z1 + 8*fA1*Z2 + e,
    Y2 =  2 - 1.5 * Z1 - 2.7*Z2 + 3*Z3 - Z4  + 0.5*A2 + 6*fA2 + A2*Z1 + 8*fA2*Z2 + e
  )

# coef_y_inter.true=c(2,-1.5,-2.7,3,-1,0.5,6,1,8)

# hist(sims_1000x_m500_n4$pA1, freq = TRUE)
# hist(sims_1000x_m500_n4$pA2, freq = TRUE)

# ipw_results1 <- sims_500x_m500_n4 %>%
#   # filter(simID <= 500) %>%
#   {plyr::dlply(., plyr::.(simID), .progress = 'text', .parallel = TRUE,
#                function(sim){
#                  this_sim <- interference(formula =  Y1 | A1 ~ Z1 + Z2 + Z3 + Z4 | group,
#                                           allocations = alphas,
#                                           data = sim,
#                                           model_method = 'oracle',
#                                           model_options = list(fixed.effects = c(0.5,-1,0.5,-0.25,-0.1),
#                                                                random.effects = 1),
#                                           causal_estimation_options = list(variance_estimation ='naive'),
#                                           method = 'simple')
#                })}

ipw_results2 <- sims_400x_m500_n4 %>%
  {plyr::dlply(., plyr::.(simID), .progress = 'text', .parallel = TRUE,
               function(sim){
                 this_sim <- interference(formula =  Y2 | A2 ~ Z1 + Z2 + Z3 + Z4 + (1|group) | group,
                                          allocations = alphas,
                                          data = sim,
#                                           model_method = 'oracle',
#                                           model_options = list(fixed.effects = c(0.5,-1,0.5,-0.25,-0.1),
#                                                                random.effects = 1),
                                          causal_estimation_options = list(variance_estimation ='naive'),
                                          method = 'simple')
               })}

truth <- data.frame(alpha1 = alphas, truth = 2+0.5*1+6*(alphas*3/4))

# ipw_estimates1 <- lapply(ipw_results1, function(x){
#   e <- x$estimates
#   e %>% filter(effect_type == 'outcome', trt1 == 1) %>%
#     left_join(truth, by = 'alpha1')
# }) %>% bind_rows()
# 
# ipw_estimates1 %>%
#   group_by(alpha1) %>%
#   summarise(bias = mean(estimate - truth))

ipw_estimates2 <- lapply(ipw_results2, function(x){
  e <- x$estimates
  e %>% filter(effect_type == 'outcome', trt1 == 1) %>%
    left_join(truth, by = 'alpha1')
}) %>% bind_rows()

ipw_estimates2 %>%
  group_by(alpha1) %>%
  summarise(bias = mean(estimate - truth, na.rm = T))


# alpha1 mean(estimate)
# (dbl)          (dbl)
# 1    0.1       2.559289
# 2    0.5       4.702558
# 3    0.9       6.546338

