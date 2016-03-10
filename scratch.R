
DRresults <- DRsims %>% 
  #   filter(simID <= 1) %>%
{plyr::dlply(., plyr::.(simID), function(sim){
  this_sim <- dr_ipw(formula_outcome = Y ~ Z1 + Z2 + Z3 + Z4 + A + fA + A*Z1 + fA*Z2,
                     formula_interference = Y | A ~ -1 + Z1 + Z2 + Z3 + Z4 + (1|group) | group,
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

#   
# DRresults %>%
#   group_by(type, alpha1) %>%
#   summarise(median = median(estimate),
#             mean   = mean(estimate))
# 
# DRresults %>%
#   filter(alpha1 == .1)
# 
# temp <- dr_ipw(formula_outcome = Y ~ X1 + A + fA + A*Z1 + fA*Z2,
#                formula_interference = Y | A ~ -1 + Z1 + Z2 + Z3 + Z4 + (1|group) | group,
#                method_outcome  = geepack::geeglm,
#                method_opts_outcome = list(id = quote(group), family = gaussian),
#                method_opts_interference = list(allocations = alphas,
#                                                method = 'simple', 
#                                                runSilent = T),
#                dr_term2_function = dr_term2,
#                data = DRsims %>% filter(simID == 450))
# summary(temp$ipw$model$propensity_model)$coefficients
# temp$ipw$estimates %>% head()

# dr_ipw_estimate(temp)
# rm(temp)
# testdt <- DRsims %>% filter(simID == 2)
# mtest <- geepack::geeglm(Y ~ Z1 + Z2 + Z3 + Z4 + A + fA + A*Z1 + fA*Z2, data = testdt,
#                 id = group)
# head(dr_term2(.1, mtest, testdt)) 
# 
# - head(mtest$y)
# mtest$residuals
