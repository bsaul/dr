
library(doMC)
registerDoMC(4)
allocations <- c(.3, .45, .6)

## Make the simulated datasets and compute oracle
library(interferenceSim)
library(plyr)
simdt <- interferenceSim::sim_interference(n = 3000, 
                                           N = 250, 
                                           nsims= 250, 
                                           alphas = allocations)
detach('package:interferenceSim')
detach('package:plyr')

#### Estimate results ####
estimates <- plyr::llply(simdt$sims[1], .parallel = FALSE, function(x){
   x <- x %>%
     group_by_(~group) %>% 
     mutate_(fA = ~ mean(A))
  
   x_splt <- split(x, x$group)
   
  ## Component models
  tmodel <- lme4::glmer(B ~ X1 + X2 + (1|group), data = x, 
                        family = binomial(link = 'logit') )
  omodel <- geepack::geeglm(y ~ A + fA + X1 + X2 + A:fA, data = x, 
                            id = group,
                            family = binomial(link = 'logit'))
  
  ## Grab component model parameter estimates
  theta_t <- unlist(lme4::getME(tmodel, c('beta', 'theta')))
  theta_o <- coef(omodel)
  
  ## Estimate parameters for each allocation
  lapply(allocations, function(allocation){
    temp_list <- append(list(eeFUN = dr_eefun, splitdt = x_splt),
                        list(ee_args = list(alpha = allocation)))

    # Point estimates
    target <- lapply(temp_list$splitdt, function(x){
      est <- dr_estimators(x, t_model = tmodel, o_model = omodel)
      est(c(theta_t, theta_o), alpha = allocation)
    }) %>% 
      list_matrix() %>% 
      apply(., 2, mean)
    
    # vcov estimate
    mats <- compute_matrices(temp_list,
                             theta   = c(theta_t, theta_o, target),
                             numDeriv_options = list(method = 'simple'),
                             t_model = tmodel,
                             o_model = omodel)
    Sigma <- compute_sigma(mats$A, mats$B) 
    
    # return
    list(estimates = target, vcov = Sigma)
  }) 
})

estimates

# ## Organize the oracle
# oracle <- simdt$truth %>%
#   tidyr::gather(a, truth, -alpha) %>%
#   mutate(a = ifelse(a == 'a0', 0, 
#                     ifelse(a == 'a1', 1, NA)))
# 
# ## Summarize results
# results <- estimates %>%
#   left_join(oracle, by = c('a', 'alpha')) %>%
#   mutate(bias      = estimate - truth,
#          conf.low  = estimate - 1.96 * std.err,
#          conf.high = estimate + 1.96 * std.err,
#          covered   = conf.low < truth & conf.high > truth) %>%
#   group_by(method, a, alpha) %>%
#   summarise(bias = mean(bias),
#             ase  = mean(std.err),
#             ese  = sd(estimate),
#             coverage = mean(covered))
# 
# 
# # parameter_index <- (p + 1):length(z$parameters)
# # data_frame(
# #   method   = rep(c('ipw', 'otc', 'dbr'), each = 3),
# #   a        = rep(c(0, 1, NA), 3),
# #   estimate = z$parameters[parameter_index],
# #   std.err  = sqrt(diag(z$vcov)[parameter_index]),
# #   alpha    = allocation
# 
# 
# save(results, file = paste0('development/simresults_', Sys.Date(), '.rda'))
