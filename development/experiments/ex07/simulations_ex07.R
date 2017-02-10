
library(inferference)
library(geex)
library(dr)
allocations <- c(0.1, 0.5, 0.9)

sims <- split(sims_500x_m300_n4_s198_ex07, sims_500x_m300_n4_s198_ex07$simID)
library(doMC)
registerDoMC(4)

estimates <- plyr::llply(sims, .parallel = TRUE, function(x){
  
  x <- x %>%
    group_by_(~group) %>% 
    mutate_(fA = ~ mean(A))
  
  ## Component models
  tmodel <- lme4::glmer(A ~ Z1 + (1|group) , data = x, 
                        family = binomial(link = 'logit') )
  omodel <- geepack::geeglm(Y ~ Z1 + A + fA + A:fA, data = x, 
                            id = group, family = gaussian(link = 'identity'))
  
  ## Grab component model parameter estimates
  theta_t <- unlist(lme4::getME(tmodel, c('beta', 'theta')))
  theta_o <- coef(omodel)
  # theta_o <- c(0.5, -0.788, -2.953, -0.098, -0.145, 0.351)
  
  temp <- list(eeFUN = dr_eefun, splitdt = split(x, x$group))
  
  ## Estimate parameters for each allocation
  
  lapply(allocations, function(allocation){
    temp <- append(temp, list(ee_args = list(alpha = allocation)))
    
    # Point estimates
    target <- lapply(temp$splitdt, function(gdt){
      est <- dr_estimators(gdt, t_model = tmodel, o_model = omodel, randomization = 1)
      est(c(theta_t, theta_o), alpha = allocation)
    }) %>% 
      list_matrix() %>% 
      apply(., 2, mean)
    
    # vcov estimate
    # mats <- compute_matrices(temp,
    #                          theta   = c(theta_t, theta_o, target),
    #                          numDeriv_options = list(method = 'simple'),
    #                          t_model = tmodel,
    #                          o_model = omodel)
    # Sigma <- compute_sigma(mats$A, mats$B) 
    
    # return
    # list(estimates = target, vcov = Sigma)
    target
  }) 
  
})


oracle <- expand.grid(alpha = allocations, 
                      a = c(0, 1, NA), 
                      truth = 2)

## Summarize results
results <- lapply(estimates, function(x){
  lapply(seq_along(x), function(k){
    data_frame(
      method   = rep(c('ipw', 'otc', 'dbr'), each = 3),
      a        = rep(c(0, 1, NA), 3),
      estimate = x[[k]],
      # std.err  = sqrt(diag(x[[k]][['vcov']]))[11:19],
      alpha    = allocations[k])
  }) %>% bind_rows()
})  %>% bind_rows() %>%
  left_join(oracle, by = c('a', 'alpha')) %>%
  mutate(bias      = estimate - truth
         # conf.low  = estimate - 1.96 * std.err,
         # conf.high = estimate + 1.96 * std.err,
         # covered   = conf.low < truth & conf.high > truth
         ) %>%
  group_by(method, a, alpha) %>%
  summarise(meanest = mean(estimate),
            bias    = mean(bias),
            medbias = median(bias),
            # ase     = mean(std.err),
            ese     = sd(estimate)
            # coverage = mean(covered)
            )



save(results, estimates, oracle, file = paste0('development/experiments/ex07/simresults_', Sys.Date(), '.rda'))


sims_500x_m300_n4_s198_ex07 %>%
  filter(simID == 1) %>%
  group_by(simID, group) %>%
  summarise(fA = mean(A)) %>%
  .$fA %>% hist()

 # IPW seems relatively unbiased, though still bias in 2nd decimnal and 