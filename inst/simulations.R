
library(doMC)
registerDoMC(4)
allocations <- c(.3, .45, .6)

## Make the simulated datasets and compute oracle
library(interferenceSim)
library(plyr)
simdt <- interferenceSim::sim_interference(n = 5000, 
                                           N = 300, 
                                           nsims= 250, 
                                           alphas = allocations)
detach('package:interferenceSim')
detach('package:plyr')

#### Estimate results ####

estimates <- plyr::llply(simdt$sims, .parallel = TRUE, function(x){
   x <- x %>%
     group_by_(~group) %>% 
     mutate_(fA = ~ mean(A))
  
  ## Component models
  tmodel <- lme4::glmer(B ~ X1 + X2 + (1|group), data = x, 
                        family = binomial(link = 'logit') )
  omodel <- geepack::geeglm(y ~ A + fA + X1 + X2 + A:fA, data = x, 
                            id = group,
                            family = binomial(link = 'logit'))
  
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
      est <- dr_estimators(gdt, t_model = tmodel, o_model = omodel)
      est(c(theta_t, theta_o), alpha = allocation)
    }) %>% 
      list_matrix() %>% 
      apply(., 2, mean)
    
    # vcov estimate
    mats <- compute_matrices(temp,
                             theta   = c(theta_t, theta_o, target),
                             numDeriv_options = list(method = 'simple'),
                             t_model = tmodel,
                             o_model = omodel)
    Sigma <- compute_sigma(mats$A, mats$B) 
    
    # return
    list(estimates = target, vcov = Sigma)
  }) 
})


## Organize the oracle
oracle <- simdt$truth %>%
  tidyr::gather(a, truth, -alpha) %>%
  mutate(a = ifelse(a == 'a0', 0,
                    ifelse(a == 'a1', 1, NA)))

## Summarize results
results <- lapply(estimates, function(x){
  lapply(seq_along(x), function(k){
    data_frame(
      method   = rep(c('ipw', 'otc', 'dbr'), each = 3),
      a        = rep(c(0, 1, NA), 3),
      estimate = x[[k]][['estimates']],
      std.err  = sqrt(diag(x[[k]][['vcov']]))[11:19],
      alpha    = allocations[k])
  }) %>% bind_rows()
})  %>% bind_rows() %>%
  left_join(oracle, by = c('a', 'alpha')) %>%
  mutate(bias      = estimate - truth,
         conf.low  = estimate - 1.96 * std.err,
         conf.high = estimate + 1.96 * std.err,
         covered   = conf.low < truth & conf.high > truth) %>%
  group_by(method, a, alpha) %>%
  summarise(meanest = mean(estimate),
            bias    = mean(bias),
            medbias = median(bias),
            ase     = mean(std.err),
            ese     = sd(estimate),
            coverage = mean(covered))

 save(results, estimates, oracle, file = paste0('development/simresults_', Sys.Date(), '.rda'))
