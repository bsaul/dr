#------------------------------------------------------------------------------#
#   Title: Doubly Robust Functions
#  Author: B. Saul
#    Date: 2016-03-07
# Purpose: functions for IPW and DR estimators
#------------------------------------------------------------------------------#

dr_ipw <- function(formula_outcome,
                   formula_interference,
                   method_outcome,
                   method_opts_outcome,
                   method_opts_interference,
                   dr_term2_function,
                   data)
{
  
  ## Outcome Model
  model_args_outcome <- append(list(formula = formula_outcome), method_opts_outcome)
  model_args_outcome <- append(model_args_outcome, list(data = data))
  model_outcome <- do.call(method_outcome, args = model_args_outcome)
  
  ## IPW Model
  model_args_ipw <- append(list(formula = formula_interference), method_opts_interference)
  model_args_ipw <- append(list(data = data), model_args_ipw)
  model_ipw <- do.call(inferference::interference, args = model_args_ipw)

  ## M-estimation
#   ee_outcome <- estfun(model_outcome)
#   ee_ipw <- model_ipw$scores
#   # Stack EEs
#   ee <- cbind(ee_outcome, ee_ipw)
  
  ## DR
  #mu <- predict(model_outcome, type = 'response')
  mu <- model.matrix(model_outcome) %*% coef(model_outcome)
  w <- apply(model_ipw$weights, 2, function(row) {matrix(rep(row, each = 4))})
  
  w_mu <- apply(w, 2, function(col) col * mu)

  ## Term 2
  term2 <- lapply(model_ipw$summary$allocations, function(alpha) {
    dr_term2_function(alpha, outcome_model = model_outcome, data = data) 
  })
  
  out <- list(outcome = model_outcome, 
              ipw = model_ipw, 
              mu = mu, 
              w = w, 
              w_mu = w_mu, 
              term2 = term2)
  return(out)
}


dr_term2 <- function(alpha, outcome_model, data)
{
  data %>%
    # this needs to be done by group. plyr splits the data frame by group
    # then applies a function to each piece
    {plyr::dlply(., plyr::.(group), function(x){
      x %>%
        # Generate all possible sum(a_i) for each subject
        merge(expand.grid(sum_a = 0:(n_distinct(.$ID) - 1)), all = T)  %>%
        group_by_(~ID) %>%
        # Compute pi and p_ij for each sum(a_i)
        mutate_(A =~ 1,
                fA = ~ sum_a/n(),
                # fAn = ~ sum_a/n(),
                pi = ~ choose(n() - 1, sum_a) * alpha^sum_a * (1 - alpha)^(n() - 1 - sum_a)) %>%
        ungroup() %>%
        # Compute mu_ij for each a_i per subject
        # For some reason predict() is throwing error with geeglm model
        #mutate_(mu_ij = ~predict(outcome_model, newdata = ., type = 'response')) %>%
        mutate_(mu =~ as.numeric(model.matrix(outcome_model$formula, .) %*% coef(outcome_model) )  ) %>%
        # Sum by individual to compute term2
        group_by_(~ID) %>%
        summarize_(term2 = ~sum(mu * pi)) 
    })} %>%
    # Put the split data back together
    bind_rows()
}

dr_ipw_estimate <- function(obj)
{

  ## IPW
  ipw_estimates <- obj$ipw$estimates %>%
    filter(effect_type == 'outcome', trt1 == 1) %>%
    select(alpha1, trt1, estimate, std.error, conf.low, conf.high) %>%
    mutate(method = 'ipw')
  
  ## DR
  dr_results <- numeric(ncol(obj$w))
  outcome_results <- numeric(ncol(obj$w))
  for(k in 1:ncol(obj$w)){
    term2 <- obj$term2[[k]]
    alphak <- as.numeric(dimnames(obj$w)[[2]])[k]
    
    outcome_results[k] <- obj$outcome$data %>%
      left_join(term2, by = 'ID') %>%
      mutate_(Yhat_ij = ~ term2) %>%
      group_by_(~group) %>%
      summarise_(Yhat_i = ~mean(Yhat_ij)) %>%
      summarise_(Yhat = ~mean(Yhat_i) ) %>%
      as.numeric()

    dr_results[k] <- obj$outcome$data %>%
      mutate_(mu =~ as.numeric(obj$mu) ) %>%
      mutate_(term1 =~ as.numeric( ((A == 1) * (Y - mu)) * obj$w[ , k]/alphak) ) %>%
      left_join(term2, by = 'ID') %>%
      mutate_(Yhat_ij = ~term1 + term2) %>%
      group_by_(~group) %>%
      summarise_(Yhat_i = ~mean(Yhat_ij)) %>%
      summarise_(Yhat = ~mean(Yhat_i) ) %>%
      as.numeric()
  }
  dr_estimates <- data.frame(alpha1 = as.numeric(dimnames(obj$w)[[2]]),
                             trt1 = 1, estimate = dr_results, method = 'dr')
  
  outcome_estimates <- data.frame(alpha1 = as.numeric(dimnames(obj$w)[[2]]),
                                  trt1 = 1, estimate = outcome_results, 
                                  method = 'outcome')
  
  bind_rows(ipw_estimates, dr_estimates, outcome_estimates)
}

dr_results <- function(sims, formula_outcome, formula_interference)
{
  sims %>% 
    {plyr::dlply(., plyr::.(simID), .progress = 'text', .parallel = TRUE,
                 function(sim){
      this_sim <- dr_ipw(formula_outcome = formula_outcome,
                         formula_interference = formula_interference,
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
}