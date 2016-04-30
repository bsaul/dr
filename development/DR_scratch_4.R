system.time(temp <- estimation(tru_propen, tru_outcome, DRsims,
                   allocations = c(.1, .5, .65, .9), target_a = 1))

temp  %>%
  mutate( psi_prime_val = list(evaluate_df_function(numDeriv::grad, func = psi_func, 
                                                    x = c(theta, mean(estimate)), method = 'simple', alpha = alpha, a = 1)),
                psi_val  = estimate - mean(estimate) )




# this works 
temp2 <- temp %>%
  group_by(estimator_type) %>%
  rowwise() %>%
  do(estimate = as.numeric(.$estimator(theta = .$theta, alpha = .5, a = 1)))


temp2 <- temp %>%
  group_by(estimator_type) %>%
  rowwise() %>%
  do(data_frame(
    estimate = as.numeric(.$estimator(theta = .$theta, alpha = .5, a = 1)),
    psi_fun  = list(psi(.$estimator))
    ) ) 

evaluate_df_function2 <- function(x, ...){
#   print(x)
#   print(list(...))
   x(...)  
}


temp2 <- temp %>%
  group_by(estimator_type) %>%
  rowwise() %>%
  mutate(estimate = evaluate_df_function2(estimator, theta = theta, alpha = alpha, a = a))
 


 mutate(estimate = .$estimator[[1]](theta = .$theta[[1]], alpha = .5, a = 1))


  do(estimate = .$estimator(theta = .$theta, alpha = .5, a = 1))

create_df_function <- function(x, f, ...){
  lapply(x, function(y) f(y, ...) )
  # lapply(seq_along(x), function(i) x[[i]](...) ) %>% unlist() 
}

# lapply(temp2$ipw, function(x) print(x))
# str(temp2[4, 1]$ipw)

temp2 <- temp %>%
  mutate(ipw_estimate = evaluate_df_function(ipw_estimator, theta = c(0, 0, 0, 0, 1), alpha = .5, a = 1),
         ipw_psi = create_df_function(ipw_estimator, psi),
         ipw_target = mean(ipw_estimate),
         # ipw_theta = lapply(ipw_target, function(x) c(0, 0, 0, 0, 1, x) ) ,
         ipw_psi_value = evaluate_df_function(ipw_psi, theta = c(0, 0, 0, 0, 1, mean(ipw_target)), alpha = .5, a = 1) ) 

temp2


expand.grid(type = c('ipw', 'otc', 'dbr'), allocation = alphas)


melt


temp3 <- data_frame(temp2$ipw)
temp3$ipw[[1]]


lapply(temp[1:2], function(x) names(x))


theta <- list(theta2, theta1, c(theta2, theta1))
lapply(temp[1:2], function(funcs) {
   lapply(seq_along(funcs), function(f){
     # print(funcs[[f]])
     funcs[[f]](theta[[f]], alpha = .5, a = 1)
  })
})
split_data <- split(DRsims, DRsims$group)

f <- function(x) x

data_frame(func = list(f, f))
