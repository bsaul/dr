#------------------------------------------------------------------------------#
#   Title: Analysis of cholera data with IPW, OTC, and DR estimators
#  Author: B. Saul
#    Date: 2016-05-20
# Purpose:
#     Log:
# 20170623 - running all but WLS estimator; will do seperately as I'm getting
#            an error from the WLS estimator when estimating the variance
#------------------------------------------------------------------------------#

library(dplyr)
library(dr)
library(doMC)
library(geex)
registerDoMC(4)
# library(geex)
# alphas <- .45
 alphas <- lapply(seq(.3, .6, by = .02), function(x) sort(c(.4, x)))
 alphas <- list(c(.4))
# alphas <- c(.4, .6)

load( pipe( 'ssh saulb@diamond.bios.unc.edu "cat /home/groups/projects/mhudgens/emch/data/R_data/emch_analysis_data.Rdata"' ))

choleradt <- analysis_c %>%
  group_by(group) %>%
  mutate(fA = mean(A)) %>%
  filter(n() < 1074) %>%
  ungroup() 
# %>%
#   filter(group < 100)

analysis_model_args <- list(
  t_model = 
    list(method = lme4::glmer,
         formula = B ~ age + rivkm + age2 + rivkm2 + (1|group),
         options = list(family = binomial(link = 'logit'))),
  o_model =
    list(method  = geepack::geeglm,
         formula = y_obs ~ A + fA + age + rivkm + age2 + rivkm2,
         options = list(
           family  = binomial(link = 'logit'),
           id      = quote(group))),
  wls_model_0 = 
    list(method = stats::glm,
         formula = as.integer(y_obs) ~ fA + age + rivkm + age2 + rivkm2,
         options = list(family = quasibinomial(link = 'logit'))),
  wls_model_1 = 
    list(method = stats::glm,
         formula = as.integer(y_obs) ~ fA + age + rivkm + age2 + rivkm2,
         options = list(family = quasibinomial(link = 'logit'))),
  pcov_model_0 = list(
    method  = stats::glm,
    formula = y_obs ~ fA + age + rivkm + age2 + rivkm2 + ipw0,
    options = list(family = binomial(link = 'logit'))
  ),
  pcov_model_1 = list(
    method  = stats::glm,
    formula = y_obs ~ fA + age + rivkm + age2 + rivkm2 + ipw1,
    options = list(family = binomial(link = 'logit'))
  )
)


models0 <- estimate_cholera_parms_step0(
  data        = choleradt,
  model_args  = analysis_model_args
)

# results <- estimate_cholera_parms_step2(
#   data        = choleradt,
#   allocations = list(alphas),
#   models      = models0,
#   model_args  = analysis_model_args,
#   compute_se  = FALSE
# )

results <- lapply(seq_along(alphas), function(i){
  attempt <- try(estimate_cholera_parms_step2(
    data        = choleradt,
    allocations = alphas[i],
    models      = models0,
    model_args  = analysis_model_args,
    randomization  = 2/3,
    compute_se  = TRUE
  ))
  if(is(attempt, 'try-error')){
    NULL
  } else{
    attempt
  }
})

results2 <- results

# Method = "Richardson" & 2 target quantities
sigma <- results2[[1]][[1]]$wls_dbr$vcov
estimate <- results2[[1]][[1]]$wls_dbr$estimate * 1000
C <- matrix(c(1, -1), byrow = TRUE, ncol = 2 )
Cvcov   <- cbind(matrix(0, nrow = nrow(C), ncol = nrow(sigma)- 2), C)
std_err <- as.numeric(sqrt(diag(Cvcov %*% sigma1 %*% t(Cvcov))) * 1000)

de <- estimate1[1] - estimate1[2]
de + 1.96 * std_err
de - 1.96 * std_err

# Method = "Simple" & 4 target quantities
sigma <- results[[6]][[1]]$wls_dbr$vcov
estimate <- results[[6]][[1]]$wls_dbr$estimate * 1000
C <- matrix(c(1, 0, -1, 0), byrow = TRUE, ncol = 4 )
Cvcov   <- cbind(matrix(0, nrow = nrow(C), ncol = nrow(sigma)- 4), C)
std_err <- as.numeric(sqrt(diag(Cvcov %*% sigma %*% t(Cvcov))) * 1000)

de <- estimate[1] - estimate[3]
de + 1.96 * std_err
de - 1.96 * std_err


save(results, file = paste0('inst/cholera_analysis/cholera_results_dr_wls_', Sys.Date(),'.rda'))
  
