#------------------------------------------------------------------------------#
#   Title: Analysis of cholera data with IPW, OTC, and DR estimators
#  Author: B. Saul
#    Date: 2016-05-20
# Purpose: 
#------------------------------------------------------------------------------#

library(dplyr)
library(dr)
library(doMC)
registerDoMC(4)
# library(geex)
# alphas <- .45
alphas <- lapply(seq(.3, .6, by = .02), function(x) sort(c(.4, x)))

load( pipe( 'ssh saulb@diamond.bios.unc.edu "cat /home/groups/projects/mhudgens/emch/data/R_data/emch_analysis_data.Rdata"' ))

choleradt <- analysis_c %>%
  group_by(group) %>%
  mutate(fA = mean(A)) %>%
  filter(n() < 1074)

analysis_model_args <- list(
  t_model = 
    list(method = lme4::glmer,
         formula = B ~ age + rivkm + (1|group),
         options = list(family = binomial(link = 'logit'))),
  o_model =
    list(method  = geepack::geeglm,
         # formula = y_obs ~ A + fA + A*fA + age + rivkm,
         formula = y_obs ~ A + fA + age + rivkm,
         options = list(
           family  = binomial(link = 'logit'),
           id      = quote(group)))
)


results <- estimate_cholera_parms(
  data        = choleradt,
  allocations = alphas,
  model_args  = analysis_model_args
)



save(results, file = 'inst/cholera_analysis/cholera_results.rda')

  