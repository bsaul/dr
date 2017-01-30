## 
library(inferference)
library(geex)
##### 

vaccinesim

mm <- interference(Y | A ~ X1 + (1|group) | group,
                   allocations = c(.3, .6),
                   data = vaccinesim)
pmodel <- mm$models$propensity_model
omodel <- glm(Y ~ X1, data = vaccinesim)
theta_t <- unlist(lme4::getME(pmodel, c('beta', 'theta')))


mylist <-  append(list(eeFUN = ipw_eefun, splitdt = split(vaccinesim, vaccinesim$group)), 
                  list(ee_args = list(alpha1 = .3, alpha2 = .6, a1 = 0, a2 = 0)))

str(mylist)
ff <- ipw_eefun(mylist$splitdt[[1]], t_model = pmodel, o_model = omodel)
ff(theta = c(theta_t, .5, .5), alpha1 = .5, alpha2 = .5, a1 = 0, a2 = 1)


mats <- compute_matrices(mylist,  
                 theta   = c(theta_t, .5, .5),
                 numDeriv_options = list(method = 'simple'),
                 t_model = pmodel,
                 o_model = omodel)

sig <- compute_sigma(mats$A, mats$B)

sqrt(diag(sig))

summary(pmodel)

ptm <- proc.time()
test <- geex::estimate_equations(
  eeFUN   = ipw_eefun,
  data    = vaccinesim,
  units   = 'group',
  roots   = c(theta_t, .5, .5),
  numDeriv_options = list(method = 'simple'),
  t_model = pmodel,
  o_model = omodel,
  ee_args = list(alpha1 = .3, alpha2 = .6, a1 = 0, a2 = 0)
)
proc.time() - ptm

test$parameters
sqrt(diag(test$vcov))
