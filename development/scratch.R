## 
library(inferference)
library(geex)
library(dplyr)
##### 

vaccinesim <- vaccinesim %>%
  group_by(group) %>%
  mutate(fA = mean(A))

mm <- interference(Y | A ~ X1 + (1|group) | group,
                   allocations = c(.3, .6),
                   data = vaccinesim)
pmodel <- mm$models$propensity_model
omodel <- geepack::geeglm(Y ~ A + fA + X1, data = vaccinesim, 
                          id = group,
                          family = binomial(link = 'logit'))


theta_t <- unlist(lme4::getME(pmodel, c('beta', 'theta')))
theta_o <- coef(omodel)

mylist <-  append(list(eeFUN = lan_eefun, 
                       splitdt = split(vaccinesim, vaccinesim$group)), 
                  list(ee_args = list(alpha1 = .3, alpha2 = .6, a1 = 0, a2 = 0)))

ff <- lan_eefun(mylist$splitdt[[1]], t_model = pmodel, o_model = omodel)
ff(theta = c(theta_t, theta_o, .5, .5, .5, .5), alpha1 = .5, alpha2 = .5, a1 = 0, a2 = 1)

mats <- compute_matrices(mylist,  
                 theta   = c(theta_t, .5, .5),
                 numDeriv_options = list(method = 'simple'),
                 t_model = pmodel,
                 o_model = omodel)

sig <- compute_sigma(mats$A, mats$B)


ptm <- proc.time()
test <- geex::estimate_equations(
  eeFUN   = lan_eefun,
  data    = vaccinesim,
  units   = 'group',
  roots   = c(theta_t, theta_o, .5, .5, .5, .5),
  numDeriv_options = list(method = 'simple'),
  t_model = pmodel,
  o_model = omodel,
  ee_args = list(alpha1 = .3, alpha2 = .6, a1 = 0, a2 = 0)
)
proc.time() - ptm

test$parameters
sqrt(diag(test$vcov))
