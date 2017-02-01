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
                  list(ee_args = list(alpha = c(.5))))

ff <- lan_eefun(mylist$splitdt[[1]], t_model = pmodel, o_model = omodel)
ff(theta = c(theta_t, theta_o, rep(.5, 9)), alpha = c(.5))
# 
# mats <- compute_matrices(mylist,  
#                  theta   = c(theta_t, .5, .5),
#                  numDeriv_options = list(method = 'simple'),
#                  t_model = pmodel,
#                  o_model = omodel)
# 
# sig <- compute_sigma(mats$A, mats$B)
wf <- weight_estimator(A = mylist$splitdt[[1]]$A, 
                 X = get_design_matrix(get_fixed_formula(pmodel), mylist$splitdt[[1]]))
wf(theta_t, c(.5, .6, .7))

test <- geex::eeroot(
  geex_list   = mylist,
  start   = c(theta_t, theta_o, rep(.5, 9)),
  t_model = pmodel,
  o_model = omodel)

test_matrices <- geex::compute_matrices(
  geex_list   = mylist,
  theta   = test$root,
  numDeriv_options = list(method = 'simple'),
  t_model = pmodel,
  o_model = omodel)

ptm <- proc.time()
test <- geex::estimate_equations(
  eeFUN   = lan_eefun,
  data    = vaccinesim,
  units   = 'group',
  roots   = c(theta_t, theta_o, rep(.5, 9)),
  numDeriv_options = list(method = 'simple'),
  t_model = pmodel,
  o_model = omodel,
  ee_args = list(alpha = .5)
)
proc.time() - ptm

test$parameters
sqrt(diag(test$vcov))
