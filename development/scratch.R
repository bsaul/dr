## 
library(inferference)
library(geex)
library(dplyr)
##### 

temp <- simdt$sims[[1]] %>%
  group_by(group) %>%
  mutate(fA = mean(A))

# mm <- interference(y | A | B ~ X1 + X2 + (1|group) | group,
#                    allocations = c(.3, .6),
#                    data = simdt$sims[[1]])
tmodel <-  tmodel <- lme4::glmer(B ~ X1 + X2 + (1|group), data = temp, 
                                  family = binomial(link = 'logit') )
omodel <- geepack::geeglm(y ~ A + fA + X1 + A:fA, data = temp, 
                          id = group,
                          family = binomial(link = 'logit'))


theta_t <- unlist(lme4::getME(pmodel, c('beta', 'theta')))
theta_o <- coef(omodel)

mylist <-  append(list(eeFUN = dr_eefun,
                       splitdt = split(temp, temp$group)),
                  list(ee_args = list(alpha = c(.5))))


f2 <- dr_estimators(mylist$splitdt[[1]], t_model = pmodel, o_model = omodel)
f2(c(theta_t, theta_o), alpha = .5)

xx <- lapply(mylist$splitdt, function(x){
  f <- dr_estimators(x, t_model = pmodel, o_model = omodel)
  f(c(theta_t, theta_o), alpha = .5)
})
targ <- xx %>% list_matrix() %>% apply(., 2, mean)
mats <- compute_matrices(mylist,
                 theta   = c(theta_t, theta_o, targ),
                 numDeriv_options = list(method = 'simple'),
                 t_model = tmodel,
                 o_model = omodel)
head(mats)
compute_sigma(mats$A, mats$B) %>% diag() %>% sqrt()

summary(pmodel)
# sig <- compute_sigma(mats$A, mats$B)
# wf <- weight_estimator(A = mylist$splitdt[[1]]$A, 
#                  X = get_design_matrix(get_fixed_formula(pmodel), mylist$splitdt[[1]]))
# wf(theta_t, c(.5, .6, .7))
# 
# test <- geex::eeroot(
#   geex_list   = mylist,
#   start   = c(theta_t, theta_o, rep(.5, 9)),
#   t_model = pmodel,
#   o_model = omodel)
# 
# test_matrices <- geex::compute_matrices(
#   geex_list   = mylist,
#   theta   = test$root,
#   numDeriv_options = list(method = 'simple'),
#   t_model = pmodel,
#   o_model = omodel)

ptm <- proc.time()
test <- geex::estimate_equations(
  eeFUN   = dr_eefun,
  data    = vaccinesim,
  units   = 'group',
  roots   = c(theta_t, theta_o, rep(.5, 9)),
  numDeriv_options = list(method = 'simple'),
  t_model = pmodel,
  o_model = omodel,
  ee_args = list(alpha = .45)
)
proc.time() - ptm

test$parameters
sqrt(diag(test$vcov))
