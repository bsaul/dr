## 
library(inferference)
library(geex)
library(dplyr)
##### 

temp1 <- simdt$sims[[1]] %>%
  group_by_(~group) %>%
  mutate_(fA = ~mean(A))

tmodel1 <-  tmodel <- lme4::glmer(B ~ X1 + X2 + (1|group), data = temp1, 
                                  family = binomial(link = 'logit') )
omodel1 <- geepack::geeglm(y ~ A + fA + X1 + X2 +  A:fA, data = temp1, 
                          id = group,
                          family = binomial(link = 'logit'))


theta_t1 <- unlist(lme4::getME(tmodel1, c('beta', 'theta')))
theta_o1 <- coef(omodel1)

mylist <-  append(list(eeFUN = dr_eefun,
                       splitdt = split(temp1, temp1$group)),
                  list(ee_args = list(alpha = c(.3))))

fitted(omodel1)


xx <- lapply(mylist$splitdt, function(x){
  f <- dr_estimators(x, t_model = tmodel1, o_model = omodel1)
  f(c(theta_t1, theta_o1), alpha = .3)
})
targ <- xx %>% list_matrix() %>% apply(., 2, mean)

f1 <- dr_eefun(data = mylist$splitdt[[1]], t_model = tmodel1, o_model = omodel1)
f1(c(theta_t1, theta_o1, targ), alpha = .3) 

mats <- compute_matrices(mylist,
                 theta   = c(theta_t1, theta_o1, targ),
                 numDeriv_options = list(method = 'simple'),
                 t_model = tmodel1,
                 o_model = omodel1)
head(mats)
compute_sigma(mats$A, mats$B) 


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
