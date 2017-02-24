## 
library(inferference)
library(geex)
library(dplyr)
library(dr)
##### 

temp1 <- simulations[[12]][[1]] 

test_margs <- list(
  t_model = 
    list(method = lme4::glmer,
         formula = A ~ Z1_abs + Z2 + (1|group),
         options = list(family = binomial(link = 'logit'))),
  o_model =
    list(method  = geepack::geeglm,
         formula = Y ~ A + fA + Z1_abs + Z2,
         options = list(
           family  = gaussian(link = 'identity'),
           id      = quote(group)))
)

mods <- make_models(model_args = test_margs, data = temp1)
extract_model_info(mods, data = temp1 %>% filter(group == 1), estimator_type = 'dbr')

theta_t1 <- unlist(lme4::getME(mods$t_model, c('beta', 'theta')))
theta_o1 <- coef(mods$o_model)




f2 <- ipw_estimator(data= temp1 %>% filter(group == 2),
                    models = mods, randomization = 1)
f2(c(theta_t1), alpha = .5)


f2 <- otc_estimator(data= temp1 %>% filter(group == 1),
                    models = mods, randomization = 1)
f2(c(theta_o1), alpha = .5)


mylist <-  append(list(eeFUN = generic_eefun,
                       splitdt = split(temp1, temp1$group)),
                  list(ee_args = list(alpha = c(.3))))

f3 <-generic_eefun(
  data = temp1 %>% filter(group == 5),
  models = mods,
  randomization =  1,
  'dbr')

f3(c(theta_t1, theta_o1, 0, 0), alpha= .5 )

xx <- lapply(mylist$splitdt, function(x){
  f <- otc_estimator(x, models = mods, randomization = 1)
  f(c(theta_o1), alpha = .5)
})
target <- xx %>% list_matrix() %>% apply(., 2, mean)
target

mats <- geex::compute_matrices(
  mylist,
   theta   = c(theta_o1, target),
   numDeriv_options = list(method = 'simple'),
   models  = mods,
   randomization = 1,
   estimator_type = 'otc')
library(Matrix)

# Sigma <- 
  chol2inv(chol(mats$A))
%*% mats$B %*% t(qr.solve(mats$A))
sqrt(diag(Sigma)[-(1:(length(theta_t1) + length(theta_o1)))])
matrix(c(1,2,3,4), nrow = 2) %>% qr.solve()

