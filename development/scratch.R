## 
library(inferference)
library(geex)
library(dplyr)
library(dr)
##### 

temp1 <- simulations[[5]][[1]] 

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

theta_t1 <- unlist(lme4::getME(mods$t_model, c('beta', 'theta')))
theta_o1 <- coef(mods$o_model)




f2 <- ipw_estimator(data= temp1 %>% filter(group == 2),
                    models = mods, randomization = 1, hajek = TRUE)
f2(c(theta_t1), alpha = .5)


f2 <- otc_estimator(data= temp1 %>% filter(group == 1),
                    models = mods, randomization = 1)
f2(c(theta_o1), alpha = .5)


mylist <-  append(list(eeFUN = generic_eefun,
                       splitdt = split(temp1, temp1$group)),
                  list(ee_args = list(alpha = c(.5))))


lapply(mylist$splitdt, function(x) {
  f3 <-generic_eefun(
    data = x,
    models = mods,
    randomization =  1,
    'ipw',
    hajek = TRUE)
  
  f3(c(theta_t1, 0, 0), alpha= .5 )
})

temp1  %>% filter(group == 133)

