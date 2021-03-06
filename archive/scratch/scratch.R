## 
library(inferference)
library(geex)
library(dplyr)
library(dr)
##### 

temp1 <- gen_sim(
  gamma = c(0.1, 0.2, 0, .2),
  theta = c(.3),
  beta  = c(2, 10, 0, -1.5, 2, -3),
  m     = 50,
  ni    = 20,
  nsims = 1
)
temp1 <- temp1[[1]]


test_margs <- list(
  t_model = 
    list(method = lme4::glmer,
         formula = A ~ Z1_abs + Z2 + Z1_abs*Z2 + (1|group),
         options = list(family = binomial(link = 'logit'))),
  o_model =
    list(method  = geepack::geeglm,
         formula = Y ~ A + fA + Z1_abs + Z2 + Z1_abs*Z2,
         options = list(
           family  = gaussian(link = 'identity'),
           id      = quote(group)))
)

mods <- make_models(model_args = test_margs, data = temp1)

theta_t1 <- unlist(lme4::getME(mods$t_model, c('beta', 'theta')))
theta_o1 <- coef(mods$o_model)


mylist <-  append(list(eeFUN = generic_eefun,
                       splitdt = split(temp1, temp1$group)),
                  list(ee_args = list(alpha = c(.5, .6))))


lapply(mylist$splitdt, function(x) {
  f3 <-generic_eefun(
    data = x,
    models = mods,
    randomization =  1,
    'dbr',
    hajek = FALSE)
  
  f3(c(theta_t1, theta_o1, 0, 0, 0, 0, 0, 0), alpha= c(.5, .6, .7))
})

temp1 %>% filter(group == 20) %>% View()

