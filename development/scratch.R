## 
library(inferference)
library(geex)
library(dplyr)
library(dr)
##### 

temp1 <- DRsims %>%
  filter(simID == 1)

tmodel1 <-   lme4::glmer(A ~ Z1 + Z2 + Z3 + Z4 + (1|group), data = temp1, 
                                  family = binomial(link = 'logit') )
omodel1 <- geepack::geeglm(Y ~ A + fA + Z1 + Z2 + Z3 + Z4 + A:fA, data = temp1, 
                          id = group)

theta_t1 <- unlist(lme4::getME(tmodel1, c('beta', 'theta')))
theta_o1 <- coef(omodel1)

mylist <-  append(list(eeFUN = dr_eefun,
                       splitdt = split(temp1, temp1$group)),
                  list(ee_args = list(alpha = c(.3))))




f2 <- dr_estimators(data= mylist$splitdt[[3]], t_model = tmodel1, o_model = omodel1, randomization = 1)
f2(c(theta_t1, theta_o1), alpha = .5)


xx <- lapply(mylist$splitdt, function(x){
  f <- dr_estimators(x, t_model = tmodel1, o_model = omodel1, randomization = 1)
  f(c(theta_t1, theta_o1), alpha = .5)
})
targ <- xx %>% list_matrix() %>% apply(., 2, mean)
targ





