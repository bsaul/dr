## 
library(inferference)
library(geex)
library(dplyr)
library(dr)
##### 

temp1 <- simdt$sims[[1]] %>%
  group_by_(~group) %>%
  mutate_(k  =~ sum(A) - A,
          fA =~ (A + k)/n(),
          fA2 =~ mean(A))
hist(temp1$fA)
tmodel1 <-   lme4::glmer(B ~ X1 + X2 + (1|group), data = temp1, 
                                  family = binomial(link = 'logit') )
omodel1 <- geepack::geeglm(y ~ A + fA + X1 + X2 +  A:fA, data = temp1, 
                          id = group,
                          family = binomial(link = 'logit'))
summary(omodel1)

theta_t1 <- unlist(lme4::getME(tmodel1, c('beta', 'theta')))
theta_o1 <- coef(omodel1)

mylist <-  append(list(eeFUN = dr_eefun,
                       splitdt = split(temp1, temp1$group)),
                  list(ee_args = list(alpha = c(.3))))

fitted(omodel1)

test <- lapply()


test <- lapply(simdt$sims[1:50], function(x){
 x <- x %>% group_by_(~group) %>%
    mutate_(fA  =~ mean(A))
  geepack::geeglm(y ~ A + fA + X1 + X2 +  A:fA, data = x, 
                  id = group,
                  family = binomial(link = 'logit')) %>%
    fitted() %>%
    cbind(plogis(model.matrix(y ~ A + fA + X1 + X2 +  A:fA, data = x) %*% c(.5, -.788, -2.953, -0.098, -0.145, 0.351))) %>%
    as.data.frame() 
})
test3 <- expand_outcome_frame(as.data.frame(model.matrix(y ~ A + fA + X1 + X2 +  A:fA, data = temp1 %>% filter(group == 1))), ~ A + fA + X1 + X2 +  A:fA)

head(test3)
test2 <- test %>%
  bind_rows() %>% 
  mutate(V3 = V2 - V1)


xx <- lapply(mylist$splitdt, function(x){
  f <- dr_estimators(x, t_model = tmodel1, o_model = omodel1)
  f(c(theta_t1, theta_o1), alpha = .3)
})
targ <- xx %>% list_matrix() %>% apply(., 2, mean)
targ
f1 <- dr_eefun(data = mylist$splitdt[[1]], t_model = tmodel1, o_model = omodel1)
f1(c(theta_t1, theta_o1, targ), alpha = .3) 

f2 <- dr_estimators(data= mylist$splitdt[[2]], t_model = tmodel1, o_model = omodel1)

mats <- compute_matrices(mylist,
                 theta   = c(theta_t1, theta_o1, targ),
                 numDeriv_options = list(method = 'simple'),
                 t_model = tmodel1,
                 o_model = omodel1)
head(mats)
compute_sigma(mats$A, mats$B) 



f2 <- dr_estimators(data= mylist$splitdt[[1]], t_model = tmodel1, o_model = omodel1)
f2(c(theta_t1, theta_o1), alpha = 0.3)

f3 <- make_otc_estimator(as.data.frame(model.matrix(y ~ A + fA + X1 + X2 +  A:fA, data = temp1 %>% filter(group == 1))), 
                         ~ A + fA + X1 + X2 +  A:fA, inv_link = plogis)
f3(theta_o1, a = 0, alpha = 0.3)

temp <- sims_250x_m300_n4_s198_X010 %>%
  filter(simID ==1) %>%
  interference(formula = Y | A ~ Z1 + Z2 + Z3 + Z4 + (1|group) | group, 
               data = ., 
               allocations = c(.1, .5, .9))

  
  temp$estimates %>%
    filter(effect == 'outcome', trt1 == 0)
  
  estimates[[1]][[1]]$estimates[1:3]
  oracle
  
  
 temp2 <- sims_250x_m300_n10_s198_X010 %>%
    group_by(simID, group) %>%
    summarise(fA = mean(A)) %>%
    ungroup() %>%
    select(fA)
  
  temp2$fA%>% hist()
  