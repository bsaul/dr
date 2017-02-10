### Comparing weight estimators from dr and inferference

testdt <- sims_100x_m300_n4_s198_ex06 %>%
  filter(simID == 1)
test_tmodel <- lme4::glmer(A ~ Z1 + (1|group) , data = testdt, 
                      family = binomial(link = 'logit') )
test_coef <- unlist(lme4::getME(test_tmodel, c('beta', 'theta')))
testdt_g1 <- testdt %>%
  filter(group == 1)

X <- geex::get_design_matrix(~ Z1, data = testdt_g1)
A <- geex::get_response(A ~ Z1 + (1|group), data = testdt_g1)


logit_integrand(b = c(0, 1), X = X, A = A, parameters = test_coef, allocation = .5)
integrand(b = c(0, 1), response = A, xmatrix = X, theta = test_coef, alpha = .5)
