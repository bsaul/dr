
## Generate scenarios ##
scenarios <- data_frame(
  sid = 1:12,
  gamma = rep(list(c(0.1, 0.2, 0, .2),
                   c(2.1, 0.4, -4.2, .2),
                   c( -1, 0.5, -1, .2), 
                   c(0  , 0.75, 0, .2)),
              times = 3),
  theta = rep(list(.3, .3, .3, 1.5), 3),
  beta  = rep(list(c(2, 0, 0, -1.5, 2, -3), 
                   c(2, 2, 0, -1.5, 2, -3),
                   c(2, 2, 1, -1.5, 2, -3)), 
              each = 4),
  m     = 100,
  ni    = 30
)

## Oracle ## 
oracle <- expand.grid(
  sid    = 1:12,
  alpha = c(0.1, 0.5, 0.9), 
  a     = c(0, 1, NA)) %>%
  mutate(truth = 2 + # Expected value of intercept 
           (-1.5 * sqrt(2/pi)) +   # sqrt(2/pi) is mean of half normal distr when sigma = 1
           2 * 0.5 + # expected value of Z2 term 
           (-3 * sqrt(2/pi) * 0.5) +# expected value of interaction term
           (I(sid > 4) * a * 2) + # Direct effect
           (I(sid > 8) * alpha * 1) # Indirect effect
         ) 

## Component model arguments ##
tcor_ocor_margs <- list(
  t_model = 
    list(method = lme4::glmer,
         formula = A ~ Z1_abs + Z2 + Z1_abs*Z2 + (1|group),
         options = list(family = binomial(link = 'logit'))),
  o_model =
    list(method  = geepack::geeglm,
         formula = Y ~ A + fA + Z1_abs + Z2 + Z1_abs*Z2,
         options = list(
           family  = gaussian(link = 'identity'),
           id      = quote(group))),
  wls_model_0 = list(
    method  = glm,
    formula = Y ~ fA + Z1_abs + Z2 + Z1_abs*Z2,
    options = list(
      family  = gaussian(link = 'identity'))
  ),
  wls_model_1 = list(
    method  = glm,
    formula = Y ~ fA + Z1_abs + Z2 + Z1_abs*Z2,
    options = list(
      family  = gaussian(link = 'identity'))
  )
)

tmis_ocor_margs <- tcor_ocor_margs
tmis_ocor_margs$t_model$formula <- A ~ Z1 + Z2 + (1|group)

tcor_omis_margs <- tcor_ocor_margs
tcor_omis_margs$o_model$formula      <- Y ~ A + fA + Z1 + Z2
tcor_omis_margs$wls_model_0$formula  <- Y ~ fA + Z1 + Z2
tcor_omis_margs$wls_model_1$formula  <- Y ~ fA + Z1 + Z2

tmis_omis_margs <- tcor_ocor_margs
tmis_omis_margs$t_model$formula      <- A ~ Z1 + Z2 + (1|group)
tmis_omis_margs$o_model$formula      <- Y ~ A + fA + Z1 + Z2
tmis_omis_margs$wls_model_0$formula  <- Y ~ fA + Z1 + Z2
tmis_omis_margs$wls_model_1$formula  <- Y ~ fA + Z1 + Z2

margs <- list(
  tcor_ocor = tcor_ocor_margs,
  tmis_ocor = tmis_ocor_margs,
  tcor_omis = tcor_omis_margs,
  tmis_omis = tmis_omis_margs
)
