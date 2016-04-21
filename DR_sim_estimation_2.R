
#### Fit nuisance models ####

model_args = list(
  model_treatment = list(
    formula = tru_propenA,
    method  = lme4::glmer,
    options = list(family = binomial)
  ) ,
  model_outcome = list(
    formula = tru_outcome,
    method  = geepack::geeglm,
    options = list(family = gaussian, id = quote(group))
  ) )

models <- make_models(model_args = model_args, data = DRsims)

theta_treatment <- unlist(lme4::getME(models$model_treatment, c('beta', 'theta')))
theta_outcome <- coef(models$model_outcome)
theta_observed <- c(theta_treatment, theta_outcome)

#### Create list of group level data  ####

groups     <- DRsims$group
xmat_treatment <- as.data.frame(model.matrix(models$model_treatment))
xmat_outcome <- as.data.frame(model.matrix(models$model_outcome))
treatment  <- model.response(model.frame(models$model_treatment))
outcome    <- DRsims$Y

frame <- list(outcome = outcome, 
              treatment = treatment,
              X_treatment = xmat_treatment,
              X_outcome   = xmat_outcome) %>%
  lapply(., function(x) split(x, groups, drop = FALSE))

split_frame <- Map(list, outcome      = frame[['outcome']],
                         treatment    = frame[['treatment']],
                         X_treatment  = frame[['X_treatment']],
                         X_outcome    = frame[['X_outcome']])

#### Create list of group level estimators ####
estimators  <- make_group_estimators(split_frame,
                              estimator = 'make_dr_estimator',
                              formula_outcome = formula(models$model_outcome))

estimators <- lapply(estimators, function(f) memoise::memoise(f))

#### Estimate target parameters #### 
target <- lapply(estimators, function(f) {
  f(theta = theta_observed,
    alpha = .5,
    a = 1)
}) %>%
  list_matrix() %>%
  mean()
target

#### Estimate variance parameters ####

psi_observed <- lapply(estimators, function(f) {
  f2 <- psi(f)
  f2(c(theta_observed, target), alpha = .5, a= 1)
}) %>%
  list_matrix()

ee_0 <- estfun_stacker(models)
ee <- cbind(ee_0, psi_observed)
V <- crossprod(ee)/nrow(ee)

psi_prime_target <- lapply(estimators, function(f) {
  f2 <- psi(f)
  numDeriv::grad(f2, x = c(theta_observed, target), alpha = .5, a = 1)
}) %>%
  list_matrix() %>%
  apply(2, mean) %>%
  matrix(nrow = 1)

bread_model <- make_bread(models, grad_method = 'Richardson')

U <- rbind_fill_zero(list(bread_model, -psi_prime_target))
U <- as.matrix(U)
sigma <- (solve(U) %*% V %*% t(solve(U)))/nrow(ee)

ase <- sigma[nrow(sigma), ncol(sigma)]
sqrt(ase)

target - 1.96 * sqrt(ase)
