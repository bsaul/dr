library(sandwich)
library(inferference)
library(lme4)


example1 <- interference(
           formula = y | A | B ~ X1 + X2 | group,
           allocations = c(.3, .45,  .6),
           model_method = 'glm',
           model_options = list(family = binomial()),
           data = vaccinesim,
           randomization = 2/3,
           method = 'simple')

s <- example1$scores

fit <- glm(B ~ X1 + X2, data = vaccinesim, family = binomial())
u <- apply(estfun(fit), 2, function(x) tapply(x, vaccinesim$group, sum)) 

head(s)
head(u)

# Close
all.equal(s, u)

s - u


##### Glmer estfunning

example2 <- interference(
  formula = y | A | B ~ X1 + X2 + (1 | group) | group,
  allocations = c(.3, .45,  .6),
  data = vaccinesim,
  randomization = 2/3,
  method = 'simple')

s2 <- example2$scores

test <- vaccinesim
test$g <- rep(1:200, 15)
gfit <- glmer(B ~ X1 + X2 + (1|group ), data = test, family = binomial)

dd <- update(gfit,devFunOnly=TRUE)
dd(unlist(getME(gfit,c("theta","beta"))))

numDeriv::grad(dd, x = unlist(getME(gfit,c("theta","beta"))) ) * -1/2


u <- apply(estfun(fit), 2, function(x) tapply(x, vaccinesim$group, sum)) 


# glmer.estfun <- function(x, ...)
# {
#   xmat  <- model.matrix(x)
#   resp  <- getME(x, 'y')
#   parms <- unlist(getME(x, c('beta', 'theta')))
#   nparms <- length(parms)
#   
#   integrand <- function(b, parms){
#     lc <- outer(xmat %*% parms[-nparms], b, '+')
#     h  <- apply(lc, 3, function(x) dbinom(resp, 1, plogis(x) ) )
#     hh <- apply(h, 2, prod)
#     hh * dnorm(b, mean = 0, sd = parms[nparms])
#   }
#   
#   objective.fun <- function(parms){
#     log(integrate(integrand, lower = -Inf, upper = Inf, parms = parms)$value)
#   }
#   
#   numDeriv::grad(objective.fun, x = parms)
# 
# }

test <- microbenchmark::microbenchmark(
  estfun.glmer(gfit, grad.method = 'simple'),
  score_matrix(logit_integrand, X = getME(gfit, 'X'), A = getME(gfit, 'y'),
               G = getME(gfit, 'flist')[[1]], fixed.effects = getME(gfit, 'beta'),
               random.effects = getME(gfit, 'theta'), method = 'simple') )


