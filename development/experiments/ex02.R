#------------------------------------------------------------------------------#
#   Title: Experiment 02: speeding up the outcome estimator
#  Author: B. Saul
#    Date: 2016-04-30
# Purpose: develop faster outcome estimators
#------------------------------------------------------------------------------#

library(lme4)
library(dplyr)
library(geepack)
library(magrittr)
library(sandwich)
library(sandwichShop)
library(lineprof)

DRsims <- sims_500x_m500_n4 %>%  filter(simID < 2)

tru_outcome <- Y ~ Z1 + Z2 + Z3 + Z4 + A + fA + A*Z1 + fA*Z2
tru_treatment  <- A ~ -1 + Z1 + Z2 + Z3 + Z4 + (1|group) 

m_Y <- geeglm(tru_outcome, data = DRsims, id = ID)
m_A <- glmer(tru_treatment, data = DRsims, family = binomial)
theta_o <- coef(m_Y)
theta_t <- unlist(getME(m_A, c('beta', 'theta')))

theta_a <- c(theta_t, theta_o)
rhs_o <- formula(m_Y)[-2]
rhs_t <- formula(m_A)[-2]

split_dt <- split(DRsims, DRsims$group)
g1_dt <- split_dt[[1]]

#------------------------------------------------------------------------------#
#### Experiment 02-00 ####
# examine performance of dbr estimator as it is in DR_functions_4.R - with the 
# exception of using otc_estimator_03
#------------------------------------------------------------------------------#

source('development/experiments/ex02_functions.R')

# estimator 
system.time(ex02_00_1())
# derivative of estimator - single group
system.time(ex02_00_2())
# derivative of estimator over all groups - not an option - taking over 5 minutes!
# system.time(ex02_00_3()
# derivative of estimator - single group; method = 'simple'
# shaves off alot of time!
system.time(ex02_00_4())

#------------------------------------------------------------------------------#
#### Experiment 02-01 ####
# examine performance of dbr estimator as it is in DR_functions_4.R - with the 
# exception of using otc_estimator_03
#------------------------------------------------------------------------------#

source('development/experiments/ex02_functions.R')

# U21 for single group 
system.time(ex02_01_1())
# U21 for all groups
system.time(ex02_01_2())

system.time(ex02_01_3())

system.time(ex02_01_4())

# compare values - close in most cases. different in terms of A, Z1:A and Z2:fA; why?
ex02_01_1()
ex02_00_2()

# Directly computing derivatives takes 1/4 of the time compared to calling numDeriv
microbenchmark::microbenchmark(ex02_01_1(), ex02_00_2(), times = 50L)


# Directly computing derivatives takes 1/4 of the time compared to calling numDeriv;
# how about using method = 'simple'; still quite a bit faster
microbenchmark::microbenchmark(ex02_01_1(), ex02_00_4(), times = 50L)
#------------------------------------------------------------------------------#
#### Experiment 02_02 ####
# create a function that directly computes the closed form of the U21 matrix;
# but passing in the Weights and derivatives precomputed
#------------------------------------------------------------------------------#

#U21 for single group
microbenchmark::microbenchmark(ex02_01_1(), ex02_03_1(), times = 50L)
# U21 for all groups 
system.time(ex02_03_2())


#------------------------------------------------------------------------------#
#### Experiment 02_03 ####
# create a function that directly computes the closed form of the U21 matrix;
# but passing in the Weights and derivatives precomputed; using method = 'simple'
#------------------------------------------------------------------------------#

system.time(ex02_03_1())
#U21 for single group
microbenchmark::microbenchmark(ex02_03_1(), ex02_02_1(), ex02_01_1(), ex02_00_4(), times = 50L)
#option3 is the current choice: 

make_dbr_U21_estimator_03 <- function(Y, A, X_outcome, X_treatment, rhs_formula_outcome){
  w <- weight_estimator(A = A, X = X_treatment)
  ws <- w(theta_t)
  wd <- numDeriv::grad(w, x = theta_t, method = 'simple')
  pi_t <- pi_term(A = A)
  dr_term1 <- make_dr_term1(X_outcome)(theta_o)
  n <- nrow(X_outcome)
  function(alpha, a = NULL){
    pi_t1 <- pi_t(alpha) / {if(!is.null(a)) dbinom(a, 1, alpha) else 1}
    # U21 corresponding to treatment parameters
    Ia <- if(is.null(a)) 1 else (A == a) * 1
    Ybar <- mean(Ia * (Y - dr_term1) ) 
    U21_t <- Ybar * wd * pi_t1
    
    # U21 corresponding to  outcome parameters
    xmat <- model.matrix(rhs_formula_outcome, data = X_outcome)
    tt <- 1 - (Ia * pi_t1 * ws)
    
    U21_o <- apply(xmat, 2, function(col) {
      sum(col * tt)/n
    })
    
    c(U21_t, U21_o)
  }
}