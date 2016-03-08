#------------------------------------------------------------------------------#
#   Title: Lan et al Doubly Robust paper simulations
#  Author: B. Saul
#    Date: 2016-03-07
# Purpose: Simulate data for Lan et al Doubly Robust paper
#------------------------------------------------------------------------------#

library(dplyr)
library(magrittr)
library(simcausal) # need version >= 0.5.0 which support multivariate nodes
library(mvtnorm)

n_i <- 4
m <- 500
nsims <- 1000
totalobs <- n_i * m * nsims
seed <- 42


#### Functions ####

group_assign <- function(n, n_i)
{
  ni <- n_i[1]
  m <- n/ni
  out <- numeric(n)
  for(i in 0:(m-1)){
    start <- i*ni + 1
    end   <- i*ni + ni 
    out[start:end] <- i + 1
  }
  return(out)
}

fA <- function(n, A, groups)
{
  unlist(tapply(A, groups, function(x) {(sum(x) - x)/length(x)} ) )
}

# fAn <- function(n, A, groups)
# {
#   unlist(tapply(A, groups, function(x) {rep((sum(x) - x)/length(x), length(x))} ) )
# }

rnorm_group <- function(n, mean, sd, groups)
{
  m <- length(unique(groups))
  hold <- rnorm(m, mean = mean, sd = sd)
  out <- numeric(n)
  for(i in 1:n){
    out[which(groups == i)] <- hold[i]
  }
  return(out)
}

#### Creating DAG ####
D <- DAG.empty()
D <- D + 
  node('group', 
       distr = 'group_assign', 
       n_i = 4) + 
  node(c("Z1","Z2","Z3","Z4"), 
       distr = "rmvnorm", 
       mean = c(0,0,0,0)) + 
  node('X1', 
       distr = "rconst", 
       const = exp(Z1/2) ) + 
  node('X2', 
       distr = "rconst", 
       const = Z2/(1 + exp(Z1)) + 10) +
  node('X3',
       distr = "rconst",
       const = ((Z1*Z3)/25 + 0.6)^3 ) + 
  node('X4',
       distr = "rconst",
       const = (Z1 + Z4 + 20)^2 ) + 
  node('b',
     distr  = 'rnorm_group',
     mean   = 0,
     sd     = 1,
     groups = group) + 
  node('A',
       distr = 'rbern',
       prob  = plogis(-Z1 + 2*Z2 - 1.25*Z3 - 0.1*Z4 + b)) +
  node('fA',
       distr = 'fA',
       A = A, 
       groups = group) + 
  node('epsilon',
       distr  = 'rnorm_group',
       mean   = 0,
       sd     = 1,
       groups = group) + 
  node('Y',
       distr  = 'rconst',
       const  = 2 - Z1 - 2.7*Z2 + 3*Z3 - Z4  + 0.5*A + 6*fA + A*Z1 + 8*fA*Z2 )

D <- set.DAG(D)

#### Simulate Data ####

DRsims <- simobs(D, n = totalobs, rndseed = seed)
DRsims$simID <- sort(rep(1:nsims, n_i * m))