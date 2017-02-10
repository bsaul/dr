#------------------------------------------------------------------------------#
#   Title: Lan et al Doubly Robust paper simulations
#  Author: B. Saul
#    Date: 2016-03-07
# Purpose: Simulate data for Lan et al Doubly Robust paper
#------------------------------------------------------------------------------#

library(dplyr)
library(magrittr)
library(simcausal) 
n_i <- 10
m <- 300
nsims <- 250
totalobs <- n_i * m * nsims
seed <- 198
experimentID <- 'X010'
#### Functions ####

group_assign <- function(n, n_i)
{
  ni <- n_i[1]
  m <- ceiling(n/ni)
  out <- rep(1:m, each = ni)[1:n]

  return(out)
}

fA <- function(n, A, groups)
{
  unlist(tapply(A, groups, function(x) {(sum(x) - x)/length(x)} ) )
}

rnorm_group <- function(n, mean, sd, groups)
{
  m <- length(unique(groups))
  g <- as.numeric(table(groups))
  hold <- rnorm(m, mean = mean, sd = sd)
  
  out <- rep(hold, times = g)
  return(out)
}

#### Creating DAG ####
D <- DAG.empty()
D <- D + 
  node('group', 
       distr = 'group_assign', 
       n_i = .(n_i) ) + 
  node('Z1',
       distr = 'rnorm',
       mean = 0, sd = 1) + 
  node('Z2',
       distr = 'rnorm',
       mean = 0, sd = 1) + 
  node('Z3',
       distr = 'rnorm',
       mean = 0, sd = 1) + 
  node('Z4',
       distr = 'rnorm',
       mean = 0, sd = 1) + 
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
#   But anyways, it is very close to what in my table 1, so I think my coefficient is 
# > coef_pi.true
# [1]  0.50 -1.00  0.50 -0.25 -0.10  
    node('A',
         distr = 'rbern',
         prob  = plogis(0.5 - Z1 + 0.5*Z2 - 0.25*Z3 - 0.1*Z4 + b)) + 
  node('fA',
       distr = 'fA',
       A = A, 
       groups = group) + 
  node('e',
       distr  = 'rnorm',
       mean   = 0,
       sd     = 1)  + 
  node('Y',
       distr  = 'rconst',
       const  = 2 - 1.5 * Z1 - 2.7 * Z2 + 3 * Z3 - Z4  + 0.5 * A + 6 * fA + A*Z1 + 8 * fA * Z2 + e)

D <- set.DAG(D)

#### Simulate Data ####

DRsims <- simobs(D, n = totalobs, rndseed = seed) 
DRsims$simID <- sort(rep(1:nsims, n_i * m))

name <- paste0('sims_', nsims, 'x_', 'm', m, '_n', n_i, '_s', seed, '_', experimentID)
assign(name, DRsims)
save(list = name, file = paste0('data/', name, '.rda'))

rm(n_i, name, m, totalobs, seed, fA, group_assign, rnorm_group, nsims, D)
