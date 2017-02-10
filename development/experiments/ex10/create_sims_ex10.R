#------------------------------------------------------------------------------#
#   Title: Lan et al Doubly Robust paper simulations
#  Author: B. Saul
#    Date: 2017-02-08
# Purpose: 
#------------------------------------------------------------------------------#

library(dplyr)
library(magrittr)
library(simcausal) 
n_i <- 20
m <- 300
nsims <- 200
totalobs <- n_i * m * nsims
seed <- 198
experimentID <- 'ex10'

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
  # unlist(tapply(A, groups, function(x) {(sum(x) - x)/length(x)} ) )
  unlist(tapply(A, groups, function(x) {rep(sum(x)/length(x), length(x)) }  ) )
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
       mean = 0, sd = .2) + 
  node('b',
     distr  = 'rnorm_group',
     mean   = 0,
     sd     = 1,
     groups = group) + 
  node('A',
         distr = 'rbern',
         prob  = plogis(0.5 - Z1 + b)) + 
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
       const  = 2 - 1.5*Z1 + 2*A + 0*fA + 0*A*fA + e)

D <- set.DAG(D)

#### Simulate Data ####

DRsims <- simobs(D, n = totalobs, rndseed = seed) 
DRsims$simID <- sort(rep(1:nsims, n_i * m))

name <- paste0('sims_', nsims, 'x_', 'm', m, '_n', n_i, '_s', seed, '_', experimentID)
assign(name, DRsims)
save(list = name, file = paste0('data/', name, '.rda'))

rm(n_i, name, m, totalobs, seed, fA, group_assign, rnorm_group, nsims, D)
