#------------------------------------------------------------------------------#
#   Title: Experiment 05: comparing different integrands for IPW
#  Author: B. Saul
#    Date: 2016-06-13
# Purpose: 
#------------------------------------------------------------------------------#

temp <- DRsims %>%
  filter(simID == 1, group == 1)

A <- temp$A
X <- cbind(1, temp$Z1, temp$Z2, temp$Z3, temp$Z4)


integrand <- function(b, response, xmatrix, theta){
  lp <- outer(xmatrix %*% theta[-length(theta)], b, '+')
  h  <- apply(lp, 3, function(x) dbinom(response, 1, plogis(x) ) )
  # print(h)
  hh <- apply(h, 2, prod)
  # print(hh)
  hh * dnorm(b, mean = 0, sd = theta[length(theta)])
}

integrand2 <- function(b, response, xmatrix, theta, alpha){
  # print(b)
  lp <- outer(xmatrix %*% theta[-length(theta)], b, '+')
  plp <- plogis(lp)
  # print(plp)
  h  <- apply(plp, 3, function(x) response * log(x/alpha) + (1 - response) * log((1 - x)/(1 - alpha)))
  # print(h)
  hh <- apply(h, 2, function(x) exp(sum(x)))
  # print(hh)
  hh * dnorm(b, mean = 0, sd = theta[length(theta)])
}

integrand3 <- function(b, response, xmatrix, theta, alpha){
  # print(b)
  lp <- outer(xmatrix %*% theta[-length(theta)], b, '+')
  plp <- plogis(lp)
  h <- apply(plp, 3, function(x) ifelse(response == 1, log(x/alpha), log((1 - x)/(1 - alpha))))
  # print(h)
  hh <- apply(h, 2, function(x) exp(sum(x)))
  # print(hh)
  hh * dnorm(b, mean = 0, sd = theta[length(theta)])
}

# Note that integrand is not divided by alpha so the values will be slightly different

integrand(b        = c(0, 20,  100, 200),
          response = A,
          xmatrix  = X,
          theta    = c(0.5, -1, 0.5, -0.25, -0.1, 1))


integrand2(b        = c(0, 20,  100, 200),
           response = A,
           xmatrix  = X,
           theta    = c(0.5, -1, 0.5, -0.25, -0.1, 1),
           alpha    = .5)

integrand3(b        = c(0, 20,  100, 200),
           response = A,
           xmatrix  = X,
           theta    = c(0.5, -1, 0.5, -0.25, -0.1, 1),
           alpha    = .5)

