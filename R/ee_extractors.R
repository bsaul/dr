#------------------------------------------------------------------------------#
#' glmer Estimating Equations
#'
#' Currently supports Logistic-Normal model with univariate random effect.
#'
#' @param x a \link[lme4]{merMod-class} object
#' @param grad_method method passed to \code{\link[numDeriv]{grad}}
#' @param grad_options list of options passed to \code{\link[numDeriv]{grad}}
#'  'method.opts' argument
#' @importFrom sandwich estfun
#' @return a (number of groups) x (number of parameters) matrix
#' @export
#------------------------------------------------------------------------------#

make_ee.merMod <- function(x, grad_method = 'simple', grad_options = NULL)
{
 psi.glmerMod(x = x, grad_method  = grad_method, grad_options = grad_options)
}


#------------------------------------------------------------------------------#
#' glmer Estimating Equations
#'
#' Used by \link{estfun.glmerMod} to create estimating equation per group
#'
#' @param x a \code{merMod} object
#' @param grad_method method passed to \code{\link[numDeriv]{grad}} or
#' \code{\link[numDeriv]{hessian}} when \code{deriv = TRUE}
#' @param grad_options method.args passed to \code{\link[numDeriv]{grad}} or
#'  \code{\link[numDeriv]{hessian}} when \code{deriv = TRUE}
#' @param deriv logical indicating whether to create derivative matrix of psi
#' @export
#------------------------------------------------------------------------------#

psi.glmerMod <- function(x, grad_method, grad_options = NULL)
{
  xmat   <- lme4::getME(x, 'X')
  resp   <- lme4::getME(x, 'y')
  frame  <- as.data.frame(cbind(resp, xmat))
  parms  <- unlist(lme4::getME(x, c('beta', 'theta')))
  clust  <- lme4::getME(x, 'flist')
  family <- x@resp$family
  linkinv <- family$linkinv
  objective_fun <- objFun_glmerMod(family$family)
  
  deriv_fun <- if(deriv == TRUE) numDeriv::hessian else numDeriv::grad
  
  ## Warnings ##
  if(length(lme4::getME(x, 'theta')) > 1){
    stop('estfun.glmer currently does not handle >1 random effect')
  }
  
  function(theta){
    objective_fun(parms = theta, response = Y, xmatrix = X, linkinv = linkinv)
  }
}

#------------------------------------------------------------------------------#
#' glmer Objective Fundtion
#'
#'@param family
#'@param ... additional arguments pass to objective function
#'@export
#------------------------------------------------------------------------------#

objFun_glmerMod <- function(family, ...){
  switch(family,
         binomial = objFun_glmerMod_binomial,
         stop('Objective function not defined'))
}

#------------------------------------------------------------------------------#
#' glmer Objective Function for Logistic-Normal Likelihood
#'
#' @param parms vector of parameters
#' @param response
#' @param xmatrix
#' @param linkinv inverse link function
#' @export
#------------------------------------------------------------------------------#

objFun_glmerMod_binomial <- function(parms, response, xmatrix, linkinv)
{
  # Logistic-Normal Model
  integrand <- function(b, response, xmatrix, parms){
    pr_response <- linkinv( drop(outer(xmatrix %*% parms[-length(parms)], b, '+') ) )
    hh <- dbinom(response, 1, prob = pr_response )
    hha <- apply(hh, 2, prod)
    hha * dnorm(b, mean = 0, sd = parms[length(parms)])
  }
  
  log(integrate(integrand, lower = -Inf, upper = Inf, parms = parms,
                response = response, xmatrix = xmatrix)$value)
  
}


# #------------------------------------------------------------------------------#
# #' glmer Objective Function for Likelihoods specified by glmerMod object
# #'
# #'@param vector of parameters
# #'@param response
# #'@param xmatrix
# #'@param inv_link_fun
# #'@export
# #------------------------------------------------------------------------------#
#
# objFun.glmerMod.generalized <- function(parms, response, xmatrix, inv_link_fun)
# {
#   # Logistic-Normal Model
#   integrand <- function(b, response, xmatrix, parms){
#     lp <- outer(xmatrix %*% parms[-length(parms)], b, '+')
#     h  <- apply(lp, 3, function(x) dbinom(response, 1, plogis(x) ) )
#     hh <- apply(h, 2, prod)
#     hh * dnorm(b, mean = 0, sd = parms[length(parms)])
#   }
#
#   inv_link_fun(integrate(integrand, lower = -Inf, upper = Inf, parms = parms,
#                 response = response, xmatrix = xmatrix)$value)
#
# }