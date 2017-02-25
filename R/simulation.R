#------------------------------------------------------------------------------#
#' Generate data
#' @export
#------------------------------------------------------------------------------#
gen_data <- function(m, ni, gamma, theta, beta){
  n <- m * ni
  data_frame(
    m   = m,
    ni  = ni,
    group = rep(1:m, each = ni),
    Z1  = rnorm(n, sd = 1),
    Z1_abs = abs(Z1),
    Z2  = rbinom(n, size = 1, prob = .5),
    b   = rep(rnorm(m, sd = theta), each = ni),
    p   = as.numeric(plogis(cbind(1, Z1_abs, Z2, Z1_abs*Z2) %*% gamma + b)),
    A   = rbinom(n, size = 1, prob = p),
    fA  = rep(tapply(A, group, mean), each = ni),
    Y   = rnorm(n, mean = (cbind(1, A, fA, Z1_abs, Z2, Z1_abs*Z2) %*% beta), sd = 1)
  )
}

#------------------------------------------------------------------------------#
#' Generate data
#' @export
#------------------------------------------------------------------------------#

gen_sim <- function(nsims, ...){
  dots <- list(...)
  args <- dots[pmatch(names(dots), names(formals((dr::gen_data))))]
  x    <- as.call(append(list(dr::gen_data), args))
  replicate(nsims, eval(x), simplify = FALSE)
}
