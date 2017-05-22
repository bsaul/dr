# Functions
integrand_picov <- function(b, response, xmatrix, theta, alpha = NULL, randomization = 1){
  if(is.null(alpha)){
    alpha <- response # so that alpha does not affect anything when NULL
  }
  lp  <- (xmatrix %*% theta[-length(theta)]) +  b
  plp <- randomization * plogis(lp)
  h   <- apply(plp, 2, function(x) ifelse(response == 1, log(x/alpha), log((1 - x)/(1 - alpha))))
  hh  <- apply(h, 2, function(x) exp(sum(x)))
  hh 
}

weight_estimator_picov <- function(A, X, b, lower = -Inf, upper = Inf, randomization = 1)
{
  X <- as.matrix(X)
  f <- function(theta, alpha){
    vapply(alpha, function(x){
      w <- try(integrand_picov(theta = theta, alpha = x, response = A, xmatrix = X, b = b,
                         randomization = randomization),
               silent = FALSE)
      if(is(w, 'try-error')){
        NA
      } else {
        1/w
      } }, numeric(1))
  }
  f
}

gen_data <- function(m, ni, gamma, theta, beta){
  n <- m * ni
  data_frame(
    m   = m,
    ni  = ni,
    group = rep(1:m, each = ni),
    Z1  = rnorm(n, sd = 1),
    Z1_abs = abs(Z1),
    Z2  = rbinom(n, size = 1, prob = .5),
    X1  = exp(Z1/2),
    X2  = plogis(Z2) * Z1,
    b   = rep(rnorm(m, sd = theta), each = ni),
    p   = as.numeric(plogis((cbind(1, Z1_abs, Z2, Z1_abs*Z2) %*% gamma) + b)),
    A   = rbinom(n, size = 1, prob = p),
    fA  = rep(tapply(A, group, mean), each = ni),
    Y   = rnorm(n, mean = (cbind(1, A, fA, Z1_abs, Z2, Z1_abs*Z2) %*% beta), sd = 1)
  )
}

gen_sim <- function(nsims, ...){
  dots <- list(...)
  args <- dots[pmatch(names(dots), names(formals((gen_data))))]
  x    <- as.call(append(list(gen_data), args))
  replicate(nsims, eval(x), simplify = FALSE)
}

get_fixed_formula <- function(model){
  formula(model, fixed.only = TRUE)[-2]
}
