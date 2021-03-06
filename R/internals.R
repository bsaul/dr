#------------------------------------------------------------------------------#
#' Binomial integrand
#' 
#' For logistic-Normal model
#' 
#' @export
#------------------------------------------------------------------------------#

integrand <- function(b, response, xmatrix, theta, alpha = NULL, randomization = 1){
  if(is.null(alpha)){
    alpha <- response # so that alpha does not affect anything when NULL
  }
  lp  <- outer(xmatrix %*% theta[-length(theta)], b, '+')
  plp <- randomization * plogis(lp)
  h   <- apply(plp, 3, function(x) ifelse(response == 1, log(x/alpha), log((1 - x)/(1 - alpha))))
  hh  <- apply(h, 2, function(x) exp(sum(x)))
  hh * dnorm(b, mean = 0, sd = theta[length(theta)])
}

#------------------------------------------------------------------------------#
#' pi term
#' @export
#------------------------------------------------------------------------------#

pi_term <- function(A){
  f <- function(alpha){
    prod(dbinom(x = A, 1, prob = alpha))
  }
  memoise::memoise(f)
}

grids <- new.env()
#------------------------------------------------------------------------------#
#' Expand a grid 
#' 
#' For use in the outcome estimator
#' 
#' @export
#------------------------------------------------------------------------------#

expand_grid_n <- function(n1, n2)
{
  lookup <- paste(n1, n2, sep = '_')
  if(!exists(lookup, envir = grids)){
    grids[[lookup]] <- expand.grid(ID = 1:n1, sum_a = 0:n2,
                                   A = 0:1) %>%
      mutate_(fA = ~ sum_a / n1,
              `A:fA` =~ A * fA)
  }
  return(grids[[lookup]])
}

#------------------------------------------------------------------------------#
#' Expand a outcome data frame
#' 
#' For use in the outcome estimator
#'
#' @export
#------------------------------------------------------------------------------#

expand_outcome_frame <- function(X_outcome, rhs_formula_outcome){
  n <- nrow(X_outcome)
  X_outcome %>%
    mutate_(ID = ~ row_number()) %>%
    # Remove treatment variables so that they are updated
    # select(-A, -fA, -`A:fA`) %>%
    select(-A, -fA) %>%
    # Generate the relevant values of A and
    # all possible sum(a_i) for each subject
    full_join(expand_grid_n(n, n - 1), by = "ID")
}
