#------------------------------------------------------------------------------#
#   Title: Doubly Robust Functions
#  Author: B. Saul
#    Date: 2016-04-20
# Purpose: functions for IPW and DR estimators
#------------------------------------------------------------------------------#

#### Functions ####

integrand <- function(b, response, xmatrix, theta){
  lp <- outer(xmatrix %*% theta[-length(theta)], b, '+')
  h  <- apply(lp, 3, function(x) dbinom(response, 1, plogis(x) ) )
  hh <- apply(h, 2, prod)
  hh * dnorm(b, mean = 0, sd = theta[length(theta)])
}

weight_estimator <- function(A, X)
{
  X <- as.matrix(X)
  function(theta){
    1/integrate(integrand, lower = -Inf, upper = Inf,
                theta = theta, response = A, xmatrix = X)$value
  }
}

pi_term <- function(A){
  f <- function(alpha){
    prod(dbinom(x = A, 1, prob = alpha))
  }
}

make_dr_term1 <- function(X){
  X <- as.matrix(X)
  function(theta){
    X %*% theta
  }
}

make_dr_term2 <- function(X, formula){
  function(theta, alpha, a){
    n <- nrow(X) - {if(is.null(a)) 0 else 1}
    
    X %>%
      mutate(ID = row_number()) %>%
      select(-A) %>%
      # Generate all possible sum(a_i) for each subject
      merge(expand.grid(sum_a = 0:n,
                        A = {if(is.null(a)) 0:1 else a }),
            all = T) %>%
      group_by_(~ID) %>%
      # Compute pi and p_ij for each sum(a_i)
      mutate_(fA = ~ sum_a/n(),
              # fAn = ~ sum_a/n(),
              pi = ~ dbinom(sum_a, n, prob = alpha) ) %>%
      ungroup()  %>%
      # Compute mu_ij for each a_i per subject
      mutate_(mu =~ as.numeric(model.matrix(formula[-2], data = .)%*% theta )  ) %>%
      # Sum by individual to compute term2
      group_by_(~ID) %>%
      summarize_(term2 = ~sum(mu * pi)) %>%
      summarize_(mean = ~mean(term2)) %>%
      as.numeric()
  }
}

make_dr_estimator <- function(Y, A, X_outcome, X_treatment, formula_outcome){
  w <- weight_estimator(A = A, X = X_treatment)
  pi_t <- pi_term(A = A)
  dr_term1 <- make_dr_term1(X_outcome)
  dr_term2 <- make_dr_term2(X_outcome, formula_outcome)
  
  q_treatment <- ncol(X_treatment) + 1

  function(theta, alpha, a = NULL){
    Ia <- if(is.null(a)) 1 else (A == a) * 1
    ( mean(Y * Ia - dr_term1(theta[(q_treatment + 1):length(theta)]) ) * w(theta[1:q_treatment]) * pi_t(alpha) /
      (dbinom(a, 1, alpha) * !is.null(a) * 1) ) +
      dr_term2(theta[(q_treatment + 1):length(theta)], alpha, a)
  }
}

psi <- function(group_estimator){
  function(theta, alpha, a){
    q <- length(theta)
    group_estimator(theta[1:(q-1)], alpha, a) - theta[q]
  }
}

make_group_estimators <- function(split_data,
                                  estimator,
                                  formula_outcome)
{
  out <- lapply(split_data, function(group_data)
  {
    X_treatment <- group_data[['X_treatment']]
    X_outcome <- group_data[['X_outcome']]
    A <- group_data[['treatment']]
    Y <- group_data[['outcome']]
    eval(call(estimator, Y = Y, A = A,
              X_treatment = X_treatment,
              X_outcome = X_outcome,
              formula_outcome = formula_outcome))
  } )
  return(out)
}


list_matrix <- function(this_list)
{
  ulist <- unlist(this_list)
  m <- length(this_list)
  p <- length(ulist)/m
  matrix(ulist, nrow = m, ncol = p, byrow = T)
}

rbind_fill_zero <- function(this_list)
{
  p <- max(unlist(lapply(this_list, ncol)))
  out <- lapply(this_list, function(x){
    xp <- ncol(x)
    if(xp < p ){
      cbind(x,matrix(0, nrow = nrow(x), ncol = p - xp))
    } else {
      x
    }
  })
  do.call('rbind', out)
}