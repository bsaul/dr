## ----filehead, echo = FALSE----------------------------------------------
#   Title: Doubly Robust Estimator with Interference
#  Author: B. Saul
#    Date: 2016-03-01
# Purpose: code the DR estimator in Liu et al

## ----libraries, echo=TRUE, message=FALSE, warning=FALSE------------------
library(inferference)
library(dplyr)
library(geepack)
library(magrittr)
#library(plyr)

## ---- echo=TRUE----------------------------------------------------------
vaccinesim <- vaccinesim %>%
  mutate(id = row_number()) %>%
  group_by(group) %>%
  # compute p_ij for each subject
  mutate(#pA = (sum(A) - A)/(n() - 1),
         pA = sum(A)/(n() ),
         # py = ifelse(sum(y)==0, .5/n(), sum(y)/n()),
         py = sum(y)/n(),
         By = log(py/(1-py))) %>%
  ungroup() 

## ----outcome_model, echo=TRUE--------------------------------------------
m_outcome <- geepack::geeglm(y ~ A*pA + X1 + X2, data = vaccinesim,
                             id = id, family = binomial)

## ----ipw_model, echo = TRUE----------------------------------------------
alphas <- c(.2, .3, .4, .5, .6, .7)

ipw_pieces <- interference(formula = y | A | B ~ X1 + X2 + (1|group) | group,
                     data = as.data.frame(vaccinesim), 
                     allocations = alphas,
                     method = 'simple')

## ----term1, echo=TRUE----------------------------------------------------
## First term ####
dr_term1 <- function(data, alpha, ipw_part, outcome_model)
{
  # Get the weights from the interference object
  weights <- ipw_part$weights[ , as.character(alpha)]
  weights <- data.frame(group = as.integer(names(weights)), w = weights)
  
  data %>%   
    # Compute mu_ij for each subject
    # mutate_(muA = ~predict(outcome_model, type = 'response')) %>%
    mutate_(muA = ~ model.matrix(outcome_model$formula, .) %*% coef(outcome_model) ) %>%
    # Add group IPW weights to the data
    left_join(weights, by = 'group') %>%
    # Compute term 1
    mutate_(term1 = ~ plogis(By - muA) * w )
}

## ----term2, echo=TRUE----------------------------------------------------
## Second term ####
dr_term2 <- function(data, alpha, outcome_model)
{
  data %>%
    select(-pA) %>%
    # this needs to be done by group. plyr splits the data frame by group
    # then applies a function to each piece
    {plyr::dlply(., plyr::.(group), function(x){
      x %>%
        # Generate all possible sum(a_i) for each subject
        merge(expand.grid(sum_a = 0:(n_distinct(.$id) - 1)), all = T)  %>%
        group_by_(~id) %>%
        # Compute pi and p_ij for each sum(a_i)
        mutate_(pA = ~ sum_a/(n() - 1),
                pi = ~ alpha^sum_a * (1 - alpha)^(n() - 1 - sum_a)) %>%
        ungroup() %>%
        # Compute mu_ij for each a_i per subject
        # For some reason predict() is throwing error with geeglm model
        #mutate_(mu_ij = ~predict(outcome_model, newdata = ., type = 'response')) %>%
        mutate_(mu =~ as.numeric(plogis(model.matrix(outcome_model$formula, .) %*% coef(outcome_model) ) ) ) %>%
        # mutate_(mu =~0) %>%
        # Sum by individual to compute term2
        group_by_(~id) %>%
        summarize_(term2 = ~sum(mu * pi)) 
    })} %>%
  # Put the split data back together
  bind_rows()
}


## ----DR, echo = TRUE-----------------------------------------------------
## All together now ####
dr <- function(data, alpha, ipw_part, outcome_model)
{
 dr_term1(data, alpha, ipw_part, outcome_model) %>%
    left_join(dr_term2(data, alpha, outcome_model), by = 'id') %>%
    mutate_(Y_hat_ij = ~ term1 + term2  ) %>%
    group_by_(~group) %>%
    summarise_(Y_hat_i = ~mean(Y_hat_ij)) %>%
    summarise_(Y_hat = ~mean(Y_hat_i) )
}

## ----DR_estimate, echo=TRUE----------------------------------------------
lapply(alphas, function(alpha) dr(vaccinesim, alpha, ipw_pieces, m_outcome)) %>%
  bind_rows()

## ----IPW_estimate, echo=TRUE---------------------------------------------
ipw_pieces$estimates %>%
  filter(effect_type == 'outcome', marginal == TRUE) %>%
  select(alpha1, estimate)

# x[2, 1] - x[5, 1]
# y
# y[2, 2] - y[5, 2]
