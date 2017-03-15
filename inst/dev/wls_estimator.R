
library(dr)
library(dplyr)
exdt <- gen_data(20, 10, c(0.1, 0.2, 0, .2), 1, c(2, 2, 1, -1.5, 2, -3))

mods <- make_models(model_args = list(
    t_model = 
      list(method = lme4::glmer,
           formula = A ~ Z1_abs + Z2 + Z1_abs*Z2 + (1|group),
           options = list(family = binomial(link = 'logit'))),
    o_model =
      list(method  = geepack::geeglm,
           formula = Y ~ A + fA + Z1_abs + Z2 + Z1_abs*Z2,
           options = list(
             family  = gaussian(link = 'identity'),
             id      = quote(group)))),
    data = exdt)


wls_preprocess <- function(data, models, randomization = 1){
  splitdt <- split(data, data$group)
  lapply(splitdt, function(x){
    comp <- extract_model_info(model = models, data = x, 'dbr')
    Y   <- comp$Y
    A   <- comp$A
    N   <- comp$N
    X_o <- comp$X_o
    f_o <- terms.formula(comp$rhs_o)
    L   <- model.matrix(drop.terms(f_o, attr(f_o, 'term.labels') == 'A'),
                        data = x)

    ## components for IPW part
    ip_fun <- weight_estimator(
      A = comp$A, 
      X = comp$X_t, 
      randomization = randomization)
    W_WLS <- function(a, alpha){
      IPW <- ip_fun(unlist(lme4::getME(models$t_model, c('beta', 'theta'))), alpha = alpha)
      IPW <- IPW / dbinom(A, 1, prob = alpha)
      (A == a) * t(L) * IPW
    }

  })
}

wls_preprocess(exdt, mods)[[1]](1, .5)




dbr_wls_estimator <- function(data, models, randomization, hajek, ...){
  
  ## component data
  comp <- extract_model_info(model = models, data = data, 'dbr')
  Y <- comp$Y
  A <- comp$A
  N <- comp$N
  ## components for IPW part
  ip_fun <- weight_estimator(
    A = comp$A, 
    X = comp$X_t, 
    randomization = randomization)
  ## components for OTC part
  dr_term1_fun <- make_dr_term1(
    comp$X_o, 
    inv_link = comp$inv_link_o)
  
  dr_term2_fun <- otc_estimator(data, models, randomization)
  
  ## indices
  p_t <- comp$p_t
  p_o <- comp$p_o
  p   <- p_t + p_o
  index_t <- 1:p_t
  index_o <- (p_t + 1):(p_t + p_o)
  
  function(theta, alpha){
    
    fY      <- dr_term1_fun(theta[index_o])
    ipw     <- ip_fun(theta[index_t], alpha)
    Ybar0   <- sum((A == 0) * (Y - fY) )
    Ybar1   <- sum((A == 1) * (Y - fY) )
    term1_0 <- Ybar0 * ipw / (1 - alpha)
    term1_1 <- Ybar1 * ipw / alpha
    term2   <- dr_term2_fun(theta[index_o], alpha)
    
    dbr_0 <- (term1_0 + term2[1:length(alpha)]*N)
    dbr_1 <- (term1_1 + term2[(length(alpha) + 1):length(term2)]*N)
    
    if(hajek){
      Nhat0 <- sum(A == 0) * ipw / (1 - alpha)
      Nhat1 <- sum(A == 1) * ipw / alpha
      out <- c(
        ifelse(sum(A == 0) > 0, term1_0/Nhat0 + term2[1], 0),
        ifelse(sum(A == 1) > 0, term1_1/Nhat1 + term2[2], 0))
      names(out) <- paste0(rep(c('dbr_hjk_Y0_', 'dbr_hjk_Y1_'), each = length(alpha)), alpha)
    } else {
      out <- c(dbr_0/N, dbr_1/N)
      names(out) <- paste0(rep(c('dbr_Y0_', 'dbr_Y1_'), each = length(alpha)), alpha)
    }
    out
  }
}
