##### 
# estfUN


wls_estFUN <- function(data, wls_formula, ipw_formula, invlnk, randomization){
  
  L <- grab_design_matrix(data, wls_formula)
  L_ex <- expand_outcome_frame(data, geex::grab_fixed_formula(wls_formula))
  L_ex_0 <- L_ex %>% filter(A == 0)
  L_ex_1 <- L_ex %>% filter(A == 1)
  MM_0 <- grab_design_matrix(L_ex_0, wls_formula)
  MM_1 <- grab_design_matrix(L_ex_1, wls_formula)
  X <- grab_design_matrix(data, ipw_formula)
  Y <- grab_response(data, wls_formula)
  A <- grab_response(data, ipw_formula)
  n <- length(Y)
  p <- ncol(L)
  
  ipwFUN <- weight_estimator(
    A = A, 
    X = X, 
    randomization = randomization)
  
  function(theta, alpha, ipw_theta){
    nalpha <- length(alpha)
    index0 <- 1:p
    index1 <- (p*nalpha + 1):(p*nalpha + p)
    
    wls0_ee <- lapply(alpha, function(x){
      ipw <- ipwFUN(ipw_theta, alpha = x)/(1 - x)
      W  <- diag(ipw * (A == 0), nrow = n, ncol = n)
      ee <- t(L) %*% W %*% (Y - invlnk(L %*% theta[index0]))

      index0 <<- index0 + p
      return(ee)
    }) %>% unlist()


    wls1_ee <- lapply(alpha, function(x){
      ipw <- ipwFUN(ipw_theta, alpha = x)/(x)
      W  <- diag(ipw * (A == 1), nrow = n, ncol = n)
      ee <- t(L) %*% W %*% (Y - invlnk(L %*% theta[index1]))
      index1 <<- index1 + p
      return(ee)
    }) %>% unlist()
    
    ce0 <- ce1 <- numeric(nalpha)
    ### Regression-based DRR estimator ###
    for(k in 1:length(alpha)){
      if(k == 1){
        index0 <- 1:p
        index1 <- (p*nalpha + 1):(p*nalpha + p)
      }
      if(k > 1){
        index0 <- index0 + p
        index1 <- index1 + p
      }
      
      # compute fitted value for expanded data.frame
      mu_0 <- as.numeric(invlnk(MM_0 %*% theta[index0]))
      mu_1 <- as.numeric(invlnk(MM_1 %*% theta[index1]))
      # compute pi term per number treated in group per subject
      pi_term_a <- dbinom(L_ex_0$sum_a, n - 1, alpha[k])
      
      # mulitply mu_ij by the pi term rowwise
      piXmu_a_0 <- mu_0 * pi_term_a
      piXmu_a_1 <- mu_1 * pi_term_a
      
      # sum within levels of A (0:1) WITHIN subjects
      piXmu_a <- tapply(
        X     = c(piXmu_a_0, piXmu_a_1), 
        INDEX = paste(rep(0:1, each = nrow(L_ex_0)), c(L_ex_0$ID, L_ex_0$ID)), 
        FUN   = sum)
      
      # sum within levels of A (0:1) ACROSS subjects
      wls_ce_a <- tapply(
        X     = piXmu_a, 
        INDEX = rep(0:1, each = n), 
        FUN   = sum)
      
      ce0[k] <- wls_ce_a[1]/n
      ce1[k] <- wls_ce_a[2]/n
    }
    
    ntheta <- length(theta)
    
    c(wls0_ee, wls1_ee, 
      ce0 - theta[(ntheta - 2*nalpha + 1):(ntheta - nalpha)], 
      ce1 - theta[(ntheta - nalpha + 1):ntheta])

  }
}



library(dplyr)
library(dr)
library(geex)
library(lme4)
alphas <- list(c(.3, .4))

load( pipe( 'ssh saulb@diamond.bios.unc.edu "cat /home/groups/projects/mhudgens/emch/data/R_data/emch_analysis_data.Rdata"' ))
source('inst/cholera_analysis/1a_define_analysis_functions.R')

choleradt <- analysis_c %>%
  group_by(group) %>%
  mutate(fA = mean(A)) %>%
  filter(n() < 100) %>%
  ungroup()

analysis_model_args1 <- list(
  t_model = 
    list(method = lme4::glmer,
         formula = B ~ age + rivkm + (1|group),
         options = list(family = binomial(link = 'logit'), nAGQ = 25L)),
  o_model =
    list(method  = geepack::geeglm,
         formula = y_obs ~ A + fA + age + rivkm,
         options = list(
           family  = binomial(link = 'logit'),
           id      = quote(group))),
  wls_model_0 =
    list(method = stats::glm,
         formula = as.integer(y_obs) ~ fA + age + rivkm,
         options = list(family = quasibinomial(link = 'logit'))),
  wls_model_1 =
    list(method = stats::glm,
         formula = as.integer(y_obs) ~ fA + age + rivkm,
         options = list(family = quasibinomial(link = 'logit')))
)

models0 <- estimate_cholera_parms_step0(
  data        = choleradt,
  model_args  = analysis_model_args1
)



theta_t   <- unlist(getME(models0$t_model, c('beta', 'theta')))
wls_start <- glm(y_obs ~ fA + age + rivkm, data = choleradt, 
                 family = binomial) %>% coef()

# testfun <- wls_estFUN(data = choleradt %>% filter(group == 10),
#                        y_obs ~ fA + age + rivkm,
#                        B ~ age + rivkm,
#                        invlnk = plogis,
#                        randomization = 2/3)
# #
# #
# # 
# testfun(theta = c(rep(wls_start, times = 3), rep(0, 8)),
#         alpha = c(.3, .4), ipw_theta = theta_t)

ptm <- proc.time()
mtest_wls <- m_estimate(
  estFUN = wls_estFUN,
  data   = choleradt,
  units  = 'group',
  compute_vcov = FALSE,
  root_control = setup_root_control(start = c(rep(wls_start, times = 4), rep(0, 4))),
  outer_args = list(wls_formula = y_obs ~ fA + age + rivkm,
                    ipw_formula = B ~ age + rivkm,
                    invlnk = plogis,
                    randomization = 2/3),
  inner_args = list(alpha = c(.3, .4), ipw_theta = theta_t)
)
proc.time() - ptm
roots(mtest_wls)