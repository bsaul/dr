library(lme4)
library(dplyr)
library(dr)

#### Generate simulated datasets #### 
sim_data <- gen_sim(
  nsims = 50,
  m = 100, #  of groups
  ni = 30, #  of subject per group  
  gamma = c(0.1, 0.2, 0.0, 0.2),
  theta = 0.3,
  beta  = c(2.0, 2.0, 1.0, -1.5, 2.0, -3.0))


#### Compute Yhat(little_a, ALPHA) per simulation ####
ALPHA <- 0.5
little_a <- 0

lapply(sim_data, function(dt){
  
  t_model <- lme4::glmer(
    A ~ Z1_abs + Z2 + Z1_abs*Z2 + (1|group),
    data   = dt,
    family = binomial(link = 'logit')
  )
  
  ## Grab component model parameter estimates
  theta_t <- unlist(lme4::getME(t_model, c('beta', 'theta')))
  dt$bhat <- lme4::getME(t_model, 'Z') %*% lme4::getME(t_model, 'b')
  ## Compute pi/f(A) 
  lapply(split(dt, f = dt$group), function(x){
    ip_fun <- weight_estimator(
      A = x$A, 
      X = model.matrix(get_fixed_formula(t_model), data = x))
    
    IPW <- ip_fun(theta_t, alpha = ALPHA)
    # IPW has pi = prod_i^n, so to remove contribution of jth subject,
    # divide by alpha^A_{ij} (1 - alpha)^(1 - A_ij)
    IPW / dbinom(x$A, 1, prob = ALPHA) 
  }) %>% unlist -> 
    dt$ipw
  
  ## fit outcome model for a == 0
  picovm <- glm(
    Y ~ fA + Z1 + Z2 + ipw, # incorrect
    # Y ~ fA + Z1_abs + Z2 + Z1_abs*Z2 + ipw, # correct,
    data = dt, 
    weights = (dt$A == little_a) * 1,
    family = gaussian(link = 'identity'))
  
  N <- nrow(dt)
  
  ## Method 1 ##
  k <- 1 # number of resamples for A_tilde
  out <- numeric(k)
  for(i in 1:k){
    # Sample an A_tilde for each k
    A_tilde     <- rbinom(N, 1, prob = ALPHA)
    
    # Replace A with A_tilde
    new_data    <- dt
    new_data$A  <- A_tilde
    
    # With A_tilde, compute pi/f PER subject after setting 
    # a given subject's a_ij to 0. This needs to happen within a group,
    # hence split the data by group, then do the computations.
    split_data <- split(new_data, f = new_data$group)

    m_i <- lapply(split_data, function(grp_data){
      hold <- numeric(nrow(grp_data))
      # compute ip weight and m_ij PER subject
      for(j in 1:nrow(grp_data)){
        A_new    <- grp_data$A
        A_new[j] <- little_a # set A_ij to a_ij = 0
        
        ip_fun <- weight_estimator(
          A = A_new,
          X = model.matrix(get_fixed_formula(t_model), data = grp_data))

        grp_data$ipw[j] <- ip_fun(theta_t, ALPHA)/( (ALPHA^little_a) * ((1 - ALPHA)^(1 - little_a)) )
        # IPW has pi = prod_i^n, so to remove contribution of jth subject,
        # divide by ALPHA^little_a * (1 - ALPHA)^(1 - little_a0)

        j_data    <- grp_data[j, ]
        j_data$fA <- mean(A_new)
        m_ij <- predict(picovm, newdata = j_data)
        hold[j] <- m_ij
      }
      hold
    })

    m_i <- unlist(m_i)

    out[i] <- mean(tapply(m_i, new_data$group, mean))
  }

  mean(out)
    
})  -> results

results <- unlist(results)
results

truth <- 3.106346^little_a * (1.106346^(1 - little_a)) # Y(0, ALPHA)

mean(results - truth)
