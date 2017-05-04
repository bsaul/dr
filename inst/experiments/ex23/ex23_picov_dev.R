library(lme4)

#### Generate data #### 
sim_data <- gen_sim(
  nsims = 10,
  m = 100, #  of groups
  ni = 30, #  of subject per group  
  gamma = c(0.1, 0.2, 0.0, 0.2),
  theta = 0.3,
  beta  = c(2.0, 2.0, 1.0, -1.5, 2.0, -3.0))


#### Compute Yhat(0, 0.5) per simulation ####

lapply(sim_data, function(dt){
  
  t_model <- lme4::glmer(
    A ~ Z1_abs + Z2 + Z1_abs*Z2 + (1|group),
    data   = dt,
    family = binomial(link = 'logit')
  )
  
  ## Grab component model parameter estimates
  theta_t <- unlist(lme4::getME(t_model, c('beta', 'theta')))
  
  ## Compute pi/f(A) 
  lapply(split(dt, f = dt$group), function(x){
    ip_fun <- weight_estimator(
      A = x$A, 
      X = model.matrix(get_fixed_formula(t_model), data = x), 
      randomization = 1)
    
    IPW <- ip_fun(theta_t, alpha = 0.5)
    IPW / dbinom(x$A, 1, prob = 0.5)
  }) %>% unlist -> 
    dt$ipw
  
  ## fit outcome model for a == 0
  picovm <- glm(
    Y ~ fA + Z1 + Z2 + ipw, # incorrect
    # Y ~ fA + Z1_abs + Z2 + Z1_abs*Z2 + ipw, # correct,
    data = dt, 
    weights = (dt$A == 0) * 1,
    family = gaussian(link = 'identity'))
  
  N <- nrow(dt)
  
  ## Method 1 ##
  k <- 10 # number of resample for A_tilde
  out <- numeric(k)
  for(i in 1:k){
    A_tilde     <- rbinom(N, 1, prob = .5)
    new_data    <- dt
    new_data$A  <- A_tilde

    split_data <- split(new_data, f = new_data$group)

    m_i <- lapply(split_data, function(grp_data){
      hold <- numeric(nrow(grp_data))
      # compute ip weight and m_ij PER subject
      for(j in 1:nrow(grp_data)){
        A_new    <- grp_data$A
        A_new[j] <- 0 # set A_ij to a_ij = 0

        ip_fun <- weight_estimator(
          A = A_new,
          X = model.matrix(get_fixed_formula(t_model), data = grp_data),
          randomization = 1)

        grp_data$ipw[j] <- ip_fun(theta_t, .5)

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
  #   ## Method 2 ## - testing
  #   resamples <- 100
  # 
  #   split_dt <- split(dt, f = dt$group)
  #   lapply(split_dt, function(x){
  #     ni   <- nrow(x)
  # 
  #     # Generate a_i vectors #
  #     aMAT <- matrix(0, nrow = ni, ncol = resamples)
  #     for(j in 1:resamples){
  #       # how many get treated?
  #       ntreated <- sample(0:ni, 1)
  #       # which ones get treated
  #       index_treated <- sample(0:ni, ntreated)
  #       # replace 0 with 1 in aMAT for those treated
  #       
  #       aMAT[index_treated, j] <- 1
  #     }
  #     
  #     # Compute m_i
  #     mpi <- numeric(ni)
  #     # compute ip weight and m_ij PER subject
  #     for(j in 1:nrow(x)){
  #       
  #       apply(aMAT, 2, function(a_i){
  #           newdata <- x
  #           newdata$fA <- mean(a_i)
  #           a_i[j] <- 0 # set A_ij to a_ij = 0
  #         
  #         ip_fun <- weight_estimator(
  #           A = a_i, 
  #           X = model.matrix(geex::get_fixed_formula(t_model), data = x), 
  #           randomization = 1)
  #         
  #         newdata$ipw[j] <- ip_fun(theta_t, .5)
  #         
  #         j_data    <- x[j, ]
  #         j_data$fA <- mean(a_i)
  #         m_ij <- predict(picovm, newdata = j_data)
  #         pi <- prod(.5^a_i * (1 - .5)^(1 - a_i))/(1 - .5)
  #         # print(m_ij)
  #         m_ij * pi
  #       }) -> hold
  #       mpi[j] <- (sum(hold) * 2^ni / resamples)
  #     }
  #     mean(mpi)
  #   }) -> grp_estimator
  # 
  # mean(unlist(grp_estimator))
    
})  -> results

results <- unlist(results)
results

truth <- 1.106346 # Y(0, 0.5)