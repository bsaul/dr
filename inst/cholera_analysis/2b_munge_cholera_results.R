#------------------------------------------------------------------------------#
#   Title: Analysis of cholera data with IPW, OTC, and DR estimators
#  Author: B. Saul
#    Date: 2017-03-02
# Purpose: 
#------------------------------------------------------------------------------#

library(dplyr)

# results2 <- results
# 
# dr_results2 <- lapply(dr_results, function(x) x$dbr[[1]])
# results$dbr <- dr_results2
# results[[3]]
# results

methods <- names(results)

cholera_results <- lapply(seq_along(results), function(m){
  x <- results[[m]]
  lapply(seq_along(x), function(i){
    xx <- x[[i]]

    if(xx$alpha[1] < 0.4){
      C <- matrix(
        c(1, 0, 0, 0,   # Y(0, alpha)
          0, 0, 1, 0,   # Y(1, alpha)
          1,  0, -1, 0, # Y(0, alpha') - Y(1, alpha')
          -1, 1, 0, 0,  # Y(0, .4) - Y(0, alpha')
          0, 1, -1, 0), # Y(0, .4) - Y(1, alpha')
        byrow = TRUE, ncol = 4
      )
      
      alpha1  <- xx$alpha[2]
      alpha2  <- xx$alpha[1]
    } else {
      C <- matrix(
        c(0, 1, 0, 0,   # Y(0, alpha)
          0, 0, 0, 1,   # Y(1, alpha)
          0, 1, 0, -1,
          1, -1, 0, 0,
          1, 0, 0, -1),
        byrow = TRUE, ncol = 4
      )
      
      alpha1  <- xx$alpha[1]
      alpha2  <- xx$alpha[2]
    }
    
    p <- ncol(xx$vcov)
    # Direct Effects

    
    est <- as.numeric(C %*% xx$estimate * 1000)
    
    Cvcov <- cbind(matrix(0, nrow = nrow(C), ncol = p- 4), C)
    
    std_err <- as.numeric(sqrt(diag(Cvcov %*% xx$vcov %*% t(Cvcov))) * 1000)
    
    data_frame(
      method    = methods[m],
      effect    = factor(c('Y0', 'Y1', 'de', 'ie', 'te'),
                         levels = c('Y0', 'Y1', 'de', 'ie', 'te'), ordered = TRUE),
      alpha1    = alpha1,
      alpha2    = alpha2,
      estimate  = est,
      std_err   = std_err,
      conf_low  = est - 1.96 * std_err,
      conf_high = est + 1.96 * std_err
    )
  }) %>%
    bind_rows
}) %>%
  bind_rows
cholera_results





