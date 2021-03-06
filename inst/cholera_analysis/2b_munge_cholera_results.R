#------------------------------------------------------------------------------#
#   Title: Analysis of cholera data with IPW, OTC, and DR estimators
#  Author: B. Saul
#    Date: 2017-03-02
# Purpose: 
#------------------------------------------------------------------------------#

library(dplyr)

load('inst/cholera_analysis/cholera_results_2017-09-20.rda')
load('inst/cholera_analysis/cholera_results_wls_2017-09-22.rda')


results <- results_nonwls
lapply(seq_along(results_nonwls), function(i){
  ## Subset estimates to target estimates in wls
  p <- length(results_wls[[i]][[1]]$wls_dbr$estimate)
  target_index <- (p -3):p
  results_wls[[i]][[1]]$wls_dbr$estimate <<- results_wls[[i]][[1]]$wls_dbr$estimate[target_index]
  results[[i]][[1]]$wls_dbr <<- results_wls[[i]][[1]]$wls_dbr[c('alpha', 'estimate', 'vcov')]
})

# results[[1]][[1]]$wls_dbr$vcov
# results_wls[[1]][[1]]$wls_dbr$
# dr_results2 <- lapply(dr_results, function(x) x$dbr[[1]])
# results$dbr <- dr_results2
# results[[3]]
# results

methods <- names(results[[1]][[1]])


cholera_results <- lapply(seq_along(results), function(i){
  x <- results[[i]]
  lapply(seq_along(x), function(j){
    xx <- x[[j]]
    lapply(seq_along(xx), function(m){
      xxx <- xx[[m]]

      if(xxx$alpha[1] < 0.4){
        alpha1  <- xxx$alpha[2]
        alpha2  <- xxx$alpha[1]
        C <- matrix(
          c(1, 0, 0, 0,   # Y(0, alpha)
            0, 1, 0, 0, 
            0, 0, 1, 0,   # Y(1, alpha)
            0, 0, 0, 1,
            1,  0, -1, 0, # Y(0, alpha') - Y(1, alpha')
            -1, 1, 0, 0,  # Y(0, .4) - Y(0, alpha')
            0, 1, -1, 0, # Y(0, .4) - Y(1, alpha')
            (1 - alpha1), -(1- alpha2), alpha1, -alpha2),
          byrow = TRUE, ncol = 4
        )
        

      } else {
        alpha1  <- xxx$alpha[1]
        alpha2  <- xxx$alpha[2]
        
        C <- matrix(
          c(0, 1, 0, 0,   # Y(0, alpha1)
            1, 0, 0, 0,   # Y(0, alpha2)
            0, 0, 0, 1,   # Y(1, alpha1)          
            0, 0, 1, 0,   # Y(1, alpha2)
            0, 1, 0, -1,
            1, -1, 0, 0,
            1, 0, 0, -1,
            (1 - alpha1), -(1- alpha2), alpha1, -alpha2),
          byrow = TRUE, ncol = 4
        )
        
      }
      
      # print(C)
      est <- as.numeric(C %*% xxx$estimate * 1000)
      
      p <- ncol(xxx$vcov)
      Cvcov <- cbind(matrix(0, nrow = nrow(C), ncol = p- 4), C)
      # print(Cvcov)
      # print(sqrt(diag(Cvcov %*% xxx$vcov %*% t(Cvcov))))
      std_err <- as.numeric(sqrt(diag(Cvcov %*% xxx$vcov %*% t(Cvcov))) * 1000)
      # std_err <- NA
      
      data_frame(
        method    = methods[m],
        effect    = factor(c('Y0_1', 'Y0_2', 'Y1_1', 'Y1_2', 'de', 'ie', 'te', 'oe'),
                           levels = c('Y0_1', 'Y0_2', 'Y1_1', 'Y1_2', 'de', 'ie', 'te', 'oe'), ordered = TRUE),
        alpha1    = alpha1,
        alpha2    = alpha2,
        estimate  = est,
        std_err   = std_err,
        conf_low  = est - 1.96 * std_err,
        conf_high = est + 1.96 * std_err)
    })  -> hold 
    bind_rows(hold)
  }) %>%
    bind_rows
}) %>%
  bind_rows



