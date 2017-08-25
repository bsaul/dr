#------------------------------------------------------------------------------#
#   Title: Analysis of cholera data with IPW, OTC, and DR estimators
#  Author: B. Saul
#    Date: 2017-03-02
# Purpose: 
#------------------------------------------------------------------------------#

library(dplyr)

load('inst/cholera_analysis/cholera_results_dr_wls_2017-06-27.rda')
results_wls <- results

load('inst/cholera_analysis/cholera_results_dr_2017-06-25.rda')
lapply(seq_along(results), function(i){
  results[[i]][[1]]$wls_dbr <<- results_wls[[i]][[1]]$wls_dbr
})

# 
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
        C <- matrix(
          c(1, 0, 0, 0,   # Y(0, alpha)
            0, 0, 1, 0,   # Y(1, alpha)
            1,  0, -1, 0, # Y(0, alpha') - Y(1, alpha')
            -1, 1, 0, 0,  # Y(0, .4) - Y(0, alpha')
            0, 1, -1, 0), # Y(0, .4) - Y(1, alpha')
          byrow = TRUE, ncol = 4
        )
        
        alpha1  <- xxx$alpha[2]
        alpha2  <- xxx$alpha[1]
      } else {
        C <- matrix(
          c(1, 0, 0, 0,   # Y(0, alpha1)
            0, 1, 0, 0,   # Y(0, alpha2)
            0, 0, 1, 0,   # Y(1, alpha1)          
            0, 0, 0, 1,   # Y(1, alpha2)
            0, 1, 0, -1,
            1, -1, 0, 0,
            1, 0, 0, -1
            (1 - alpha1), -(1- alpha2), alpha1, -alpha2),
          byrow = TRUE, ncol = 4
        )
        
        alpha1  <- xxx$alpha[1]
        alpha2  <- xxx$alpha[2]
      }
      est <- as.numeric(C %*% xxx$estimate * 1000)
      
      p <- ncol(xxx$vcov)
      Cvcov <- cbind(matrix(0, nrow = nrow(C), ncol = p- 4), C)
      print(sqrt(diag(Cvcov %*% xxx$vcov %*% t(Cvcov))))
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



