#!/usr/bin/env Rscript 
 library(dplyr) 
 library(dr, lib.loc='/Rlibs/') 
 source() 
 source() 
 nsims <- 1 
 which_scenarios <- 9 
 allocations <- list(c(0.1, 0.5, 0.9)) 
 estimates <- do_scenarios(nsims, which_scenarios, allocations, 
                                all_model_args = margs,
                                compute_se = TRUE,
                                verbose    = FALSE,
                                .parallel = FALSE) %>% bind_rows() 
 save(estimates, file =  'results1.rda' ) 
 ###endfile
