
##  Originally created 2016-01-27 by Brian Barkley bribark@live.unc.edu

### Step 1A: First, set up some stuff for bookkeeping ####---------------

## This section helps me to save all my functions to 
##  the same directory.
# base_dir <- 'pine/scr/s/b/saulb/'
Today <- Sys.Date()
ScriptDumpToday <- paste(Today, "/", sep="")
dir.create(ScriptDumpToday, showWarnings=T)
ScriptDump <- function(file_name){
  ##This function date stamps the file_name
  paste(ScriptDumpToday, Today, "-", file_name, sep="")
}

### Step 1B: Set up whatever variables you want to input #### ------

basedir <- "/nas/longleaf/home/saulb"

num_scripts <- 3
outdir <- paste0("results/", Today, '/')
# outdir <- paste0(getwd(), '/inst/experiments/ex28/temp/')
dir.create(outdir, showWarnings=T)
libloc <- paste0("'", basedir, '/Rlibs', "'")
funsR  <- paste0("'", basedir, "/drsims/ex28_funs.R", "'")
settingsR  <- paste0("'", basedir, "/drsims/ex28_settings.R", "'")
### Step 2: Create an R script for each Simulation #####---------
for (script_num in 1:(num_scripts)) {
  ## Create a new .R file to write to
  R_script_out <- ScriptDump(paste('Rscript-',script_num, ".R", sep=""))
  L_file_out   <- ScriptDump(paste('BashScript-',script_num, ".sh", sep=""))
  R_file_out   <- ScriptDump(paste('Rscript-',script_num, ".Rout", sep=""))
  results_out  <- paste0("'results", script_num, ".rda'")
  seed         <- paste0("set.seed(", script_num, ")")
  cat("#!/usr/bin/env Rscript \n",
      "library(dplyr) \n",
      seed, "\n",
      "library(dr, lib.loc=", libloc, ") \n",
      "source(", funsR, ") \n",
      "source(", settingsR, ") \n",
      "nsims <- 1 \n",
      if(num_scripts > 700){ "margs <- margs[1] \n" },
      "which_scenarios <- 9 \n",
      "allocations <- list(c(0.1, 0.5, 0.9)) \n",
      "estimates <- do_scenarios(nsims, which_scenarios, allocations, 
                                 all_model_args = margs,
                                 compute_se = TRUE,
                                 verbose    = FALSE,
                                 .parallel = FALSE) %>% bind_rows() \n",
      # "estimates <- data_frame(Y = rnorm(100)) \n",
      "save(estimates, file = ", results_out, ") \n",
      "###endfile\n",
      file = R_script_out)
  
  cat("#!/bin/bash -l\n",
      "#SBATCH --mail-type=FAIL \n",
      "#SBATCH --mail-user=saulb@live.unc.edu \n",
      "#SBATCH --job-name=NAME-", script_num, '\n',
      "srun R CMD BATCH ", R_script_out, " ", R_file_out, "\n",
      file=L_file_out)
}

### Step 6: Create an Overall Bash Script that Runs all Simulations #####

## Name yet another bash script
B_file_out <- ScriptDump('ScriptAll.sh')

## Open it
cat("#!/bin/bash -l \n",
    "echo \"This script can make the other scripts run\"\n",
    "#SBATCH --mail-type=ALL \n",
    "#SBATCH --mail-user=saulb@live.unc.edu \n",
    "#SBATCH --job-name=BIGBASH \n\n", file=B_file_out)


## Now we make one line for each Bash script we want to run.
for (script_num in 1:(num_scripts)) {
  #Write sbatch command
  Lscriptfile <-  ScriptDump(paste('BashScript-',script_num, ".sh", sep=""))
  cat("sbatch ", Lscriptfile, "\n",sep="",file=B_file_out,append=TRUE)
}



