################################ Program Description ###
# This program will show you how to 'automate' 
# many R scripts at once onto the BIOS computing clusters 
# (Diamond or Sapphire)
################################ Program Description ###


## Let's say you have a function that you want to run, 
##  and it is called myfun(). Let's say myfun() is saved 
##  in the script myprog.R. Let's say myfun() takes a few
##  objects as input. 
##
## For the sake of argument, we will say myfun() has three
##  major steps:
##    a) myfun() first generates a simulated dataset
##    b) myfun() executes a statistical method on the dataset
##    c) myfun() saves output from the method to a file
##
## Note that if you run X=15 scripts then you will have
##  X=15 many sets of output. You will likely need to go back
##  and aggregate all those outputs together which I'll show.

## Created 2016-01-27 by Brian Barkley bribark@live.unc.edu
## 
## There are many ways to do this, and this is just one that seems to work.
## If you have better ideas, I would love to hear them! Please email me.
##
## Sorry for any typos or mistakes

# The last section is left uncommented so you have to read it!

#
##
### Step 1A: First, set up some stuff for bookkeeping ####---------------

## This section helps me to save all my functions to 
##  the same directory.
base_dir <- 'home/users/saulb/'
# base_dir <- 'temp/'
Today <- Sys.Date()
ScriptDumpToday <- paste(base_dir, 'drsims/', Today, "/", sep="")
# ScriptDumpToday <- paste(getwd(), '/inst/experiments/ex28/', Today, "/", sep="")
dir.create(ScriptDumpToday, showWarnings=T)
ScriptDump <- function(file_name){
  ##This function date stamps the file_name
  paste(ScriptDumpToday, Today, "-", file_name, sep="")
}

### Step 1B: Set up whatever variables you want to input #### ------

num_scripts <- 1
outdir <- paste0("home/users/saulb/Results/", Today, '/')
# outdir <- paste0(getwd(), '/inst/experiments/ex28/temp/')
dir.create(outdir, showWarnings=T)

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
      # "library(dr, lib.loc='/Rlibs/') \n",
      # "source() \n",
      # "source() \n",
      # "nsims <- 1 \n",
      # "which_scenarios <- 9 \n",
      # "allocations <- list(c(0.1, 0.5, 0.9)) \n",
      # "estimates <- do_scenarios(nsims, which_scenarios, allocations, 
      #                           all_model_args = margs,
      #                           compute_se = TRUE,
      #                           verbose    = FALSE,
      #                           .parallel = FALSE) %>% bind_rows() \n",
      "estimates <- data_frame(Y = rnorm(100)) \n",
      "save(estimates, file = ", results_out, ") \n",
      "###endfile\n",
      file = R_script_out)
  
  cat("#!/bin/bash -l\n",
      "#SBATCH --mail-type=FAIL \n",
      "#SBATCH --mail-user=saulb@live.unc.edu \n",
      "#SBATCH --job-name=NAME-", script_num, '\n',
      "srun R312 CMD BATCH ", R_script_out, " ", R_file_out, "\n",
      file=L_file_out)
}

### Step 6: Create an Overall Bash Script that Runs all Simulations #####

## Name yet another bash script
B_file_out <- ScriptDump('ScriptAll.sh')

## Open it
cat("#!/bin/bash -l \n",
    "echo \"This script can make the other scripts run\"",
    "#SBATCH --mail-type=ALL \n",
    "#SBATCH --mail-user=saulb@live.unc.edu \n",
    "#SBATCH --job-name=BIGBASH \n\n", file=B_file_out)


## Now we make one line for each Bash script we want to run.
for (script_num in 1:(num_scripts)) {
  #Write sbatch command
  Lscriptfile <-  ScriptDump(paste('BashScript-',script_num, ".sh", sep=""))
  cat("sbatch ", Lscriptfile, "\n",sep="",file=B_file_out,append=TRUE)
}



