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
base_dir <- 'netscr/saulb/'
Today <- Sys.Date()
ScriptDumpToday <- paste(base_dir, 'drsims/', Today, "/", sep="")
dir.create(ScriptDumpToday, showWarnings=T)
ScriptDump <- function(file_name){
  ##This function date stamps the file_name
  paste(ScriptDumpToday, Today, "-", file_name, sep="")
}

### Step 1B: Set up whatever variables you want to input #### ------

num_scripts <- 1
my_output_folder <- paste0("pine/scr/s/b/saulb/MyOutputFolder/", Today, '/')
dir.create(my_output_folder, showWarnings=T)

### Step 2: Create an R script for each Simulation #####---------
for (script_num in 1:(num_scripts)) {
  ## Create a new .R file to write to
  R_script_out <- ScriptDump(paste('Rscript-',script_num, ".R", sep=""))
  L_file_out   <- ScriptDump(paste('BashScript-',script_num, ".sh", sep=""))
  R_file_out   <- ScriptDump(paste('Rscript-',script_num, ".Rout", sep=""))
  
  cat("#!/usr/bin/env Rscript \n",
      "library(dr, lib.loc='/netscr/yourOnyen/Rlibs/') \n",
      "do_simulation(nsims = 1, m = c(seq(10, 30, by = 5), 150), seeds = c(",
      paste(sample(1e6, 6, replace = FALSE), collapse = ','), "))) \n",
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

# ## Name yet another bash script
# B_file_out <- ScriptDump('ScriptAll.sh')
# 
# ## Open it
# cat("#!/bin/bash -l \n",
#     "echo \"This script can make the other scripts run\"",
#     "#SBATCH --mail-type=ALL \n",
#     "#SBATCH --mail-user=saulb@live.unc.edu \n",
#     "#SBATCH --job-name=BIGBASH \n\n", file=B_file_out)
# 
# 
# ## Now we make one line for each Bash script we want to run.
# for (script_num in 1:(num_scripts)) {
#   #Write sbatch command
#   Lscriptfile <-  ScriptDump(paste('BashScript-',script_num, ".sh", sep=""))
#   cat("sbatch ", Lscriptfile, "\n",sep="",file=B_file_out,append=TRUE)
# }



