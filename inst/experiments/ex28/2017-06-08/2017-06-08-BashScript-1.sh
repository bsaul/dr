#!/bin/bash -l
 #SBATCH --mail-type=FAIL 
 #SBATCH --mail-user=saulb@live.unc.edu 
 #SBATCH --job-name=NAME- 1 
 srun R312 CMD BATCH  /Users/bradley/Dropbox/Research/projects/landr/inst/experiments/ex28/2017-06-08/2017-06-08-Rscript-1.R   /Users/bradley/Dropbox/Research/projects/landr/inst/experiments/ex28/2017-06-08/2017-06-08-Rscript-1.Rout 
