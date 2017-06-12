#!/bin/bash

#SBATCH --output=attempt.out
#SBATCH --mail-type=ALL    
#SBATCH --mail-user=saulb@live.unc.edu
#SBATCH --job-name=Autom8

srun R312 CMD BATCH --vanilla diamond_example2.R diamond_example2.Rout