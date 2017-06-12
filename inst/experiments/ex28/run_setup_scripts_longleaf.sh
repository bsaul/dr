#!/bin/bash

#SBATCH --output=attempt.out
#SBATCH --mail-type=ALL    
#SBATCH --mail-user=saulb@live.unc.edu
#SBATCH --job-name=Autom8
module add r/3.3.1
srun R CMD BATCH setup_scripts_longleaf.R setup_scripts_longleaf.Rout