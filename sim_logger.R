#------------------------------------------------------------------------------#
#   TITLE: Simulation Logger
#  AUTHOR: Bradley Saul
#    DATE: 5/02/16
# PURPOSE:
#------------------------------------------------------------------------------#

ts <- Sys.time()
# estimation_file <- 'DR_sim_estimation_6'
experimentID <- 'X006'

simdir <- 'logs/'
filename <- paste0(ts, '_', experimentID, '.log')
fn <- paste0(simdir, filename)

sink(fn, append = TRUE, type = 'output')
# sink(fn, append = TRUE, type = 'message')

timestamp(prefix = "##------ Begin: ")

# Notes

# Create similuated dataset
source('R/create_sims.R', echo = T, max.deparse.length = 5000)

# Load necessary functions
source('R/functions_v5.R', echo = T, max.deparse.length = 5000)

# Run esimtations
source('development/DR_sim_estimation_9.R', echo = T, max.deparse.length = 5000)

timestamp(prefix = "##------ End: ")

sink()

rm(ts, simdir, filename, fn)
