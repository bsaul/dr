#------------------------------------------------------------------------------#
#   TITLE: Determine sample size for simulations
#    DATE: 2017-06-05
#  AUTHOR: Bradley Saul
#   NOTES: targets coverage level, which should be close to 95%
#------------------------------------------------------------------------------#

sim_power <- function(target_p, accuracy, alpha){
  (target_p * (1 - target_p) * qnorm(1 - alpha/2)^2 / accuracy^2)
}


# to second decimal accuracy for coverage of 95%, .01 level
sim_power(.95, 0.005, .01) ## chosen sample size (rounded to 13000)

# to first decimal accuracy for coverage of 5%, .01 level
sim_power(.5, 0.05, .01) ## chosen sample size (rounded to 700)
