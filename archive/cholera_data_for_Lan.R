# Make sample analysis dataset for Lan

library(dplyr)

load( pipe( 'ssh saulb@diamond.bios.unc.edu "cat /home/groups/projects/mhudgens/emch/data/R_data/emch_analysis_data.Rdata"' ))




lan_data <- analysis_c %>%
  mutate_(
    ## Add noise to age, age2, rivkm, rivkm2
    age    = ~ age + rnorm(n(), 0, 1),
    age2   = ~ age^2,
    rivkm  = ~ pmax(0, rivkm + rnorm(n(), 2, .5)),
    rivkm2 = ~ rivkm^2,
    ## Permute y, A, B
    y_obs = ~ y_obs[sample(1:n(), n(), replace = FALSE)],
    A     = ~ A[sample(1:n(), n(), replace = FALSE)],
    B     = ~ B[sample(1:n(), n(), replace = FALSE)]
  ) %>%
  group_by(group) %>%
  sample_frac(.9)

save(lan_data, file = 'cholera_analysis/lan_data.rda' )





