#------------------------------------------------------------------------------#
#   Title: Numbers for manuscript 
#  Author: B. Saul
#    Date: 2017-07-17
# Purpose: 
#------------------------------------------------------------------------------#

results %>%
  filter(model_spec == "tcor_ocor", a == 1, alpha == 0.5)

results %>%
  filter(model_spec == "tcor_omis", a == 1, alpha == 0.5)

results %>%
  filter(model_spec == "tmis_ocor", a == 1, alpha == 0.5)

results %>%
  filter(model_spec == "tmis_omis", a == 1, alpha == 0.5)



results %>%
  filter(model_spec == "tmis_omis", a == 1, alpha == 0.9)
