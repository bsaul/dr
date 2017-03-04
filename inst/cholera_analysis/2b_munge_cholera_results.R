#------------------------------------------------------------------------------#
#   Title: Analysis of cholera data with IPW, OTC, and DR estimators
#  Author: B. Saul
#    Date: 2017-03-02
# Purpose: 
#------------------------------------------------------------------------------#

de <- results %>%
  group_by(method, hajek, alpha) %>%
  summarise(estimate = -diff(estimate) * 1000)


ie <- results %>%
  filter(a == 0) %>%
  group_by(method, hajek, a) %>%
  mutate(estimate = -(estimate - estimate[alpha == 0.4]) * 1000)
