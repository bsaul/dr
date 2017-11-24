nrow(choleradt)
length(unique(choleradt$group))


cholera_results %>%
  filter(effect == "oe", alpha2 == 0.3)

cholera_results %>%
  filter(effect == "oe", alpha2 == 0.6)
