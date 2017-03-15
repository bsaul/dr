#------------------------------------------------------------------------------#
#   Title: Plots of cholera results
#  Author: B. Saul
#    Date: 2017-03-06
# Purpose: 
#------------------------------------------------------------------------------#

library(ggplot2)


p <- ggplot(
  data = cholera_results %>% filter(effect %in% c('de', 'ie', 'te')),
  aes(x = alpha2, 
      y = estimate,
      group = method,
      linetype = method,
      color = method)
) + 
  geom_hline(
    yintercept = 0
  ) + 

  geom_ribbon(
    aes(ymin = conf_low,
        ymax = conf_high,
        # alpha = method,
        fill  = method),
    alpha = .2,
    size = .2
  ) + 
  geom_line() + 
  scale_color_manual(
    values = c('ipw' = "#658b83", 'otc' = "#a4044d", 'dbr'= "#359721")
  ) +
  scale_fill_manual(
    values = c('ipw' = "#658b83", 'otc' = "#a4044d", 'dbr'= "#359721")
  ) +
  scale_x_continuous(
    name = expression(alpha)
  ) + 
  facet_wrap(
    ~ effect, 
    nrow = 1, 
    labeller = labeller(
      effect = c(de = 'Direct Effect',
                 ie = 'Indirect Effect',
                 te = 'Total Effect'))
    ) + 
  theme_light() +
  theme(
    # legend.position = c(.75, .25),
    strip.text.x = element_text(color = 'black')
  )
p

ggsave(p, filename = 'inst/cholera_analysis/cholera_results.pdf',
       width = 8, height = 3 )
