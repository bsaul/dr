#------------------------------------------------------------------------------#
#   Title: Plots of cholera results
#  Author: B. Saul
#    Date: 2017-03-06
# Purpose: 
#------------------------------------------------------------------------------#

library(ggplot2)

## Color values
color_vals <- c("ipw" = "#EFC583",
                "otc" = rgb(86, 180, 233, max = 255),
                "dbr" = rgb(204, 121, 167, max = 255),
                "wls_dbr" = rgb(213, 94, 0, max = 255))


p0 <- ggplot(
  data = cholera_results %>% filter(effect %in% c('Y0_1', 'Y1_1')),
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
    values = color_vals
    # values = c('ipw' = "#658b83", 'otc' = "#a4044d", 'dbr'= "#359721", 'wls_dbr' = 'black')
  ) +
  scale_fill_manual(
    values = color_vals
    # values = c('ipw' = "#658b83", 'otc' = "#a4044d", 'dbr'= "#359721", 'wls_dbr' = 'black')
  ) +
  scale_x_continuous(
    name = expression(alpha)
  ) +
  facet_wrap(
    ~ effect, 
    nrow = 1, 
    labeller = labeller(
      effect = c(Y0 = 'Y(0, alpha)',
                 Y1 = 'Y(1, alpha)'))
  ) + 
  theme_light() +
  theme(
    # legend.position = c(.75, .25),
    strip.text.x = element_text(color = 'black')
  )
p0

labs <- labeller(
  effect = c(de = "DE(~alpha~)", ie = 'Indirect Effect', te = 'Total Effect', oe = 'Overall Effect'),
  label_parsed)

cholera_results %>% filter(effect %in% c('de', 'ie', 'te', 'oe')) %>%
  mutate(plot_label = case_when(
    .$effect == 'de' ~ "'DE('*alpha*')'",
    .$effect == 'ie' ~ "'IE(0.4,'~alpha*')'",
    .$effect == 'te' ~ "'TE(0.4,'~alpha*')'",
    .$effect == 'oe' ~ "'OE(0.4,'~alpha*')'")) -> plot_data


p <- ggplot(
  data = plot_data,
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
    alpha = .5,
    size = .2
  ) +
  geom_line(size = 1) + 
  scale_color_manual(
    name = '',
    values = color_vals,
    # values = c('ipw' = "#658b83", 'otc' = "#a4044d", 'dbr'= "#359721", 'wls_dbr' = 'black'),
    labels = c("ipw" = "IPW", "otc" = "OTC", "dbr" = "DBR", "wls_dbr" = "DBR(WLS)")
  ) +
  scale_linetype_discrete(
    name = '',
    # values = c('ipw' = "#658b83", 'otc' = "#a4044d", 'dbr'= "#359721", 'wls_dbr' = 'black'),
    labels = c("ipw" = "IPW", "otc" = "OTC", "dbr" = "DBR", "wls_dbr" = "DBR(WLS)")
  ) +
  scale_fill_manual(
    values = color_vals
    # values = c('ipw' = "#658b83", 'otc' = "#a4044d", 'dbr'= "#359721", 'wls_dbr' = 'black')
  ) +
  scale_x_continuous(
    name = expression(alpha)
  ) + 
  facet_wrap(
    ~ plot_label, 
    nrow = 2,
    ncol = 2,
    labeller = label_parsed,
    strip.position = "top"
  ) + 
  theme_light() +
  theme(
    # legend.position = c(.75, .25),
    # strip.text.x = element_text(color = 'black'),
    strip.text = element_text(color = 'black'),
    strip.background = element_blank(),
    axis.title.y = element_blank()
  )
p

ggsave(p, filename = 'inst/cholera_analysis/cholera_results.pdf',
       width = 6, height = 6 )
