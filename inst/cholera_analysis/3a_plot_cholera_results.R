#------------------------------------------------------------------------------#
#   Title: Plots of cholera results
#  Author: B. Saul
#    Date: 2017-03-06
# Purpose: 
#------------------------------------------------------------------------------#

library(ggplot2)

plot_vers <- "V003"

## Color values
color_vals <- c("1ipw" = "#EFC583",
                "2otc" = rgb(86, 180, 233, max = 255),
                "3dbr" = rgb(204, 121, 167, max = 255),
                "4wls_dbr" = rgb(213, 94, 0, max = 255))


# p0 <- ggplot(
#   data = cholera_results %>% filter(effect %in% c('Y0_1', 'Y1_1')),
#   aes(x = alpha2, 
#       y = estimate,
#       group = method,
#       linetype = method,
#       color = method)
# ) + 
#   geom_hline(
#     yintercept = 0
#   ) +  
#   geom_ribbon(
#     aes(ymin = conf_low,
#         ymax = conf_high,
#         # alpha = method,
#         fill  = method),
#     alpha = .2,
#     size = .2
#   ) + 
#   geom_line() + 
#   scale_color_manual(
#     values = color_vals
#     # values = c('ipw' = "#658b83", 'otc' = "#a4044d", 'dbr'= "#359721", 'wls_dbr' = 'black')
#   ) +
#   scale_fill_manual(
#     values = color_vals
#     # values = c('ipw' = "#658b83", 'otc' = "#a4044d", 'dbr'= "#359721", 'wls_dbr' = 'black')
#   ) +
#   scale_x_continuous(
#     name = expression(alpha)
#   ) +
#   facet_wrap(
#     ~ effect, 
#     nrow = 1, 
#     labeller = labeller(
#       effect = c(Y0 = 'Y(0, alpha)',
#                  Y1 = 'Y(1, alpha)'))
#   ) + 
#   theme_light() +
#   theme(
#     # legend.position = c(.75, .25),
#     strip.text.x = element_text(color = 'black')
#   )
# p0

labs <- labeller(
  effect = c(de = "DE(~alpha~)", ie = 'Indirect Effect', te = 'Total Effect', oe = 'Overall Effect'),
  label_parsed)

cholera_results %>% filter(effect %in% c('de', 'ie', 'te', 'oe')) %>%
  mutate(plot_label = case_when(
    .$effect == 'de' ~ "'DE('*alpha*')'",
    .$effect == 'ie' ~ "'IE(0.4,'~alpha*')'",
    .$effect == 'te' ~ "'TE(0.4,'~alpha*')'",
    .$effect == 'oe' ~ "'OE(0.4,'~alpha*')'")) %>%
  # change order of methods
  mutate(method = case_when(
    method == "ipw" ~ "1ipw",
    method == "otc" ~ "2otc",
    method == "dbr" ~ "3dbr",
    method == "wls_dbr" ~ "4wls_dbr"
  )) -> plot_data

plotter <- function(effects){
  p <- ggplot(
    data = plot_data %>% filter(effect %in% !!effects),
    aes(x = alpha2, 
        y = estimate,
        group = method,
        linetype = method,
        color = method)
  ) + 
    geom_hline(
      yintercept = 0
    ) + 
    
    # geom_ribbon(
    #   aes(ymin = conf_low,
    #       ymax = conf_high,
    #       # alpha = method,
    #       fill  = method),
    #   alpha = .5,
    #   size = .2
    # ) +
    geom_line(size = 1) + 
    scale_color_manual(
      name = '',
      values = color_vals,
      # values = c('ipw' = "#658b83", 'otc' = "#a4044d", 'dbr'= "#359721", 'wls_dbr' = 'black'),
      labels = c("1ipw" = "IPW", "2otc" = "REG", "3dbr" = "DR (BC)", "4wls_dbr" = "DR (WLS)")
    ) +
    scale_linetype_discrete(
      name = '',
      # values = c('ipw' = "#658b83", 'otc' = "#a4044d", 'dbr'= "#359721", 'wls_dbr' = 'black'),
      labels = c("1ipw" = "IPW", "2otc" = "REG", "3dbr" = "DR (BC)", "4wls_dbr" = "DR (WLS)")
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
}

p1 <- plotter(c("de",  "ie"))
p2 <- plotter(c("te",  "oe"))

ggsave(p1, filename = paste0('inst/cholera_analysis/figures/cholera_results1_', plot_vers, '.pdf'),
       width = 6, height = 3 )

ggsave(p2, filename = paste0('inst/cholera_analysis/figures/cholera_results2_', plot_vers, '.pdf'),
       width = 6, height = 3 )
