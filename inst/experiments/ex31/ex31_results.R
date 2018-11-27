#------------------------------------------------------------------------------#
#   Title: Results for ex30 and ex29
#  Author: B. Saul
#    Date: 2017-03-08
# Purpose: 
#------------------------------------------------------------------------------#

library(ggplot2)
library(ggbeeswarm)
library(grid)
library(gridExtra)

vers <- "V002"

# load(file = 'inst/experiments/ex31/ex31_results.rda')

# estimates <- estimates %>%
#   ungroup() %>%
#   mutate_(method =~ paste0(method, ifelse(regression_type == 'none', '', regression_type)))
# 
# results <- results %>%
#   ungroup() %>%
#   mutate_(method =~ paste0(method, ifelse(regression_type == 'none', '', regression_type)))

### Forest plots

# ggplot(
#   data = estimates %>% filter(alpha == 0.5, a == 0, model_spec == 'tcor_ocor') %>%
#     group_by(method) %>% arrange(estimate) %>% mutate(plotid = 1:n()),
#   aes(x = plotid,
#       y = estimate,
#       color = covered)) +
#   geom_hline(
#     yintercept = 1.106346
#   ) +
#   geom_point() +
#   geom_segment(
#     aes(xend = plotid,
#         y    = conf.low,
#         yend = conf.high)
#   ) + facet_grid(
#     method ~ .,
#     scale = 'free_y'
#   )


### Level one
plot_one <- function(.sid, .model_spec, .method, .alpha, .a, .axis_labels = FALSE, .method_labels = FALSE){
  
  estimates <- estimates %>%
    ungroup() %>%
    filter_(~ sid == .sid, ~model_spec == .model_spec, ~method == .method, ~alpha == .alpha, ~a == .a) %>%
    mutate(log10absbias = -log10(abs(bias)))
  
  summary_vals <- results %>% ungroup() %>%
    filter_(~ sid == .sid, ~ model_spec == .model_spec, ~method == .method, ~alpha == .alpha, ~a == .a) %>%
    mutate(log10absbias = -mean_lbias)
    # mutate(log10absbias = -log10(abs(mean_bias))) 
  
  # Plot prep  
  y_coverage <- round(summary_vals$coverage, 2)
  
  if(.axis_labels){
    yaxis_text <- element_text(size = 7, hjust = 1)
  } else {
    yaxis_text <- element_blank()
  }
  
  if(.method == 'ipw'){
    # color <- "#658b83"
    # color <- rgb(230, 159, 0, .6*255, max = 255)
    color <- "#EFC583"
    method_lab <- 'IPW'
  } else if(.method == 'otc'){
    # color <- "#a4044d"     
    color <- rgb(86, 180, 233, max = 255)
    method_lab <- 'REG'
  } else if(.method == 'dbr'){
    # color <- "#359721"
    # color <- rgb(0, 114, 178, max = 255)
    color <- rgb(204, 121, 167, max = 255)
    method_lab <- 'DR\n(BC)'
  } else if(.method =='wls_dbr'){
    # color <- "#94d3bc"
    color <- rgb(213, 94, 0, max = 255)
    method_lab <- 'DR\n(WLS)'
  } 
  # else if(.method =='pcov_dbrpcov'){
  #   color <- "#0b5313"
  #   method_lab <- 'DBR\n(PCov)'
  # }
  
  ## Bias plot
  bias_plot <- ggplot(estimates, 
                      aes(x = 1, y =log10absbias )) +
    geom_hline(
      yintercept = c(0, 1, 2, 3), 
      linetype = 'solid',
      color    = 'grey85'
    ) +
    geom_violin(
      scale = 'width',
      fill = color,
      size = .25,
      color = 'grey50'
    ) +
    geom_point(data = summary_vals) +
    scale_x_continuous(
      name = '',
      expand = c(0,.1)
    ) + 
    scale_y_continuous(
      name   = '',
      breaks = -c(0, -1, -2, -3),
      labels = c(1, 0.1, 0.01, 0.001),
      expand = c(0, 0)
    ) + 
    coord_cartesian(ylim = -c(2, -6)) + 
    annotate(geom = 'text', label = '0.01 bias') + 
    theme_void() + 
    theme(
      plot.title   = element_blank(),
      panel.spacing = unit(c(0, 0, 0, 0), 'cm'),
      plot.margin  = unit(c(0, 0, 0, 0), 'cm'),
      axis.text.y  = yaxis_text,
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    )
  
  if(.method_labels){
    bias_plot <- bias_plot +
      geom_text(
        data = summary_vals,
        aes(x = 1, y = -1, label = method_lab),
        color = color,
        size = 2)
  }
  
  coverage_plot <- ggplot(
    summary_vals,
    aes(x = 1, y = coverage)) +
    geom_bar(
      aes(x = 1, y = 1),
      stat = 'identity', width = .1,
      fill = NA,
      color = 'grey95',
      inherit.aes = FALSE
    ) + 
    geom_bar(
      stat = 'identity', 
      width = .1,
      fill = color
    ) +
    geom_hline(
      yintercept = 0.95,
      color      = 'grey50'
    ) + 
    geom_text(
      aes(label = round(coverage, 2), y = coverage - .025),
      size = 2,
      vjust = 1,
      color = 'grey25'
    ) +     
    # geom_text(
    #   aes(label = round(coverage, 2), y = coverage - .025),
    #   size = 2.01,
    #   vjust = 1,
    #   color = 'grey25'
    # ) + geom_text(
    #   aes(label = round(coverage, 2), y = coverage - .025),
    #   size = 2.01,
    #   vjust = 1,
    #   color = 'grey25'
    # ) +
    scale_x_continuous(
      name = '',
      expand = c(0, 0)
    ) + 
    scale_y_continuous(
      name ='',
      limits = c(0, 1),
      breaks = y_coverage,
      labels = y_coverage,
      expand = c(0, 0)
    ) + 
    theme_void() + 
    theme(
      plot.margin  = unit(c(0, 0, 0, 0), 'cm'),
      axis.text.y  = yaxis_text,
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    )
  
  
  p3 <- arrangeGrob(
    ggplotGrob(bias_plot),
    ggplotGrob(coverage_plot),
    ncol = 1,
    widths = unit(c(.3), 'in'),
    heights = unit(c(2, 1), 'in'),
    padding = unit(0, 'lines'))
  
  p3
  
}

grid.newpage()

plot_one(9, 'tcor_ocor', .method = 'ipw', .5, 0) %>% grid.draw()


### Level 2: 
plot_two <- function(.sid, .model_spec, .method_labels, .alpha, .a){
  arrangeGrob(
    grobs = lapply(c('ipw', 'otc', 'dbr', 'wls_dbr'), function(x) {
      plot_one(.sid         = .sid, 
               .model_spec  = .model_spec, 
               .method      = x, 
               .alpha       = .alpha,
               .a           = .a,
               .axis_labels = FALSE, 
               .method_labels = .method_labels)
      }),
    widths = unit(rep(.30, 4), 'in'),
    ncol = 4,
    padding = unit(0, 'lines'))
}

grid.newpage()
plot_two(9, 'tcor_ocor', TRUE, .5, 0) %>% grid.draw()

### Level 3
plot_three <- function(.sid, .model_spec, .method_labels, .alpha, .a, .spec_label){
  plots <- arrangeGrob(
    grobs = list(
      textGrob(label = .spec_label, 
               gp    = gpar(fontsize = 10),
               just = 'left',
               x = 0),
      plot_two(.sid, .model_spec, .method_labels, .alpha, .a)),
    nrow = 2,
    widths = unit(1.2, 'in'),
    heights = unit(c(.2, 3), 'in'))
  
  rect <- rectGrob(
    # height = unit(c(0.15, row_height*rows), 'in'), 
    # width  = unit(8/3 * cols, 'in'), 
    height = unit(1, 'npc'), 
    width  = unit(1, 'npc'), 
    gp     = gpar(lwd = 1, col = "grey80", fill = NA)) 
  
  gTree(children = gList(plots, rect))
}

### AXES
bias_plot_axes <- ggplot(
  data_frame(log10absbias = c(10, .000001)), 
  aes(x = 1, y =log10absbias )) +
  geom_blank() + 
  scale_x_continuous(
    name = '',
    expand = c(0,0)
  ) + 
  scale_y_continuous(
    # name   = expression('-log'[10]*'(|bias|)'),
    name   = '|Bias|',
    # name   = '',
    breaks = -c(0, -1, -2, -3),
    labels = c(1, 0.1, 0.01, 0.001),
    expand = c(0, 0)
  ) + 
  coord_cartesian(ylim = -c(2, -6)) + 
  theme_void() + 
  theme(
    plot.title   = element_blank(),
    panel.spacing = unit(c(0, 0, 0, 0), 'cm'),
    plot.margin  = unit(c(0, 0, 0, 0.9), 'cm'),
    axis.text.y  = element_text(size = 7, hjust = 1, angle = 0, debug = FALSE),
    axis.ticks.y = element_blank(),
    axis.title.y = element_text(size = 7, hjust = 0.5, angle = 90, debug = FALSE),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()
  )

coverage_plot_axes <- ggplot(
  data_frame(coverage = c(0, 1)),
  aes(x = 1, y = coverage)) +
  geom_blank() + 
  scale_x_continuous(
    name = '',
    expand = c(0, 0)
  ) + 
  scale_y_continuous(
    name ='Coverage',
    limits = c(0, 1),
    breaks = c(0, 0.95),
    labels = c(0, 0.95),
    expand = c(0, 0)
  ) + 
  theme_void() + 
  theme(
    plot.margin  = unit(c(0, 0, 0, 0.95), 'cm'),
    axis.text.y  = element_text(size = 7, hjust = 1, angle = 0, debug = FALSE),
    axis.ticks.y = element_blank(),
    axis.title.y = element_text(size = 7, hjust = .5, angle = 90, debug = FALSE),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()
  )

axes <- arrangeGrob(
  ggplotGrob(bias_plot_axes),
  ggplotGrob(coverage_plot_axes),
  ncol = 1,
  widths = unit(c(.75), 'in'),
  heights = unit(c(2, 1), 'in'),
  padding = unit(0, 'lines'))

axes_aligned <- arrangeGrob(
  grobs = list(
    textGrob(label = '', 
             just = 'left',
             x = 0),
    axes),
  nrow = 2,
  widths = unit(.75, 'in'),
  heights = unit(c(.2, 3), 'in'))

grid.newpage()
grid.draw(axes_aligned)
### Level 4: Put it all together

plot_four <- function(.sid, .alpha, .a){
  p1 <- plot_three(
    .sid           = .sid,
    .alpha         = .alpha,
    .a             = .a,
    .model_spec   = 'tcor_ocor', 
    .method_labels = TRUE, 
    .spec_label    = expression("(i) "*"f"*": true"*" "*"m"*": true"))
    # .spec_label    =  'T: correct\nO: correct')
  p2 <- plot_three(
    .sid           = .sid,
    .alpha         = .alpha,
    .a             = .a,
    .model_spec   = 'tmis_ocor',
    .method_labels = FALSE, 
    .spec_label    = expression("(ii) "*"f"*": false"*" "*"m"*": true"))
  p3 <- plot_three(
    .sid           = .sid,
    .alpha         = .alpha,
    .a             = .a,
    .model_spec   = 'tcor_omis', 
    .method_labels = FALSE, 
    .spec_label    = expression("(iii) "*"f"*": true"*" "*"m"*": false"))
  p4 <- plot_three(
    .sid           = .sid,
    .alpha         = .alpha,
    .a             = .a,
    .model_spec   = 'tmis_omis', 
    .method_labels = FALSE, 
    .spec_label    = expression("(iv) "*"f"*": false"*" "*"m"*": false"))
  arrangeGrob(
    axes_aligned, p1, p2, p3, p4,
    widths = unit(c(.8, rep(1.2 + .1, 4)), 'in'),
    ncol = 5)
}


### Print results


# y_1_5 <- arrangeGrob(
#   textGrob(label = expression('Y(1, 0.5)')),
#   plot_four(.sid = 9, .alpha = .5, .a  = 1),
#   nrow = 2,
#   heights = unit(c(.5, 3.5), 'in')
# ) 

y_1_5 <- arrangeGrob(
  plot_four(.sid = 9, .alpha = .5, .a  = 1),
  nrow = 1, heights = unit(3.5, 'in'))

grid.newpage()
grid.draw(y_1_5)

ggsave(filename = paste0('inst/experiments/ex31/figures/', 'y_1_5_', vers, '.pdf'),
       plot = y_1_5,
       width = 6.5, height = 3.6)

y_0_5 <- arrangeGrob(
  textGrob(label = expression('Y(0, 0.5)')),
  plot_four(.sid = 9, .alpha = .5, .a  = 0),
  nrow = 2,
  heights = unit(c(.5, 3.5), 'in')
) 

ggsave(filename = paste0('inst/experiments/ex31/figures/', 'y_0_5_', vers, '.pdf'),
       plot = y_0_5,
       width = 8.5, height = 5)




#------------------------------------------------------------------------------#
# data generating story ####
#------------------------------------------------------------------------------#

# Make plots to demonstrate example of scenario
make_example_plots <- function(...){
  dots <- list(...)
  args <- dots[pmatch(names(dots), names(formals((gen_data))))]
  sampledt    <- do.call(gen_data, args)
  
  p1 <- ggplot(data = sampledt, aes(x = p)) + 
    geom_histogram(bins = 30 ) + 
    scale_x_continuous(
      name = '',
      limits = c(0, 1)) +
    ggtitle('Distribution of Pr(A|X)') + 
    theme(
      plot.title = element_text(size = 9)
    )
  
  p2 <- ggplot(data = distinct(sampledt, group, fA), aes(x = fA)) + 
    geom_histogram(bins = 20) + 
    scale_x_continuous(
      name = '',
      limits = c(0, 1)) + 
    ggtitle('Distribution of f(A)') + 
    theme(
      plot.title = element_text(size = 9)
    )
  arrangeGrob(
    grobs = list(ggplotGrob(p1), ggplotGrob(p2)), 
    nrow = 2,
    heights = unit(c(1.5, 1.5), 'in'),
    widths = unit(c(2), 'in'))
}


plot_four_b <- function(.sid, .a){
  p1 <- plot_three(
    .sid           = .sid, 
    .model_spec    = 'tcor_ocor',
    .method_labels = TRUE, 
    .alpha         = .1, 
    .a             = .a, 
    .spec_label    = paste0('Y(', .a, ', 0.1)')) 
  p2 <- plot_three(
    .sid           = .sid, 
    .model_spec    = 'tcor_ocor',
    .method_labels = TRUE, 
    .alpha         = .5, 
    .a             = .a, 
    .spec_label    = paste0('Y(', .a, ', 0.5)')) 
  p3 <- plot_three(
    .sid           = .sid, 
    .model_spec    = 'tcor_ocor',
    .method_labels = TRUE, 
    .alpha         = .9, 
    .a             = .a, 
    .spec_label    = paste0('Y(', .a, ', 0.9)')) 
  arrangeGrob(
    axes_aligned, p1, p2, p3,
    widths = unit(c(.8, rep(.5*3 + .25, 3)), 'in'),
    ncol = 4)
}

plot_five_b <- function(.sid, .a){
  xx <- scenarios %>% filter(sid == .sid)
  yy <- make_example_plots(
    beta = xx$beta[[1]], 
    gamma = xx$gamma[[1]], 
    theta = xx$theta[[1]],
    m = xx$m, ni = xx$ni) 
  
  arrangeGrob(
    arrangeGrob(
      textGrob(label = 'Data Generating Scenario'),
      yy,
      nrow = 2,
      heights = unit(c(.5, 3), 'in')),
    plot_four_b(.sid = .sid, .a = .a),
    nrow = 1,
    widths = unit(c(2, sum(c(.8, rep(.5*3 + .25, 3)))), 'in')
  )
}
grid.newpage()

lapply(c(9), function(sid){
  p <- plot_five_b(sid, 1) 
  
  ggsave(filename = paste0('inst/experiments/ex31/figures/scn', sid, '.pdf'),
         plot = p,
         width = 8.5, height = 4)
})

