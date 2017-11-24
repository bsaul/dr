# Trying to debug/understand what's going on with WLS estimator

library(dplyr)
library(dr)
library(geex)
alphas <- lapply(seq(.3, .4, by = .02), function(x) sort(c(.4, x)))

load( pipe( 'ssh saulb@diamond.bios.unc.edu "cat /home/groups/projects/mhudgens/emch/data/R_data/emch_analysis_data.Rdata"' ))
source('inst/cholera_analysis/1a_define_analysis_functions.R')

choleradt <- analysis_c %>%
  group_by(group) %>%
  mutate(fA = mean(A)) %>%
  filter(n() < 100) %>%
  ungroup()

analysis_model_args1 <- list(
  t_model = 
    list(method = lme4::glmer,
         formula = B ~ age + rivkm + (1|group),
         options = list(family = binomial(link = 'logit'))),
  o_model =
    list(method  = geepack::geeglm,
         formula = y_obs ~ A + fA + age + rivkm,
         options = list(
           family  = binomial(link = 'logit'),
           id      = quote(group))),
  wls_model_0 =
    list(method = stats::glm,
         formula = as.integer(y_obs) ~ fA + age + rivkm,
         options = list(family = quasibinomial(link = 'logit'))),
  wls_model_1 =
    list(method = stats::glm,
         formula = as.integer(y_obs) ~ fA + age + rivkm,
         options = list(family = quasibinomial(link = 'logit')))
)

models0 <- estimate_cholera_parms_step0(
  data        = choleradt,
  model_args  = analysis_model_args1
)

results1 <- lapply(seq_along(alphas), function(i){
  estimate_cholera_parms_step2(
    data        = choleradt,
    allocations = alphas[i],
    models      = models0,
    model_args  = analysis_model_args1,
    randomization  = 2/3,
    compute_se  = TRUE
  )
})


# save(results1, file = paste0('inst/cholera_analysis/wls_testing1_', Sys.Date(),'.rda'))
# save(results2, file = paste0('inst/cholera_analysis/wls_testing2_', Sys.Date(),'.rda'))

lapply(seq_along(results2), function(i) {
  results2[[i]][[1]]$wls_dbr <<- results1[[i]][[1]]$wls_dbr
})

methods <- names(results1[[1]][[1]])

munger <- function(results){
  lapply(seq_along(results), function(i){
    x <- results[[i]]
    lapply(seq_along(x), function(j){
      xx <- x[[j]]
      lapply(seq_along(xx), function(m){
        xxx <- xx[[m]]
        
        if(xxx$alpha[1] < 0.4){
          alpha1  <- xxx$alpha[2]
          alpha2  <- xxx$alpha[1]
          C <- matrix(
            c(1, 0, 0, 0,   # Y(0, alpha)
              0, 1, 0, 0, 
              0, 0, 1, 0,   # Y(1, alpha)
              0, 0, 0, 1,
              1,  0, -1, 0, # Y(0, alpha') - Y(1, alpha')
              -1, 1, 0, 0,  # Y(0, .4) - Y(0, alpha')
              0, 1, -1, 0, # Y(0, .4) - Y(1, alpha')
              (1 - alpha1), -(1- alpha2), alpha1, -alpha2),
            byrow = TRUE, ncol = 4
          )
          
          
        } else {
          alpha1  <- xxx$alpha[1]
          alpha2  <- xxx$alpha[2]
          
          C <- matrix(
            c(1, 0, 0, 0,   # Y(0, alpha1)
              0, 1, 0, 0,   # Y(0, alpha2)
              0, 0, 1, 0,   # Y(1, alpha1)          
              0, 0, 0, 1,   # Y(1, alpha2)
              0, 1, 0, -1,
              1, -1, 0, 0,
              1, 0, 0, -1,
              (1 - alpha1), -(1- alpha2), alpha1, -alpha2),
            byrow = TRUE, ncol = 4
          )
          
        }
        
        print(xxx$estimate)
        est <- as.numeric(C %*% xxx$estimate * 1000)
        
        p <- ncol(xxx$vcov)
        Cvcov <- cbind(matrix(0, nrow = nrow(C), ncol = p- 4), C)
        print(sqrt(diag(Cvcov %*% xxx$vcov %*% t(Cvcov))))
        std_err <- as.numeric(sqrt(diag(Cvcov %*% xxx$vcov %*% t(Cvcov))) * 1000)
        # std_err <- NA
        
        data_frame(
          method    = methods[m],
          effect    = factor(c('Y0_1', 'Y0_2', 'Y1_1', 'Y1_2', 'de', 'ie', 'te', 'oe'),
                             levels = c('Y0_1', 'Y0_2', 'Y1_1', 'Y1_2', 'de', 'ie', 'te', 'oe'), ordered = TRUE),
          alpha1    = alpha1,
          alpha2    = alpha2,
          estimate  = est,
          std_err   = std_err,
          conf_low  = est - 1.96 * std_err,
          conf_high = est + 1.96 * std_err)
      })  -> hold 
      bind_rows(hold)
    }) %>%
      bind_rows
  }) %>%
    bind_rows
}

dt1 <- munger(results1)
dt2 <- munger(results2)
#### Plots ####
library(ggplot2)

## Color values
color_vals <- c("ipw" = "#EFC583",
                "otc" = rgb(86, 180, 233, max = 255),
                "dbr" = rgb(204, 121, 167, max = 255),
                "wls_dbr" = rgb(213, 94, 0, max = 255))

plotter <- function(dt){
  p0 <- ggplot(
    data = dt %>% filter(effect %in% c('Y0_1', 'Y1_1')),
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
}

plotter(dt1)
plotter(dt2)

## Look at A matrices