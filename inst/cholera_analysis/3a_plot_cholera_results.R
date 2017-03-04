library(ggplot2)


ggplot(
  data = de,
  aes(x = alpha, 
      y = estimate,
      group = method,
      linetype = method)
) + 
  geom_line()

ggplot(
  data = ie,
  aes(x = alpha, 
      y = estimate,
      group = method,
      linetype = method)
) + 
  geom_line()



plot_effect <- function(var, ylim, xlab, ylab){
  with(target, {
    dbr <- target[estimator_type == 'dbr', ] 
    otc <- target[estimator_type == 'otc', ] 
    ipw <- target[estimator_type == 'ipw', ] 
    
    plot(dbr$alpha, dbr[[var]], 
         type = 'n', bty = 'l', 
         ylim = ylim,
         xlab = xlab,
         ylab = ylab,
         mgp  = c(2, .5, 0))
    lines(dbr$alpha, dbr[[var]])
    lines(otc$alpha, otc[[var]], lty = 2)
    lines(ipw$alpha, ipw[[var]], lty = 4)
  })
}

pdf(file = 'cholera_analysis/cholera_ce_diff_plots_V001.pdf', width = 7, height = 3)
par(mfrow=c(1,3))
plot_effect('de_diff', c(-2, 8), expression(alpha), expression(widehat("DE")*"("*alpha*")"))
legend(.45, 6.5, c('dbr', 'otc', 'ipw'),  lty = c(1, 2, 4), bty = 'n')
plot_effect('ie_diff', c(-2, 8), expression(alpha), expression(widehat("IE")*"("*0.4*","*alpha*"'"*")"))
plot_effect('te_diff', c(-2, 8), expression(alpha), expression(widehat("TE")*'('*0.4*","*alpha*"'"*")"))
dev.off()


pdf(file = 'cholera_analysis/cholera_ce_diff_plots_V001.pdf', width = 7, height = 3)
par(mfrow=c(1,3))
plot_effect('de_diff', c(-2, 8), expression(alpha), expression(widehat("DE")*"("*alpha*")"))
legend(.45, 6.5, c('dbr', 'otc', 'ipw'),  lty = c(1, 2, 4), bty = 'n')
plot_effect('ie_diff', c(-2, 8), expression(alpha), expression(widehat("IE")*"("*0.4*","*alpha*"'"*")"))
plot_effect('te_diff', c(-2, 8), expression(alpha), expression(widehat("TE")*'('*0.4*","*alpha*"'"*")"))
dev.off()

pdf(file = 'cholera_analysis/cholera_y_plots_V001.pdf', width = 7, height = 3)
par(mfrow=c(1,3))
plot_effect('y0', c(-2, 8), expression(alpha), expression(widehat("Y")*"(0,"*alpha*")"))
legend(.45, 8, c('dbr', 'otc', 'ipw'),  lty = c(1, 2, 4), bty = 'n')
plot_effect('y1', c(-2, 8), expression(alpha), expression(widehat("Y")*"(1,"*alpha*")"))
dev.off()






# test_glmer <- lme4::glmer(A ~ age + rivkm + (1|group), data = choleradt,
#                           family = binomial)
# test_geeglm <- geepack::geeglm(y_obs ~ A + fA + A*fA + age + rivkm, data = choleradt, 
#                                family = binomial, id = group)
# 
# test_glm <- glm(y_obs ~ A + fA + A*fA + age + rivkm, data = choleradt, 
#                                family = binomial)
