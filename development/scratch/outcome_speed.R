
m_A <- lme4::glmer(tru_propen, family = binomial, data = DRsims)
m_Y <- geepack::geeglm(tru_outcome, family = gaussian, data = DRsims, id = ID)
theta <- coef(m_Y)



xmat_outcome <- as.data.frame(model.matrix(m_Y))
treatment      <- model.response(model.frame(m_Y))
rhs_outcome <- get_fixed_formula(m_Y)

frame <- list(outcome   = DRsims$Y,
              X_outcome = xmat_outcome) %>%
  lapply(., function(x) split(x, DRsims$group, drop = FALSE))

split_frame <- Map(list,
                   outcome      = frame[['outcome']],
                   X_outcome  = frame[['X_outcome']])

#### new ####
expand_funcs <- lapply(split_frame, function(x){
  expand_outcome_frame2(X_outcome = x$X_outcome, rhs_outcome)
})

system.time({
test_frame <-  lapply(expand_funcs, function(f) 
  {
    f(theta) 
  } ) } )
test_frame[[1000]]
as.matrix(get_design_frame(rhs_outcome, test_frame[[1000]])) %*% theta

outcome_funcs <- lapply(split_frame, function(x){
  make_otc_estimator(X_outcome = x$X_outcome, rhs_outcome) 
})

system.time({test1 <- lapply(outcome_funcs, function(f) f(theta, alpha =.5, a = 1))})
system.time({test2 <- lapply(outcome_funcs, function(f) f(theta, alpha =.5, a = 1))})
system.time({test3 <- lapply(outcome_funcs, function(f) f(theta, alpha =.5, a = 0))})
system.time({test4 <- lapply(outcome_funcs, function(f) f(theta, alpha =.25, a = 0))})



outcome_funcs2 <- lapply(split_frame, function(x){
  make_otc_estimator2(X_outcome = x$X_outcome, rhs_outcome) 
})

system.time({test1 <- lapply(outcome_funcs2, function(f) f(theta, alpha =.5, a = NULL))})
system.time({test2 <- lapply(outcome_funcs2, function(f) f(theta, alpha =.5, a = 1))})
system.time({test3 <- lapply(outcome_funcs2, function(f) f(theta, alpha =.5, a = 0))})
system.time({test4 <- lapply(outcome_funcs2, function(f) f(theta, alpha =.25, a = 0))})

#### old #### 
outcome_funcs_old <- lapply(split_frame, function(x){
  make_outcome_estimator(X_outcome = x$X_outcome, formula(m_Y), theta = theta)
})

system.time({test1_old <- lapply(outcome_funcs_old, function(f) f(alpha =.5, a = NULL))})
system.time({test <- lapply(outcome_funcs_old, function(f) f(alpha =.1, a = 1))})

# %>%
#   group_by(A) %>%
#       summarize_(mean = ~mean(term2)) %>%
#       summarize_(mean = ~mean(mean)) 
# 
# 
#       as.numeric()
# %>% mutate_(term2 = ~ fitted_value * pi_value) %>%
#   group_by_(~ID) %>% summarise(term2 = mean(term2)) %>% summarise(mean(term2))
# test1_old[[1000]] %>% arrange(sum_a)

# expand.grid.alt <- function(seq1,seq2) {
#   cbind(rep.int(seq1, length(seq2)),
#         c(t(matrix(rep.int(seq2, length(seq1)), nrow=length(seq2)))))
# }

