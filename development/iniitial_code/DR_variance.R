#------------------------------------------------------------------------------#
#   Title: Doubly Robust variance estimation
#  Author: B. Saul
#    Date: 2016-03-09
# Purpose: variance estimation for Lan's DR
#------------------------------------------------------------------------------#

library(geepack)
library(sandwich)


test_sim <- dr_ipw(formula_outcome = Y ~ Z1 + Z2 + Z3 + Z4 + A + fA + A*Z1 + fA*Z2,
                   formula_interference = Y | A ~ -1 + Z1 + Z2 + Z3 + Z4 +(1|group) | group,
                   method_outcome  = geepack::geeglm,
                   method_opts_outcome = list(id = quote(group), family = gaussian),
                   method_opts_interference = list(allocations = alphas,
                                                   method = 'simple',
                                                   runSilent = T),
                   dr_term2_function = dr_term2,
                   data = filter(DRsims, simID == 1))

ipw_scores <- test_sim$ipw$scores
estfun_outcome <- estfun(test_sim$outcome)

ee <- cbind(ipw_scores, estfun_outcome)
test_sim$ipw$point_estimates
V <- V_matrix(ee, test_sim$ipw$point_estimates, allocation1 = .1, trt.lvl1 = 1,
                       effect_type = 'outcome', marginal = F)
vdim <- dim(V)[1]
V21 <- V[vdim, 1:(vdim - 1)] # Last row, up to last column
V11 <- V[1:(vdim - 1), 1:(vdim - 1)] # up to last row, up to last column
V22 <- V[vdim, vdim] # bottom right element

U21 <- (t(as.matrix(apply(-U_pe_grp, 2, sum, na.rm = T))))/N

((U21 - 2*V21) %*% solve(V11) %*% t(U21) + V22)/N

#   ee_outcome <- estfun(model_outcome)
#   ee_ipw <- model_ipw$scores
#   # Stack EEs
#   ee <- cbind(ee_outcome, ee_ipw)
