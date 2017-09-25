library(geex)
library(geepack, quietly = TRUE)
data("ohio")

test_binomial <- geeglm(resp ~ age, data = ohio, id = id,
                     family = binomial(link = "logit"))

test_eefun2 <- function(data, model){
  f <- grab_psiFUN(object = model, data = data)
  function(theta){
    f(theta)
  }
}

basis_t <- create_basis(
  estFUN = test_eefun2,
  data   = ohio,
  units  = 'id',
  outer_args = list(model = test_binomial)
)

basis_t@.GFUN(coef(test_binomial))
