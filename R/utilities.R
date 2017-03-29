#' #------------------------------------------------------------------------------#
#' #' Get the formula from a model object
#' #' @export
#' #------------------------------------------------------------------------------#
#' get_fixed_formula <- function(model_object){
#'   formula(model_object, fixed.only = TRUE)[-2]
#' }
#' 
#' #------------------------------------------------------------------------------#
#' #' Get a matrix of fixed effects from a model object
#' #' @export
#' #------------------------------------------------------------------------------#
#' get_design_frame <- function(rhs_formula, data){
#'   as.data.frame(model.matrix(rhs_formula, data))
#' }
#' 
#' #------------------------------------------------------------------------------#
#' #' Get a vector of responses from a model object
#' #' @export
#' #------------------------------------------------------------------------------#
#' 
#' get_response <- function(formula, data){
#'   model.response(model.frame(formula, data = data))
#' }

#------------------------------------------------------------------------------#
#' Converts a matrix to a list
#' @export
#------------------------------------------------------------------------------#

list_matrix <- function(this_list)
{
  ulist <- unlist(this_list)
  m <- length(this_list)
  p <- length(ulist)/m
  matrix(ulist, nrow = m, ncol = p, byrow = T)
}

#' #'-----------------------------------------------------------------------------#
#' #' Fills a matrix with zeros
#' #' @export
#' #'-----------------------------------------------------------------------------#
#' 
#' rbind_fill_zero <- function(this_list)
#' {
#'   if(!is.list(this_list)){
#'     this_list <- list(this_list)
#'   }
#'   p <- max(unlist(lapply(this_list, ncol)))
#'   out <- lapply(this_list, function(x){
#'     xp <- ncol(x)
#'     if(xp < p ){
#'       cbind(x,matrix(0, nrow = nrow(x), ncol = p - xp))
#'     } else {
#'       x
#'     }
#'   })
#'   do.call('rbind', out)
#' }

#------------------------------------------------------------------------------#
#' Evaluate a function on a dataframe
#' 
#' For use in dplyr
#' 
#' @export
#------------------------------------------------------------------------------#

evaluate_df_function <- function(f, ...){
  f(...)
}

#------------------------------------------------------------------------------#
#' Capture the model warnings from an object
#' @export
#------------------------------------------------------------------------------#

model_warnings <- function(object){
  if(class(object) %in% 'merMod'){
    paste(object@optinfo$conv$lme4$messages, collapse = ' ')
  } else {
    paste(warnings(objects), collapse = ' ')
  }
}

#------------------------------------------------------------------------------#
#' Creates a list of models from model_args
#' @export
#------------------------------------------------------------------------------#
make_models <- function(model_args, data)
{
  out <- lapply(model_args, function(x){
    method <- match.fun(x$method)
    if(is.null(x$user)) x$user <- FALSE
    
    if(x$user == TRUE){
      NULL # For now
    } else {
      args <- append(x$options, list(formula = x$formula, data = data) )
      do.call(method, args = args)
    }
  })
 out
}

#------------------------------------------------------------------------------#
#' Extract relevant information from list of models
#' @export
#------------------------------------------------------------------------------#

extract_model_info <- function(models, data, estimator_type, regression_type = 'none'){

  stopifnot(estimator_type %in% c('ipw', 'otc', 'dbr', 'reg_dbr'))
  t_model <- models$t_model
  stopifnot(class(t_model) == 'glmerMod') # currently designed to only work with merMod objects with random intercept
  o_model <- models$o_model
  if(estimator_type == 'reg_dbr' & regression_type == 'wls'){
    model_0 <- models$wls_model_0
    model_1 <- models$wls_model_1
  }
  
  if(estimator_type == 'reg_dbr' & regression_type == 'pcov'){
    model_0 <- models$pcov_model_0
    model_1 <- models$pcov_model_1
  }
  
  ## component data
  out <- list()
  if(estimator_type %in% c('ipw', 'dbr', 'reg_dbr')){
    out$X_t <- geex::get_design_matrix(geex::get_fixed_formula(t_model), data = data)
    out$Y   <- geex::get_response(formula(o_model), data = data)
    out$A   <- geex::get_response(A ~ 1, data = data)
    out$N   <- length(out$Y)  
    out$p_t <- ncol(out$X_t) + 1
    out$p   <- out$p_t
  } 
  if (estimator_type %in% c('otc', 'dbr')){
    out$X_o    <- as.data.frame(geex::get_design_matrix(geex::get_fixed_formula(o_model), data = data))
    out$X_o_ex <- expand_outcome_frame(out$X_o, geex::get_fixed_formula(o_model))
    out$N      <- nrow(out$X_o)
    out$p_o    <- ncol(out$X_o)
    out$p      <- out$p_o
    out$rhs_o  <- geex::get_fixed_formula(o_model)
    out$inv_link_o <- family(o_model)$linkinv
  } 
  if (estimator_type == 'reg_dbr'){
    out$X_o_reg_1 <- as.data.frame(geex::get_design_matrix(geex::get_fixed_formula(model_1), data = data))
    out$X_o_reg_0 <- as.data.frame(geex::get_design_matrix(geex::get_fixed_formula(model_0), data = data))
    out$rhs_o_0   <- update.formula(geex::get_fixed_formula(model_0), ~ A + .)
    out$rhs_o_1   <- update.formula(geex::get_fixed_formula(model_1), ~ A + .)
    out$X_o    <- as.data.frame(geex::get_design_matrix(out$rhs_o_0, data = data))
    out$X_o_ex <- expand_outcome_frame(out$X_o, out$rhs_o)
    out$N      <- nrow(out$X_o)
    out$p_o_1  <- ncol(out$X_o_reg_1)
    out$p_o_0  <- ncol(out$X_o_reg_0)
    out$p_o    <- out$p_o_1 + out$p_o_0
    out$p      <- out$p_t + out$p_o
    out$rhs_o_reg_1  <- geex::get_fixed_formula(model_1)
    out$rhs_o_reg_0  <- geex::get_fixed_formula(model_0)
    out$inv_link_o <- family(model_0)$linkinv
  } 
  if (estimator_type == 'dbr'){
    out$p <- out$p_t + out$p_o
  }
  out
}
