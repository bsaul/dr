#------------------------------------------------------------------------------#
#' Get the formula from a model object
#' @export
#------------------------------------------------------------------------------#
get_fixed_formula <- function(model_object){
  formula(model_object, fixed.only = TRUE)[-2]
}

#------------------------------------------------------------------------------#
#' Get a matrix of fixed effects from a model object
#' @export
#------------------------------------------------------------------------------#
get_design_frame <- function(rhs_formula, data){
  as.data.frame(model.matrix(rhs_formula, data))
}

#------------------------------------------------------------------------------#
#' Get a vector of responses from a model object
#' @export
#------------------------------------------------------------------------------#

get_response <- function(formula, data){
  model.response(model.frame(formula, data = data))
}

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