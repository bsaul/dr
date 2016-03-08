
get_args <- function(FUN, args_list = NULL, ...){
  dots <- append(args_list, list(...))
  arg_names <- names(formals(match.fun(FUN)))
  
  args <- dots[arg_names]
  args[sapply(args, is.null)] <- NULL
  
  return(args)
}

log_likelihood <- function(x, 
                           pos, 
                           integrand,
                           ...)
{
  ## Necessary pieces ##
  integrand <- match.fun(integrand)
  dots      <- list(...)
  dot.names <- names(dots)
  
  ## Integrate() arguments ##
  if(!'lower' %in% dot.names){
    dots$lower <- -Inf
  }
  
  if(!'upper' %in% dot.names){
    dots$upper <- Inf
  }
  
  int.args <- append(get_args(integrate, dots),
                     get_args(integrand, dots))
  args <- append(int.args, list(f = integrand, x = x, pos = pos))
  
  ## Calculuation ##
  attempt <- try(do.call(integrate, args = args))
  val <- ifelse(is(attempt, 'try-error'), NA, attempt$value)
  
  return(log(val))
}

score_calc <- function(integrand,
                       hide.errors = TRUE,
                       fixed.effects,
                       random.effects,
                       ...)
{
  ## Necessary bits ##
  params <- c(fixed.effects, random.effects)
  integrand <- match.fun(integrand)
  dots <- list(...)
  
  ## Function arguments ##
  int.args <- append(get_args(integrand, dots),
                     get_args(integrate, dots))
  fargs    <- append(int.args, get_args(numDeriv::grad, dots))
  
  ## Compute the derivative of the log likelihood for each parameter ##
  scores <- sapply(1:length(params), function(i){
    args <- append(fargs,
                   list(func = log_likelihood,
                        integrand = integrand, 
                        fixed.effects = fixed.effects,
                        random.effects = random.effects,
                        x = params[i], 
                        pos = i))
    
    attempt <- try(do.call(numDeriv::grad, args = args), silent = hide.errors)
    return(ifelse(is(attempt, 'try-error'), NA, attempt))
  })
  
  return(scores)
}

score_matrix <- function(integrand,
                         X, A, G, 
                         fixed.effects,
                         random.effects,
                         runSilent = F, #BB 2015-06-23 #Pass in from ipw_interference()
                         ...)
{
  ## Warnings ##
  if(length(fixed.effects) != ncol(X)){
    stop("The length of params is not equal to the number of predictors")
  }
  
  ## Necessary bits ##
  integrand <- match.fun(integrand)
  dots <- list(...)
  XX <- cbind(X, A)
  pp <- length(fixed.effects)
  gg <- sort(unique(G))
  
  ## Compute score for each group and parameter ##
  int.args <- append(get_args(integrand, dots),
                     get_args(integrate, dots))
  fargs <- append(int.args, get_args(numDeriv::grad, dots))
  
  if(runSilent != T){print("Calculating matrix of scores...")} #BB 2015-06-23
  s.list <- by(XX, INDICES = G, simplify = TRUE, 
               FUN = function(xx) {
                 args <- append(fargs, 
                                list(integrand = integrand, 
                                     fixed.effects = fixed.effects,
                                     random.effects = random.effects,
                                     A = xx[ , (pp + 1)],
                                     X = xx[ , 1:pp]))
                 return(do.call(score_calc, args = args))})
  
  ## Reshape list into matrix ##
  out <- matrix(unlist(s.list, use.names = FALSE), 
                ncol = pp + length(random.effects), 
                byrow = TRUE,
                dimnames = list(gg, names(c(fixed.effects, random.effects))))
  
  return(out)
}
