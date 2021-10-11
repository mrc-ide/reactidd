# Save this file as `R/time_lag_model_param_new.R

#' MCMC model fitting - return new parameters
#'
#' @export
#' @param pars two item list of values for beta (time lag) and alpha (scaling factor)
#' @return New parameter values

time_lag_model_param_new <- function(pars){
  pars_new <- rnorm(2, mean = pars, sd = c(1,0.05))
  pars_new
}
