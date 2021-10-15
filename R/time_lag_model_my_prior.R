# Save this file as `R/time_lag_model_my_prior.R

#' MCMC model fitting - return likelihood from prior
#'
#' @export
#' @param pars list of 2 parameter values for beta and alpha respectively
#' @return likelihood of parameter values from prior - Uniform from 0-40 for the time lag parameter

time_lag_model_my_prior <- function(pars){
  beta <- dunif(pars[1], 0, 40)
  sum(log(beta))
}
