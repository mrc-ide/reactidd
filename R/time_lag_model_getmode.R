# Save this file as `R/time_lag_model_getmode.R

#' MCMC model fitting - get mode of time lag vaector
#'
#' @export
#' @param v vector of values
#' @return mode of values


time_lag_model_getmode <- function(v){
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v,uniqv)))]
}
