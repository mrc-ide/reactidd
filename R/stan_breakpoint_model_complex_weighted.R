# Save this file as `R/stan_exp_model.R`

#' Bayesian exponential model using stan
#'
#' @export
#' @param X date vector.
#' @param Y Numeric vector of number of positive samples
#' @param N Numeric vector of total number of samples
#' @param restriction_dates vector of the dates of the key restriction changes
#' @param Vax1 proportion protected by single dose of vaccine only on days covering round 8-13
#' @param Vax2 proportion protected by second dose of vaccine on days covering round 8-13
#' @param Delta proportion of cases due to Delta variant on days covering round 8-13
#' @param ... Arguments passed to `rstan::sampling` (iter, warmup).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
stan_breakpoint_model_complex_weighted <- function(X, Y, N,restriction_dates,Vax1, Vax2, Delta, iter = 5000, warmup =1000, chains=4, cores = 1){

  X <- as.numeric(X)
  num_data <- length(X)

  Y <- as.numeric(Y)
  N <- as.numeric(N)
  min_X <- min(X)
  X <- X - min(X)
  max_X <- max(X)

  X_all <- seq(min(X), max(X))
  num_all <- length(X_all)
  indexes = as.integer(X+1)


  restriction_dates <- as.numeric(rest_dates) - min_X

  data_new = list(num_data = num_data,
                  num_all = num_all,
                  Y = Y,
                  N =N,
                  Vax1 = Vax1,
                  Vax2 = Vax2,
                  Delta = Delta,
                  X = X,
                  X_all = X_all,
                  max_X = max_X,
                  indexes = indexes,
                  restriction_dates = restriction_dates,
                  gen = 1)


  #' Load and run stan model
  rstan::rstan_options(auto_write = TRUE)
  options(mc.cores = cores)

  fit_spline <- rstan::sampling(stanmodels$breakpoint_complex,
                         iter=iter,
                         warmup =warmup,
                         chains=chains,
                         control = list(adapt_delta=0.95,
                                        max_treedepth = 10),
                         data = data_new)



  return(fit_spline)

}
