# Save this file as `R/stan_exp_model.R`

#' Bayesian exponential model using stan
#'
#' @export
#' @param X date vector.
#' @param Y Numeric vector of number of positive samples
#' @param N Numeric vector of total number of samples
#' @param ... Arguments passed to `rstan::sampling` (iter, warmup).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
stan_exp_model <- function(X, Y, N, iter = 5000, warmup =1000, cores = 1){

  X <- as.numeric(X)
  num_data <- length(X)

  Y <- as.integer(Y)
  N <- as.integer(N)

  X_adj <- min(X)
  X <- X - X_adj


  #' Load and run stan model
  rstan::rstan_options(auto_write = TRUE)
  options(mc.cores = cores)
  fit_spline <- rstan::sampling(stanmodels$linear,
                         iter=iter,
                         warmup = warmup,
                         chains=4,
                         control = list(adapt_delta=0.95,
                                        max_treedepth = 10),
                         data = list(num_data = num_data,
                                     Y = Y,
                                     N =N,
                                     X = X))


  return(fit_spline)

}
