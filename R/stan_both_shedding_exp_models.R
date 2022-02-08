# Save this file as `R/stan_p_spline.R`

#' Fit exponetial decay models for shedding data
#'
#' @export
#' @param dat list of input data created by `reactidd::shedding_prepare_data_for_binary_model`
#' @param ... Arguments passed to `rstan::sampling` (iter, warmup).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'

stan_both_shedding_exp_models <- function(dat, iter = 20000, warmup = 1000, cores = 1, chains = 4){

  #' Load and run stan model
  rstan::rstan_options(auto_write = TRUE)
  options(mc.cores = cores)
  fit_spline1 <- rstan::sampling(stanmodels$shedding_exp,
                                 iter=iter,
                                 warmup =warmup,
                                 chains=chains,
                                 control = list(adapt_delta=0.95,
                                                max_treedepth = 10),
                                 data = dat)
  fit_spline2 <- rstan::sampling(stanmodels$shedding_exp_tdelay,
                                 iter=iter,
                                 warmup =warmup,
                                 chains=chains,
                                 control = list(adapt_delta=0.95,
                                                max_treedepth = 10),
                                 data = dat)



  return(list(fit_spline1, fit_spline2))
}
