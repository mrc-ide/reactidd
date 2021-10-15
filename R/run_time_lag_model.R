# Save this file as `R/run_time_lag_model.R

#' MCMC model fitting a time lag and multiplicative factor to death/hosp data to REACT data
#'
#' @export
#' @param N number of iterations.
#' @param init_pars Set initial value of parameters
#' @param ensemble
#' @param dates
#' @param ENG
#' @param burnin number of samples to remove from the MCMC posterior as burnin
#' @param name label of the model fit
#' @param rounds label descirbing the number of rounds fit to
#' @param data_source label descirbing the public data set fit to
#' @return An object of class `list` with the first item being a data frame with best fit parameters/95% CI and second item being the whole posterior of the MCMC

run_time_lag_model <- function(N, init_pars, ensemble, dates, ENG ,burnin, name,
                      rounds ="Rounds 1-12", data_source = "Deaths"){
  chain_pre_burn<-reactidd::time_lag_model_runMCMC(init_pars, N, ensemble, dates, ENG)
  chain <- chain_pre_burn[burnin:nrow(chain_pre_burn),]
  df<- data.frame(name =name,
                  rounds = rounds,
                  data_source = data_source,
                  alpha = quantile(chain$alpha, c(0.5, 0.025,0.975))[1][[1]],
                  alpha_lb = quantile(chain$alpha, c(0.5, 0.025,0.975))[2][[1]],
                  alpha_ub = quantile(chain$alpha, c(0.5, 0.025,0.975))[3][[1]],
                  beta = time_lag_model_getmode(round(chain$beta)),
                  beta_lb = quantile(round(chain$beta), c(0.5, 0.025,0.975))[2][[1]],
                  beta_ub = quantile(round(chain$beta), c(0.5, 0.025,0.975))[3][[1]])
  list(df, chain)
}
