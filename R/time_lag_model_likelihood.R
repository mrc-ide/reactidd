# Save this file as `R/time_lag_model_likelihood.R

#' MCMC model fitting - return likelihood for current parameter values
#'
#' @export
#' @param ensemble
#' @param dates
#' @param ENG
#' @param beta time lag parameter
#' @param alpha scaling factor parameter
#' @return Average log-likelihood value over spline posterior specific values of beta and alpha

time_lag_model_likelihood <- function(ensemble,dates,ENG, beta, alpha){
  time_lag <- round(beta)
  day = as.numeric(dates) - time_lag
  log_likelihood <- 0

  hosp_df <- data.frame(day = as.numeric(dates)-time_lag,
                        hosp =seq(1,ncol(ENG)))

  temp_ensemble <- merge(ensemble, hosp_df)
  ENG_temp <- ENG[,temp_ensemble$hosp]
  ENG_temp <- ENG_temp * exp(alpha)
  log_likelihood <- temp_ensemble$pos*log(t(ENG_temp)) +
    (temp_ensemble$obs-temp_ensemble$pos)*log(1-t(ENG_temp))

  log_likelihood<- sum(log_likelihood)/nrow(ENG)
  return(log_likelihood)

}
