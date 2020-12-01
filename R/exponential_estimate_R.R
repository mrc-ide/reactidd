# Save this file as `R/exponential_estimate_R.R`

#' Bayesian exponential model using stan
#'
#' @export
#' @param fit_exponential a model fit using reacttemporal2::stan_exp_model() or reactemporal2::stan_exp_model_phe()
#' @param n_mean shape parameter for inverse gamma serial interval (default from Bi et al)
#' @param b_mean rate parameter for inverse gamma serial interval (default from Bi et al)
#' @return An object of class `dataframe` with calculations of growth rate, R, doubling/halving time with 95% CIs and the probability R>1
#'

exponential_estimate_R <- function(fit_exponential, n_mean = 2.29, b_mean = 0.36, label = "label"){
  ff <- rstan::extract(fit_exponential)

  growth_rate <- ff$beta

  #' 95% CI for r and median
  r<-quantile(growth_rate, c(0.025,0.5,0.975))
  prob_r_pos <- length(growth_rate[growth_rate>0])*100/ length(growth_rate)
  R <- (1+r/b_mean)**n_mean
  t_half <- log(2)/r

  return(data.frame(label = label,
                    growth_rate = r[[2]],
                    growth_rate_lb = r[[1]],
                    growth_rate_ub = r[[3]],
                    R = R[[2]],
                    R_lb = R[[1]],
                    R_ub = R[[3]],
                    t_half_double = t_half[[2]],
                    t_half_double_lb = t_half[[1]],
                    t_half_double_ub = t_half[[3]]))


}



