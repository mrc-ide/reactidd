# Save this file as `R/stan_p_spline.R`

#' Bayesian p-spline model using stan
#'
#' @export
#' @param X date vector.
#' @param Y Numeric vector of number of positive samples
#' @param N Numeric vector of total number of samples
#' @param target_distance_between_knots sets the number of days between adjacent knots (default = 5)
#' @param spline_degree sets the degree of the splines (default = 3)
#' @param ... Arguments passed to `rstan::sampling` (iter, warmup).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'


shedding_get_params <- function(mod, label){
  ff <- rstan::extract(mod)
  sens <- quantile(ff$sens, c(0.5,0.025,0.975))
  k_par <- quantile(ff$k, c(0.5,0.025,0.975))
  decay_mean <- 1/k_par
  decay_median <- log(2)/k_par
  td<- quantile(ff$t_d, c(0.5,0.025,0.975))
  pos_mean <- quantile((1/ff$k)+ff$t_d, c(0.5, 0.025, 0.975))
  pos_median <- quantile((log(2)/ff$k)+ff$t_d, c(0.5, 0.025, 0.975))
  N= ncol(ff$log_lik)
  row <- data.frame(label=label,
                    N=N,
                    k = k_par[1],
                    k_lb = k_par[2],
                    k_ub = k_par[3],
                    decay_mean = decay_mean[1],
                    decay_mean_lb = decay_mean[2],
                    decay_mean_ub = decay_mean[3],
                    decay_median = decay_median[1],
                    decay_median_lb = decay_median[2],
                    decay_median_ub = decay_median[3],
                    sens = sens[1],
                    sens_lb = sens[2],
                    sens_ub = sens[3],
                    td = td[1],
                    td_lb = td[2],
                    td_ub = td[3],
                    pos_mean = pos_mean[1],
                    pos_mean_lb = pos_mean[2],
                    pos_mean_ub = pos_mean[3],
                    pos_median = pos_median[1],
                    pos_median_lb = pos_median[2],
                    pos_median_ub = pos_median[3])
  return(row)

}

