# Save this file as `R/stan_p_spline_phe.R`

#' Bayesian p-spline model using stan
#'
#' @export
#' @param X date vector.
#' @param Y Numeric vector of number of cases at date X
#' @param target_distance_between_knots sets the number of days between adjacent knots (default = 5)
#' @param spline_degree sets the degree of the splines (default = 3)
#' @param ... Arguments passed to `rstan::sampling` (iter, warmup, chains, cores).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'

stan_p_spline_phe <- function(X, Y, target_dist_between_knots = 5, spline_degree = 3,
                          iter = 20000, warmup = 1000, cores = 1, chains = 4){


  #' Convert date to numeric
  X <- as.numeric(X)
  min_date_numeric <- min(X)
  max_date_numeric <- max(X)

  num_knots <- ceiling((max_date_numeric- min_date_numeric)/target_dist_between_knots) + 7
  days_per_knot <- (max_date_numeric - min_date_numeric)/(num_knots - 7)
  num_basis <- num_knots + spline_degree - 1


  num_data <- length(X)
  knots <- unname(seq(min(X)-3*days_per_knot, max(X)+3*days_per_knot, length.out = num_knots))

  Y <- as.integer(Y)

  #' Load and run stan model
  rstan::rstan_options(auto_write = TRUE)
  options(mc.cores = cores)
  fit_spline <- rstan::sampling(stanmodels$b_splines_actual_phe,
                                iter=iter,
                                warmup =warmup,
                                chains=chains,
                                control = list(adapt_delta=0.95,
                                               max_treedepth = 10),
                                data = list(num_data = num_data,
                                            num_knots = num_knots,
                                            knots = knots,
                                            Y = Y,
                                            X = X,
                                            spline_degree = spline_degree))



  return(fit_spline)
}
