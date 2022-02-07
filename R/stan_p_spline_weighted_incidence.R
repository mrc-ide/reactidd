# Save this file as `R/stan_p_spline_weighted.R`

#' Bayesian p-spline model using stan
#'
#' @export
#' @param X date vector.
#' @param Y Numeric vector of weighted number of positive samples
#' @param N Numeric vector of weighted total number of samples
#' @param target_distance_between_knots sets the number of days between adjacent knots (default = 5)
#' @param spline_degree sets the degree of the splines (default = 3)
#' @param ... Arguments passed to `rstan::sampling` (iter, warmup).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'

stan_p_spline_weighted_incidence <- function(X, Y, N, target_dist_between_knots = 5, spline_degree = 3,
                                   k=0.07126, sens=0.79, iter = 20000, warmup = 1000, cores = 1, chains = 4){


  #' Convert date to numeric
  X <- as.numeric(X)
  X_data <- X
  min_date_numeric <- min(X)
  max_date_numeric <- max(X)

  incidence<- data.frame(day=seq(min_date_numeric-50, max_date_numeric))
  min_date_numeric <- min(incidence$day)
  max_date_numeric <- max(incidence$day)

  #' Set initial conditions and variable names to pass to stan
  # Knots are now based on incidence
  num_knots <- ceiling((max_date_numeric- min_date_numeric)/target_dist_between_knots)+7
  days_per_knot <- (max_date_numeric - min_date_numeric)/(num_knots -7)

  num_basis <- num_knots + spline_degree - 1

  X <- incidence$day
  num_days <- length(X)
  knots <- unname(seq(min(X)-3*days_per_knot, max(X)+3*days_per_knot, length.out = num_knots))

  Y <- as.numeric(Y)
  N <- as.numeric(N)
  num_data <- length(X_data)


  conv_matrix <- matrix(data=NA, nrow=num_days, ncol = num_data)


  for(i in seq_len(num_data)){
    for(j in seq_len(num_days)){
      t_dif <- X_data[i] - X[j]
      if(t_dif>49){
        conv_matrix[j,i]<-0
      } else if(t_dif<0){
        conv_matrix[j,i]<-0
      } else{
        conv_matrix[j,i] <- exp(-k*t_dif)
      }
    }

  }
  conv_matrix <- conv_matrix*(1/k)/sum(conv_matrix[,1])


  #' Load and run stan model
  rstan::rstan_options(auto_write = TRUE)
  options(mc.cores = cores)
  fit_spline <- rstan::sampling(stanmodels$b_splines_actual_weighted_incidence,
                                iter=iter,
                                warmup =warmup,
                                chains=chains,
                                control = list(adapt_delta=0.95,
                                               max_treedepth = 10),
                                data = list(num_data = num_data,
                                            num_days = num_days,
                                            num_knots = num_knots,
                                            knots = knots,
                                            Y = Y,
                                            N = N,
                                            X = X,
                                            X_data = X_data,
                                            spline_degree = spline_degree,
                                            sens = sens,
                                            conversion_matrix = conv_matrix))



  return(fit_spline)
}
