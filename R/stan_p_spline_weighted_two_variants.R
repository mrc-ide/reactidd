# Save this file as `R/stan_p_spline_weighted.R`

#' Bayesian p-spline model using stan
#'
#' @export
#' @param X date vector.
#' @param Y Numeric vector of weighted number of positive samples
#' @param N Numeric vector of weighted total number of samples
#' @param V data.frame with columns of date, variant 1, variant 2
#' @param target_distance_between_knots sets the number of days between adjacent knots (default = 5)
#' @param spline_degree sets the degree of the splines (default = 3)
#' @param ... Arguments passed to `rstan::sampling` (iter, warmup).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'

stan_p_spline_weighted_two_variants <- function(X, Y, N, V, target_dist_between_knots = 5, spline_degree = 3,
                                   iter = 20000, warmup = 1000, cores = 1, chains = 4){



  V1 <- c()
  V2 <- c()

  for(i in seq_len(length(X))){

    if(X[i] %in% V[,1]){
      index <- which(X[i] == V[,1])
      V1 <- c(V1,V[index,2])
      V2 <- c(V2,V[index,3])

    } else{

      V1 <- c(V1,as.numeric(0))
      V2 <- c(V2,as.numeric(0))

    }
  }
  #' Convert date to numeric
  #'
  #'
  X <- as.numeric(X)
  min_date_numeric <- min(X)
  max_date_numeric <- max(X)

  num_knots <- ceiling((max_date_numeric- min_date_numeric)/target_dist_between_knots) + 7
  days_per_knot <- (max_date_numeric - min_date_numeric)/(num_knots - 7)
  num_basis <- num_knots + spline_degree - 1


  num_data <- length(X)
  knots <- unname(seq(min(X)-3*days_per_knot, max(X)+3*days_per_knot, length.out = num_knots))

  Y <- as.numeric(Y)
  N <- as.numeric(N)

  YL <- as.numeric(V1)  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  NL <- as.numeric(V1) + as.numeric(V2)

  #' Load and run stan model
  rstan::rstan_options(auto_write = TRUE)
  options(mc.cores = cores)
  fit_spline <- rstan::sampling(stanmodels$b_splines_actual_weighted_two_variants,
                                iter=iter,
                                warmup =warmup,
                                chains=chains,
                                control = list(adapt_delta=0.95,
                                               max_treedepth = 10),
                                data = list(num_data = num_data,
                                            num_knots = num_knots,
                                            knots = knots,
                                            Y = Y,
                                            N =N,
                                            YL = YL,
                                            NL = NL,
                                            X = X,
                                            spline_degree = spline_degree))



  return(fit_spline)
}
