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

shedding_simulate_data <- function(N_3tests=662, N_2tests=212,min_T=6, max_T=17,sens=0.79, k=0.071, tau=0){

  df <- data.frame()
  dif <- (max_T - min_T)/2

  for(i in seq_len(N_3tests)){

    # Random time for test to occur
    t1 <- round(runif(1,min_T,min_T+dif))
    t2 <- round(runif(1,min_T+dif, max_T))
    # True positive at those days?
    p1 <- rbinom(1, 1, min(exp(-k*(t1-tau)), 1) )
    p2 <- p1 * rbinom(1, 1, min(exp(-k*(t2-t1)), 1) )
    # Test positive at those days?
    e1 <- rbinom(1, 1, p1*sens)
    e2 <- rbinom(1, 1, p2*sens)

    row_df <- data.frame(passcode=i,
                         time_shed0 = 0,
                         time_shed1 = t1,
                         time_shed2 = t2,
                         estbinres_shed0 = 1,
                         estbinres_shed1 = e1,
                         estbinres_shed2 = e2)
    df <- rbind(df, row_df)
  }

  for(i in seq_len(N_2tests)){

    # Random time for test to occur
    t1 <- round(runif(1,min_T,max_T))
    # True positive at those days?
    p1 <- rbinom(1, 1, min(exp(-k*(t1-tau)), 1) )
    # Test positive at those days?
    e1 <- rbinom(1, 1, p1*sens)

    row_df <- data.frame(passcode=i+N_3tests,
                         time_shed0 = 0,
                         time_shed1 = t1,
                         time_shed2 = NA,
                         estbinres_shed0 = 1,
                         estbinres_shed1 = e1,
                         estbinres_shed2 = NA)
    df <- rbind(df, row_df)
  }


  return(df)
}
