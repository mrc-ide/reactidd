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


shedding_prepare_data_for_binary_model <- function(shed_dat){
  IDS_missing <- shed_dat[is.na(shed_dat$estbinres_shed2)==TRUE,]$passcode

  time1_11 <- shed_dat[shed_dat$passcode %in% IDS_missing& shed_dat$estbinres_shed1==1,]$time_shed1
  time1_10 <- shed_dat[shed_dat$passcode %in% IDS_missing& shed_dat$estbinres_shed1==0,]$time_shed1

  time1_111 <- shed_dat[!(shed_dat$passcode %in% IDS_missing)& shed_dat$estbinres_shed1==1& shed_dat$estbinres_shed2==1,]$time_shed1
  time1_110 <- shed_dat[!(shed_dat$passcode %in% IDS_missing)& shed_dat$estbinres_shed1==1& shed_dat$estbinres_shed2==0,]$time_shed1
  time1_100 <- shed_dat[!(shed_dat$passcode %in% IDS_missing)& shed_dat$estbinres_shed1==0& shed_dat$estbinres_shed2==0,]$time_shed1
  time1_101 <- shed_dat[!(shed_dat$passcode %in% IDS_missing)& shed_dat$estbinres_shed1==0& shed_dat$estbinres_shed2==1,]$time_shed1
  time2_111 <- shed_dat[!(shed_dat$passcode %in% IDS_missing)& shed_dat$estbinres_shed1==1& shed_dat$estbinres_shed2==1,]$time_shed2
  time2_110 <- shed_dat[!(shed_dat$passcode %in% IDS_missing)& shed_dat$estbinres_shed1==1& shed_dat$estbinres_shed2==0,]$time_shed2
  time2_100 <- shed_dat[!(shed_dat$passcode %in% IDS_missing)& shed_dat$estbinres_shed1==0& shed_dat$estbinres_shed2==0,]$time_shed2
  time2_101 <- shed_dat[!(shed_dat$passcode %in% IDS_missing)& shed_dat$estbinres_shed1==0& shed_dat$estbinres_shed2==1,]$time_shed2

  N_11 <- length(time1_11)
  N_10 <- length(time1_10)
  N_111 <- length(time1_111)
  N_101 <- length(time1_101)
  N_110 <- length(time1_110)
  N_100 <- length(time1_100)

  data_for_binary_model<-list(N_11 = N_11,
                              time1_11 = time1_11,
                              N_10 = N_10,
                              time1_10 = time1_10,
                              N_111 = N_111,
                              time1_111 = time1_111,
                              time2_111 = time2_111,
                              N_101 = N_101,
                              time1_101 = time1_101,
                              time2_101 = time2_101,
                              N_110 = N_110,
                              time1_110 = time1_110,
                              time2_110 = time2_110,
                              N_100 = N_100,
                              time1_100 = time1_100,
                              time2_100 = time2_100)

  return(data_for_binary_model)

}
