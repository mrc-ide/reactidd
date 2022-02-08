#' Formatting simulated individual level data to be fed to stan model
#'
#' @export
#' @param she_dat simulated indiviudal level data obtained from `reactidd::shedding_simultate_data`function.
#' @return A list of input data that can be given to the function `reactidd::stan_both_shedding_exp_models`
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
