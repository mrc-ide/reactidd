# Save this file as `R/clean_hosp_data.R`

#' Bayesian p-spline model using stan
#'
#' @export
#' @param death_dat public data set of hosp data by age as downloaded from ...
#' @param min_date minimum date to use
#' @param max_date maximum date to use
#' @return A dataframe with columns for all hospitalizations, hospitalizations in those aged 64 and under and for those aged 65 and over
#'

clean_hosp_data <- function(hosp_dat, min_date = as.Date("2020-05-01"), max_date=as.Date("2021-08-31")){
  hosp_dat <- maditr::dcast(hosp_dat, date~age, value.var = "value")
  hosp_dat[2:nrow(hosp_dat),2:6] <- hosp_dat[2:nrow(hosp_dat),2:6] - hosp_dat[1:(nrow(hosp_dat)-1),2:6]
  hosp_dat$age_64_under <- hosp_dat[,2]+hosp_dat[,3]+hosp_dat[,4]
  hosp_dat$age_65_over <- hosp_dat[,5]+hosp_dat[,6]
  hosp_dat$all <- hosp_dat$age_64_under + hosp_dat$age_65_over
  hosp_dat$date <- as.Date(hosp_dat$date)
  hosp_dat <- hosp_dat[hosp_dat$date>= min_date,]
  hosp_dat <- hosp_dat[hosp_dat$date<= max_date,]
}
