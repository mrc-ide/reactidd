# Save this file as `R/clean_death_data.R`

#' Bayesian p-spline model using stan
#'
#' @export
#' @param death_dat public data set of death data by age as downloaded from ...
#' @param min_date minimum date to use
#' @param max_date maximum date to use
#' @return A dataframe with columns for all deaths, deaths in those aged 64 and under and for those aged 65 and over
#'

clean_death_data <- function(death_dat, min_date = as.Date("2020-05-01"), max_date=as.Date("2021-08-31")){
  death_dat <- maditr::dcast(death_dat, date~age, value.var = "deaths")
  death_dat$age_64_under <- death_dat[,3]+death_dat[,16]
  death_dat$age_65_over <- death_dat[,17]+death_dat[,18]+death_dat[,19]+death_dat[,20]+death_dat[,21]+death_dat[,22]
  death_dat$all <- death_dat$age_64_under + death_dat$age_65_over
  death_dat$date <- as.Date(death_dat$date)
  death_dat <- death_dat[death_dat$date>= min_date,]
  death_dat <- death_dat[death_dat$date<= max_date,]

}
