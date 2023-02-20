# Save this file as `R/clean_death_data.R`

#' Bayesian p-spline model using stan
#'
#' @export
#' @param X1  dates of first array
#' @param Y_array1 full posterior of the log of a p-spline with ncol = number of dates
#' @param X2 dates of second array
#' @param Y_array2 full posterior of the log of a p-spline with ncol = number of dates
#' @param time_delay time lag between Y_array1 and Y_array2
#' @return A dataframe with columns for the log difference between the responses of two different p-spline model fits
#'

plot_spline_difference <- function(X1, Y_array1,
                                   X2, Y_array2,
                                   time_delay){
  minrow<-min(nrow(Y_array1), nrow(Y_array2))
  Y_array1 <- Y_array1[1:minrow,]
  Y_array2 <- Y_array2[1:minrow,]

  X2 <- X2 - time_delay
  common_dates <- X1[X1 %in% X2]

  diff_df <- data.frame(date = common_dates)

  for(i in seq_len(nrow(diff_df))){
    temp_date <- diff_df$date[i]
    index1 <- match(temp_date, X1)
    index2 <- match(temp_date, X2)

    diff_vec = Y_array1[,index1] - Y_array2[,index2]
    diff_df$mean[i] <- mean(diff_vec)
    diff_df$lb[i] <- quantile(diff_vec, c(0.025,0.975))[[1]]
    diff_df$ub[i] <- quantile(diff_vec, c(0.025,0.975))[[2]]
  }

  diff_df
}


