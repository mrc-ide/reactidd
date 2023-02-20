# Save this file as `R/clean_death_data.R`

#' Bayesian p-spline model using stan
#'
#' @export
#' @param ff p-spline model posterior obtained using rstan::extract() on a p-spline model stanfit object
#' @param X all dates over which the p-spline is to be determined
#' @param days_per_knot parameter used in p-spline model fit
#' @param spline_degree parameter used in p-spline model fit
#' @return An array of the posterior P-spline response on all dates set by X

get_response_posterior <- function(ff, X,
                                   days_per_knot = 5,
                                   spline_degree = 3){
  min_date_numeric <- min(X)
  max_date_numeric <- max(X)
  num_knots <- ceiling((max_date_numeric- min_date_numeric)/days_per_knot)+7
  days_per_knot <- (max_date_numeric - min_date_numeric)/(num_knots -7)
  num_basis <- num_knots + spline_degree - 1
  num_data <- length(X)
  knots <- unname(seq(min(X)-3*days_per_knot, max(X)+3*days_per_knot, length.out = num_knots))



  X_new <- seq(min(X)-3*days_per_knot, max(X)+3*days_per_knot, 0.1)
  B_true <- splines::bs(X_new, df=num_basis, degree=spline_degree, intercept = TRUE)
  B_true <- t(predict(B_true, X))
  Y_array <- array(data=NA, dim=c(nrow(ff$a), length(X)))

  #a0<-mean(ff$a0)
  for(i in seq_len(nrow(ff$a))){
    a <- array(NA, num_basis)
    #a0 <- ff$a0[i]
    for(j in seq_len(length(a))){
      a[j] <- ff$a[i,j]
      #a[j] <- mean(ff$a[,j])
    }
    Y_array[i,] <- as.vector(a%*%B_true)
  }


  Y_array
}
