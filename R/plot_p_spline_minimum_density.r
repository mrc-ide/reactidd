
# Save this file as `R/plot_p_spline_minimum_density.R`

#' Plotting function for the P_spline model
#'
#' @export
#' @param X date vector.
#' @param p_splinefit fit of the model to the same set of X values using reactidd::stan_p_spline() or reactidd::stan_p_spline_phe()
#' @param target_dist_between_knots sets the number of days between adjacent knots (default = 5)
#' @param spline_degree sets the degree of the splines (default = 3)
#' @return A list of the created plot, and the list of values used in the histogram.
#'
plot_p_spline_minimum_density <- function(X, p_spline_fit, target_dist_between_knots = 5, spline_degree = 3){

  inv_logit <- function(Num){
    1/(1+exp(-Num))
  }

  ff <- rstan::extract(p_spline_fit)
  X <- as.numeric(X)
  X <- seq(min(X),max(X),by=1)
  deriv <- function(x, y) diff(y)/diff(x)

  days_per_knot<-5
  min_date_numeric <- min(X)
  max_date_numeric <- max(X)
  num_knots <- ceiling((max_date_numeric- min_date_numeric)/target_dist_between_knots)+7
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
    Y_array[i,] <- as.vector(a%*%B_true) #as.vector(a0*X+a%*%B_true)
  }

  index_list<-c()
  for(i in seq_len(nrow(Y_array))){
    min_index<- which.min(Y_array[i,(1):(ncol(Y_array))])
    index_list<-c(index_list, min_index)
  }


  index_list <- as.Date(index_list - 1, origin = as.Date("2020-05-01"))
  df<- as.data.frame(index_list)
  plot1 <- ggplot2::ggplot(data = df,ggplot2::aes(x=index_list, y=..density..))+
    ggplot2::coord_cartesian(xlim = c(as.Date(min_date_numeric - 18383, origin = as.Date("2020-05-01")),
                                      as.Date(max_date_numeric - 18383, origin = as.Date("2020-05-01"))))+
    ggplot2::geom_histogram(binwidth = 1)+
    ggplot2::theme_bw()+
    ggplot2::xlab("Date of minimum prevalence")+
    ggplot2::ylab("Probability Density")



  return(list(plot1, index_list))
}

