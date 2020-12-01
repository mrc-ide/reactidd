
# Save this file as `R/plot_p_spline_phe.R`

#' Plotting function for the P_spline model
#'
#' @export
#' @param X date vector.
#' @param Y Numeric vector of number of positive samples
#' @param p_splinefit fit of the model to the same set o data using reacttemporal2::stan_p_spline_phe()
#' @param target_dist_between_knots sets the number of days between adjacent knots (default = 5)
#' @param spline_degree sets the degree of the splines (default = 3)
#' @param ylim sets the ylimit of the plot created
#' @return A list of the created plot, the raw data and CI's used in the plot, the raw data for the model fit in the plot.
#'
plot_p_spline_phe <- function(X, Y, p_spline_fit, target_dist_between_knots = 5, spline_degree = 3, ylim=5000.0){

  ff <- rstan::extract(p_spline_fit)

  X <- as.numeric(X)
  min_date_numeric <- min(X)
  max_date_numeric <- max(X)

  num_knots <- ceiling((max_date_numeric- min_date_numeric)/target_dist_between_knots) + 7
  days_per_knot <- (max_date_numeric - min_date_numeric)/(num_knots - 7)
  num_basis <- num_knots + spline_degree - 1

  #' Plot of the model fit
  X_new <- seq(min(X)-3*days_per_knot, max(X)+3*days_per_knot, 0.1)

  Y_array <- array(data=NA, dim=c(nrow(ff$a), length(X_new)))
  B_true <- t(splines::bs(X_new, df=num_basis, degree = spline_degree, intercept=TRUE))

  for(i in seq_len(nrow(ff$a))){
    a <- array(NA, num_basis)

    for(j in seq_len(length(a))){
      a[j] <- ff$a[i,j]
    }

    Y_array[i,] <- as.vector(a%*%B_true)
  }


  dfY <- data.frame(x = X_new)
  for(i in seq_len(length(X_new))){
    dfY$p[i] <-median(Y_array[,i])
    dfY$lb_2.5[i] <- quantile(Y_array[,i], probs=0.025)
    dfY$lb_25[i] <- quantile(Y_array[,i], probs=0.25)
    dfY$ub_97.5[i] <- quantile(Y_array[,i], probs=0.975)
    dfY$ub_75[i] <- quantile(Y_array[,i], probs=0.75)

  }

  df_plot_model <- dfY
  df_plot_model$d_comb <- as.Date(df_plot_model$x-18383, origin=as.Date("2020-05-01"))

  df_plot <- data.frame(X=X, p=Y)
  df_plot$d_comb <- as.Date(df_plot$X-18383, origin=as.Date("2020-05-01"))

  max_date<-max(df_plot$d_comb)
  min_date<-min(df_plot$d_comb)

  df_plot_model<-df_plot_model[df_plot_model$d_comb>=min_date & df_plot_model$d_comb<=max_date,]

  plot1 <- ggplot2::ggplot(data = df_plot, ggplot2::aes(x= d_comb, y =p))+
    ggplot2::geom_point()+
    ggplot2::coord_cartesian(ylim=c(0,ylim), xlim=c(min(df_plot$d_comb), max(df_plot$d_comb)))+
    ggplot2::theme_bw(base_size = 18)+
    ggplot2::xlab("Day of swab")+
    ggplot2::ylab("Number of cases")+
    ggplot2::scale_x_date(date_breaks = "1 month", date_labels = "%b")+
    ggplot2::geom_line(data= df_plot_model, ggplot2::aes(y=exp(p)))+
    ggplot2::geom_ribbon(data = df_plot_model,
                         ggplot2::aes(y=exp(p),
                                      ymin=exp(lb_2.5),
                                      ymax=exp(ub_97.5)),
                         alpha=0.2)+
    ggplot2::geom_ribbon(data = df_plot_model,
                         ggplot2::aes(y=exp(p),
                                      ymin=exp(lb_25),
                                      ymax=exp(ub_75)),
                         alpha=0.2)

  return(list(plot1, df_plot, df_plot_model))
}

