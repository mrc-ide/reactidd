
# Save this file as `R/plot_p_spline_prev.R`

#' Plotting function for the P_spline model
#'
#' @export
#' @param X date vector.
#' @param Y Numeric vector of number of positive samples
#' @param N Numeric vector of total number of samples
#' @param p_splinefit fit of the model to the same set of data using reactidd::stan_p_spline()
#' @param target_dist_between_knots sets the number of days between adjacent knots (default = 5)
#' @param spline_degree sets the degree of the splines (default = 3)
#' @param ylim sets the ylimit of the plot
#' @return A list of the created plot, the raw data and CI's used in the plot, the raw data for the model fit in the plot.
#'
plot_p_spline_prev <- function(X, Y, N, p_spline_fit, target_dist_between_knots = 5, spline_degree = 3, ylim=1.0){

  inv_logit <- function(Num){
    1/(1+exp(-Num))
  }

  X_og <- X
  ff <- rstan::extract(p_spline_fit)
  X <- as.numeric(X)
  X <- seq(min(X),max(X),by=1)

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


  dfY <- data.frame(x = X)
  for(i in seq_len(length(X))){
    dfY$p[i] <-median(Y_array[,i])
    dfY$lb_2.5[i] <- quantile(Y_array[,i], probs=0.025)
    dfY$lb_25[i] <- quantile(Y_array[,i], probs=0.25)
    dfY$ub_97.5[i] <- quantile(Y_array[,i], probs=0.975)
    dfY$ub_75[i] <- quantile(Y_array[,i], probs=0.75)

  }

  df_plot_model <- dfY
  df_plot_model$d_comb <- as.Date(df_plot_model$x-18383, origin=as.Date("2020-05-01"))


  CI <- prevalence::propCI(Y, N, level=0.95, method="wilson")
  df_plot <- data.frame(X=X_og, p = CI$p, lb= CI$lower, ub = CI$upper)
  df_plot$d_comb <- as.Date(df_plot$X)

  max_date<-max(df_plot$d_comb)
  min_date<-min(df_plot$d_comb)
  #df_plot_model<-df_plot_model[df_plot_model$d_comb>=min_date & df_plot_model$d_comb<=max_date,]

  plot1 <- ggplot2::ggplot(data = df_plot, ggplot2::aes(x= d_comb, y =p*100))+
    ggplot2::geom_point()+
    ggplot2::geom_errorbar(ggplot2::aes(ymin=lb*100, ymax=ub*100))+
    ggplot2::coord_cartesian(ylim=c(0,ylim), xlim=c(min(df_plot$d_comb), max(df_plot$d_comb)))+
    ggplot2::theme_bw(base_size = 18)+
    ggplot2::xlab("Day of swab")+
    ggplot2::ylab("Prevalence (%)")+
    ggplot2::scale_x_date(date_breaks = "1 month", date_labels = "%b")+
    ggplot2::geom_line(data= df_plot_model, ggplot2::aes(y=inv_logit(p)*100))+
    ggplot2::geom_ribbon(data = df_plot_model,
                         ggplot2::aes(y=inv_logit(p)*100,
                    ymin=inv_logit(lb_2.5)*100,
                    ymax=inv_logit(ub_97.5)*100),
                alpha=0.2)+
    ggplot2::geom_ribbon(data = df_plot_model,
                         ggplot2::aes(y=inv_logit(p)*100,
                    ymin=inv_logit(lb_25)*100,
                    ymax=inv_logit(ub_75)*100),
                alpha=0.2)

  return(list(plot1, df_plot, df_plot_model))
}

