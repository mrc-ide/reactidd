
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
#' @param link_function the link function used in the original fit (logit for REACT-1, log for PHE case numbers)
#' @return A list of the created plot, the raw data and CI's used in the plot, the raw data for the model fit in the plot.
#'
plot_p_spline_igr <- function(X, p_spline_fit, target_dist_between_knots = 5, spline_degree = 3, ylim=1.0, link_function="logit"){

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
    Y_array[i,] <- as.vector(a%*%B_true)
  }
  dfY <- data.frame(x = X)

  dfGrowth <- matrix(NA, nrow=nrow(Y_array), ncol=length(X)-1)

  for(i in seq_len(nrow(dfGrowth))){
    first_d <- deriv(X, Y_array[i,])
    exp_term <- exp(Y_array[i,])
    first_d <- first_d[1:length(first_d)]
    exp_term <- exp_term[2:length(exp_term)]

    if(link_function =="logit"){
      growth_rate <- (first_d/(exp_term+1))
    }
    else if(link_function =="log"){
      growth_rate <- first_d
    }

    dfGrowth[i,] <- growth_rate

  }

  for(i in seq_len(length(X)-1)){
    grad <- dfGrowth[,i]
    dfY$r[i] <-median(grad)
    dfY$lb_2.5[i] <- quantile(grad, probs=0.025)
    dfY$lb_25[i] <- quantile(grad, probs=0.25)
    dfY$ub_97.5[i] <- quantile(grad, probs=0.975)
    dfY$ub_75[i] <- quantile(grad, probs=0.75)
    dfY$prob[i] <- length(grad[grad>0.0])/length(grad)

  }


  df_plot_model <- dfY[1:nrow(dfY)-1 ,]
  df_plot_model$d_comb <- as.Date(df_plot_model$x-18383, origin=as.Date("2020-05-01"))

  max_date<-max(df_plot_model$d_comb)
  min_date<-min(df_plot_model$d_comb)

  plot1 <- ggplot2::ggplot(data = df_plot_model, ggplot2::aes(x= d_comb, y =r))+
    ggplot2::coord_cartesian(ylim=c(-ylim,ylim), xlim=c(min(df_plot_model$d_comb), max(df_plot_model$d_comb)))+
    ggplot2::theme_bw(base_size = 18)+
    ggplot2::xlab("Day of swab")+
    ggplot2::ylab("Prevalence (%)")+
    ggplot2::scale_x_date(date_breaks = "1 month", date_labels = "%b")+
    ggplot2::geom_line(data= df_plot_model, ggplot2::aes(y=r))+
    ggplot2::geom_ribbon(data = df_plot_model,
                         ggplot2::aes(y=r,
                                      ymin=lb_2.5,
                                      ymax=ub_97.5),
                         alpha=0.2)+
    ggplot2::geom_ribbon(data = df_plot_model,
                         ggplot2::aes(y=r,
                                      ymin=lb_25,
                                      ymax=ub_75),
                         alpha=0.2)+
    ggplot2::geom_hline(yintercept = 0.0, linetype="dashed")

  return(list(plot1, df_plot_model))
}




