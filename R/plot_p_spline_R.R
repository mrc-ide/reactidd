
# Save this file as `R/plot_p_spline_prev.R`

#' Plotting function for the P_spline model
#'
#' @export
#' @param X date vector.
#' @param p_splinefit fit of the model to the same set of data using reactidd::stan_p_spline()
#' @param target_dist_between_knots sets the number of days between adjacent knots (default = 5)
#' @param spline_degree sets the degree of the splines (default = 3)
#' @param ylim sets the ylimit of the plot
#' @param link_function sets the link function used for the model (logit for REACT-1, log for phe case data)
#' @param n shape parameter of the gamma distribution used for the generation time
#' @param b rate parameter of the gamma distribution used for the generation time
#' @return A list of the created plot, the raw data and CI's used in the plot, the raw data for the model fit in the plot.
#'
plot_p_spline_R <- function(X, p_spline_fit, target_dist_between_knots = 5, spline_degree = 3, ylim=2.0,
                            link_function="logit",n=2.29, b=0.36, tau_max=14){

  inv_logit <- function(Num){
    1/(1+exp(-Num))
  }

  ff <- rstan::extract(p_spline_fit)
  X <- as.numeric(X)
  X <- seq(min(X),max(X),by=1)

  gammaDist <- function(b, n, a) (b**n) * (a**(n-1)) * exp(-b*a) / gamma(n)

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

  dfGrowth <- matrix(NA, nrow=nrow(Y_array), ncol=length(X))

  calc_R <- function(b,n,tau_max,Y_vals,X, link_function){
    if(link_function =="logit"){
      Y_vals <- inv_logit(Y_vals)
    }
    else if(link_function =="log"){
      Y_vals <- exp(Y_vals)
    }

    R_list <- matrix(NA, nrow=1,ncol=length(X))
    for(i in seq_len(length(X))){
      if(i>tau_max){
        integral <- 0
        g_a <- 0
        for(j in seq_len(tau_max)){
          integral <- integral + Y_vals[i-j]*gammaDist(b,n,j)
          g_a <- g_a +gammaDist(b,n,j)
        }
        R_list[i] <- Y_vals[i]/ (integral/g_a)
      }
    }
    R_list
  }


  for(i in seq_len(nrow(dfGrowth))){
    dfGrowth[i,] <- calc_R(b,n,tau_max,Y_array[i,],X, link_function)
  }

  for(i in seq_len(length(X))){
    grad <- dfGrowth[,i]
    dfY$r[i] <-median(grad)
    dfY$lb_2.5[i] <- quantile(grad, probs=0.025, na.rm = TRUE)
    dfY$lb_25[i] <- quantile(grad, probs=0.25, na.rm = TRUE)
    dfY$ub_97.5[i] <- quantile(grad, probs=0.975, na.rm = TRUE)
    dfY$ub_75[i] <- quantile(grad, probs=0.75, na.rm = TRUE)
    dfY$prob[i] <- length(grad[grad>1.0])/length(grad)

  }

  dfY



  df_plot_model <- dfY
  df_plot_model$d_comb <- as.Date(df_plot_model$x-18383, origin=as.Date("2020-05-01"))


  max_date<-max(df_plot_model$d_comb)
  min_date<-min(df_plot_model$d_comb)
  df_plot_model<-df_plot_model[df_plot_model$d_comb>=min_date+tau_max & df_plot_model$d_comb<=max_date,]

  plot1 <- ggplot2::ggplot(data = df_plot_model, ggplot2::aes(x= d_comb, y =r))+
    ggplot2::coord_cartesian(ylim=c(0,ylim), xlim=c(min(df_plot_model$d_comb)-tau_max, max(df_plot_model$d_comb)))+
    ggplot2::theme_bw(base_size = 18)+
    ggplot2::xlab("Day of swab")+
    ggplot2::ylab("Reprodcution number")+
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
    ggplot2::geom_line(data = df_plot_model,
                       ggplot2::aes(y=prob),
                       color = 'red')

  return(list(plot1, df_plot_model))
}
