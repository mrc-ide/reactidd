
# Save this file as `R/plot_p_spline_prev.R`

#' Plotting function for the P_spline model
#'
#' @export
#' @param X date vector.
#' @param p_splinefit fit of the model to the same set of data using reactidd::stan_p_spline()
#' @param target_dist_between_knots sets the number of days between adjacent knots (default = 5)
#' @param spline_degree sets the degree of the splines (default = 3)
#' @param ylim sets the ylimit of the plot
#' @param link_function the link function used in the original fit (logit for REACT-1, log for PHE case numbers)
#' @param labs names of the the two variants compared
#' @param colors1 sets the three colours used for growth rate advantage and both variant's growth rates
#' @return A list of the created plot for instantaneous growth rate and the posterior estimates for growth rate used in the plot.
#'
plot_p_spline_igr_two_variants <- function(X, p_spline_fit, target_dist_between_knots = 5, spline_degree = 3, ylim=1.0, link_function="logit",
                                           labs = c("Delta","Omicron"), colors1=c("black","blue","red"),
                                           mindateVar1 = as.Date("2021-09-09"), maxdateVar1 = as.Date("2022-02-14"),
                                           mindateVar2 = as.Date("2021-12-03"), maxdateVar2 = as.Date("2022-03-03")){

  inv_logit <- function(Num){
    1/(1+exp(-Num))
  }

  ff <- rstan::extract(p_spline_fit)
  X <- as.numeric(X)
  X <- seq(min(X),max(X),by=1)
  deriv <- function(x, y) diff(y)/diff(x)


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
  Y_array1 <- array(data=NA, dim=c(nrow(ff$a1), length(X)))
  Y_array2 <- array(data=NA, dim=c(nrow(ff$a2), length(X)))

  #a0<-mean(ff$a0)
  for(i in seq_len(nrow(ff$a1))){
    a1 <- array(NA, num_basis)
    a2 <- array(NA, num_basis)
    #a0 <- ff$a0[i]
    for(j in seq_len(length(a1))){
      a1[j] <- ff$a1[i,j]
      a2[j] <- ff$a2[i,j]
    }
    Y_array1[i,] <- as.vector(a1%*%B_true)
    Y_array2[i,] <- as.vector(a2%*%B_true)
  }
  dfY  <- data.frame(x = X[1:length(X)-1])
  dfY1 <- data.frame(x = X[1:length(X)-1])
  dfY2 <- data.frame(x = X[1:length(X)-1])

  dfGrowth1 <- matrix(NA, nrow=nrow(Y_array1), ncol=length(X)-1)
  dfGrowth2 <- matrix(NA, nrow=nrow(Y_array2), ncol=length(X)-1)
  dfGrowthC <- matrix(NA, nrow=nrow(Y_array1), ncol=length(X)-1)

  for(i in seq_len(nrow(dfGrowth1))){
    first_d1 <- deriv(X, Y_array1[i,])
    exp_term1 <- exp(Y_array1[i,])
    first_d1 <- first_d1[1:length(first_d1)]
    exp_term1 <- exp_term1[2:length(exp_term1)]

    first_d2 <- deriv(X, Y_array2[i,])
    exp_term2 <- exp(Y_array2[i,])
    first_d2 <- first_d2[1:length(first_d2)]
    exp_term2 <- exp_term2[2:length(exp_term2)]

    if(link_function =="logit"){
      growth_rate1 <- (first_d1/(exp_term1+1))
      growth_rate2 <- (first_d2/(exp_term2+1))
    }
    else if(link_function =="log"){
      growth_rate1 <- first_d1
      growth_rate2 <- first_d2
    }

    dfGrowth1[i,] <- growth_rate1
    dfGrowth2[i,] <- growth_rate2
    dfGrowthC[i,] <- growth_rate1 - growth_rate2

  }

  for(i in seq_len(length(X)-1)){
    grad <- dfGrowthC[,i]
    dfY$r[i] <-median(grad)
    dfY$lb_2.5[i] <- quantile(grad, probs=0.025)
    dfY$lb_25[i] <- quantile(grad, probs=0.25)
    dfY$ub_97.5[i] <- quantile(grad, probs=0.975)
    dfY$ub_75[i] <- quantile(grad, probs=0.75)
    dfY$prob[i] <- length(grad[grad>0.0])/length(grad)

    grad <- dfGrowth1[,i]
    dfY1$r[i] <-median(grad)
    dfY1$lb_2.5[i] <- quantile(grad, probs=0.025)
    dfY1$lb_25[i] <- quantile(grad, probs=0.25)
    dfY1$ub_97.5[i] <- quantile(grad, probs=0.975)
    dfY1$ub_75[i] <- quantile(grad, probs=0.75)
    dfY1$prob[i] <- length(grad[grad>0.0])/length(grad)

    grad <- dfGrowth2[,i]
    dfY2$r[i] <-median(grad)
    dfY2$lb_2.5[i] <- quantile(grad, probs=0.025)
    dfY2$lb_25[i] <- quantile(grad, probs=0.25)
    dfY2$ub_97.5[i] <- quantile(grad, probs=0.975)
    dfY2$ub_75[i] <- quantile(grad, probs=0.75)
    dfY2$prob[i] <- length(grad[grad>0.0])/length(grad)

  }


  ###############################################################################################################

  df_plot_model_r1 <- dfY1#[1:nrow(dfY1)-1 ,]
  df_plot_model_r2 <- dfY2#[1:nrow(dfY2)-1 ,]
  df_plot_model_rC <- dfY#[1:nrow(dfY)-1 ,]

  df_plot_model_r1$d_comb <- as.Date(df_plot_model_r1$x-18383, origin=as.Date("2020-05-01"))
  df_plot_model_r2$d_comb <- as.Date(df_plot_model_r2$x-18383, origin=as.Date("2020-05-01"))
  df_plot_model_rC$d_comb <- as.Date(df_plot_model_rC$x-18383, origin=as.Date("2020-05-01"))

  min_omi_date <- as.Date("2021-12-03")
  max_delta_date <- as.Date("2022-02-14")

  max_date<-max(df_plot_model_r1$d_comb)
  min_date<-min(df_plot_model_r1$d_comb)

  df_plot_model_r1 <- df_plot_model_r1[df_plot_model_r1$d_comb>=mindateVar1 &df_plot_model_r1$d_comb<=maxdateVar1,]
  df_plot_model_r2 <- df_plot_model_r2[df_plot_model_r2$d_comb>=mindateVar2 &df_plot_model_r2$d_comb<=maxdateVar2,]
  df_plot_model_rC <- df_plot_model_rC[df_plot_model_rC$d_comb>=mindateVar2 &df_plot_model_rC$d_comb<=maxdateVar1,]




  plot1 <- ggplot2::ggplot()+
    ggplot2::geom_line(data=df_plot_model_rC, ggplot2::aes(x = d_comb, y=-r, ymin=-lb_2.5 , ymax=-ub_97.5, col='black'))+
    ggplot2::geom_ribbon(data=df_plot_model_rC, ggplot2::aes(x = d_comb, y=-r , ymin=-lb_2.5 , ymax=-ub_97.5, fill='black' ),alpha=0.3)+
    ggplot2::geom_line(data=df_plot_model_r1, ggplot2::aes(x = d_comb, y=r , ymin=lb_2.5 , ymax=ub_97.5 , col='blue'))+
    ggplot2::geom_ribbon(data=df_plot_model_r1, ggplot2::aes(x = d_comb, y=r , ymin=lb_2.5 , ymax=ub_97.5 , fill='blue'),alpha=0.3)+
    ggplot2::geom_line(data=df_plot_model_r2, ggplot2::aes(x = d_comb, y=r , ymin=lb_2.5 , ymax=ub_97.5 , col='red'))+
    ggplot2::geom_ribbon(data=df_plot_model_r2, ggplot2::aes(x = d_comb, y=r , ymin=lb_2.5 , ymax=ub_97.5 , fill='red'),alpha=0.3)+
    ggplot2::coord_cartesian(ylim=c(-ylim,ylim), xlim=c(min_date, max_date))+
    ggplot2::theme_bw(base_size = 10)+
    ggplot2::xlab("Date (2021-2022)")+
    ggplot2::ylab("Growth rate (1/day)")+
    ggplot2::scale_x_date(date_breaks = "1 month", date_labels = "%b\n%Y")+
    ggplot2::scale_color_manual(values=colors1,
                       name= "Variant",
                       labels=c("Advantage",labs))+
    ggplot2::scale_fill_manual(values=colors1,
                      name= "Variant",
                      labels=c("Advantage",labs))+
    ggplot2::theme(legend.position = "bottom")+
    ggplot2::geom_hline(yintercept = 0.0, linetype="dashed")




  return(list(plot1, df_plot_model_rC,df_plot_model_r1,df_plot_model_r2 ))
}
