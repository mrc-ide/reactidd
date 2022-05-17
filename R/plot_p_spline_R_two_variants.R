
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
#' @param n1 shape parameter of the gamma distribution used for the generation time of Variant 1
#' @param b1 rate parameter of the gamma distribution used for the generation time of Variant 1
#' @param n2 shape parameter of the gamma distribution used for the generation time of Variant 2
#' @param b2 rate parameter of the gamma distribution used for the generation time of Variant 2
#' @param tau_max width of window to perform R calculations over (default is two weeks)
#' @param mindateVar1 first date that variant 1 was detected
#' @param mindateVar2 first date that variant 2 was detected
#' @param maxdateVar1 last date that variant 1 was detected
#' @param maxdateVar2 last date that variant 2 was detected
#' @return A list of the created plot of R over time and posterior estimates of R_t used in the plots.
#'
plot_p_spline_R_two_variants <- function(X, p_spline_fit, target_dist_between_knots = 5, spline_degree = 3, ylim=2.5,
                            link_function="logit",n1=2.29, b1=0.36,n2=2.29, b2=0.36, tau_max=14,
                            labs = c("Delta","Omicron"), colors1=c("blue","red"),
                            mindateVar1 = as.Date("2021-09-09"), maxdateVar1 = as.Date("2022-02-14"),
                            mindateVar2 = as.Date("2021-12-03"), maxdateVar2 = as.Date("2022-03-03")){

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
  Y_array1 <- array(data=NA, dim=c(nrow(ff$a1), length(X)))
  Y_array2 <- array(data=NA, dim=c(nrow(ff$a2), length(X)))

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
  dfY1 <- data.frame(x = X)
  dfY2 <- data.frame(x = X)

  dfGrowth1 <- matrix(NA, nrow=nrow(Y_array1), ncol=length(X))
  dfGrowth2 <- matrix(NA, nrow=nrow(Y_array2), ncol=length(X))

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


  for(i in seq_len(nrow(dfGrowth1))){
    dfGrowth1[i,] <- calc_R(b1,n1,tau_max,Y_array1[i,],X, link_function)
    dfGrowth2[i,] <- calc_R(b2,n2,tau_max,Y_array2[i,],X, link_function)
  }

  for(i in seq_len(length(X))){
    grad <- dfGrowth1[,i]
    dfY1$r[i] <-median(grad)
    dfY1$lb_2.5[i] <- quantile(grad, probs=0.025, na.rm = TRUE)
    dfY1$lb_25[i] <- quantile(grad, probs=0.25, na.rm = TRUE)
    dfY1$ub_97.5[i] <- quantile(grad, probs=0.975, na.rm = TRUE)
    dfY1$ub_75[i] <- quantile(grad, probs=0.75, na.rm = TRUE)
    dfY1$prob[i] <- length(grad[grad>1.0])/length(grad)

    grad <- dfGrowth2[,i]
    dfY2$r[i] <-median(grad)
    dfY2$lb_2.5[i] <- quantile(grad, probs=0.025, na.rm = TRUE)
    dfY2$lb_25[i] <- quantile(grad, probs=0.25, na.rm = TRUE)
    dfY2$ub_97.5[i] <- quantile(grad, probs=0.975, na.rm = TRUE)
    dfY2$ub_75[i] <- quantile(grad, probs=0.75, na.rm = TRUE)
    dfY2$prob[i] <- length(grad[grad>1.0])/length(grad)

  }

  #########


  nat_plot_1 <- dfY1
  nat_plot_1$d_comb <- as.Date(nat_plot_1$x-18383, origin=as.Date("2020-05-01"))


  nat_plot_2 <- dfY2
  nat_plot_2$d_comb <- as.Date(nat_plot_2$x-18383, origin=as.Date("2020-05-01"))


  max_date<-max(nat_plot_1$d_comb)
  min_date<-min(nat_plot_1$d_comb)

  nat_plot_2<-nat_plot_2[nat_plot_2$d_comb>=mindateVar2+tau_max & nat_plot_2$d_comb<=maxdateVar2,]
  #nat_plot_2[is.na(nat_plot_2$r)==TRUE,]$prob <- NA
  nat_plot_1<-nat_plot_1[nat_plot_1$d_comb>=mindateVar1+tau_max & nat_plot_1$d_comb<=maxdateVar1,]
  #nat_plot_1[is.na(nat_plot_1$r)==TRUE,]$prob <- NA



  plot9<-ggplot2::ggplot(data = nat_plot_2, ggplot2::aes(x=d_comb,y=r))+
    ggplot2::geom_line(ggplot2::aes(color="Omicron"))+
    ggplot2::geom_ribbon(ggplot2::aes(ymin=lb_2.5, ymax=ub_97.5, fill="Omicron"), alpha=0.2)+
    ggplot2::geom_ribbon(ggplot2::aes(ymin=lb_25, ymax=ub_75, fill="Omicron"), alpha=0.2)+
    ggplot2::geom_line(data=nat_plot_1, ggplot2::aes(color="Delta"))+
    ggplot2::geom_ribbon(data=nat_plot_1,ggplot2::aes(ymin=lb_2.5, ymax=ub_97.5,fill="Delta"), alpha=0.2)+
    ggplot2::geom_ribbon(data=nat_plot_1,ggplot2::aes(ymin=lb_25, ymax=ub_75,fill="Delta"), alpha=0.2)+
    ggplot2::coord_cartesian(ylim=c(0,ylim), xlim=c(min_date, max_date))+
    ggplot2::theme_bw(base_size = 10)+
    ggplot2::geom_hline(yintercept = 1, linetype='dashed')+
    ggplot2::xlab("Date")+
    ggplot2::ylab("Reproduction number\n(2 week average)")+
    ggplot2:: scale_x_date(date_breaks = "1 month", date_labels = "%b\n%Y")+
    ggplot2::theme(legend.position = "bottom")+
    ggplot2::scale_color_manual(values=colors1,
                                name= "Variant",
                                labels=labs)+
    ggplot2::scale_fill_manual(values=colors1,
                               name= "Variant",
                               labels=labs)



  return(list(plot9,nat_plot_1, nat_plot_2))
}


