
# Save this file as `R/plot_p_spline_prev.R`

#' Plotting function for the P_spline model
#'
#' @export
#' @param X date vector.
#' @param stan_fit fit of the model to the same set of data using reactidd::stan_p_spline()
#' @param ylim sets the ylimit of the plot
#' @param n shape parameter of the gamma distribution used for the generation time
#' @param b rate parameter of the gamma distribution used for the generation time
#' @return A list of the created plot of R over time and posterior estimates of R_t used in the plots.
#'
plot_breakpoint_R <- function(X, stan_fit, ylim=2.0, n=2.29, b=0.36, name="England"){

  ff <- rstan::extract(stan_fit)

  X <- seq(min(X),max(X),by=1)


  p_delay <- quantile(ff$delay, c(0.5,0.025,0.975))

  beta1   <- (1+(quantile(ff$beta[,1], c(0.5,0.025,0.975))/b))**n
  beta2   <- (1+(quantile(ff$beta[,2], c(0.5,0.025,0.975))/b))**n
  beta3   <- (1+(quantile(ff$beta[,3], c(0.5,0.025,0.975))/b))**n
  beta4   <- (1+(quantile(ff$beta[,4], c(0.5,0.025,0.975))/b))**n
  beta5   <- (1+(quantile(ff$beta[,5], c(0.5,0.025,0.975))/b))**n
  beta6   <- (1+(quantile(ff$beta[,6], c(0.5,0.025,0.975))/b))**n

  df1 <- data.frame(name = name,period = "0", beta = beta1[[1]],beta_lb=beta1[[2]],beta_ub=beta1[[3]])
  df2 <- data.frame(name = name,period = "1", beta = beta2[[1]],beta_lb=beta2[[2]],beta_ub=beta2[[3]])
  df3 <- data.frame(name = name,period = "2", beta = beta3[[1]],beta_lb=beta3[[2]],beta_ub=beta3[[3]])
  df4 <- data.frame(name = name,period = "3", beta = beta4[[1]],beta_lb=beta4[[2]],beta_ub=beta4[[3]])
  df5 <- data.frame(name = name,period = "4", beta = beta5[[1]],beta_lb=beta5[[2]],beta_ub=beta5[[3]])
  df6 <- data.frame(name = name,period = "5", beta = beta6[[1]],beta_lb=beta6[[2]],beta_ub=beta6[[3]])

  df_R <- rbind(df1, df2, df3, df4, df5, df6)
  df_delay <- data.frame(name=name, delay = p_delay[[1]], delay_lb = p_delay[[2]], delay_ub = p_delay[[3]])

  test<-(1+(ff$beta[,2]/b))**n/(1+(ff$beta[,1]/b))**n
  test[ff$beta[,2]< -b] <- 0.0
  test[ff$beta[,1]< -b] <- 10000
  beta1   <- quantile( test, c(0.5,0.025,0.975))

  test<-(1+(ff$beta[,3]/b))**n/(1+(ff$beta[,2]/b))**n
  test[ff$beta[,3]< -b] <- 0.0
  test[ff$beta[,2]< -b] <- 10000
  beta2   <- quantile( test, c(0.5,0.025,0.975))

  test<-(1+(ff$beta[,4]/b))**n/(1+(ff$beta[,3]/b))**n
  test[ff$beta[,4]< -b] <- 0.0
  test[ff$beta[,3]< -b] <- 10000
  beta3   <- quantile( test, c(0.5,0.025,0.975))

  test<-(1+(ff$beta[,5]/b))**n/(1+(ff$beta[,4]/b))**n
  test[ff$beta[,5]< -b] <- 0.0
  test[ff$beta[,4]< -b] <- 10000
  beta4   <- quantile( test, c(0.5,0.025,0.975))

  test<-(1+(ff$beta[,6]/b))**n/(1+(ff$beta[,5]/b))**n
  test[ff$beta[,6]< -b] <- 0.0
  test[ff$beta[,5]< -b] <- 10000
  beta5   <- quantile( test, c(0.5,0.025,0.975))

  df1 <- data.frame(name = name,period = "1", beta = beta1[[1]],beta_lb=beta1[[2]],beta_ub=beta1[[3]])
  df2 <- data.frame(name = name,period = "2", beta = beta2[[1]],beta_lb=beta2[[2]],beta_ub=beta2[[3]])
  df3 <- data.frame(name = name,period = "3", beta = beta3[[1]],beta_lb=beta3[[2]],beta_ub=beta3[[3]])
  df4 <- data.frame(name = name,period = "4", beta = beta4[[1]],beta_lb=beta4[[2]],beta_ub=beta4[[3]])
  df5 <- data.frame(name = name,period = "5", beta = beta5[[1]],beta_lb=beta5[[2]],beta_ub=beta5[[3]])

  df_mult <- rbind(df1, df2, df3, df4, df5)

  eng <- list(df_R, df_delay, df_mult)



  #########################################################################
  X_all <- seq(0,194,1)
  min_X <- 18626
  df_plot_model<-data.frame()
  df_beta<-data.frame()
  df_betat<-data.frame()
  for(i in seq_len(ncol(ff$prev_t))){

    row_temp <- data.frame(X=X_all[i]+min_X,
                           p = 100*quantile(ff$prev_t[,i], c(0.5,0.025,0.975))[[1]],
                           lb1= 100*quantile(ff$prev_t[,i], c(0.5,0.025,0.975))[[2]],
                           ub1= 100*quantile(ff$prev_t[,i], c(0.5,0.025,0.975))[[3]],
                           lb2= 100*quantile(ff$prev_t[,i], c(0.5,0.25,0.75))[[2]],
                           ub2= 100*quantile(ff$prev_t[,i], c(0.5,0.25,0.75))[[3]])

    Rlist <- (1+(log(ff$beta_t[,i])/b))^n
    row_tempbeta <- data.frame(X=X_all[i]+min_X,
                               beta = quantile(Rlist, c(0.5,0.025,0.975))[[1]],
                               lb1= quantile(Rlist, c(0.5,0.025,0.975))[[2]],
                               ub1= quantile(Rlist, c(0.5,0.025,0.975))[[3]],
                               lb2= quantile(Rlist, c(0.5,0.25,0.75))[[2]],
                               ub2= quantile(Rlist, c(0.5,0.25,0.75))[[3]])

    Rlist <- (1+(log(ff$beta_t[,i])/b))^n
    row_tempbetat <- data.frame(X=X_all[i]+min_X,
                                beta = quantile(Rlist, c(0.5,0.025,0.975))[[1]],
                                lb1= quantile(Rlist, c(0.5,0.025,0.975))[[2]],
                                ub1= quantile(Rlist, c(0.5,0.025,0.975))[[3]],
                                lb2= quantile(Rlist, c(0.5,0.25,0.75))[[2]],
                                ub2= quantile(Rlist, c(0.5,0.25,0.75))[[3]])

    df_plot_model <- rbind(df_plot_model, row_temp)
    df_beta <- rbind(df_beta, row_tempbeta)
    df_betat <- rbind(df_betat, row_tempbetat)
  }

  df_plot_model$d_comb <- as.Date(df_plot_model$X-18383, origin=as.Date("2020-05-01"))
  df_beta$d_comb <- as.Date(df_beta$X-18383, origin=as.Date("2020-05-01"))
  df_betat$d_comb <- as.Date(df_betat$X-18383, origin=as.Date("2020-05-01"))


  plot2<-ggplot2::ggplot(data=df_beta, ggplot2::aes(x= d_comb, y =p*100))+
    ggplot2::coord_cartesian(ylim=c(0.0,3))+
    ggplot2::theme_bw(base_size = 18)+
    ggplot2::xlab("Date")+
    ggplot2::ylab("Reproduction number")+
    ggplot2::scale_x_date(date_breaks = "1 month", date_labels = "%b\n%Y")+
    ggplot2::geom_line(data= df_beta, ggplot2::aes(x=d_comb,y=beta))+
    ggplot2::geom_ribbon(data= df_beta, ggplot2::aes(x=d_comb,y=beta, ymin=lb1, ymax=ub1), alpha=0.2)+
    ggplot2::geom_ribbon(data= df_beta, ggplot2::aes(x=d_comb,y=beta, ymin=lb2, ymax=ub2), alpha=0.2)+
    ggplot2::geom_hline(yintercept = 1, linetype='dashed',color='black')


  eng[[3]]$d_comb <- rest_dates + mean(ff$delay)

  plot2 <- plot2+
    ggplot2::geom_vline(xintercept = rest_dates[[1]], linetype="dashed")+
    ggplot2::geom_vline(xintercept = rest_dates[[2]], linetype="dashed")+
    ggplot2::geom_vline(xintercept = rest_dates[[3]], linetype="dashed")+
    ggplot2::geom_vline(xintercept = rest_dates[[4]], linetype="dashed")+
    ggplot2::geom_vline(xintercept = rest_dates[[5]], linetype="dashed")+
    ggplot2::geom_point(data=eng[[3]], ggplot2::aes(y=beta, ymin=beta_lb, ymax=beta_ub, color=period))+
    ggplot2::geom_errorbar(data=eng[[3]], ggplot2::aes(y=beta, ymin=beta_lb, ymax=beta_ub, color=period), width=0.1)+
    ggplot2::scale_color_brewer(palette="Set1",labels=c("Lockdown","Step 1a","Step 1b","Step 2","Step 3"))+
    ggplot2::theme(legend.position = c(0.83,0.73),
          legend.title = ggplot2::element_blank(),
          panel.grid = ggplot2::element_blank())+
    ggplot2::scale_y_continuous(sec.axis = ggplot2::sec_axis(~., name = "Multiplicative growth"))+
    ggplot2::coord_cartesian(ylim=c(0,4),xlim=c(as.Date("2020-12-30"), as.Date("2021-07-12")))+
    ggplot2::scale_x_date(date_breaks = "1 month", date_labels = "%b\n%Y")




  return(list(plot2, df_beta, eng[[3]]))
}
