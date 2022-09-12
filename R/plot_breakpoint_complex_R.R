
# Save this file as `R/plot_p_spline_prev.R`

#' Plotting function for the P_spline model
#'
#' @export
#' @param X date vector.
#' @param stan_fit fit of the model to the same set of data using reactidd::stan_p_spline()
#' @param rest_dates date of key restriction changes
#' @param ylim sets the ylimit of the plot
#' @param n shape parameter of the gamma distribution used for the generation time
#' @param b rate parameter of the gamma distribution used for the generation time
#' @return A list of the created plot of R over time and posterior estimates of R_t used in the plots.
#'
plot_breakpoint_complex_R <- function(X, stan_fit,rest_dates, ylim=4.0, n=2.29, b=0.36){

  ff <- rstan::extract(stan_fit)

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

    Rlist <- (1+(log(ff$beta_int[,i])/b))^n
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


  plot4<-ggplot2::ggplot()+
    ggplot2::coord_cartesian(ylim=c(0.0,ylim))+
    ggplot2::theme_bw(base_size = 18)+
    ggplot2::xlab("Date")+
    ggplot2::ylab("Reproduction number")+
    ggplot2::scale_x_date(date_breaks = "1 month", date_labels = "%b")+
    ggplot2::geom_line(data= df_betat, ggplot2::aes(x=d_comb,y=beta), color='red')+
    ggplot2::geom_ribbon(data= df_betat, ggplot2::aes(x=d_comb,y=beta, ymin=lb1, ymax=ub1), alpha=0.2, fill='red')+
    ggplot2::geom_ribbon(data= df_betat, ggplot2::aes(x=d_comb,y=beta, ymin=lb2, ymax=ub2), alpha=0.2, fill='red')+
    ggplot2::geom_line(data= df_beta, ggplot2::aes(x=d_comb,y=beta), color='blue')+
    ggplot2::geom_ribbon(data= df_beta, ggplot2::aes(x=d_comb,y=beta, ymin=lb1, ymax=ub1), alpha=0.2, fill='blue')+
    ggplot2::geom_ribbon(data= df_beta, ggplot2::aes(x=d_comb,y=beta, ymin=lb2, ymax=ub2), alpha=0.2, fill='blue')+
    ggplot2::geom_hline(yintercept = 1, linetype='dashed',color='black')+
    ggplot2::scale_y_continuous(sec.axis = ggplot2::sec_axis(~., name = "Intrinsic reproduction number"))+
    ggplot2::geom_vline(xintercept = rest_dates[[1]], linetype="dashed")+
    ggplot2::geom_vline(xintercept = rest_dates[[2]], linetype="dashed")+
    ggplot2::geom_vline(xintercept = rest_dates[[3]], linetype="dashed")+
    ggplot2::geom_vline(xintercept = rest_dates[[4]], linetype="dashed")+
    ggplot2::geom_vline(xintercept = rest_dates[[5]], linetype="dashed")


  return(list(plot4, df_beta, df_betat))
}




