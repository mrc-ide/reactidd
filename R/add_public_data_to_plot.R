# Save this file as `R/add_public_data_to_plot.R

#' MCMC model fitting - get mode of time lag vaector
#'
#' @export
#' @param p_spline_plot initial plot of REACT data with Bayesian P-spline fit to it
#' @param Y_array array with posterior of response variables for death/hosp bayesian P-spline fit
#' @param table_eng table with best fit parameter values for scaling parameter and time lag parameter
#' @param X list of dates for which death/hosp bayesian P-spline model was fit
#' @param dat raw data with date given by a 'date' column and outcome variable (death/ hosp) given by 'all' column
#' @return plot of react data with bayesian P-spline with another P-spline scaled and translated added on top

add_public_data_to_plot <- function(p_spline_plot, Y_array, table_eng, X, dat){

  scale_par <- exp(table_eng$alpha)*100
  lag_par <- table_eng$beta


  dfY <- data.frame(x = X)

  for(i in seq_len(length(X))){
    dfY$p[i] <-median(Y_array[,i])
    dfY$lb_2.5[i] <- quantile(Y_array[,i], probs=0.025)
    dfY$lb_25[i] <- quantile(Y_array[,i], probs=0.25)
    dfY$ub_97.5[i] <- quantile(Y_array[,i], probs=0.975)
    dfY$ub_75[i] <- quantile(Y_array[,i], probs=0.75)
  }


  plot1<-p_spline_plot+
    ggplot2::coord_cartesian(ylim=c(0.008,10))+
    ggplot2::scale_y_log10(name = "Swab positivity (%)",breaks = c(seq(0.01,0.09,0.01),seq(0.1,0.9,0.1),seq(1,9,1),seq(10,90,10),seq(100,900,100),seq(1000,9000,1000)), labels = replace(character(54), c(1,10,19,28,37,46), c(0.01,0.10,1.00,10.00,100,1000)), sec.axis = ggplot2::sec_axis(~.*1/scale_par, name ="Daily deaths", breaks = c(seq(0.01,0.09,0.01),seq(0.1,0.9,0.1),seq(1,9,1),seq(10,90,10),seq(100,900,100),seq(1000,9000,1000)), labels = replace(character(54), c(19,28,37,46), c(1,10,100,1000))))+
    ggplot2::theme_bw(base_size = 18)+
    ggplot2::xlab("Day of swab")+
    ggplot2::ylab("Swab positivity (%)")+
    ggplot2::scale_x_date(date_breaks = "1 month", date_labels = "%b")+
    ggplot2::geom_point(data=dat, ggplot2::aes(x=date-lag_par, y=all*scale_par),color='red', alpha=0.4, cex=0.5)+
    ggplot2::geom_line(data = dfY,ggplot2::aes(x=x-lag_par, y=p*scale_par), color='red')+
    ggplot2::geom_ribbon(data = dfY,ggplot2::aes(x=x-lag_par, y=p*scale_par, ymin=lb_2.5*scale_par, ymax=ub_97.5*scale_par), fill='red', alpha=0.2)+
    ggplot2::theme(panel.grid.minor.y = ggplot2::element_blank())

  return(list(plot1, dfY))

}
