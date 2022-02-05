# Save this file as `R/stan_p_spline.R`

#' Bayesian p-spline model using stan
#'
#' @export
#' @param X date vector.
#' @param Y Numeric vector of number of positive samples
#' @param N Numeric vector of total number of samples
#' @param target_distance_between_knots sets the number of days between adjacent knots (default = 5)
#' @param spline_degree sets the degree of the splines (default = 3)
#' @param ... Arguments passed to `rstan::sampling` (iter, warmup).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'


shedding_exp_model_plot_density <- function(mod,label_name, ylims = c(0,1), xlims=c(0,40)){

  ff <- rstan::extract(mod)
  dfplot <- data.frame(sens=ff$sens,
                       k   =ff$k)

  plot1<- ggplot2::ggplot(data = dfplot, ggplot2::aes(x=1/k, y=sens))+
    ggplot2::geom_bin2d()+
    ggplot2::scale_fill_continuous(low=ggplot2::alpha("red",0.1),high=ggplot2::alpha("red",1))+
    ggplot2::theme_bw()+
    ggplot2::ylab("Sensitivity")+
    ggplot2::xlab("Mean duration positive")+
    ggplot2::theme(legend.position = "None")+
    ggplot2::coord_cartesian(xlim=xlims,ylim=ylims)+
    ggplot2::ggtitle(label_name)

  print(plot1)
}
