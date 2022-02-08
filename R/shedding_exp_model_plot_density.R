# Save this file as `R/stan_p_spline.R`

#' Bayesian p-spline model using stan
#'
#' @export
#' @param mod object of class `rstan::sampling` from exponetial model fit to shedding data
#' @param label_name name to label the figue produced with.
#' @param ... edit xlims or ylims to change bounds of the plot
#' @return A plot returned by `ggplot2::ggplot`
#'


shedding_exp_model_plot_density <- function(mod,label_name, ylims = c(0.5,1), xlims=c(0,30)){

  ff <- rstan::extract(mod)

  if(length(ff$t_d) == 0){
    dfplot <- data.frame(sens=ff$sens,
                         k   =ff$k,
                         td  = 0)

  } else{
    dfplot <- data.frame(sens=ff$sens,
                         k   =ff$k,
                         td  = ff$t_d)

  }


  plot1<- ggplot2::ggplot(data = dfplot, ggplot2::aes(x=td+1/k, y=sens))+
    ggplot2::geom_bin2d()+
    ggplot2::scale_fill_continuous(low=ggplot2::alpha("red",0.1),high=ggplot2::alpha("red",1))+
    ggplot2::theme_bw()+
    ggplot2::ylab("Sensitivity")+
    ggplot2::xlab("Mean duration positive")+
    ggplot2::theme(legend.position = "None")+
    ggplot2::coord_cartesian(xlim=xlims,ylim=ylims)+
    ggplot2::ggtitle(label_name)

  return(plot1)
}
