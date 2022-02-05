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


shedding_plot_exp_models <- function(dat, mod1, mod2, label_name, tmax =20){

  # Define functions for simulating models for plotting
  get_response1 <- function(tmax, ff){
    time <- seq(0, tmax, by=0.1)

    Y_array <- array(data=NA, dim=c(nrow(ff$k), length(time)))

    for(i in seq_len(nrow(Y_array))){
      k <- ff$k[i]
      sens <- ff$sens[i]

      Y_array[i,] <- sens*exp(-k*time)
    }

    dfY <- data.frame(time = time)
    for(i in seq_len(length(time))){
      dfY$p[i] <-median(Y_array[,i])
      dfY$lb_2.5[i] <- quantile(Y_array[,i], probs=0.025)
      dfY$ub_97.5[i] <- quantile(Y_array[,i], probs=0.975)
    }

    dfY

  }


  get_response2 <- function(tmax, ff){
    timer <- seq(0, tmax, by=0.1)

    Y_array <- array(data=NA, dim=c(nrow(ff$k), length(timer)))

    for(i in seq_len(nrow(Y_array))){
      k <- ff$k[i]
      sens <- ff$sens[i]
      t_d <- ff$t_d[i]

      Y_array[i,] <- sens*exp(-k*(timer-t_d))
      Y_array[i,Y_array[i,]>sens] <- sens
    }

    dfY <- data.frame(time = timer)
    for(i in seq_len(length(timer))){
      dfY$p[i] <-median(Y_array[,i])
      dfY$lb_2.5[i] <- quantile(Y_array[,i], probs=0.025)
      dfY$ub_97.5[i] <- quantile(Y_array[,i], probs=0.975)
    }

    dfY

  }

  ff1 <- rstan::extract(mod1)
  ff2 <- rstan::extract(mod2)
  dat1<-dat[c("time_shed1","estbinres_shed1")]
  dat2<-dat[c("time_shed2","estbinres_shed2")]
  colnames(dat1) <- c("time", "estbinres")
  colnames(dat2) <- c("time", "estbinres")
  new_dat <- rbind(dat1, dat2)
  tab <- table(new_dat$estbinres, new_dat$time)
  prop<-prevalence::propCI(tab[2,], tab[1,]+tab[2,], method = "wilson")
  prop <- data.frame(prop)
  prop$time <- as.numeric(rownames(prop))


  model1 <- get_response1(tmax, ff1)
  model2 <- get_response2(tmax, ff2)

  plot1<-ggplot2::ggplot()+
    ggplot2::geom_point(data=prop, ggplot2::aes(x=time,y=p))+
    ggplot2::geom_errorbar(data=prop, ggplot2::aes(x=time, ymin=lower,ymax=upper))+
    ggplot2::geom_ribbon(data=model2, ggplot2::aes(x=time, y=p, ymin=lb_2.5, ymax=ub_97.5, fill = 'red'), alpha=0.4)+
    ggplot2::geom_ribbon(data=model1, ggplot2::aes(x=time, y=p, ymin=lb_2.5, ymax=ub_97.5, fill = 'blue'), alpha=0.4)+
    ggplot2::geom_line(data=model2, ggplot2::aes(x=time, y=p ,color = 'red'))+
    ggplot2::geom_line(data=model1, ggplot2::aes(x=time, y=p ,color = 'blue'))+
    ggplot2::theme_bw(base_size = 18)+
    ggplot2::coord_cartesian(xlim = c(0,tmax), ylim=c(0,1))+
    ggplot2::xlab("Time since first test")+
    ggplot2::ylab("Proportion positive")+
    ggplot2::theme(legend.position = "none")+
    ggplot2::geom_point(ggplot2::aes(x=0, y=1), shape=4)+
    ggplot2::ggtitle(label_name)+
    ggplot2::scale_color_brewer(palette = "Dark2")+
    ggplot2::scale_fill_brewer(palette = "Dark2")

  return(plot1)
}

