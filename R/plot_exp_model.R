# Save this file as `R/plot_exp_model.R`

#' Bayesian exponential model using stan
#'
#' @export
#' @param X date vector of data you want on the graph
#' @param Y Numeric vector of number of positive samples
#' @param N Numeric vector of total number of samples
#' @param fit_exp either fitted exponential model using reacttemporal2::stan_exp_model() or list of such models
#' @param X_model the X data for the fittd model or a list() of the X data for the fitted models (length must match number of models given for 'fit_exp')
#' @param color_list list() of color characters defaults to list("red") of same length as number of models.
#' @param ylim sets the ylimit of the plot
#' @return A ggplot of the exp model/models, the raw data and binomial CI's and the raw data of the model fits used in plotting (all in list format)
#'
plot_exp_model <- function(X, Y, N, fit_exp, X_model, color_list = as.list(rep('red', length(X_model) )), ylim = 2.0){

  lappend <- function (lst, ...){
    lst <- c(lst, list(...))
    return(lst)
  }

  plot_lin_fit <- function(ff, X){
    X = seq(min(X), max(X), by=0.2)
    dfY <- data.frame(x = X)
    Y_array <- array(data=NA, dim=c(nrow(ff$alpha), length(X)))
    for(i in seq_len(nrow(ff$alpha))){
      alpha<-ff$alpha[i]
      beta <- ff$beta[i]
      Y_array[i,] <- exp(alpha+beta*X)

    }

    alpha<-median(ff$alpha)
    beta <- median(ff$beta)

    for(i in seq_len(length(X))){
      dfY$p[i] <- exp(alpha+beta*X[i])
      dfY$ub[i] <- quantile(Y_array[,i], c(0.975))
      dfY$lb[i] <- quantile(Y_array[,i], c(0.025))
    }

    dfY
  }


  CI <- prevalence::propCI(Y,N, level=0.95, method="wilson")
  df_plot <- data.frame(X=X, p = CI$p, lb= CI$lower, ub = CI$upper)
  df_plot$d_comb <- as.Date(df_plot$X)

  if(class(X_model) == "list"){
    fit <- list()
    for(i in seq_len(length(X_model))){
      X_model[[i]] <- as.numeric(X_model[[i]])
      X_adj_temp <- min(X_model[[i]])
      X_model[[i]] <- X_model[[i]] - X_adj_temp


      ff_temp <- rstan::extract(fit_exp[[i]])

      fit_temp <- plot_lin_fit(ff_temp, X_model[[i]])
      fit_temp$x <- as.Date(fit_temp$x - 18383 + X_adj_temp, origin = as.Date("2020-05-01"))

      fit <- lappend(fit, fit_temp)
    }

    plt <- ggplot2::ggplot(data = df_plot, ggplot2::aes(x=d_comb, y=p*100))+
      ggplot2::theme_bw()+
      ggplot2::geom_point(data=df_plot, ggplot2::aes(x=d_comb, y=p*100))+
      ggplot2::geom_errorbar(data = df_plot, ggplot2::aes(x=d_comb, y=p*100, ymin=lb*100, ymax=ub*100), alpha=0.4, width=0.001)+
      ggplot2::coord_cartesian( ylim = c(0,ylim))+
      ggplot2::xlab("Date")+
      ggplot2::ylab("Prevalence (%)")+
      ggplot2::scale_x_date(date_breaks = "1 month", date_labels = "%b")

    for(i in seq_len(length(X_model))){
      plt<- plt + ggplot2::geom_line(data= fit[[i]],ggplot2::aes(y=p*100, x=x), col=color_list[[i]])+
        ggplot2::geom_ribbon(data= fit[[i]],ggplot2::aes(y=p*100, x=x,ymin=lb*100, ymax=ub*100), fill=color_list[[i]], alpha=0.2)
    }


  } else{
    X_model <- as.numeric(X_model)
    X_adj <- min(X_model)
    X_model <- X_model - X_adj


    ff <- rstan::extract(fit_exp)

    fit <- plot_lin_fit(ff, X_model)
    fit$x <- as.Date(fit$x - 18383 + X_adj, origin = as.Date("2020-05-01"))

    plt <- ggplot2::ggplot(data = df_plot, ggplot2::aes(x=d_comb, y=p*100))+
      ggplot2::theme_bw()+
      ggplot2::geom_point(data=df_plot, ggplot2::aes(x=d_comb, y=p*100))+
      ggplot2::geom_errorbar(data = df_plot, ggplot2::aes(x=d_comb, y=p*100, ymin=lb*100, ymax=ub*100), alpha=0.4, width=0.001)+
      ggplot2::coord_cartesian( ylim = c(0,ylim))+
      ggplot2::xlab("Date")+
      ggplot2::ylab("Prevalence (%)")+
      ggplot2::scale_x_date(date_breaks = "1 month", date_labels = "%b")+
      ggplot2::geom_line(data= fit, ggplot2::aes(y=p*100, x=x), col=color_list[[1]])+
      ggplot2::geom_ribbon(data= fit, ggplot2::aes(y=p*100, x=x,ymin=lb*100, ymax=ub*100), fill=color_list[[1]], alpha=0.2)
  }


  return(list(plt, df_plot, fit))

}
