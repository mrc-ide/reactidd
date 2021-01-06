# Save this file as `R/plot_exp_model_phe.R`

#' Bayesian exponential model using stan
#'
#' @export
#' @param X date vector of data you want on the graph
#' @param Y Numeric vector of number of positive samples
#' @param fit_exp either fitted exponential model using reactidd::stan_exp_model_phe() or list of such models
#' @param X_model the X data for the fittd model or a list() of the X data for the fitted models (length must match number of models)
#' @param color_list list() of color characters defaults to list("red") of same length as number of models.
#' @param ylim sets the y limit of the plot
#' @return A ggplot of the exp model/models, the raw data and binomial CI's and the raw data of the model fits used in plotting (all in list format)
#'
plot_exp_model_phe <- function(X, Y, fit_exp, X_model, color_list = as.list(rep('red', length(X_model) )), ylim = 5000.0){
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


  df_plot <- data.frame(X=X, p = Y)
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

    plt <- ggplot2::ggplot(data = df_plot, ggplot2::aes(x=d_comb, y=p))+
      ggplot2::theme_bw()+
      ggplot2::geom_point(data=df_plot, ggplot2::aes(x=d_comb, y=p))+
      ggplot2::coord_cartesian( ylim = c(0,ylim))+
      ggplot2::xlab("Date")+
      ggplot2::ylab("Number of cases")+
      ggplot2::scale_x_date(date_breaks = "1 month", date_labels = "%b")

    for(i in seq_len(length(X_model))){
      plt<- plt + ggplot2::geom_line(data= fit[[i]],ggplot2::aes(y=p, x=x), col=color_list[[i]])+
        ggplot2::geom_ribbon(data= fit[[i]],ggplot2::aes(y=p, x=x,ymin=lb, ymax=ub), fill=color_list[[i]], alpha=0.2)
    }


  } else{
    X_model <- as.numeric(X_model)
    X_adj <- min(X_model)
    X_model <- X_model - X_adj


    ff <- rstan::extract(fit_exp)

    fit <- plot_lin_fit(ff, X_model)
    fit$x <- as.Date(fit$x - 18383 + X_adj, origin = as.Date("2020-05-01"))

    plt <- ggplot2::ggplot(data = df_plot, ggplot2::aes(x=d_comb, y=p))+
      ggplot2::theme_bw()+
      ggplot2::geom_point(data=df_plot, ggplot2::aes(x=d_comb, y=p))+
      ggplot2::coord_cartesian( ylim = c(0,ylim))+
      ggplot2::xlab("Date")+
      ggplot2::ylab("Number of cases")+
      ggplot2::scale_x_date(date_breaks = "1 month", date_labels = "%b")+
      ggplot2::geom_line(data= fit, ggplot2::aes(y=p, x=x), col=color_list[[1]])+
      ggplot2::geom_ribbon(data= fit, ggplot2::aes(y=p, x=x,ymin=lb, ymax=ub), fill=color_list[[1]], alpha=0.2)
  }


  return(list(plt, df_plot, fit))

}
