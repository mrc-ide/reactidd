
# Save this file as `R/plot_p_spline_prev.R`

#' Plotting function for the P_spline model
#'
#' @export
#' @param X date vector.
#' @param Y Numeric vector of number of positive samples
#' @param N Numeric vector of total number of samples
#' @param p_splinefit fit of the model to the same set of data using reactidd::stan_p_spline()
#' @param target_dist_between_knots sets the number of days between adjacent knots (default = 5)
#' @param spline_degree sets the degree of the splines (default = 3)
#' @param ylim sets the ylimit of the plot
#' @return A list of the created plot, the raw data and CI's used in the plot, the raw data for the model fit in the plot.
#'
plot_p_spline_prev <- function(X, Y, N, p_spline_fit, target_dist_between_knots = 5, spline_degree = 3, ylim=1.0){

  inv_logit <- function(Num){
    1/(1+exp(-Num))
  }

  X_og <- X
  ff <- rstan::extract(p_spline_fit)
  X <- as.numeric(X)
  X <- seq(min(X),max(X),by=1)

  days_per_knot<-5
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
  Y_array <- array(data=NA, dim=c(nrow(ff$a), length(X)))

  #a0<-mean(ff$a0)
  for(i in seq_len(nrow(ff$a))){
    a <- array(NA, num_basis)
    #a0 <- ff$a0[i]
    for(j in seq_len(length(a))){
      a[j] <- ff$a[i,j]
      #a[j] <- mean(ff$a[,j])
    }
    Y_array[i,] <- as.vector(a%*%B_true) #as.vector(a0*X+a%*%B_true)
  }


  dfY <- data.frame(x = X)
  for(i in seq_len(length(X))){
    dfY$p[i] <-median(Y_array[,i])
    dfY$lb_2.5[i] <- quantile(Y_array[,i], probs=0.025)
    dfY$lb_25[i] <- quantile(Y_array[,i], probs=0.25)
    dfY$ub_97.5[i] <- quantile(Y_array[,i], probs=0.975)
    dfY$ub_75[i] <- quantile(Y_array[,i], probs=0.75)

  }

  df_plot_model <- dfY
  df_plot_model$d_comb <- as.Date(df_plot_model$x-18383, origin=as.Date("2020-05-01"))


  CI <- prevalence::propCI(Y, N, level=0.95, method="wilson")
  df_plot <- data.frame(X=X_og, p = CI$p, lb= CI$lower, ub = CI$upper)
  df_plot$d_comb <- as.Date(df_plot$X)

  max_date<-max(df_plot$d_comb)
  min_date<-min(df_plot$d_comb)
  #df_plot_model<-df_plot_model[df_plot_model$d_comb>=min_date & df_plot_model$d_comb<=max_date,]

  plot1 <- ggplot2::ggplot(data = df_plot, ggplot2::aes(x= d_comb, y =p*100))+
    ggplot2::geom_point()+
    ggplot2::geom_errorbar(ggplot2::aes(ymin=lb*100, ymax=ub*100))+
    ggplot2::coord_cartesian(ylim=c(0,ylim), xlim=c(min(df_plot$d_comb), max(df_plot$d_comb)))+
    ggplot2::theme_bw(base_size = 18)+
    ggplot2::xlab("Day of swab")+
    ggplot2::ylab("Prevalence (%)")+
    ggplot2::scale_x_date(date_breaks = "1 month", date_labels = "%b")+
    ggplot2::geom_line(data= df_plot_model, ggplot2::aes(y=inv_logit(p)*100))+
    ggplot2::geom_ribbon(data = df_plot_model,
                         ggplot2::aes(y=inv_logit(p)*100,
                                      ymin=inv_logit(lb_2.5)*100,
                                      ymax=inv_logit(ub_97.5)*100),
                         alpha=0.2)+
    ggplot2::geom_ribbon(data = df_plot_model,
                         ggplot2::aes(y=inv_logit(p)*100,
                                      ymin=inv_logit(lb_25)*100,
                                      ymax=inv_logit(ub_75)*100),
                         alpha=0.2)

  return(list(plot1, df_plot, df_plot_model))
}





plot_gam_fit_lineage <- function(ff, X,
                                 days_per_knot = 5,
                                 spline_degree = 3){
  min_date_numeric <- min(X)
  max_date_numeric <- max(X)
  num_knots <- ceiling((max_date_numeric- min_date_numeric)/days_per_knot)+7
  days_per_knot <- (max_date_numeric - min_date_numeric)/(num_knots -7)
  num_basis <- num_knots + spline_degree - 1
  num_data <- length(X)
  knots <- unname(seq(min(X)-3*days_per_knot, max(X)+3*days_per_knot, length.out = num_knots))

  X_new <- seq(min(X)-3*days_per_knot, max(X)+3*days_per_knot, 0.1)
  B_true <- splines::bs(X_new, df=num_basis, degree=spline_degree, intercept = TRUE)
  B_true <- t(predict(B_true, X))
  Y_array1 <- array(data=NA, dim=c(nrow(ff$a1), length(X)))
  Y_array2 <- array(data=NA, dim=c(nrow(ff$a2), length(X)))

  #a0<-mean
  for(i in seq_len(nrow(ff$a1))){
    a1 <- array(NA, num_basis)
    a2 <- array(NA, num_basis)
    #a0 <- ff$a0[i]
    for(j in seq_len(length(a1))){
      a1[j] <- ff$a1[i,j]
      a2[j] <- ff$a2[i,j]
      #a[j] <- mean(ff$a[,j])
    }
    Y_array1[i,] <- as.vector(a1%*%B_true)
    Y_array2[i,] <- as.vector(a2%*%B_true)
  }
  dfY <- data.frame(x = X)
  dfY1 <- data.frame(x = X)
  dfY2 <- data.frame(x = X)
  dfP <- data.frame(x = X)

  for(i in seq_len(length(X))){
    dfY1$p[i]   <-  median(inv_logit(Y_array1[,i]))
    dfY1$lb_2.5[i] <- quantile(inv_logit(Y_array1[,i]), probs=0.025)
    dfY1$lb_25[i] <- quantile(inv_logit(Y_array1[,i]), probs=0.25)
    dfY1$ub_97.5[i] <- quantile(inv_logit(Y_array1[,i]), probs=0.975)
    dfY1$ub_75[i] <- quantile(inv_logit(Y_array1[,i]), probs=0.75)

    dfY2$p[i]   <-  median(inv_logit(Y_array2[,i]))
    dfY2$lb_2.5[i] <- quantile(inv_logit(Y_array2[,i]), probs=0.025)
    dfY2$lb_25[i] <- quantile(inv_logit(Y_array2[,i]), probs=0.25)
    dfY2$ub_97.5[i] <- quantile(inv_logit(Y_array2[,i]), probs=0.975)
    dfY2$ub_75[i] <- quantile(inv_logit(Y_array2[,i]), probs=0.75)

    dfY$p[i]   <-  median(inv_logit(Y_array1[,i])+inv_logit(Y_array2[,i]))
    dfY$lb_2.5[i] <- quantile(inv_logit(Y_array1[,i])+inv_logit(Y_array2[,i]), probs=0.025)
    dfY$lb_25[i] <- quantile(inv_logit(Y_array1[,i])+inv_logit(Y_array2[,i]), probs=0.25)
    dfY$ub_97.5[i] <- quantile(inv_logit(Y_array1[,i])+inv_logit(Y_array2[,i]), probs=0.975)
    dfY$ub_75[i] <- quantile(inv_logit(Y_array1[,i])+inv_logit(Y_array2[,i]), probs=0.75)

    dfP$p[i]   <-  median(inv_logit(Y_array2[,i])/(inv_logit(Y_array1[,i])+inv_logit(Y_array2[,i])))
    dfP$lb_2.5[i] <- quantile(inv_logit(Y_array2[,i])/(inv_logit(Y_array1[,i])+inv_logit(Y_array2[,i])), probs=0.025)
    dfP$lb_25[i] <- quantile(inv_logit(Y_array2[,i])/(inv_logit(Y_array1[,i])+inv_logit(Y_array2[,i])), probs=0.25)
    dfP$ub_97.5[i] <- quantile(inv_logit(Y_array2[,i])/(inv_logit(Y_array1[,i])+inv_logit(Y_array2[,i])), probs=0.975)
    dfP$ub_75[i] <- quantile(inv_logit(Y_array2[,i])/(inv_logit(Y_array1[,i])+inv_logit(Y_array2[,i])), probs=0.75)




  }


  list(dfY1, dfY2, dfY, dfP)

}




plot_p_spline_prev <- function(X, Y, N,V, p_spline_fit, target_dist_between_knots = 5, spline_degree = 3, ylim=1.0,
                               labs = c("Delta","Omicron"), colors1=c("black","blue","red")){

  inv_logit <- function(Num){
    1/(1+exp(-Num))
  }

  X_og <- X
  ff <- rstan::extract(p_spline_fit)
  X <- as.numeric(X)
  X <- seq(min(X),max(X),by=1)

  min_date_numeric <- min(X)
  max_date_numeric <- max(X)
  num_knots <- ceiling((max_date_numeric- min_date_numeric)/days_per_knot)+7
  days_per_knot <- (max_date_numeric - min_date_numeric)/(num_knots -7)
  num_basis <- num_knots + spline_degree - 1
  num_data <- length(X)
  knots <- unname(seq(min(X)-3*days_per_knot, max(X)+3*days_per_knot, length.out = num_knots))

  X_new <- seq(min(X)-3*days_per_knot, max(X)+3*days_per_knot, 0.1)
  B_true <- splines::bs(X_new, df=num_basis, degree=spline_degree, intercept = TRUE)
  B_true <- t(predict(B_true, X))
  Y_array1 <- array(data=NA, dim=c(nrow(ff$a1), length(X)))
  Y_array2 <- array(data=NA, dim=c(nrow(ff$a2), length(X)))

  #a0<-mean
  for(i in seq_len(nrow(ff$a1))){
    a1 <- array(NA, num_basis)
    a2 <- array(NA, num_basis)
    #a0 <- ff$a0[i]
    for(j in seq_len(length(a1))){
      a1[j] <- ff$a1[i,j]
      a2[j] <- ff$a2[i,j]
      #a[j] <- mean(ff$a[,j])
    }
    Y_array1[i,] <- as.vector(a1%*%B_true)
    Y_array2[i,] <- as.vector(a2%*%B_true)
  }
  dfY <- data.frame(x = X)
  dfY1 <- data.frame(x = X)
  dfY2 <- data.frame(x = X)
  dfP <- data.frame(x = X)

  for(i in seq_len(length(X))){
    dfY1$p[i]   <-  median(inv_logit(Y_array1[,i]))
    dfY1$lb_2.5[i] <- quantile(inv_logit(Y_array1[,i]), probs=0.025)
    dfY1$lb_25[i] <- quantile(inv_logit(Y_array1[,i]), probs=0.25)
    dfY1$ub_97.5[i] <- quantile(inv_logit(Y_array1[,i]), probs=0.975)
    dfY1$ub_75[i] <- quantile(inv_logit(Y_array1[,i]), probs=0.75)

    dfY2$p[i]   <-  median(inv_logit(Y_array2[,i]))
    dfY2$lb_2.5[i] <- quantile(inv_logit(Y_array2[,i]), probs=0.025)
    dfY2$lb_25[i] <- quantile(inv_logit(Y_array2[,i]), probs=0.25)
    dfY2$ub_97.5[i] <- quantile(inv_logit(Y_array2[,i]), probs=0.975)
    dfY2$ub_75[i] <- quantile(inv_logit(Y_array2[,i]), probs=0.75)

    dfY$p[i]   <-  median(inv_logit(Y_array1[,i])+inv_logit(Y_array2[,i]))
    dfY$lb_2.5[i] <- quantile(inv_logit(Y_array1[,i])+inv_logit(Y_array2[,i]), probs=0.025)
    dfY$lb_25[i] <- quantile(inv_logit(Y_array1[,i])+inv_logit(Y_array2[,i]), probs=0.25)
    dfY$ub_97.5[i] <- quantile(inv_logit(Y_array1[,i])+inv_logit(Y_array2[,i]), probs=0.975)
    dfY$ub_75[i] <- quantile(inv_logit(Y_array1[,i])+inv_logit(Y_array2[,i]), probs=0.75)

    dfP$p[i]   <-  median(inv_logit(Y_array2[,i])/(inv_logit(Y_array1[,i])+inv_logit(Y_array2[,i])))
    dfP$lb_2.5[i] <- quantile(inv_logit(Y_array2[,i])/(inv_logit(Y_array1[,i])+inv_logit(Y_array2[,i])), probs=0.025)
    dfP$lb_25[i] <- quantile(inv_logit(Y_array2[,i])/(inv_logit(Y_array1[,i])+inv_logit(Y_array2[,i])), probs=0.25)
    dfP$ub_97.5[i] <- quantile(inv_logit(Y_array2[,i])/(inv_logit(Y_array1[,i])+inv_logit(Y_array2[,i])), probs=0.975)
    dfP$ub_75[i] <- quantile(inv_logit(Y_array2[,i])/(inv_logit(Y_array1[,i])+inv_logit(Y_array2[,i])), probs=0.75)




  }



  ########## XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  CI <- prevalence::propCI(Y, N, level=0.95, method="wilson")
  df_plot <- data.frame(X=X_og, p = CI$p, lb= CI$lower, ub = CI$upper)
  df_plot$d_comb <- as.Date(df_plot$X)

  max_date<-max(df_plot$d_comb)
  min_date<-min(df_plot$d_comb)


  df_plot_modelT <- dfY
  df_plot_model1 <- dfY1
  df_plot_model2 <- dfY2
  df_plot_modelP <- dfP

  df_plot_modelT$d_comb <- as.Date(df_plot_modelT$x-18383, origin=as.Date("2020-05-01"))
  df_plot_model1$d_comb <- as.Date(df_plot_model1$x-18383, origin=as.Date("2020-05-01"))
  df_plot_model2$d_comb <- as.Date(df_plot_model2$x-18383, origin=as.Date("2020-05-01"))
  df_plot_modelP$d_comb <- as.Date(df_plot_modelP$x-18383, origin=as.Date("2020-05-01"))


  CI <- prevalence::propCI(V[,3], V[,2]+V[,3], level=0.95, method="wilson")
  df_plotL <- data.frame(X=V[,1], p = CI$p, lb= CI$lower, ub = CI$upper)
  df_plotL$d_comb <- as.Date(df_plotL$X)


  plot1 <- ggplot2::ggplot()+
    ggplot2::geom_point(data = df_plot, ggplot2::aes(x=d_comb,y=p*100))+
    ggplot2::geom_errorbar(data = df_plot, ggplot2::aes(x=d_comb,y=p*100, ymin=lb*100,ymax=ub*100),width=0)+
    ggplot2::geom_line(data=df_plot_modelT, ggplot2::aes(x = d_comb, y=p*100, ymin=lb_2.5*100, ymax=ub_97.5*100,col='black'))+
    ggplot2::geom_ribbon(data=df_plot_modelT, ggplot2::aes(x = d_comb, y=p*100, ymin=lb_2.5*100, ymax=ub_97.5*100, fill='black'),alpha=0.3)+
    ggplot2::geom_line(data=df_plot_model1, ggplot2::aes(x = d_comb, y=p*100, ymin=lb_2.5*100, ymax=ub_97.5*100, col='red'))+
    ggplot2::geom_ribbon(data=df_plot_model1, ggplot2::aes(x = d_comb, y=p*100, ymin=lb_2.5*100, ymax=ub_97.5*100, fill='red'),alpha=0.3)+
    ggplot2::geom_line(data=df_plot_model2, ggplot2::aes(x = d_comb, y=p*100, ymin=lb_2.5*100, ymax=ub_97.5*100, col='blue'))+
    ggplot2::geom_ribbon(data=df_plot_model2, ggplot2::aes(x = d_comb, y=p*100, ymin=lb_2.5*100, ymax=ub_97.5*100, fill='blue'),alpha=0.3)+
    ggplot2::coord_cartesian(ylim=c(0.1,ylim), xlim=c(min_date, max_date))+
    ggplot2::theme_bw(base_size = 10)+
    ggplot2::xlab("Date (2021-2022)")+
    ggplot2::ylab("Prevalence (%)")+
    ggplot2::scale_color_manual(values=colors1,
                       name= "Variant",
                       labels=c("Total",labs))+
    ggplot2::scale_fill_manual(values=colors1,
                      name= "Variant",
                      labels=c("Total",labs))+
    ggplot2::scale_x_date(date_breaks = "1 month", date_labels = "%b\n%Y")+
    ggplot2::scale_y_log10()+
    ggplot2::theme(legend.position = "bottom")


  ##### XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  plot2 <- ggplot2::ggplot()+
    ggplot2::geom_line(data=df_plot_modelP, ggplot2::aes(x = d_comb, y=p, ymin=lb_2.5 , ymax=ub_97.5 ))+
    ggplot2::geom_ribbon(data=df_plot_modelP, ggplot2::aes(x = d_comb, y=p , ymin=lb_2.5 , ymax=ub_97.5 ),alpha=0.3)+
    ggplot2::geom_point(data=df_plotL, ggplot2::aes(x = d_comb, y=p, ymin=lb , ymax=ub ))+
    ggplot2::geom_errorbar(data=df_plotL, ggplot2::aes(x = d_comb, y=p , ymin=lb , ymax=ub ), width=0)+
    ggplot2::coord_cartesian(ylim=c(0,1), xlim=c(min_date, max_date))+
    ggplot2::theme_bw(base_size = 10)+
    ggplot2::xlab("Day of swab")+
    ggplot2::ylab(paste("Proportion ",labs[2]))+
    ggplot2::scale_x_date(date_breaks = "1 month", date_labels = "%b\n%Y")

  #df_plot_model<-df_plot_model[df_plot_model$d_comb>=min_date & df_plot_model$d_comb<=max_date,]



  return(list(plot1, plot2, df_plot_modelT,df_plot_model1, df_plot_model2, df_plot_modelP, df_plot, df_plotL))
}


