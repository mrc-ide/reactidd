
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
#' @param n shape parameter of the gamma distribution used for the generation time
#' @param b rate parameter of the gamma distribution used for the generation time
#' @return A list of the created plot of R over time and posterior estimates of R_t used in the plots.
#'
plot_p_spline_R <- function(X, stan_fit, target_dist_between_knots = 5, spline_degree = 3, ylim=2.0,
                            link_function="logit",n=2.29, b=0.36, tau_max=14){

  ff <- rstan::extract(stan_fit)

  X <- seq(min(X),max(X),by=1)

  gammaDist <- function(b, n, a) (b**n) * (a**(n-1)) * exp(-b*a) / gamma(n)


  get_quantities <- function(ff, name){
    p_delay <- quantile(ff$delay, c(0.5,0.025,0.975))
    beta1   <- (1+(quantile(ff$beta[,1], c(0.5,0.025,0.975))/0.36))**2.29
    beta2   <- (1+(quantile(ff$beta[,2], c(0.5,0.025,0.975))/0.36))**2.29
    beta3   <- (1+(quantile(ff$beta[,3], c(0.5,0.025,0.975))/0.36))**2.29
    beta4   <- (1+(quantile(ff$beta[,4], c(0.5,0.025,0.975))/0.36))**2.29
    beta5   <- (1+(quantile(ff$beta[,5], c(0.5,0.025,0.975))/0.36))**2.29
    beta6   <- (1+(quantile(ff$beta[,6], c(0.5,0.025,0.975))/0.36))**2.29

    df1 <- data.frame(name = name,period = "0", beta = beta1[[1]],beta_lb=beta1[[2]],beta_ub=beta1[[3]])
    df2 <- data.frame(name = name,period = "1", beta = beta2[[1]],beta_lb=beta2[[2]],beta_ub=beta2[[3]])
    df3 <- data.frame(name = name,period = "2", beta = beta3[[1]],beta_lb=beta3[[2]],beta_ub=beta3[[3]])
    df4 <- data.frame(name = name,period = "3", beta = beta4[[1]],beta_lb=beta4[[2]],beta_ub=beta4[[3]])
    df5 <- data.frame(name = name,period = "4", beta = beta5[[1]],beta_lb=beta5[[2]],beta_ub=beta5[[3]])
    df6 <- data.frame(name = name,period = "5", beta = beta6[[1]],beta_lb=beta6[[2]],beta_ub=beta6[[3]])

    df_R <- rbind(df1, df2, df3, df4, df5, df6)
    df_delay <- data.frame(name=name, delay = p_delay[[1]], delay_lb = p_delay[[2]], delay_ub = p_delay[[3]])

    test<-(1+(ff$beta[,2]/0.36))**2.29/(1+(ff$beta[,1]/0.36))**2.29
    test[ff$beta[,2]< -0.36] <- 0.0
    test[ff$beta[,1]< -0.36] <- 10000
    beta1   <- quantile( test, c(0.5,0.025,0.975))

    test<-(1+(ff$beta[,3]/0.36))**2.29/(1+(ff$beta[,2]/0.36))**2.29
    test[ff$beta[,3]< -0.36] <- 0.0
    test[ff$beta[,2]< -0.36] <- 10000
    beta2   <- quantile( test, c(0.5,0.025,0.975))

    test<-(1+(ff$beta[,4]/0.36))**2.29/(1+(ff$beta[,3]/0.36))**2.29
    test[ff$beta[,4]< -0.36] <- 0.0
    test[ff$beta[,3]< -0.36] <- 10000
    beta3   <- quantile( test, c(0.5,0.025,0.975))

    test<-(1+(ff$beta[,5]/0.36))**2.29/(1+(ff$beta[,4]/0.36))**2.29
    test[ff$beta[,5]< -0.36] <- 0.0
    test[ff$beta[,4]< -0.36] <- 10000
    beta4   <- quantile( test, c(0.5,0.025,0.975))

    test<-(1+(ff$beta[,6]/0.36))**2.29/(1+(ff$beta[,5]/0.36))**2.29
    test[ff$beta[,6]< -0.36] <- 0.0
    test[ff$beta[,5]< -0.36] <- 10000
    beta5   <- quantile( test, c(0.5,0.025,0.975))

    df1 <- data.frame(name = name,period = "1", beta = beta1[[1]],beta_lb=beta1[[2]],beta_ub=beta1[[3]])
    df2 <- data.frame(name = name,period = "2", beta = beta2[[1]],beta_lb=beta2[[2]],beta_ub=beta2[[3]])
    df3 <- data.frame(name = name,period = "3", beta = beta3[[1]],beta_lb=beta3[[2]],beta_ub=beta3[[3]])
    df4 <- data.frame(name = name,period = "4", beta = beta4[[1]],beta_lb=beta4[[2]],beta_ub=beta4[[3]])
    df5 <- data.frame(name = name,period = "5", beta = beta5[[1]],beta_lb=beta5[[2]],beta_ub=beta5[[3]])

    df_mult <- rbind(df1, df2, df3, df4, df5)

    list(df_R, df_delay, df_mult)

  }

  eng <- get_quantities(ff, name = "England")

  NW <- get_quantities(ffNW, name = "North West")
  NE <- get_quantities(ffNE, name = "North East")
  SW <- get_quantities(ffSW, name = "South West")
  SE <- get_quantities(ffSE, name = "South East")
  LN <- get_quantities(ffLN, name = "London")
  EM <- get_quantities(ffEM, name = "East Midlands")
  WM <- get_quantities(ffWM, name = "West Midlands")
  EE <- get_quantities(ffEE, name = "East of England")
  YH <- get_quantities(ffYH, name = "Yorkshire and The Humber")

  age1 <- get_quantities(ff1, name = "5-17")
  age2 <- get_quantities(ff2, name = "18-34")
  age3 <- get_quantities(ff3, name = "35-54")
  age4 <- get_quantities(ff4, name = "55+")

  df_R <- rbind(eng[[1]],
                NW[[1]], NE[[1]], SE[[1]],
                SW[[1]], LN[[1]], EE[[1]],
                EM[[1]], WM[[1]], YH[[1]],
                age1[[1]], age2[[1]], age3[[1]], age4[[1]])

  df_delay <- rbind(eng[[2]],
                    NW[[2]], NE[[2]], SE[[2]],
                    SW[[2]], LN[[2]], EE[[2]],
                    EM[[2]], WM[[2]], YH[[2]],
                    age1[[2]], age2[[2]], age3[[2]], age4[[2]])
  df_R[is.nan(df_R$beta_lb),]$beta_lb <- 0


  df_mult <- rbind(eng[[3]],
                   NW[[3]], NE[[3]], SE[[3]],
                   SW[[3]], LN[[3]], EE[[3]],
                   EM[[3]], WM[[3]], YH[[3]],
                   age1[[3]], age2[[3]], age3[[3]], age4[[3]])

  df_mult

  ggplot(data = df_R, aes(x=period, y=beta, ymin=beta_lb, ymax=beta_ub))+
    geom_point(position = position_dodge(0.9), aes(color=name))+
    geom_errorbar(width=0, position = position_dodge(0.9), aes(color=name))+
    coord_cartesian(ylim=c(0,2))+
    geom_hline(yintercept = 1.0, linetype='dashed',color='black')+
    theme_bw()



  ggplot(data = df_mult, aes(x=period, y=beta, ymin=beta_lb, ymax=beta_ub))+
    geom_point(position = position_dodge(0.9), aes(color=name))+
    geom_errorbar(width=0, position = position_dodge(0.9), aes(color=name))+
    coord_cartesian(ylim=c(0,5))+
    geom_hline(yintercept = 1.0, linetype='dashed',color='black')+
    theme_bw()


  test<-(1+(ffNE$beta[,5]/0.36))**2.29/(1+(ffNE$beta[,4]/0.36))**2.29
  test[ffNE$beta[,5]< -0.36]
  test[ffNE$beta[,4]< -0.36]

  beta4   <- quantile( (1+(ffNE$beta[,5]/0.36))**2.29/(1+(ffNE$beta[,4]/0.36))**2.29, c(0.5,0.025,0.975))





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

    Rlist <- (1+(log(ff$beta_t[,i])/0.36))^2.29
    row_tempbeta <- data.frame(X=X_all[i]+min_X,
                               beta = quantile(Rlist, c(0.5,0.025,0.975))[[1]],
                               lb1= quantile(Rlist, c(0.5,0.025,0.975))[[2]],
                               ub1= quantile(Rlist, c(0.5,0.025,0.975))[[3]],
                               lb2= quantile(Rlist, c(0.5,0.25,0.75))[[2]],
                               ub2= quantile(Rlist, c(0.5,0.25,0.75))[[3]])

    Rlist <- (1+(log(ff$beta_t[,i])/0.36))^2.29
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

  library(ggplot2)
  plot1<-ggplot(data = df_plot, aes(x= d_comb, y =p*100))+
    geom_point()+
    geom_errorbar(aes(ymin=lb*100, ymax=ub*100))+
    coord_cartesian(ylim=c(0.01,3))+
    theme_bw(base_size = 18)+
    xlab("Day of swab")+
    ylab("Prevalence (%)")+
    scale_x_date(date_breaks = "1 month", date_labels = "%b\n%Y")+
    geom_line(data= df_plot_model, aes(x=d_comb,y=p))+
    geom_ribbon(data= df_plot_model, aes(x=d_comb,y=p, ymin=lb1, ymax=ub1), alpha=0.2)+
    geom_ribbon(data= df_plot_model, aes(x=d_comb,y=p, ymin=lb2, ymax=ub2), alpha=0.2)+
    scale_y_log10()
  plot1


  plot2<-ggplot(data=df_beta, aes(x= d_comb, y =p*100))+
    coord_cartesian(ylim=c(0.0,3))+
    theme_bw(base_size = 18)+
    xlab("Date")+
    ylab("Reproduction number")+
    scale_x_date(date_breaks = "1 month", date_labels = "%b\n%Y")+
    geom_line(data= df_beta, aes(x=d_comb,y=beta))+
    geom_ribbon(data= df_beta, aes(x=d_comb,y=beta, ymin=lb1, ymax=ub1), alpha=0.2)+
    geom_ribbon(data= df_beta, aes(x=d_comb,y=beta, ymin=lb2, ymax=ub2), alpha=0.2)+
    geom_hline(yintercept = 1, linetype='dashed',color='black')
  plot2

  res0 <- as.Date("2021-01-06")
  res1 <- as.Date("2021-03-08")
  res2 <- as.Date("2021-03-29")
  res3 <- as.Date("2021-04-12")
  res4 <- as.Date("2021-05-17")
  eng[[3]]$d_comb <- c(res0, res1, res2, res3, res4)+mean(ff$delay)

  plot2 <- plot2+
    geom_vline(xintercept = rest_dates[[1]], linetype="dashed")+
    geom_vline(xintercept = rest_dates[[2]], linetype="dashed")+
    geom_vline(xintercept = rest_dates[[3]], linetype="dashed")+
    geom_vline(xintercept = rest_dates[[4]], linetype="dashed")+
    geom_vline(xintercept = rest_dates[[5]], linetype="dashed")+
    geom_point(data=eng[[3]], aes(y=beta, ymin=beta_lb, ymax=beta_ub, color=period))+
    geom_errorbar(data=eng[[3]], aes(y=beta, ymin=beta_lb, ymax=beta_ub, color=period), width=0.1)+
    scale_color_brewer(palette="Set1",labels=c("Lockdown","Step 1a","Step 1b","Step 2","Step 3"))+
    theme(legend.position = c(0.83,0.73),
          legend.title = element_blank(),
          panel.grid = element_blank())+
    scale_y_continuous(sec.axis = sec_axis(~., name = "Multiplicative growth"))+
    coord_cartesian(ylim=c(0,4),xlim=c(as.Date("2020-12-30"), as.Date("2021-07-12")))+
    scale_x_date(date_breaks = "1 month", date_labels = "%b\n%Y")

  plot2


  return(list(plot1, df_plot_model))
}
