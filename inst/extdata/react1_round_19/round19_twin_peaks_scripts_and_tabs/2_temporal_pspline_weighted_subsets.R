#!/usr/bin/env Rscript
Sys.time()
rm(list=ls(all=TRUE))
setwd("E:/Group/Ollies_workspace/")

path_name <- "E://Group/figure/temporal/MCMC_model_gam/National/England/Age/5-12/"
dir.create(path_name)
#' Then load required libraries and locally defined functions
library(prevalence)
library(mgcv)
library(stats)
library(cowplot)
library(rstan)
library(scales)
library(splines)
library(ggplot2)
library(dplyr)
source("E:/group/add_conf_ints.R")
source("./bayesian_gam_functions.R")

#' ## Read in the data and organise a little 
dfRes1<-readRDS("E:/group/saved_objects/rep1.rds")
dfRes2<-readRDS("E:/group/saved_objects/rep2.rds")
dfRes3<-readRDS("E:/group/saved_objects/rep3.rds")
dfRes4<-readRDS("E:/group/saved_objects/rep4.rds")
dfRes5<-readRDS("E:/group/saved_objects/rep5.rds")
dfRes6<-readRDS("E:/group/saved_objects/rep6.rds")
dfRes7<-readRDS("E:/group/saved_objects/rep7.rds")
dfRes8<-readRDS("E:/group/saved_objects/rep8.rds")
dfRes9<-readRDS("E:/group/saved_objects/rep9.rds")
dfRes10<-readRDS("E:/group/saved_objects/rep10.rds")
dfRes11<-readRDS("E:/group/saved_objects/rep11.rds")
dfRes12<-readRDS("E:/group/saved_objects/rep12.rds")
dfRes13<-readRDS("E:/group/saved_objects/rep13.rds")
dfRes14<-readRDS("E:/group/saved_objects/rep14.rds")
dfRes15<-readRDS("E:/group/saved_objects/rep15.rds")
dfRes16<-readRDS("E:/group/saved_objects/rep16.rds")
dfRes17<-readRDS("E:/group/saved_objects/rep17.rds")
dfRes18<-readRDS("E:/group/saved_objects/rep18.rds")
dfRes19<-readRDS("E:/group/saved_objects/rep19.rds")


dfRes1 <- dfRes1[colnames(dfRes1) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","age_group_char","age","hh_size","wt_antigen","nchild")]
dfRes2 <- dfRes2[colnames(dfRes2) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","age_group_char","age","hh_size","wt_antigen","nchild")]
dfRes3 <- dfRes3[colnames(dfRes3) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","age_group_char","age","hh_size","wt_antigen","nchild")]
dfRes4 <- dfRes4[colnames(dfRes4) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","age_group_char","age","hh_size","wt_antigen","nchild")]
dfRes5 <- dfRes5[colnames(dfRes5) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","age_group_char","age","hh_size","wt_antigen","nchild")]
dfRes6 <- dfRes6[colnames(dfRes6) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","age_group_char","age","hh_size","wt_antigen","nchild")]
dfRes7 <- dfRes7[colnames(dfRes7) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","age_group_char","age","hh_size","wt_antigen","nchild")]
dfRes8 <- dfRes8[colnames(dfRes8) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","age_group_char","age","hh_size","wt_antigen","nchild")]
dfRes9 <- dfRes9[colnames(dfRes9) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","age_group_char","age","hh_size","wt_antigen","nchild")]
dfRes10<- dfRes10[colnames(dfRes10) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","age_group_char","age","hh_size","wt_antigen","nchild")]
dfRes11<- dfRes11[colnames(dfRes11) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","age_group_char","age","hh_size","wt_antigen","nchild")]
dfRes12<- dfRes12[colnames(dfRes12) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","age_group_char","age","hh_size","wt_antigen","nchild")]
dfRes13<- dfRes13[colnames(dfRes13) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","age_group_char","age","hh_size","wt_antigen","nchild")]
dfRes14<- dfRes14[colnames(dfRes14) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","age_group_char","age","hh_size","wt_antigen","nchild")]
dfRes15<- dfRes15[colnames(dfRes15) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","age_group_char","age","hh_size","wt_antigen","nchild")]
dfRes16<- dfRes16[colnames(dfRes16) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","age_group_char","age","hh_size","wt_antigen","nchild")]
dfRes17<- dfRes17[colnames(dfRes17) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","age_group_char","age","hh_size","wt_antigen","nchild")]
dfRes18<- dfRes18[colnames(dfRes18) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","age_group_char","age","hh_size","wt_antigen","nchild")]
dfRes19<- dfRes19[colnames(dfRes19) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","age_group_char","age","hh_size","wt_antigen","nchild")]


dfRes1$round <- 1
dfRes2$round <- 2
dfRes3$round <- 3
dfRes4$round <- 4
dfRes5$round <- 5
dfRes6$round <- 6
dfRes7$round <- 7
dfRes8$round <- 8
dfRes9$round <- 9
dfRes10$round <- 10
dfRes11$round <- 11
dfRes12$round <- 12
dfRes13$round <- 13
dfRes14$round <- 14
dfRes15$round <- 15
dfRes16$round <- 16
dfRes17$round <- 17
dfRes18$round <- 18
dfRes19$round <- 19


#' Until we sort out the ensemble need this
dfRes1$d_swab <- dfRes1$d_comb
dfRes2$d_swab <- dfRes2$d_comb
dfRes3$d_swab <- dfRes3$d_comb
dfRes4$d_swab <- dfRes4$d_comb
dfRes5$d_swab <- dfRes5$d_comb
dfRes6$d_swab <- dfRes6$d_comb
dfRes7$d_swab <- dfRes7$d_comb
dfRes8$d_swab <- dfRes8$d_comb
dfRes9$d_swab <- dfRes9$d_comb
dfRes10$d_swab <- dfRes10$d_comb
dfRes11$d_swab <- dfRes11$d_comb
dfRes12$d_swab <- dfRes12$d_comb
dfRes13$d_swab <- dfRes13$d_comb
dfRes14$d_swab <- dfRes14$d_comb
dfRes15$d_swab <- dfRes15$d_comb
dfRes16$d_swab <- dfRes16$d_comb
dfRes17$d_swab <- dfRes17$d_comb
dfRes18$d_swab <- dfRes18$d_comb
dfRes19$d_swab <- dfRes19$d_comb


dfRes <- rbind(dfRes1, dfRes2, dfRes3, dfRes4, dfRes5, 
               dfRes6, dfRes7, dfRes8, dfRes9, dfRes10,
               dfRes11, dfRes12, dfRes13, dfRes14, dfRes15,
               dfRes16, dfRes17, dfRes18, dfRes19)

p_spline_function <- function(dfRes,
                              region_name = "England",
                              age_category = "All",
                              positive_variable = "Estbinres",
                              iterations = 20000,
                              warmup = 2000,
                              plot_R = FALSE,
                              path_name,
                              name){
  #dir.create(paste("E://Group/figure/temporal/MCMC_model_gam/National/",positive_variable,"/", sep=""))
  #dir.create(paste("E://Group/figure/temporal/MCMC_model_gam/National/",positive_variable,"/",region_name,"/", sep=""))
  #dir.create(paste("E://Group/figure/temporal/MCMC_model_gam/National/",positive_variable,"/",region_name,"/",age_category,"/", sep=""))
  #path_name <- paste("E://Group/figure/temporal/MCMC_model_gam/National/",positive_variable,"/",region_name,"/",age_category,"/", sep="")

  
  if(positive_variable!= "Estbinres"){
    if(positive_variable=="Asympt"){
      dfRes$estbinres <- dfRes$res_asympt
    }
    if(positive_variable=="double_pos"){
      dfRes$estbinres <- dfRes$res_ct12_gt0
    }
    if(positive_variable=="Estbinres_35"){
      dfRes$estbinres <- dfRes$estbinres35
    }
  }
  
  if(region_name != "England"){
    dfRes <- dfRes[dfRes$region == region_name,]
  }
  
  if(age_category != "All"){
    dfRes <- dfRes[dfRes$age_group_char == age_category,]
  }
  
  dfRes<-dfRes[is.na(dfRes$estbinres)==FALSE &is.na(dfRes$d_comb)==FALSE,]
  ensemble <- data.frame(day=unique(dfRes$d_comb))
  
  for(i in seq_len(nrow(ensemble))){
    ensemble$obs[i] <- sum(dfRes[dfRes$d_comb ==ensemble$day[i],]$wt_antigen, na.rm=TRUE)
    ensemble$pos[i] <- sum(dfRes[dfRes$d_comb ==ensemble$day[i] &dfRes$estbinres==1,]$wt_antigen, na.rm=TRUE)
  }
  ensemble <- ensemble[is.na(ensemble$day)==FALSE,]
  ensemble <-ensemble[order(ensemble$day),]
  
  #' Convert date to numeric 
  ensemble$day <- as.numeric(ensemble$day)
  min_date_numeric <- min(ensemble$day)
  max_date_numeric <- max(ensemble$day)
  
  #' Set initial conditions and variable names to pass to stan
  days_per_knot <- 5
  spline_degree <- 3
  
  num_knots <- ceiling((max_date_numeric- min_date_numeric)/days_per_knot)+7
  days_per_knot <- (max_date_numeric - min_date_numeric)/(num_knots -7)

  num_basis <- num_knots + spline_degree - 1
  
  X <- ensemble$day
  num_data <- length(X)
  knots <- unname(seq(min(X)-3*days_per_knot, max(X)+3*days_per_knot, length.out = num_knots))
  
  Y <- as.numeric(ensemble$pos)
  N <- as.numeric(ensemble$obs)
  
  #' Load and run stan model
  rstan_options(auto_write = TRUE)
  options(mc.cores = 4)
  spline_model <- stan_model("stan_models/b_splines_weighted.stan")
  fit_spline_test <- sampling(spline_model,
                              iter=iterations,
                              warmup =warmup,
                              chains=4,
                              control = list(adapt_delta=0.95,
                                             max_treedepth = 10),
                              data = list(num_data = num_data,
                                          num_knots = num_knots,
                                          knots = knots,
                                          Y = Y,
                                          N =N,
                                          X = X,
                                          spline_degree = spline_degree))
  
  
  saveRDS(fit_spline_test, paste('E:/Group/Ollies_workspace/stan_fits/',region_name,age_category,'.rds',sep='_'))
  fit_spline <- fit_spline_test
  ff<-rstan::extract(fit_spline)
  #' Quick check of mcmc

  ################# Plotting #######################################
  ##################Plot Prevalence ################################
  X_new <- seq(min(X)-3*days_per_knot, max(X)+3*days_per_knot, 1)
  df_plot_model<-plot_gam_fit(ff, X_new, num_basis, spline_degree)
  df_plot_model$d_comb <- as.Date(df_plot_model$x-18383, origin=as.Date("2020-05-01"))
  
  tabDaySwab <- table(dfRes$estbinres,dfRes$d_comb)
  tabDaySwab[1,] <- N-Y
  tabDaySwab[2,] <- Y
  dayPCIs <- add_conf_ints(tabDaySwab,poscol="1",negcol="0")
  df_plot<-as.data.frame(dayPCIs)
  df_plot <- df_plot[is.na(df_plot$p)==FALSE,]
  df_plot$d_comb <- as.Date(rownames(df_plot))
  max_date<-max(df_plot$d_comb)
  min_date<-min(df_plot$d_comb)
  
  df_plot_model<-df_plot_model[df_plot_model$d_comb>=min_date & df_plot_model$d_comb<=max_date,]
  plot1 <- ggplot(data = df_plot, aes(x= d_comb, y =p*100))+
    geom_point()+
    geom_errorbar(aes(ymin=lb*100, ymax=ub*100))+
    coord_cartesian(ylim=c(0,2.5), xlim=c(min(df_plot$d_comb), max(df_plot$d_comb)))+
    theme_bw(base_size = 18)+ 
    xlab("Day of swab")+
    ylab("Prevalence (%)")+
    scale_x_date(date_breaks = "1 month", date_labels = "%b")+
    geom_line(data= df_plot_model, aes(y=inv_logit(p)*100))+
    geom_ribbon(data = df_plot_model,
                aes(y=inv_logit(p)*100,
                    ymin=inv_logit(lb_2.5)*100,
                    ymax=inv_logit(ub_97.5)*100),
                alpha=0.2)+
    geom_ribbon(data = df_plot_model,
                aes(y=inv_logit(p)*100,
                    ymin=inv_logit(lb_25)*100,
                    ymax=inv_logit(ub_75)*100),
                alpha=0.2)
  
  png(paste(path_name, "Prevalence.png", sep=""), width = 2000, height =500)
  print(plot1)
  dev.off()
  pdf(paste(path_name, "Prevalence.pdf", sep=""), width = 12, height =6)
  print(plot1)
  dev.off()
  
  ####################### Plot Log of Prevalence ########################
  plot2 <- ggplot(data = df_plot, aes(x= d_comb, y =log10(p )))+
    geom_point()+
    geom_errorbar(aes(ymin=log10(lb ), ymax=log10(ub )))+
    coord_cartesian(ylim=c(-4,-1), xlim=c(min(df_plot$d_comb), max(df_plot$d_comb)))+
    theme_bw(base_size = 18)+ 
    xlab("Day of swab")+
    ylab("Log10(Prevalence)")+
    scale_x_date(date_breaks = "1 month", date_labels = "%b")+
    geom_line(data= df_plot_model, aes(y=log10(inv_logit(p) )))+
    geom_ribbon(data = df_plot_model,
                aes(y=log10(inv_logit(p) ),
                    ymin=log10(inv_logit(lb_2.5) ),
                    ymax=log10(inv_logit(ub_97.5) )),
                alpha=0.2)+
    geom_ribbon(data = df_plot_model,
                aes(y=log10(inv_logit(p) ),
                    ymin=log10(inv_logit(lb_25) ),
                    ymax=log10(inv_logit(ub_75) )),
                alpha=0.2)
  
  png(paste(path_name, "Prevalence_log10.png", sep=""), width = 2000, height =500)
  print(plot2)
  dev.off()
  pdf(paste(path_name, "Prevalence_log10.pdf", sep=""), width = 16, height =6)
  print(plot2)
  dev.off()
  
  ################### Growth Rate Plot ############################
  #################################################################
  #Calculate growth rate 
  df_plot_r <- plot_growth_rate(ff,X_new,num_basis, spline_degree, link_function = "logit")
  df_plot_r$d_comb <- as.Date(df_plot_r$x-18383, origin=as.Date("2020-05-01"))
  
  
  df_plot_r<-df_plot_r[df_plot_r$d_comb>=min_date & df_plot_r$d_comb<=max_date,]
  df_plot_r$r1 <- df_plot_r$r
  df_plot_r$r2 <- df_plot_r$r
  df_plot_r[df_plot_r$r1>0.0,]$r1 <- NA
  df_plot_r[df_plot_r$r2<0.0,]$r2 <- NA
  df_plot_r$middle <- 0.0
  second_axis_values <- c(log(2)/5, log(2)/10,log(2)/15, 0, log(2)/-15, log(2)/-10, log(2)/-5)
  second_axis_labels <- c("5", "10", "15",expression(infinity / -infinity),"-15","-10","-5")
  
  plot3 <- ggplot(data = df_plot_r, aes(x=d_comb,y=r))+
    geom_line(aes(y=r1), col = '#009900')+
    geom_line(aes(y=r2), col = 'red')+
    geom_ribbon(aes(ymin=pmin(lb_2.5,middle), ymax=pmin(lb_25,middle)),alpha=0.2, fill='#009900')+
    geom_ribbon(aes(ymin=pmin(lb_25,middle), ymax=pmin(ub_75,middle)),alpha=0.4, fill='#009900')+
    geom_ribbon(aes(ymin=pmin(ub_75,middle), ymax=pmin(ub_97.5,middle)),alpha=0.2, fill='#009900')+
    geom_ribbon(aes(ymin=pmax(lb_2.5,middle), ymax=pmax(lb_25,middle)),alpha=0.2, fill='red')+
    geom_ribbon(aes(ymin=pmax(lb_25,middle), ymax=pmax(ub_75,middle)),alpha=0.4, fill="red")+
    geom_ribbon(aes(ymin=pmax(ub_75,middle), ymax=pmax(ub_97.5,middle)),alpha=0.2, fill='red')+
    theme_bw(base_size = 18)+
    coord_cartesian(ylim=c(-0.15,0.15),xlim=c(min(df_plot$d_comb), max(df_plot$d_comb)))+
    geom_line(y=0, linetype='dashed')+
    xlab("Day of swab test")+
    ylab("Instantaenous growth rate")+
    scale_x_date(date_breaks = "1 month", date_labels = "%b")+
    scale_y_continuous(sec.axis = sec_axis(~.,labels = second_axis_labels, breaks= second_axis_values, name = "Doubling(+) / Halving(-) Time (days)"))
  
  
  png(paste(path_name, "growth_rate_cool.png", sep=""), width = 2000, height =500)
  print(plot3)
  dev.off()
  pdf(paste(path_name, "growth_rate_cool.pdf", sep=""), width = 12, height =6)
  print(plot3)
  dev.off()
  
  #################Plotting Reproduction Number ##############################################
  if(plot_R == TRUE){
    
    df_plot_R <- plot_R(ff,X_new,num_basis,spline_degree,link_function = "logit",n=2.29, b=0.36, tau_max=14)
    df_plot_R$d_comb <- as.Date(df_plot_R$x-18383, origin=as.Date("2020-05-01"))
    
    tau_max=14
    df_plot_R<-df_plot_R[df_plot_R$d_comb>=min_date+tau_max & df_plot_R$d_comb<=max_date,]
    df_plot_R[is.na(df_plot_R$r)==TRUE,]$prob <- NA
    
    plot4 <- ggplot(data = df_plot_R, aes(x=d_comb,y=r))+
      geom_line()+
      geom_ribbon(aes(ymin=lb_2.5, ymax=ub_97.5), alpha=0.2)+
      geom_ribbon(aes(ymin=lb_25, ymax=ub_75), alpha=0.2)+
      theme_bw()+
      coord_cartesian(ylim=c(0,2), xlim=c(min(df_plot$d_comb), max(df_plot$d_comb)))+
      geom_line(y=1, linetype='dashed')+
      geom_line(aes(y=prob), color='red')+
      xlab("Day of swab test")+
      scale_x_date(date_breaks = "1 month", date_labels = "%b")+
      ylab("Reproductive number R_t")
    
    png(paste(path_name, "ReproductionNumber.png", sep=""), width = 2000, height =500)
    plot4
    dev.off()
    pdf(paste(path_name, "ReproductionNumber.pdf", sep=""), width = 12, height =6)
    plot4
    dev.off()
    
    return_list<-list(ff, X_new, num_basis, spline_degree, df_plot, df_plot_model, df_plot_r, df_plot_R)
    
  } else{
    return_list<-list(ff, X_new, num_basis, spline_degree, df_plot, df_plot_model, df_plot_r)
  }
  return_list
}

#################################################################################
#path_name <- paste("E://Group/figure/temporal/MCMC_model_gam/National/",positive_variable,"/",region_name,"/",age_category,"/", sep="")
reg_NW <- p_spline_function(dfRes,
                           region_name = "North West",
                           age_category = "All",
                           positive_variable = "Estbinres",
                           iterations = 20000,
                           warmup = 5000,
                           path_name = "E://Group/figure/temporal/MCMC_model_gam/National/Estbinres/North West/All/",
                           name = "over_55")
saveRDS(reg_NW, file = "E:/Group/Ollies_workspace/pspline_func_runs/reg_NW.rds")

reg_NE <- p_spline_function(dfRes,
                            region_name = "North East",
                            age_category = "All",
                            positive_variable = "Estbinres",
                            iterations = 20000,
                            warmup = 5000,
                            path_name = "E://Group/figure/temporal/MCMC_model_gam/National/Estbinres/North East/All/",
                            name = "over_55")
saveRDS(reg_NE, file = "E:/Group/Ollies_workspace/pspline_func_runs/reg_NE.rds")

reg_SW <- p_spline_function(dfRes,
                            region_name = "South West",
                            age_category = "All",
                            positive_variable = "Estbinres",
                            iterations = 20000,
                            warmup = 5000,
                            path_name = "E://Group/figure/temporal/MCMC_model_gam/National/Estbinres/South West/All/",
                            name = "over_55")
saveRDS(reg_SW, file = "E:/Group/Ollies_workspace/pspline_func_runs/reg_SW.rds")

reg_SE <- p_spline_function(dfRes,
                            region_name = "South East",
                            age_category = "All",
                            positive_variable = "Estbinres",
                            iterations = 20000,
                            warmup = 5000,
                            path_name = "E://Group/figure/temporal/MCMC_model_gam/National/Estbinres/South East/All/",
                            name = "over_55")
saveRDS(reg_SE, file = "E:/Group/Ollies_workspace/pspline_func_runs/reg_SE.rds")

reg_EE <- p_spline_function(dfRes,
                            region_name = "East of England",
                            age_category = "All",
                            positive_variable = "Estbinres",
                            iterations = 20000,
                            warmup = 5000,
                            path_name = "E://Group/figure/temporal/MCMC_model_gam/National/Estbinres/East of England/All/",
                            name = "over_55")
saveRDS(reg_EE, file = "E:/Group/Ollies_workspace/pspline_func_runs/reg_EE.rds")

reg_LN <- p_spline_function(dfRes,
                            region_name = "London",
                            age_category = "All",
                            positive_variable = "Estbinres",
                            iterations = 20000,
                            warmup = 5000,
                            path_name = "E://Group/figure/temporal/MCMC_model_gam/National/Estbinres/London/All/",
                            name = "over_55")
saveRDS(reg_LN, file = "E:/Group/Ollies_workspace/pspline_func_runs/reg_LN.rds")

reg_EM <- p_spline_function(dfRes,
                            region_name = "East Midlands",
                            age_category = "All",
                            positive_variable = "Estbinres",
                            iterations = 20000,
                            warmup = 5000,
                            path_name = "E://Group/figure/temporal/MCMC_model_gam/National/Estbinres/East Midlands/All/",
                            name = "over_55")
saveRDS(reg_EM, file = "E:/Group/Ollies_workspace/pspline_func_runs/reg_EM.rds")

reg_WM <- p_spline_function(dfRes,
                            region_name = "West Midlands",
                            age_category = "All",
                            positive_variable = "Estbinres",
                            iterations = 20000,
                            warmup = 5000,
                            path_name = "E://Group/figure/temporal/MCMC_model_gam/National/Estbinres/West Midlands/All/",
                            name = "over_55")
saveRDS(reg_WM, file = "E:/Group/Ollies_workspace/pspline_func_runs/reg_WM.rds")

reg_YH <- p_spline_function(dfRes,
                            region_name = "Yorkshire and The Humber",
                            age_category = "All",
                            positive_variable = "Estbinres",
                            iterations = 20000,
                            warmup = 5000,
                            path_name = "E://Group/figure/temporal/MCMC_model_gam/National/Estbinres/Yorkshire and The Humber/All/",
                            name = "over_55")
saveRDS(reg_YH, file = "E:/Group/Ollies_workspace/pspline_func_runs/reg_YH.rds")

####################################################################################################

reg_north <- p_spline_function(dfRes,
                               region_name = c("North West","North East","East Midlands","West Midlands","Yorkshire and The Humber"),
                               age_category = "All",
                               positive_variable = "Estbinres",
                               iterations = 20000,
                               warmup = 5000,
                               path_name = "E://Group/figure/temporal/MCMC_model_gam/National/North_wp/",
                               name = "north_wp")
saveRDS(reg_north, file = "E:/Group/Ollies_workspace/pspline_func_runs/north.rds")

reg_south <- p_spline_function(dfRes,
                               region_name = c("South West","South East","East of England","London"),
                               age_category = "All",
                               positive_variable = "Estbinres",
                               iterations = 20000,
                               warmup = 5000,
                               path_name = "E://Group/figure/temporal/MCMC_model_gam/National/South_wp/",
                               name = "south_wp")
saveRDS(reg_south, file = "E:/Group/Ollies_workspace/pspline_func_runs/south.rds")



#############################################################################

dfRest <- dfRes[dfRes$age <=17,]
hh2<- p_spline_function(dfRest,
                        region_name = "England",
                        age_category = "All",
                        positive_variable = "Estbinres",
                        iterations = 30000,
                        warmup = 5000,
                        path_name = "E://Group/figure/temporal/MCMC_model_gam/National/over_55/",
                        name = "over_55")
saveRDS(hh2, file = "E:/Group/Ollies_workspace/pspline_func_runs/age_new1.rds")

dfRest <- dfRes[dfRes$age>17 & dfRes$age<35,]
hh2<- p_spline_function(dfRest,
                        region_name = "England",
                        age_category = "All",
                        positive_variable = "Estbinres",
                        iterations = 30000,
                        warmup = 5000,
                        path_name = "E://Group/figure/temporal/MCMC_model_gam/National/over_55/",
                        name = "over_55")
saveRDS(hh2, file = "E:/Group/Ollies_workspace/pspline_func_runs/age_new2.rds")

dfRest <- dfRes[dfRes$age>34 & dfRes$age<55,]
hh2<- p_spline_function(dfRest,
                        region_name = "England",
                        age_category = "All",
                        positive_variable = "Estbinres",
                        iterations = 30000,
                        warmup = 5000,
                        path_name = "E://Group/figure/temporal/MCMC_model_gam/National/over_55/",
                        name = "over_55")
saveRDS(hh2, file = "E:/Group/Ollies_workspace/pspline_func_runs/age_new4.rds")

dfRest <- dfRes[dfRes$age >=55,]
hh2<- p_spline_function(dfRest,
                        region_name = "England",
                        age_category = "All",
                        positive_variable = "Estbinres",
                        iterations = 30000,
                        warmup = 5000,
                        path_name = "E://Group/figure/temporal/MCMC_model_gam/National/over_55/",
                        name = "over_55")
saveRDS(hh2, file = "E:/Group/Ollies_workspace/pspline_func_runs/age_new3.rds")


##########################################################################################



dfRest <- dfRes[dfRes$age <=11,]
hh2<- p_spline_function(dfRest,
                        region_name = "England",
                        age_category = "All",
                        positive_variable = "Estbinres",
                        iterations = 25000,
                        warmup = 5000,
                        path_name = "E://Group/figure/temporal/MCMC_model_gam/National/over_55/",
                        name = "over_55")
saveRDS(hh2, file = "E:/Group/Ollies_workspace/pspline_func_runs/age_new1a.rds")

dfRest <- dfRes[dfRes$age>11 & dfRes$age<18,]
hh2<- p_spline_function(dfRest,
                        region_name = "England",
                        age_category = "All",
                        positive_variable = "Estbinres",
                        iterations = 25000,
                        warmup = 5000,
                        path_name = "E://Group/figure/temporal/MCMC_model_gam/National/over_55/",
                        name = "over_55")
saveRDS(hh2, file = "E:/Group/Ollies_workspace/pspline_func_runs/age_new1b.rds")

dfRest <- dfRes[dfRes$age>17 & dfRes$age<25,]
hh2<- p_spline_function(dfRest,
                        region_name = "England",
                        age_category = "All",
                        positive_variable = "Estbinres",
                        iterations = 25000,
                        warmup = 5000,
                        path_name = "E://Group/figure/temporal/MCMC_model_gam/National/over_55/",
                        name = "over_55")
saveRDS(hh2, file = "E:/Group/Ollies_workspace/pspline_func_runs/age_new1c.rds")

dfRest <- dfRes[dfRes$age>24 & dfRes$age<35,]
hh2<- p_spline_function(dfRest,
                        region_name = "England",
                        age_category = "All",
                        positive_variable = "Estbinres",
                        iterations = 25000,
                        warmup = 5000,
                        path_name = "E://Group/figure/temporal/MCMC_model_gam/National/over_55/",
                        name = "over_55")
saveRDS(hh2, file = "E:/Group/Ollies_workspace/pspline_func_runs/age_new1d.rds")

dfRest <- dfRes[dfRes$age>34 & dfRes$age<45,]
hh2<- p_spline_function(dfRest,
                        region_name = "England",
                        age_category = "All",
                        positive_variable = "Estbinres",
                        iterations = 25000,
                        warmup = 5000,
                        path_name = "E://Group/figure/temporal/MCMC_model_gam/National/over_55/",
                        name = "over_55")
saveRDS(hh2, file = "E:/Group/Ollies_workspace/pspline_func_runs/age_new1e.rds")

dfRest <- dfRes[dfRes$age>44 & dfRes$age<55,]
hh2<- p_spline_function(dfRest,
                        region_name = "England",
                        age_category = "All",
                        positive_variable = "Estbinres",
                        iterations = 25000,
                        warmup = 5000,
                        path_name = "E://Group/figure/temporal/MCMC_model_gam/National/over_55/",
                        name = "over_55")
saveRDS(hh2, file = "E:/Group/Ollies_workspace/pspline_func_runs/age_new1f.rds")

dfRest <- dfRes[dfRes$age>54 & dfRes$age<65,]
hh2<- p_spline_function(dfRest,
                        region_name = "England",
                        age_category = "All",
                        positive_variable = "Estbinres",
                        iterations = 25000,
                        warmup = 5000,
                        path_name = "E://Group/figure/temporal/MCMC_model_gam/National/over_55/",
                        name = "over_55")
saveRDS(hh2, file = "E:/Group/Ollies_workspace/pspline_func_runs/age_new1g.rds")

dfRest <- dfRes[dfRes$age>64 & dfRes$age<75,]
hh2<- p_spline_function(dfRest,
                        region_name = "England",
                        age_category = "All",
                        positive_variable = "Estbinres",
                        iterations = 25000,
                        warmup = 5000,
                        path_name = "E://Group/figure/temporal/MCMC_model_gam/National/over_55/",
                        name = "over_55")
saveRDS(hh2, file = "E:/Group/Ollies_workspace/pspline_func_runs/age_new1h.rds")

dfRest <- dfRes[dfRes$age>74,]
hh2<- p_spline_function(dfRest,
                        region_name = "England",
                        age_category = "All",
                        positive_variable = "Estbinres",
                        iterations = 25000,
                        warmup = 5000,
                        path_name = "E://Group/figure/temporal/MCMC_model_gam/National/over_55/",
                        name = "over_55")
saveRDS(hh2, file = "E:/Group/Ollies_workspace/pspline_func_runs/age_new1i.rds")
