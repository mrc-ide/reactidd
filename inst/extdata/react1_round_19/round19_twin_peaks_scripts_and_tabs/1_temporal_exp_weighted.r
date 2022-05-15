#!/usr/bin/env Rscript
Sys.time()
rm(list=ls(all=TRUE))
setwd("E:/Group/Ollies_workspace/")

path_name <- "E://Group/figure/temporal/MCMC_model_exp/REACT/"
dir.create(path_name)
#' Then load required libraries and locally defined functions
library(prevalence)
library(mgcv)
library(stats)
library(cowplot)
library(rstan)
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


dfRes1 <- dfRes1[colnames(dfRes1) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","wt_antigen","nchild")]
dfRes2 <- dfRes2[colnames(dfRes2) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","wt_antigen","nchild")]
dfRes3 <- dfRes3[colnames(dfRes3) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","wt_antigen","nchild")]
dfRes4 <- dfRes4[colnames(dfRes4) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","wt_antigen","nchild")]
dfRes5 <- dfRes5[colnames(dfRes5) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","wt_antigen","nchild")]
dfRes6 <- dfRes6[colnames(dfRes6) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","wt_antigen","nchild")]
dfRes7 <- dfRes7[colnames(dfRes7) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","wt_antigen","nchild")]
dfRes8 <- dfRes8[colnames(dfRes8) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","wt_antigen","nchild")]
dfRes9 <- dfRes9[colnames(dfRes9) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","wt_antigen","nchild")]
dfRes10 <- dfRes10[colnames(dfRes10) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","wt_antigen","nchild")]
dfRes11 <- dfRes11[colnames(dfRes11) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","wt_antigen","nchild")]
dfRes12 <- dfRes12[colnames(dfRes12) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","wt_antigen","nchild")]
dfRes13 <- dfRes13[colnames(dfRes13) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","wt_antigen","age", "hh_size_cat","nchild","sympt_cat")]
dfRes14 <- dfRes14[colnames(dfRes14) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","wt_antigen","age", "hh_size_cat","nchild","sympt_cat")]
dfRes15 <- dfRes15[colnames(dfRes15) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","wt_antigen","age", "hh_size_cat","nchild","sympt_cat")]
dfRes16 <- dfRes16[colnames(dfRes16) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","wt_antigen","age", "hh_size_cat","nchild","sympt_cat")]
dfRes17 <- dfRes17[colnames(dfRes17) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","wt_antigen","age", "hh_size_cat","nchild","sympt_cat")]
dfRes18 <- dfRes18[colnames(dfRes18) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","wt_antigen","age", "hh_size_cat","nchild","sympt_cat")]
dfRes19 <- dfRes19[colnames(dfRes19) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","wt_antigen","age", "hh_size_cat","nchild","sympt_cat")]

#dfRes15$postcourier <- 0

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

dfRes1$d_swab <- dfRes1$d_comb
dfRes2$d_swab <- dfRes2$d_comb
dfRes3$d_swab <- dfRes3$d_comb
dfRes4$d_swab <- dfRes4$d_comb
dfRes5$d_swab <- dfRes5$d_comb
dfRes6$d_swab <- dfRes6$d_comb
dfRes7$d_swab <- dfRes7$d_comb
dfRes8$d_swab <- dfRes8$d_comb
dfRes9$d_swab <- dfRes9$d_comb
dfRes11$d_swab <- dfRes11$d_comb
dfRes12$d_swab <- dfRes12$d_comb
dfRes13$d_swab <- dfRes13$d_comb
dfRes14$d_swab <- dfRes14$d_comb
dfRes15$d_swab <- dfRes15$d_comb
dfRes16$d_swab <- dfRes16$d_comb
dfRes17$d_swab <- dfRes17$d_comb
dfRes18$d_swab <- dfRes18$d_comb
dfRes19$d_swab <- dfRes19$d_comb




rstan_options(auto_write = TRUE)
options(mc.cores = 4)
spline_model <- stan_model("stan_models/linear_weighted.stan")


exp_r_calculator <- function(dfRes9, dfRes10, round_name, result_variable, region_variable, rounds_variable){
  if(result_variable == "estbinres"){
    dfRes9$estbinres <- dfRes9$estbinres
    dfRes10$estbinres <- dfRes10$estbinres
  } else if(result_variable == "estbinres35"){
    dfRes9$estbinres <- dfRes9$estbinres35
    dfRes10$estbinres <- dfRes10$estbinres35
  } else if(result_variable == "asympt"){
    dfRes9$estbinres <- dfRes9$res_asympt
    dfRes10$estbinres <- dfRes10$res_asympt
  } else if(result_variable == "double_pos"){
    dfRes9$estbinres <- dfRes9$res_ct12_gt0
    dfRes10$estbinres <- dfRes10$res_ct12_gt0
  }
  
  if(region_variable == "London"){
    dfRes9 <- dfRes9[dfRes9$region=="London",]
    dfRes10 <- dfRes10[dfRes10$region=="London",]
  } else if(region_variable == "North West"){
    dfRes9 <- dfRes9[dfRes9$region=="North West",]
    dfRes10 <- dfRes10[dfRes10$region=="North West",]
  } else if(region_variable == "East of England"){
    dfRes9 <- dfRes9[dfRes9$region=="East of England",]
    dfRes10 <- dfRes10[dfRes10$region=="East of England",]
  } else if(region_variable == "West Midlands"){
    dfRes9 <- dfRes9[dfRes9$region=="West Midlands",]
    dfRes10 <- dfRes10[dfRes10$region=="West Midlands",]
  } else if(region_variable == "East Midlands"){
    dfRes9 <- dfRes9[dfRes9$region=="East Midlands",]
    dfRes10 <- dfRes10[dfRes10$region=="East Midlands",]
  } else if(region_variable == "North East"){
    dfRes9 <- dfRes9[dfRes9$region=="North East",]
    dfRes10 <- dfRes10[dfRes10$region=="North East",]
  } else if(region_variable == "South East"){
    dfRes9 <- dfRes9[dfRes9$region=="South East",]
    dfRes10 <- dfRes10[dfRes10$region=="South East",]
  } else if(region_variable == "Yorkshire"){
    dfRes9 <- dfRes9[dfRes9$region=="Yorkshire and The Humber",]
    dfRes10 <- dfRes10[dfRes10$region=="Yorkshire and The Humber",]
  } else if(region_variable == "South West"){
    dfRes9 <- dfRes9[dfRes9$region=="South West",]
    dfRes10 <- dfRes10[dfRes10$region=="South West",]
  }
  
  if(rounds_variable == "double"){
    dfRes <-rbind(dfRes9, dfRes10)
    dfRes<-dfRes[is.na(dfRes$estbinres)==FALSE &is.na(dfRes$d_comb)==FALSE,]
    ensemble <- data.frame(day=unique(dfRes$d_comb))
    
    for(i in seq_len(nrow(ensemble))){
      ensemble$obs[i] <- sum(dfRes[dfRes$d_comb ==ensemble$day[i],]$wt_antigen, na.rm=TRUE)
      ensemble$pos[i] <- sum(dfRes[dfRes$d_comb ==ensemble$day[i] &dfRes$estbinres==1,]$wt_antigen, na.rm=TRUE)
    }
    ensemble <- ensemble[is.na(ensemble$day)==FALSE,]
    ensemble <-ensemble[order(ensemble$day),]

  } else if(rounds_variable == "single"){
    dfRes <-rbind(dfRes10)
    dfRes<-dfRes[is.na(dfRes$estbinres)==FALSE &is.na(dfRes$d_comb)==FALSE,]
    ensemble <- data.frame(day=unique(dfRes$d_comb))
    
    for(i in seq_len(nrow(ensemble))){
      ensemble$obs[i] <- sum(dfRes[dfRes$d_comb ==ensemble$day[i],]$wt_antigen, na.rm=TRUE)
      ensemble$pos[i] <- sum(dfRes[dfRes$d_comb ==ensemble$day[i] &dfRes$estbinres==1,]$wt_antigen, na.rm=TRUE)
    }
    ensemble <- ensemble[is.na(ensemble$day)==FALSE,]
    ensemble <-ensemble[order(ensemble$day),]
    
  }
  #' Convert date to numeric
  ensemble$day <- as.numeric(ensemble$day)
  min_date_numeric <- min(ensemble$day)
  max_date_numeric <- max(ensemble$day)
  
  #' Set initial conditions and variable names to pass to stan
  X <- ensemble$day
  num_data <- length(X)
  Y <- as.numeric(ensemble$pos)
  N <- as.numeric(ensemble$obs)
  X_adj <- min(X)
  X <- X - X_adj
  
  #' Load and run stan model
  fit_spline_test <- sampling(spline_model,
                              iter=20000,
                              warmup =5000,
                              chains=4,
                              control = list(adapt_delta=0.95,
                                             max_treedepth = 10),
                              data = list(num_data = num_data,
                                          Y = Y,
                                          N =N,
                                          X = X))
  
  
  saveRDS(fit_spline_test, paste('stan_fits/', round_name,'_react_weighted.rds', sep=""))
  fit_spline <- fit_spline_test
  ff<-rstan::extract(fit_spline)
  
  #' Plot bivariate distribution
  chain <- data.frame(littler = ff$beta,
                      I0 = exp(ff$alpha)*100)
  
  #' 95% CI for r and median
  r<-quantile(chain$littler, c(0.5,0.025,0.975))
  
  #' 95% CI for R and median
  n_mean <- 2.2 #n_upper <- 3.34 n_lower <- 1.77
  b_mean <- 0.48 #b_upper <- 0.57 b_lower <- 0.26
  
  prob <- nrow(chain[chain$littler>0,])*100/nrow(chain)
  R <- (1+r/b_mean)**n_mean
  t <- log(2)/r
  
  table_row <- data.frame(round = round_name,
                          r= r[[1]],
                          r_lb = r[[2]],
                          r_ub = r[[3]],
                          R = R[[1]],
                          R_lb = R[[2]],
                          R_ub = R[[3]],
                          prob = prob,
                          t = t[[1]],
                          t_lb = t[[2]],
                          t_ub = t[[3]])
  table_row
  
}

table1_row1 <- exp_r_calculator(dfRes18, dfRes19, round_name = "Round1819",
                                result_variable = "estbinres",
                                region_variable = "England",
                                rounds_variable = "double")

table1_row2 <- exp_r_calculator(dfRes18, dfRes19, round_name = "Round19",
                                result_variable = "estbinres",
                                region_variable = "England",
                                rounds_variable = "single")

table1_row3 <- exp_r_calculator(dfRes18, dfRes19, round_name = "Round1819_asympt",
                                result_variable = "asympt",
                                region_variable = "England",
                                rounds_variable = "double")
table1_row4 <- exp_r_calculator(dfRes18, dfRes19, round_name = "Round19_asympt",
                                result_variable = "asympt",
                                region_variable = "England",
                                rounds_variable = "single")

table1_row5 <- exp_r_calculator(dfRes18, dfRes19, round_name = "Round1819_double_pos",
                                result_variable = "double_pos",
                                region_variable = "England",
                                rounds_variable = "double")
table1_row6 <- exp_r_calculator(dfRes18, dfRes19, round_name = "Round19_double_pos",
                                result_variable = "double_pos",
                                region_variable = "England",
                                rounds_variable = "single")
table1_row7 <- exp_r_calculator(dfRes18, dfRes19, round_name = "Round1819_estbinres35",
                                result_variable = "estbinres35",
                                region_variable = "England",
                                rounds_variable = "double")
table1_row8 <- exp_r_calculator(dfRes18, dfRes19, round_name = "Round19_estbinres35",
                                result_variable = "estbinres35",
                                region_variable = "England",
                                rounds_variable = "single")

dfRes18_sympt <- dfRes18[!(dfRes18$estbinres%in%1 & !(dfRes18$sympt_cat %in%c("Other symptoms","Classic COVID symptoms"))),]
dfRes19_sympt <- dfRes19[!(dfRes19$estbinres%in%1 & !(dfRes19$sympt_cat %in%c("Other symptoms","Classic COVID symptoms"))),]
table1_row9 <- exp_r_calculator(dfRes18_sympt, dfRes19_sympt, round_name = "Round19_sympt",
                                result_variable = "estbinres",
                                region_variable = "England",
                                rounds_variable = "single")

table1_row10 <- exp_r_calculator(dfRes18_sympt, dfRes19_sympt, round_name = "Round1819_sympt",
                                result_variable = "estbinres",
                                region_variable = "England",
                                rounds_variable = "double")


table1_row11 <- exp_r_calculator(dfRes18[dfRes18$age<=18,], dfRes19[dfRes19$age<=18,], round_name = "Round19_agel",
                                 result_variable = "estbinres",
                                 region_variable = "England",
                                 rounds_variable = "single")
table1_row12a <- exp_r_calculator(dfRes18[dfRes18$age>18 & dfRes18$age<35,], dfRes19[dfRes19$age>18 & dfRes19$age<35,], round_name = "Round19_agema",
                                 result_variable = "estbinres",
                                 region_variable = "England",
                                 rounds_variable = "single")

table1_row12b <- exp_r_calculator(dfRes18[dfRes18$age>34 & dfRes18$age<55,], dfRes19[dfRes19$age>34 & dfRes19$age<55,], round_name = "Round19_agemb",
                                 result_variable = "estbinres",
                                 region_variable = "England",
                                 rounds_variable = "single")
table1_row13 <- exp_r_calculator(dfRes18[dfRes18$age>=55,], dfRes19[dfRes19$age>=55,], round_name = "Round19_ageh",
                                 result_variable = "estbinres",
                                 region_variable = "England",
                                 rounds_variable = "single")

table1_row14 <- exp_r_calculator(dfRes18[dfRes18$age<=18,], dfRes19[dfRes19$age<=18,], round_name = "Round1819_agel",
                                 result_variable = "estbinres",
                                 region_variable = "England",
                                 rounds_variable = "double")
table1_row15a <- exp_r_calculator(dfRes18[dfRes18$age>18 & dfRes18$age<35,], dfRes19[dfRes19$age>18 & dfRes19$age<35,], round_name = "Round1819_agema",
                                 result_variable = "estbinres",
                                 region_variable = "England",
                                 rounds_variable = "double")
table1_row15b <- exp_r_calculator(dfRes18[dfRes18$age>34 & dfRes18$age<55,], dfRes19[dfRes19$age>34 & dfRes19$age<55,], round_name = "Round1819_agemb",
                                 result_variable = "estbinres",
                                 region_variable = "England",
                                 rounds_variable = "double")
table1_row16 <- exp_r_calculator(dfRes18[dfRes18$age>=55,], dfRes19[dfRes19$age>=55,], round_name = "Round1819_ageh",
                                 result_variable = "estbinres",
                                 region_variable = "England",
                                 rounds_variable = "double")


table1 <- rbind(table1_row1,
                table1_row2,
                table1_row3,
                table1_row4,
                table1_row5,
                table1_row6,
                table1_row7,
                table1_row8,
                table1_row9,
                table1_row10,
                table1_row11,
                table1_row12,
                table1_row13,
                table1_row14,
                table1_row15a,
                table1_row16)

write.table(table1, "E:/group/report/round19/rate_table.csv",sep=",,", row.names = FALSE)

table2_row1 <- exp_r_calculator(dfRes18, dfRes19, round_name = "Round1819_EastMidlands",
                                result_variable = "estbinres",
                                region_variable = "East Midlands",
                                rounds_variable = "double")

table2_row2 <- exp_r_calculator(dfRes18, dfRes19, round_name = "Round1819_WestMidlands",
                                result_variable = "estbinres",
                                region_variable = "West Midlands",
                                rounds_variable = "double")
table2_row3 <- exp_r_calculator(dfRes18, dfRes19, round_name = "Round1819_EastofEngland",
                                result_variable = "estbinres",
                                region_variable = "East of England",
                                rounds_variable = "double")
table2_row4 <- exp_r_calculator(dfRes18, dfRes19, round_name = "Round1819_London",
                                result_variable = "estbinres",
                                region_variable = "London",
                                rounds_variable = "double")
table2_row5 <- exp_r_calculator(dfRes18, dfRes19, round_name = "Round1819_NorthWest",
                                result_variable = "estbinres",
                                region_variable = "North West",
                                rounds_variable = "double")
table2_row6 <- exp_r_calculator(dfRes18, dfRes19, round_name = "Round1819_NorthEast",
                                result_variable = "estbinres",
                                region_variable = "North East",
                                rounds_variable = "double")
table2_row7 <- exp_r_calculator(dfRes18, dfRes19, round_name = "Round1819_SouthEast",
                                result_variable = "estbinres",
                                region_variable = "South East",
                                rounds_variable = "double")
table2_row8 <- exp_r_calculator(dfRes18, dfRes19, round_name = "Round1819_SouthWest",
                                result_variable = "estbinres",
                                region_variable = "South West",
                                rounds_variable = "double")
table2_row9 <- exp_r_calculator(dfRes18, dfRes19, round_name = "Round1819_Yorkshire",
                                result_variable = "estbinres",
                                region_variable = "Yorkshire",
                                rounds_variable = "double")

table2_row10 <- exp_r_calculator(dfRes18, dfRes19, round_name = "Round19_EastMidlands",
                                 result_variable = "estbinres",
                                 region_variable = "East Midlands",
                                 rounds_variable = "single")
table2_row11 <- exp_r_calculator(dfRes18, dfRes19, round_name = "Round19_WestMidlands",
                                 result_variable = "estbinres",
                                 region_variable = "West Midlands",
                                 rounds_variable = "single")
table2_row12 <- exp_r_calculator(dfRes18, dfRes19, round_name = "Round19_EastofEngland",
                                 result_variable = "estbinres",
                                 region_variable = "East of England",
                                 rounds_variable = "single")
table2_row13 <- exp_r_calculator(dfRes18, dfRes19, round_name = "Round19_London",
                                 result_variable = "estbinres",
                                 region_variable = "London",
                                 rounds_variable = "single")
table2_row14 <- exp_r_calculator(dfRes18, dfRes19, round_name = "Round19_NorthWest",
                                 result_variable = "estbinres",
                                 region_variable = "North West",
                                 rounds_variable = "single")
table2_row15 <- exp_r_calculator(dfRes18, dfRes19, round_name = "Round19_NorthEast",
                                 result_variable = "estbinres",
                                 region_variable = "North East",
                                 rounds_variable = "single")
table2_row16 <- exp_r_calculator(dfRes18, dfRes19, round_name = "Round19_SouthEast",
                                 result_variable = "estbinres",
                                 region_variable = "South East",
                                 rounds_variable = "single")
table2_row17 <- exp_r_calculator(dfRes18, dfRes19, round_name = "Round19_SouthWest",
                                 result_variable = "estbinres",
                                 region_variable = "South West",
                                 rounds_variable = "single")
table2_row18 <- exp_r_calculator(dfRes18, dfRes19, round_name = "Round19_Yorkshire",
                                 result_variable = "estbinres",
                                 region_variable = "Yorkshire",
                                 rounds_variable = "single")

table2<- rbind(table2_row1, table2_row2,
               table2_row3, table2_row4,
               table2_row5, table2_row6,
               table2_row7, table2_row8,
               table2_row9, table2_row10,
               table2_row11, table2_row12,
               table2_row13, table2_row14,
               table2_row15, table2_row16,
               table2_row17, table2_row18)

write.table(table2, "E:/group/report/round19/Region_rate_table.csv",sep=",,", row.names = FALSE)


table1 <- rbind(table1_row2,
                table1_row4,
                table1_row9,
                table1_row6,
                table1_row8,
                table1_row11,
                table1_row12a,
                table1_row12b,
                table1_row13,
                table2_row10,
                table2_row11,
                table2_row12,
                table2_row13,
                table2_row14,
                table2_row15,
                table2_row15,
                table2_row16,
                table2_row17,
                table2_row18)

table2 <- rbind(table1_row1,
                table1_row3,
                table1_row10,
                table1_row5,
                table1_row7,
                table1_row14,
                table1_row15a,
                table1_row15b,
                table1_row16,
                table2_row1,
                table2_row2,
                table2_row3,
                table2_row4,
                table2_row5,
                table2_row6,
                table2_row7,
                table2_row8,
                table2_row9)

write.table(rbind(table1,table2), "E:/group/report/round19/rate_table_all.csv",sep=",,", row.names = FALSE)




table3_row1 <- exp_r_calculator(dfRes18[dfRes18$age>4 & dfRes18$age<12,], dfRes19[dfRes19$age>4 & dfRes19$age<12,], round_name = "Round19_age_a",
                                result_variable = "estbinres",
                                region_variable = "England",
                                rounds_variable = "single")
table3_row2 <- exp_r_calculator(dfRes18[dfRes18$age>11 & dfRes18$age<12,], dfRes19[dfRes19$age>11 & dfRes19$age<18,], round_name = "Round19_age_b",
                                result_variable = "estbinres",
                                region_variable = "England",
                                rounds_variable = "single")
table3_row3 <- exp_r_calculator(dfRes18[dfRes18$age>4 & dfRes18$age<12,], dfRes19[dfRes19$age>17 & dfRes19$age<25,], round_name = "Round19_age_c",
                                result_variable = "estbinres",
                                region_variable = "England",
                                rounds_variable = "single")
table3_row4 <- exp_r_calculator(dfRes18[dfRes18$age>4 & dfRes18$age<12,], dfRes19[dfRes19$age>24 & dfRes19$age<35,], round_name = "Round19_age_d",
                                result_variable = "estbinres",
                                region_variable = "England",
                                rounds_variable = "single")
table3_row5 <- exp_r_calculator(dfRes18[dfRes18$age>4 & dfRes18$age<12,], dfRes19[dfRes19$age>34 & dfRes19$age<45,], round_name = "Round19_age_e",
                                result_variable = "estbinres",
                                region_variable = "England",
                                rounds_variable = "single")
table3_row6 <- exp_r_calculator(dfRes18[dfRes18$age>4 & dfRes18$age<12,], dfRes19[dfRes19$age>44 & dfRes19$age<55,], round_name = "Round19_age_f",
                                result_variable = "estbinres",
                                region_variable = "England",
                                rounds_variable = "single")
table3_row7 <- exp_r_calculator(dfRes18[dfRes18$age>4 & dfRes18$age<12,], dfRes19[dfRes19$age>54 & dfRes19$age<65,], round_name = "Round19_age_g",
                                result_variable = "estbinres",
                                region_variable = "England",
                                rounds_variable = "single")
table3_row8 <- exp_r_calculator(dfRes18[dfRes18$age>4 & dfRes18$age<12,], dfRes19[dfRes19$age>64 & dfRes19$age<75,], round_name = "Round19_age_h",
                                result_variable = "estbinres",
                                region_variable = "England",
                                rounds_variable = "single")
table3_row9 <- exp_r_calculator(dfRes18[dfRes18$age>4 & dfRes18$age<12,], dfRes19[dfRes19$age>74 & dfRes19$age<1200,], round_name = "Round19_age_i",
                                result_variable = "estbinres",
                                region_variable = "England",
                                rounds_variable = "single")

table3 <- rbind(table3_row1,
                table3_row2,
                table3_row3,
                table3_row4,
                table3_row5,
                table3_row6,
                table3_row7,
                table3_row8,
                table3_row9)

write.table(rbind(table3), "E:/group/report/round19/rate_table_ages.csv",sep=",,", row.names = FALSE)

