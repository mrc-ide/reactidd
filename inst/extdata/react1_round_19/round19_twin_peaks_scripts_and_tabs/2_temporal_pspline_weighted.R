#!/usr/bin/env Rscript
Sys.time()
rm(list=ls(all=TRUE))
setwd("E:/Group/Ollies_workspace/")

path_name <- "E://Group/figure/temporal/MCMC_model_gam/National/England/R16Weighted/"
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

# Artificially give weights to round 14 15 16 17
#dfRes17$wt_antigen <- NA
#dfRes17<-dfRes17[is.na(dfRes17$estbinres)==FALSE &is.na(dfRes17$d_comb)==FALSE,]
#dfRes17[is.na(dfRes17$estbinres)==FALSE,]$wt_antigen <- 1.0

#dfRes18$wt_antigen <- NA
#dfRes18<-dfRes18[is.na(dfRes18$estbinres)==FALSE &is.na(dfRes18$d_comb)==FALSE,]
#dfRes18[is.na(dfRes18$estbinres)==FALSE,]$wt_antigen <- 1.0

#dfRes19$wt_antigen <- NA
#dfRes19<-dfRes19[is.na(dfRes19$estbinres)==FALSE &is.na(dfRes19$d_comb)==FALSE,]
#dfRes19[is.na(dfRes19$estbinres)==FALSE,]$wt_antigen <- 1.0



dfRes1 <- dfRes1[colnames(dfRes1) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","age_group_char","age","wt_antigen")]
dfRes2 <- dfRes2[colnames(dfRes2) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","age_group_char","age","wt_antigen")]
dfRes3 <- dfRes3[colnames(dfRes3) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","age_group_char","age","wt_antigen")]
dfRes4 <- dfRes4[colnames(dfRes4) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","age_group_char","age","wt_antigen")]
dfRes5 <- dfRes5[colnames(dfRes5) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","age_group_char","age","wt_antigen")]
dfRes6 <- dfRes6[colnames(dfRes6) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","age_group_char","age","wt_antigen")]
dfRes7 <- dfRes7[colnames(dfRes7) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","age_group_char","age","wt_antigen")]
dfRes8 <- dfRes8[colnames(dfRes8) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","age_group_char","age","wt_antigen")]
dfRes9 <- dfRes9[colnames(dfRes9) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","age_group_char","age","wt_antigen")]
dfRes10 <- dfRes10[colnames(dfRes10) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","age_group_char","age","wt_antigen")]
dfRes11 <- dfRes11[colnames(dfRes11) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","age_group_char","age","wt_antigen")]
dfRes12 <- dfRes12[colnames(dfRes12) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","age_group_char","age","wt_antigen")]
dfRes13 <- dfRes13[colnames(dfRes13) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","age_group_char","age","wt_antigen")]
dfRes14 <- dfRes14[colnames(dfRes14) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","age_group_char","age","wt_antigen")]
dfRes15 <- dfRes15[colnames(dfRes15) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","age_group_char","age","wt_antigen")]
dfRes16 <- dfRes16[colnames(dfRes16) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","age_group_char","age","wt_antigen")]
dfRes17 <- dfRes17[colnames(dfRes17) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","age_group_char","age","wt_antigen")]
dfRes18 <- dfRes18[colnames(dfRes18) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","age_group_char","age","wt_antigen")]
dfRes19 <- dfRes19[colnames(dfRes19) %in% c("estbinres","estbinres35","res_ct12_gt0","res_asympt","d_swab","d_comb","region","lacode","age_group_char","age","wt_antigen")]


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


dfRes <- rbind(dfRes1, dfRes2, dfRes3, dfRes4, dfRes5, dfRes6, dfRes7,
               dfRes8, dfRes9, dfRes10, dfRes11, dfRes12, dfRes13, dfRes14,
               dfRes15, dfRes16, dfRes17, dfRes18, dfRes19)

# Save out data for github
#pos<-table(dfRes[dfRes$estbinres==1,]$region,dfRes[dfRes$estbinres==1,]$d_comb)
#neg<-table(dfRes[dfRes$estbinres%in% c(0,1),]$region,dfRes[dfRes$estbinres %in% c(0,1),]$d_comb)
#write.csv(t(neg), 'total.csv')
#write.csv(t(pos), 'positive.csv')

############################### 
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
days_per_knot <- 5 #target value
num_knots <- ceiling((max_date_numeric- min_date_numeric)/days_per_knot)+7
days_per_knot <- (max_date_numeric - min_date_numeric)/(num_knots -7)

spline_degree <- 3
num_basis <- num_knots + spline_degree - 1

X <- ensemble$day
num_data <- length(X)
knots <- unname(seq(min(X)-3*days_per_knot, max(X)+3*days_per_knot, length.out = num_knots))

Y <- as.numeric(ensemble$pos)
N <- as.numeric(ensemble$obs)

#' Load and run stan model
rstan_options(auto_write = TRUE)
options(mc.cores = 8)
spline_model <- stan_model("stan_models/b_splines_actual_weighted.stan")
fit_spline_test <- sampling(spline_model,
                       iter=50000,
                       warmup =5000,
                       chains=8,
                       control = list(adapt_delta=0.95,
                                      max_treedepth = 10),
                       data = list(num_data = num_data,
                                   num_knots = num_knots,
                                   knots = knots,
                                   Y = Y,
                                   N =N,
                                   X = X,
                                   spline_degree = spline_degree))


saveRDS(fit_spline_test, paste('E:/Group/Ollies_workspace/stan_fits/weighted19_2','.rds',sep=''))
#fit_spline <- readRDS(paste('E:/Group/Ollies_workspace/stan_fits/weighted19','.rds',sep=''))
fit_spline <- fit_spline_test
ff<-rstan::extract(fit_spline)
#' Quick check of mcmc
ff$Y_hat[1,]
traceplot(fit_spline,pars=c('a[1]', 'a[4]', 'tau', 'a[5]'), inc_warmup=TRUE, nrow=2)

#' ## Analysis and plotting 
#' Plot of the model fit
X_new <- seq(min(X)-3*days_per_knot, max(X)+3*days_per_knot, 1)
df_plot_model<-plot_gam_fit(ff, X_new, num_basis, spline_degree)
df_plot_model$d_comb <- as.Date(df_plot_model$x-18383, origin=as.Date("2020-05-01"))

tabDaySwab <- table(dfRes$estbinres,dfRes$d_comb)
tabDaySwab[1,] <- N-Y
tabDaySwab[2,] <- Y
dayPCIs <- add_conf_ints(tabDaySwab,poscol="1",negcol="0")
df_plot <- as.data.frame(dayPCIs)
df_plot <- df_plot[is.na(df_plot$p)==FALSE,]
df_plot$d_comb <- as.Date(rownames(df_plot))

max_date<-max(df_plot$d_comb)
min_date<-min(df_plot$d_comb)

df_plot_model<-df_plot_model[df_plot_model$d_comb>=min_date & df_plot_model$d_comb<=max_date,]

mindate1 <- as.Date("2020-05-01")
mindate8<-as.Date("2021-1-06")
mindate14<-as.Date("2021-09-09")
maxdate7<-as.Date("2020-12-03")
maxdate13 <- as.Date("2021-07-12")



df_plot_model1<-df_plot_model[df_plot_model$d_comb <=maxdate7,]
df_plot_model2<-df_plot_model[df_plot_model$d_comb >=mindate8 &df_plot_model$d_comb<=maxdate13,]
df_plot_model3<-df_plot_model[df_plot_model$d_comb >=mindate14,]


plot3 <- ggplot(data = df_plot, aes(x= d_comb, y =p*100))+
  geom_point()+
  geom_errorbar(aes(ymin=lb*100, ymax=ub*100))+
  coord_cartesian(ylim=c(0.01,10), xlim=c(mindate1, max(df_plot$d_comb)))+
  theme_bw(base_size = 18)+ 
  xlab("Day of swab")+
  ylab("Prevalence (%)")+
  scale_x_date(date_breaks = "1 month", date_labels = "%b\n%Y")+
  geom_line(data= df_plot_model1, aes(y=inv_logit(p)*100))+
  geom_ribbon(data = df_plot_model1,
              aes(y=inv_logit(p)*100,
                  ymin=inv_logit(lb_2.5)*100,
                  ymax=inv_logit(ub_97.5)*100),
              alpha=0.2)+
  geom_ribbon(data = df_plot_model1,
              aes(y=inv_logit(p)*100,
                  ymin=inv_logit(lb_25)*100,
                  ymax=inv_logit(ub_75)*100),
              alpha=0.2)+
  geom_line(data= df_plot_model2, aes(y=inv_logit(p)*100))+
  geom_ribbon(data = df_plot_model2,
              aes(y=inv_logit(p)*100,
                  ymin=inv_logit(lb_2.5)*100,
                  ymax=inv_logit(ub_97.5)*100),
              alpha=0.2)+
  geom_ribbon(data = df_plot_model2,
              aes(y=inv_logit(p)*100,
                  ymin=inv_logit(lb_25)*100,
                  ymax=inv_logit(ub_75)*100),
              alpha=0.2)+
  geom_line(data= df_plot_model3, aes(y=inv_logit(p)*100))+
  geom_ribbon(data = df_plot_model3,
              aes(y=inv_logit(p)*100,
                  ymin=inv_logit(lb_2.5)*100,
                  ymax=inv_logit(ub_97.5)*100),
              alpha=0.2)+
  geom_ribbon(data = df_plot_model3,
              aes(y=inv_logit(p)*100,
                  ymin=inv_logit(lb_25)*100,
                  ymax=inv_logit(ub_75)*100),
              alpha=0.2)+
  scale_y_log10()



pdf(paste( "E:/group/report/round19/Temporal_prevlog.pdf", sep=""), width = 24, height =6)
plot3
dev.off()


plot3 <- ggplot(data = df_plot, aes(x= d_comb, y =p*100))+
  geom_point()+
  geom_errorbar(aes(ymin=lb*100, ymax=ub*100))+
  coord_cartesian(ylim=c(0.0,10), xlim=c(mindate1, max(df_plot$d_comb)))+
  theme_bw(base_size = 18)+ 
  xlab("Day of swab")+
  ylab("Prevalence (%)")+
  scale_x_date(date_breaks = "1 month", date_labels = "%b\n%Y")+
  geom_line(data= df_plot_model1, aes(y=inv_logit(p)*100))+
  geom_ribbon(data = df_plot_model1,
              aes(y=inv_logit(p)*100,
                  ymin=inv_logit(lb_2.5)*100,
                  ymax=inv_logit(ub_97.5)*100),
              alpha=0.2)+
  geom_ribbon(data = df_plot_model1,
              aes(y=inv_logit(p)*100,
                  ymin=inv_logit(lb_25)*100,
                  ymax=inv_logit(ub_75)*100),
              alpha=0.2)+
  geom_line(data= df_plot_model2, aes(y=inv_logit(p)*100))+
  geom_ribbon(data = df_plot_model2,
              aes(y=inv_logit(p)*100,
                  ymin=inv_logit(lb_2.5)*100,
                  ymax=inv_logit(ub_97.5)*100),
              alpha=0.2)+
  geom_ribbon(data = df_plot_model2,
              aes(y=inv_logit(p)*100,
                  ymin=inv_logit(lb_25)*100,
                  ymax=inv_logit(ub_75)*100),
              alpha=0.2)+
  geom_line(data= df_plot_model3, aes(y=inv_logit(p)*100))+
  geom_ribbon(data = df_plot_model3,
              aes(y=inv_logit(p)*100,
                  ymin=inv_logit(lb_2.5)*100,
                  ymax=inv_logit(ub_97.5)*100),
              alpha=0.2)+
  geom_ribbon(data = df_plot_model3,
              aes(y=inv_logit(p)*100,
                  ymin=inv_logit(lb_25)*100,
                  ymax=inv_logit(ub_75)*100),
              alpha=0.2)



pdf(paste( "E:/group/report/round19/Temporal_prev.pdf", sep=""), width = 24, height =6)
plot3
dev.off()

pdf("Prevalence_full_log.pdf", width = 16, height =4)
plot3+
  geom_vline(xintercept=as.Date("2021-05-19"), color='red',linetype="dashed")+
  geom_vline(xintercept=as.Date("2021-09-08"), color='blue',linetype="dashed")+
  geom_vline(xintercept=as.Date("2021-10-18"), color='dark green',linetype="dashed")
dev.off()


# Calculate rolling two week R number
#df_plot_R <- plot_R(ff,X_new,num_basis,spline_degree,link_function = "logit",n=2.29, b=0.36, tau_max=14)
#df_plot_R$d_comb <- as.Date(df_plot_R$x-18383, origin=as.Date("2020-05-01"))

#tau_max=14
#df_plot_R<-df_plot_R[df_plot_R$d_comb>=min_date+tau_max & df_plot_R$d_comb<=max_date,]
#df_plot_R[is.na(df_plot_R$r)==TRUE,]$prob <- NA

#saveRDS(df_plot_R, "df_plot_R.rds")

#plot4 <- ggplot(data = df_plot_R, aes(x=d_comb,y=r))+
#  geom_line()+
#  geom_ribbon(aes(ymin=lb_2.5, ymax=ub_97.5), alpha=0.2)+
#  geom_ribbon(aes(ymin=lb_25, ymax=ub_75), alpha=0.2)+
#  theme_bw()+
#  coord_cartesian(ylim=c(0,2), xlim=c(min(df_plot$d_comb), max(df_plot$d_comb)))+
#  geom_line(y=1, linetype='dashed')+
#  geom_line(aes(y=prob), color='red')+
#  scale_x_date(date_breaks = "1 month", date_labels = "%b")+
#  xlab("Day of swab test")+
#  ylab("Reproductive number R_t")

#png(paste(path_name, "ReproductionNumber.png", sep=""), width = 1000, height =500)
#plot4
#dev.off()
#pdf(paste(path_name, "ReproductionNumber.pdf", sep=""), width = 12, height =6)
#plot4
#dev.off()


#Calculate growth rate 
X_REACT <- seq(min(X), max(X), by=1)
df_plot_r <- plot_growth_rate2(ff,X_REACT,num_basis, spline_degree, link_function = "logit")
df_plot_r$d_comb <- as.Date(df_plot_r$x-18383, origin=as.Date("2020-05-01"))

saveRDS(df_plot_r, "abcdefg.rds")

maxdate19 <- as.Date("2022-04-01")
df_plot_r1<-df_plot_r[df_plot_r$d_comb <=maxdate7,]
df_plot_r2<-df_plot_r[df_plot_r$d_comb >=mindate8 &df_plot_r$d_comb<=maxdate13,]
df_plot_r3<-df_plot_r[df_plot_r$d_comb >=mindate14 &df_plot_r$d_comb<=maxdate19-1,]


df_plot_r1<-df_plot_r1[df_plot_r1$d_comb>=min_date & df_plot_r1$d_comb<=max_date,]
df_plot_r1$r1 <- df_plot_r1$r
df_plot_r1$r2 <- df_plot_r1$r
df_plot_r1[df_plot_r1$r1>0.0,]$r1 <- NA
df_plot_r1[df_plot_r1$r2<0.0,]$r2 <- NA
df_plot_r1$middle <- 0.0

df_plot_r2<-df_plot_r2[df_plot_r2$d_comb>=min_date & df_plot_r2$d_comb<=max_date,]
df_plot_r2$r1 <- df_plot_r2$r
df_plot_r2$r2 <- df_plot_r2$r
df_plot_r2[df_plot_r2$r1>0.0,]$r1 <- NA
df_plot_r2[df_plot_r2$r2<0.0,]$r2 <- NA
df_plot_r2$middle <- 0.0

df_plot_r3<-df_plot_r3[df_plot_r3$d_comb>=min_date & df_plot_r3$d_comb<=max_date,]
df_plot_r3$r1 <- df_plot_r3$r
df_plot_r3$r2 <- df_plot_r3$r
df_plot_r3[df_plot_r3$r1>0.0,]$r1 <- NA
df_plot_r3[df_plot_r3$r2<0.0,]$r2 <- NA
df_plot_r3$middle <- 0.0


second_axis_values <- c(log(2)/5, log(2)/10,log(2)/15, 0, log(2)/-15, log(2)/-10, log(2)/-5)
second_axis_labels <- c("5", "10", "15",expression(infinity / -infinity),"-15","-10","-5")

plot2<-ggplot(data = df_plot_r1, aes(x=d_comb,y=r))+
  geom_line(data = df_plot_r1,aes(y=r1), col = '#009900')+
  geom_line(data = df_plot_r1,aes(y=r2), col = 'red')+
  geom_ribbon(data = df_plot_r1,aes(ymin=pmin(lb_2.5,middle), ymax=pmin(lb_25,middle)),alpha=0.2, fill='#009900')+
  geom_ribbon(data = df_plot_r1,aes(ymin=pmin(lb_25,middle), ymax=pmin(ub_75,middle)),alpha=0.4, fill='#009900')+
  geom_ribbon(data = df_plot_r1,aes(ymin=pmin(ub_75,middle), ymax=pmin(ub_97.5,middle)),alpha=0.2, fill='#009900')+
  geom_ribbon(data = df_plot_r1,aes(ymin=pmax(lb_2.5,middle), ymax=pmax(lb_25,middle)),alpha=0.2, fill='red')+
  geom_ribbon(data = df_plot_r1,aes(ymin=pmax(lb_25,middle), ymax=pmax(ub_75,middle)),alpha=0.4, fill="red")+
  geom_ribbon(data = df_plot_r1,aes(ymin=pmax(ub_75,middle), ymax=pmax(ub_97.5,middle)),alpha=0.2, fill='red')+
  geom_line(data = df_plot_r2,aes(y=r1), col = '#009900')+
  geom_line(data = df_plot_r2,aes(y=r2), col = 'red')+
  geom_ribbon(data = df_plot_r2,aes(ymin=pmin(lb_2.5,middle), ymax=pmin(lb_25,middle)),alpha=0.2, fill='#009900')+
  geom_ribbon(data = df_plot_r2,aes(ymin=pmin(lb_25,middle), ymax=pmin(ub_75,middle)),alpha=0.4, fill='#009900')+
  geom_ribbon(data = df_plot_r2,aes(ymin=pmin(ub_75,middle), ymax=pmin(ub_97.5,middle)),alpha=0.2, fill='#009900')+
  geom_ribbon(data = df_plot_r2,aes(ymin=pmax(lb_2.5,middle), ymax=pmax(lb_25,middle)),alpha=0.2, fill='red')+
  geom_ribbon(data = df_plot_r2,aes(ymin=pmax(lb_25,middle), ymax=pmax(ub_75,middle)),alpha=0.4, fill="red")+
  geom_ribbon(data = df_plot_r2,aes(ymin=pmax(ub_75,middle), ymax=pmax(ub_97.5,middle)),alpha=0.2, fill='red')+
  geom_line(data = df_plot_r3,aes(y=r1), col = '#009900')+
  geom_line(data = df_plot_r3,aes(y=r2), col = 'red')+
  geom_ribbon(data = df_plot_r3,aes(ymin=pmin(lb_2.5,middle), ymax=pmin(lb_25,middle)),alpha=0.2, fill='#009900')+
  geom_ribbon(data = df_plot_r3,aes(ymin=pmin(lb_25,middle), ymax=pmin(ub_75,middle)),alpha=0.4, fill='#009900')+
  geom_ribbon(data = df_plot_r3,aes(ymin=pmin(ub_75,middle), ymax=pmin(ub_97.5,middle)),alpha=0.2, fill='#009900')+
  geom_ribbon(data = df_plot_r3,aes(ymin=pmax(lb_2.5,middle), ymax=pmax(lb_25,middle)),alpha=0.2, fill='red')+
  geom_ribbon(data = df_plot_r3,aes(ymin=pmax(lb_25,middle), ymax=pmax(ub_75,middle)),alpha=0.4, fill="red")+
  geom_ribbon(data = df_plot_r3,aes(ymin=pmax(ub_75,middle), ymax=pmax(ub_97.5,middle)),alpha=0.2, fill='red')+
  theme_bw(base_size = 18)+
  coord_cartesian(ylim=c(-0.15,0.15),xlim=c(min(df_plot$d_comb), max(df_plot$d_comb)))+
  geom_line(y=0, linetype='dashed')+
  xlab("Date")+
  scale_x_date(date_breaks = "1 month", date_labels = "%b\n%Y")+
  ylab("Instantaenous growth rate")+
  scale_y_continuous(sec.axis = sec_axis(~.,labels = second_axis_labels, breaks= second_axis_values, name = "Doubling(+) / Halving(-) Time (days)"))



pdf(paste( "E:/group/report/round19/Temporal_gr.pdf", sep=""), width = 24, height =6)
plot2
dev.off()



# Get csv of spline to send to gov
X_REACT <- seq(min(X), max(X), by=1)
react_array <- get_response_posterior(ff,
                                      X_REACT,
                                      X_REACT)


plot_response <- function(Y_array, X){
  dfY <- data.frame(x = X)
  
  for(i in seq_len(length(X))){
    dfY$p[i] <-median(Y_array[,i])
    dfY$lb_2.5[i] <- quantile(Y_array[,i], probs=0.025)
    dfY$lb_25[i] <- quantile(Y_array[,i], probs=0.25)
    dfY$ub_97.5[i] <- quantile(Y_array[,i], probs=0.975)
    dfY$ub_75[i] <- quantile(Y_array[,i], probs=0.75)
    
  }
  
  dfY
  
}

react_plot <- plot_response(react_array, as.Date(X_REACT-18383, origin=as.Date("2020-05-01")))
react_plot1 <- react_plot[react_plot$x<=as.Date("2020-12-03"),]
react_plot2 <- react_plot[react_plot$x>=as.Date("2021-01-06") &react_plot$x<=as.Date("2021-07-12"),]
react_plot3 <- react_plot[react_plot$x>=as.Date("2021-09-09"),]

write.csv(rbind(react_plot1, react_plot2, react_plot3), "p_spline.csv")

## Combined figure #####
fit_linear<-readRDS(paste('stan_fits/', "Round19",'_react_weighted.rds', sep=""))
fflin<-rstan::extract(fit_linear)
X1 <- as.numeric(ensemble$day[ensemble$day>=as.Date("2022-03-08")])
N1 <- as.numeric(ensemble$obs[ensemble$day>=as.Date("2022-03-08")])
X_adj1 <- min(X1)
X1<-X1 - X_adj1
fitlin<-plot_lin_fit_3(fflin, X1, X1)
fitlin$x <- as.Date(fitlin$x -18383+X_adj1, origin=as.Date("2020-05-01"))



fit_linear<-readRDS(paste('stan_fits/', "Round1819",'_react_weighted.rds', sep=""))
fflin2<-rstan::extract(fit_linear)
X1 <- as.numeric(unique(dfRes$d_comb)[unique(dfRes$d_comb)>=as.Date("2022-02-08")])
N1 <- as.numeric(unique(dfRes$d_comb)[unique(dfRes$d_comb)>=as.Date("2022-02-08")])
X_adj1 <- min(X1)
X1<-X1 - X_adj1
fitlin2<-plot_lin_fit_3(fflin2, X1, X1)
fitlin2$x <- as.Date(fitlin2$x -18383+X_adj1, origin=as.Date("2020-05-01"))

fit_linear<-readRDS(paste('stan_fits/', "Round16_early",'_react_weighted.rds', sep=""))
fflin3<-rstan::extract(fit_linear)
X1 <- as.numeric(unique(dfRes$d_comb)[unique(dfRes$d_comb)>=as.Date("2021-11-23")&unique(dfRes$d_comb)<as.Date("2021-12-01")])
N1 <- as.numeric(unique(dfRes$d_comb)[unique(dfRes$d_comb)>=as.Date("2021-11-23")&unique(dfRes$d_comb)<as.Date("2021-12-01")])
X_adj1 <- min(X1)
X1<-X1 - X_adj1
fitlin3<-plot_lin_fit_3(fflin3, X1, X1)
fitlin3$x <- as.Date(fitlin3$x -18383+X_adj1, origin=as.Date("2020-05-01"))

table(fflin2$beta-fflin3$beta<0)
59793/60000

plote<-plot3+
  geom_point(data = df_plot, aes(x= d_comb, y =p*100))+
  geom_errorbar(data = df_plot, aes(x= d_comb, y =p*100, ymin=lb*100, ymax=ub*100))+
  coord_cartesian(ylim=c(0.9,20), xlim=c(as.Date("2021-11-23"), max(df_plot$d_comb)))+
  theme_bw(base_size = 18)+
  xlab("Day of swab")+
  ylab("Prevalence (%)")+
  scale_x_date(date_breaks = "1 month", date_labels = "%b")+
  geom_line(data= fitlin,aes(y=p*100, x=x, color='red'))+
  geom_ribbon(data= fitlin,aes(y=p*100, x=x,ymin=lb*100, ymax=ub*100, fill='red'), alpha=0.2)+
  geom_line(data= fitlin2,aes(y=p*100, x=x, color='blue'))+
  geom_ribbon(data= fitlin2,aes(y=p*100, x=x,ymin=lb*100, ymax=ub*100,fill='blue'), alpha=0.2)+
  #geom_line(data= fitlin3,aes(y=p*100, x=x, color='dark green'))+
  #geom_ribbon(data= fitlin3,aes(y=p*100, x=x,ymin=lb*100, ymax=ub*100,fill='dark green'), alpha=0.2)+
  scale_x_date(date_breaks = "1 month", date_labels = "%b")+
  scale_y_log10()+
  scale_color_brewer(palette = "Set1", name= "Rounds",
                     labels=c("18-19 partial","19 partial"))+
  scale_fill_brewer(palette = "Set1", name= "Rounds",
                   labels=c("18-19 partial","19 partial"))+
  theme(panel.grid.minor = element_blank(),
        legend.position = c(0.15,0.85))

pdf(paste( "E:/group/report/round19/Temporal3.pdf", sep=""), width = 8, height =8, useDingbats = FALSE)
plote
dev.off()
inv_logit(-3.74)

react_plot$d_comb <- as.Date(react_plot$x)

plot3 <- ggplot(data = df_plot, aes(x= d_comb, y =p*100))+
  geom_point()+
  geom_errorbar(aes(ymin=lb*100, ymax=ub*100))+
  coord_cartesian(ylim=c(0.01,3), xlim=c(mindate1, max(df_plot$d_comb)))+
  theme_bw(base_size = 18)+ 
  xlab("Day of swab")+
  ylab("Prevalence (%)")+
  scale_x_date(date_breaks = "1 month", date_labels = "%b")+
  geom_line(data= react_plot, aes(y=inv_logit(p)*100))+
  geom_ribbon(data = react_plot,
              aes(y=inv_logit(p)*100,
                  ymin=inv_logit(lb_2.5)*100,
                  ymax=inv_logit(ub_97.5)*100),
              alpha=0.2)+
  geom_ribbon(data = react_plot,
              aes(y=inv_logit(p)*100,
                  ymin=inv_logit(lb_25)*100,
                  ymax=inv_logit(ub_75)*100),
              alpha=0.2)+
  scale_y_log10()
plot3





## Combined figure #####
fit_linear<-readRDS(paste('stan_fits/', "Round16",'_react_weighted.rds', sep=""))
fflin<-rstan::extract(fit_linear)
X1 <- as.numeric(ensemble$day[ensemble$day>=as.Date("2021-11-23")])
N1 <- as.numeric(ensemble$obs[ensemble$day>=as.Date("2021-11-23")])
X_adj1 <- min(X1)
X1<-X1 - X_adj1
fitlin<-plot_lin_fit_3(fflin, X1, X1)
fitlin$x <- as.Date(fitlin$x -18383+X_adj1, origin=as.Date("2020-05-01"))



fit_linear<-readRDS(paste('stan_fits/', "Round16_latter",'_react_weighted.rds', sep=""))
fflin2<-rstan::extract(fit_linear)

start<- as.numeric(as.Date("2021-12-01"))
end<- as.numeric(as.Date("2021-12-01"))+30

X1 <- seq(start, end, 1)
N1 <- X1
X_adj1 <- min(X1)
X1<-X1 - X_adj1
fitlin2<-plot_lin_fit_3(fflin2, X1, X1)
fitlin2$x <- as.Date(fitlin2$x -18383+X_adj1, origin=as.Date("2020-05-01"))




ggplot(data = df_plot, aes(x= d_comb, y =p*100))+
  geom_point()+
  geom_errorbar(aes(ymin=lb*100, ymax=ub*100))+
  coord_cartesian(ylim=c(0.5,5), xlim=c(mindate14+80, max(df_plot$d_comb)+10))+
  theme_bw(base_size = 18)+ 
  xlab("Day of swab")+
  ylab("Prevalence (%)")+
  scale_x_date(date_breaks = "1 week", date_labels = "%b%d")+
  geom_line(data=fitlin2, aes(y=p*100, x=x))+
  geom_ribbon(data=fitlin2, aes(y=p*100, x=x, ymin=lb*100, ymax=ub*100),alpha=0.2)+
  geom_line(data=fitlin2, aes(y=p*100, x=x-10),color='red')+
  geom_ribbon(data=fitlin2, aes(y=p*100, x=x-10, ymin=lb*100, ymax=ub*100),alpha=0.2, fill='red')+
  geom_line(data=fitlin2, aes(y=p*100, x=x-14),color='blue')+
  geom_ribbon(data=fitlin2, aes(y=p*100, x=x-14, ymin=lb*100, ymax=ub*100),alpha=0.2, fill='blue')+
  geom_line(data=fitlin2, aes(y=p*100, x=x-6),color='dark green')+
  geom_ribbon(data=fitlin2, aes(y=p*100, x=x-6, ymin=lb*100, ymax=ub*100),alpha=0.2, fill='dark green')+
  scale_y_continuous(sec.axis = sec_axis(~.*560000/11.2, name="Incidence"))


################################################################



find_minimum_point <- function(ff, X, num_basis, spline_degree, link_function, number_extra){
  X <- seq(min(X), max(X), 1)
  Y_array <- array(data=NA, dim=c(nrow(ff$a), length(X)))
  B_true <- t(bs(X, df=num_basis, degree = spline_degree, intercept=TRUE))
  for(i in seq_len(nrow(ff$a))){
    a <- array(NA, num_basis)
    
    for(j in seq_len(length(a))){
      a[j] <- ff$a[i,j]
    }
    Y_array[i,] <- as.vector(a%*%B_true)
  }
  
  index_list<-c()
  max_list <- c()
  for(i in seq_len(nrow(Y_array))){
    min_index<- which.max(Y_array[i,(number_extra+1):(ncol(Y_array)-number_extra)])
    index_list<-c(index_list, min_index)
    max_list <- c(max_list, max(Y_array[i,(number_extra+1):(ncol(Y_array)-number_extra)]))
  }
  
  index_list
}

index_list <- find_maximum_point2(ff,X_REACT,num_basis, spline_degree, link_function = "logit", number_extra = 0)

index_list <- as.Date(index_list-1, origin = as.Date("2020-05-01"))

hist(index_list, breaks="days")
check<-as.data.frame(index_list)
plot5<- ggplot(check,aes(x=index_list, y=..density..))+
  coord_cartesian(xlim=c(as.Date("2021-12-15"),as.Date("2022-01-30")))+
  geom_histogram(binwidth=1,alpha=0.5, color='black')+
  theme_bw(base_size = 18)+
  xlab("Date of maximum prevalence")+
  ylab("Probability density")
plot5
median(index_list)-7
index_list2 <- as.numeric(index_list)
quantile(index_list2, c(0.5,0.025,0.975))
#2022-01-02
#2021-12-30
#2022-01-11


death_eng_all<-readRDS("E:/Group/report/round17/death_eng_all.rds")
hosp_eng_all<-readRDS("E:/Group/report/round17/hosp_eng_all.rds")

pdf(paste("E:/group/report/round17/Date_of_max_prev.pdf", sep=""), width = 6, height =6)
plot5
dev.off()


dfRes<-dfRes[is.na(dfRes$estbinres)==FALSE,]

tab<-table(dfRes$round, dfRes$estbinres)
df1 <- data.frame(tab)

dfRes<-dfRes[is.na(dfRes$estbinres)==FALSE &is.na(dfRes$d_comb)==FALSE,]

tab<-table(dfRes$round, dfRes$estbinres)
df2 <- data.frame(tab)
df2

df<-data.frame(round = seq(1,19,by=1),
           total = df1$Freq[1:19]+df1$Freq[20:38],
           pos = df1$Freq[20:38],
           d_total =df2$Freq[1:19]+df2$Freq[20:38],
           d_pos = df2$Freq[20:38])
df$p_total = df$d_total*100/df$total
df$p_pos = df$d_pos*100/df$pos

df
write.csv(df, "thesis_intro.csv")










##################################################################################################
#
predB1 <- readRDS("E:/Group/report/round_19_lineage/lineage_prop_ba2.rds")
dfY <- readRDS("E:/Group/report/round_19_lineage/fgfgdhf2.rds")

plot3 <- ggplot(data = df_plot, aes(x= d_comb, y =p*100))+
  coord_cartesian(ylim=c(0.0,10), xlim=c(as.Date("2021-11-15"), max(df_plot$d_comb)))+
  theme_bw(base_size = 18)+ 
  xlab("Date")+
  ylab("Prevalence (%)")+
  scale_x_date(date_breaks = "1 month", date_labels = "%b\n%Y")+
  geom_line(data= df_plot_model3, aes(y=inv_logit(p)*100))+
  geom_ribbon(data = df_plot_model3,
              aes(y=inv_logit(p)*100,
                  ymin=inv_logit(lb_2.5)*100,
                  ymax=inv_logit(ub_97.5)*100),
              alpha=0.2)+
  geom_ribbon(data = df_plot_model3,
              aes(y=inv_logit(p)*100,
                  ymin=inv_logit(lb_25)*100,
                  ymax=inv_logit(ub_75)*100),
              alpha=0.2)+
  geom_line(data = predB1, aes(x= d_comb, y = inv_logit(p)*10, ymin=inv_logit(lb)*10, ymax=inv_logit(ub)*10, color='blue'))+
  geom_ribbon(data = predB1,  aes(x= d_comb, y =inv_logit(p)*10, ymin=inv_logit(lb)*10, ymax=inv_logit(ub)*10,fill='blue'),alpha = 0.2)+
  geom_ribbon(data = dfY[dfY$lineage=="Omicron",],  aes(x= X, y = p*10, ymin=lb_2.5*10, ymax=ub_97.5*10,fill='red'),alpha = 0.2)+
  geom_line(data = dfY[dfY$lineage=="Omicron",], aes(x= X, y = p*10, ymin=lb_2.5*10, ymax=ub_97.5*10, color='red'))+
  scale_y_continuous(sec.axis = sec_axis(~./10, breaks= c(0,0.2,0.4,0.6,0.8,1.0), name = "Proportion of variant"))+
  scale_color_manual(values=c("red","blue"),
                     name= "Variant",
                     labels=c("BA.2","All Omicron"))+
  scale_fill_manual(values=c("red","blue"),
                     name= "Variant",
                     labels=c("BA.2","All Omicron"))+
  theme(legend.position = c(0.8,0.1),
        legend.background = element_rect(color='black'),
        panel.grid = element_blank())
plot3


pdf(paste( "E:/group/report/round19/Fig0.pdf", sep=""), width = 8, height = 8)
print(plot3)
dev.off()

