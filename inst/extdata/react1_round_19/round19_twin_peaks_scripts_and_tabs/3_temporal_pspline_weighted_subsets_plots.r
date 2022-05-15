#!/usr/bin/env Rscript
Sys.time()
rm(list=ls(all=TRUE))
setwd("E:/Group/Ollies_workspace/")

path_name <- "E://Group/figure/temporal/MCMC_model_gam/National/Estbinres/England/Joined/"
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

## READ IN THE DATA
reg_YH<-readRDS("E:/Group/Ollies_workspace/pspline_func_runs/reg_YH.rds")
reg_NW<-readRDS("E:/Group/Ollies_workspace/pspline_func_runs/reg_NW.rds")
reg_NE<-readRDS("E:/Group/Ollies_workspace/pspline_func_runs/reg_NE.rds")
reg_SW<-readRDS("E:/Group/Ollies_workspace/pspline_func_runs/reg_SW.rds")
reg_SE<-readRDS("E:/Group/Ollies_workspace/pspline_func_runs/reg_SE.rds")
reg_EE<-readRDS("E:/Group/Ollies_workspace/pspline_func_runs/reg_EE.rds")
reg_LN<-readRDS("E:/Group/Ollies_workspace/pspline_func_runs/reg_LN.rds")
reg_EM<-readRDS("E:/Group/Ollies_workspace/pspline_func_runs/reg_EM.rds")
reg_WM<-readRDS("E:/Group/Ollies_workspace/pspline_func_runs/reg_WM.rds")


###############################################
reg_NW[[6]]$reg_char <- "North West"
reg_NE[[6]]$reg_char <- "North East"
reg_SE[[6]]$reg_char <- "South East"
reg_SW[[6]]$reg_char <- "South West"
reg_YH[[6]]$reg_char <- "Yorkshire and The Humber"
reg_LN[[6]]$reg_char <- "London"
reg_EE[[6]]$reg_char <- "East of England"
reg_EM[[6]]$reg_char <- "East Midlands"
reg_WM[[6]]$reg_char <- "West Midlands"

reg_plot <- rbind(reg_NW[[6]],
                  reg_NE[[6]],
                  reg_SE[[6]],
                  reg_SW[[6]],
                  reg_LN[[6]],
                  reg_EE[[6]],
                  reg_EM[[6]],
                  reg_WM[[6]],
                  reg_YH[[6]])


reg_NW[[5]]$reg_char <- "North West"
reg_NE[[5]]$reg_char <- "North East"
reg_SE[[5]]$reg_char <- "South East"
reg_SW[[5]]$reg_char <- "South West"
reg_YH[[5]]$reg_char <- "Yorkshire and The Humber"
reg_LN[[5]]$reg_char <- "London"
reg_EE[[5]]$reg_char <- "East of England"
reg_EM[[5]]$reg_char <- "East Midlands"
reg_WM[[5]]$reg_char <- "West Midlands"

reg_data <- rbind(reg_NW[[5]],
                  reg_NE[[5]],
                  reg_SE[[5]],
                  reg_SW[[5]],
                  reg_LN[[5]],
                  reg_EE[[5]],
                  reg_EM[[5]],
                  reg_WM[[5]],
                  reg_YH[[5]])

reg_plot <- as.data.frame(reg_plot)

reg_plot$reg_char <- factor(reg_plot$reg_char, levels= c("North West", "Yorkshire and The Humber", "North East", "West Midlands",
                                                         "East Midlands", "East of England","South West", "London", "South East" ))
reg_data$reg_char <- factor(reg_data$reg_char, levels= c("North West", "Yorkshire and The Humber", "North East", "West Midlands",
                                                         "East Midlands", "East of England","South West", "London", "South East" ))


mindate1<-as.Date("2020-05-01")
mindate8<-as.Date("2021-1-06")
mindate14 <- as.Date("2021-09-09")


maxdate7<-as.Date("2020-12-03")
maxdate13<-as.Date("2021-07-12")


reg_plot1<-reg_plot[reg_plot$d_comb >=mindate1 &reg_plot$d_comb<=maxdate7,]
reg_plot2<-reg_plot[reg_plot$d_comb >=mindate8 &reg_plot$d_comb<=maxdate13,]
reg_plot3<-reg_plot[reg_plot$d_comb >=mindate14,]


reg_plot_1 <- ggplot(data = reg_plot, aes(x = d_comb, y=inv_logit(p)*100, ymin=inv_logit(lb_2.5)*100, ymax=inv_logit(ub_97.5)*100))+
  geom_line(aes( col = reg_char))+
  geom_ribbon(aes(fill = reg_char), alpha=0.1)+
  coord_cartesian(ylim=c(0,3.0), xlim=c(min(reg_plot$d_comb), max(reg_plot$d_comb)))+
  theme_bw(base_size = 18)+ 
  xlab("Day of swab")+
  ylab("Prevalence (%)")+
  scale_x_date(date_breaks = "1 month", date_labels = "%b")+
  labs(fill="Region", col="Region")

reg_plot_2 <- ggplot(data = reg_plot, aes(x = d_comb, y=log10(inv_logit(p)), ymin=log10(inv_logit(lb_2.5)), ymax=log10(inv_logit(ub_97.5))))+
  geom_line(aes( col = reg_char))+
  geom_ribbon(aes(fill = reg_char), alpha=0.1)+
  coord_cartesian(ylim=c(-4,-1.5), xlim=c(min(reg_plot$d_comb), max(reg_plot$d_comb)))+
  theme_bw(base_size = 18)+ 
  xlab("Day of swab")+
  ylab("Log(Prevalence)")+
  scale_x_date(date_breaks = "1 month", date_labels = "%b")+
  labs(fill="Region", col="Region")


png(paste(path_name, "reg_plot_1.png", sep=""), width = 1000, height =500)
reg_plot_1
dev.off()
pdf(paste(path_name, "reg_plot_1.pdf", sep=""), width = 12, height =6)
print(reg_plot_1)
dev.off()

png(paste(path_name, "reg_plot_2.png", sep=""), width = 1000, height =500)
print(reg_plot_2)
dev.off()
pdf(paste(path_name, "reg_plot_2.pdf", sep=""), width = 12, height =6)
print(reg_plot_2)
dev.off()

reg_plot_3 <- ggplot(data = reg_plot, aes(x = d_comb, y=inv_logit(p)*100, ymin=inv_logit(lb_2.5)*100, ymax=inv_logit(ub_97.5)*100))+
  geom_line(data=reg_plot1,aes())+
  geom_ribbon(data=reg_plot1,aes(), alpha=0.3)+
  geom_line(data=reg_plot2,aes())+
  geom_ribbon(data=reg_plot2,aes(), alpha=0.3)+
  geom_line(data=reg_plot3,aes())+
  geom_ribbon(data=reg_plot3,aes(), alpha=0.3)+
  coord_cartesian(ylim=c(0.3,5.0), xlim=c(mindate14, max(reg_plot$d_comb)))+
  theme_bw(base_size = 18)+ 
  xlab("Day of swab")+
  ylab("Prevalence (%)")+
  scale_x_date(date_breaks = "1 month", date_labels = "%b")+
  scale_y_log10()+
  labs(fill="Region", col="Region")+
  facet_wrap(vars(reg_char), ncol=3)+
  geom_point(data = reg_data, aes(x=d_comb, y=p*100, ymin=lb*100, ymax=ub*100), alpha=0.5)+
  geom_errorbar(data = reg_data, aes(x=d_comb, y=p*100, ymin=lb*100, ymax=ub*100), alpha=0.5)+
  theme(legend.position = "none")

mindate14 <- as.Date("2021-09-09")
reg_plot3
reg_plot_3<-ggplot(data = reg_plot3, aes(x = d_comb, y=inv_logit(p)*100, ymin=inv_logit(lb_2.5)*100, ymax=inv_logit(ub_97.5)*100))+
  geom_line(data=reg_plot3)+
  geom_ribbon(alpha=0.3)+
  geom_ribbon(aes(x = d_comb, y=inv_logit(p)*100, ymin=inv_logit(lb_25)*100, ymax=inv_logit(ub_75)*100),alpha=0.3)+  
  coord_cartesian(ylim=c(0.1,10.0), xlim=c(as.Date("2021-11-20"), max(reg_plot$d_comb)))+
  theme_bw(base_size = 18)+ 
  xlab("Day of swab")+
  ylab("Prevalence (%)")+
  scale_y_log10()+
  scale_x_date(date_breaks = "1 month", date_labels = "%b")+
  facet_wrap(vars(reg_char), ncol=3)+
  geom_point(data = reg_data, aes(x=d_comb, y=p*100, ymin=lb*100, ymax=ub*100), alpha=0.5)+
  geom_errorbar(data = reg_data, aes(x=d_comb, y=p*100, ymin=lb*100, ymax=ub*100), alpha=0.5)+
  theme(legend.position = "none",
        panel.grid.minor.y = element_blank())
reg_plot_3


png(paste(path_name, "reg_plot_3.png", sep=""), width = 1000, height =600)
reg_plot_3
dev.off()
pdf(paste("E:/group/report/round17/", "reg_plot_3.pdf", sep=""), width = 10, height =12)
reg_plot_3
dev.off()



##############################################################################
#Age tertile spline plots
age_1<-readRDS("E:/Group/Ollies_workspace/pspline_func_runs/age_new1.rds")
age_2<-readRDS("E:/Group/Ollies_workspace/pspline_func_runs/age_new2.rds")
age_3<-readRDS("E:/Group/Ollies_workspace/pspline_func_runs/age_new3.rds")
age_4<-readRDS("E:/Group/Ollies_workspace/pspline_func_runs/age_new4.rds")


###############################################
age_1[[6]]$age_char <- "17 and under"
age_2[[6]]$age_char <- "18 to 34"
age_3[[6]]$age_char <- "55 and over"
age_4[[6]]$age_char <- "35 to 54"


age_plot <- rbind(age_1[[6]],
                  age_2[[6]],
                  age_3[[6]],
                  age_4[[6]])

age_1[[5]]$age_char <- "17 and under"
age_2[[5]]$age_char <- "18 to 34"
age_3[[5]]$age_char <- "55 and over"
age_4[[5]]$age_char <- "35 to 54"

age_data <- rbind(age_1[[5]],
                  age_2[[5]],
                  age_3[[5]],
                  age_4[[5]])

age_plot <- as.data.frame(age_plot)

age_plot$age_char <- factor(age_plot$age_char, levels= c("17 and under","18 to 34","35 to 54","55 and over"))
age_data$age_char <- factor(age_data$age_char, levels= c("17 and under","18 to 34","35 to 54","55 and over" ))


mindate1<-as.Date("2020-05-01")
mindate8<-as.Date("2021-1-06")



maxdate7<-as.Date("2020-12-03")
maxdate13 <- as.Date("2021-07-12")

mindate14 <- as.Date("2021-09-09")

#age_plot3<-age_plot[age_plot$d_comb >=mindate1 &age_plot$d_comb<=maxdate7,]
#age_plot2<-age_plot[age_plot$d_comb >=mindate8 &age_plot$d_comb<=maxdate13,]
age_plot1<-age_plot[age_plot$d_comb >=mindate14 &age_plot$d_comb<=max(age_plot$d_comb),]

age_plot_4<-ggplot(data = age_plot1, aes(x = d_comb, y=inv_logit(p)*100, ymin=inv_logit(lb_2.5)*100, ymax=inv_logit(ub_97.5)*100))+
  geom_line(aes(color= age_char))+
  geom_ribbon(alpha=0.3, aes( fill = age_char))+
  geom_ribbon(aes(x = d_comb, y=inv_logit(p)*100, ymin=inv_logit(lb_25)*100, ymax=inv_logit(ub_75)*100, fill = age_char),alpha=0.3)+  
  coord_cartesian(ylim=c(0.1,30), xlim=c(as.Date("2021-11-20"), max(age_plot1$d_comb)))+
  theme_bw(base_size = 18)+ 
  xlab("Day of swab")+
  ylab("Prevalence (%)")+
  scale_x_date(date_breaks = "1 month", date_labels = "%b")+
  scale_y_log10()+
  scale_color_brewer(type = "q", palette = 6)+
  scale_fill_brewer(type = "q", palette = 6)+
  labs(color = "Age", fill = "Age")+
  geom_point(data = age_data, aes(x=d_comb, y=p*100, ymin=lb*100, ymax=ub*100, color=age_char), alpha=0.5)+
  geom_errorbar(data = age_data, aes(x=d_comb, y=p*100, ymin=lb*100, ymax=ub*100, color=age_char), alpha=0.5)+
  facet_wrap(vars(age_char), ncol=2)+
  theme(legend.position = "bottom",
        strip.text = element_blank(),
        panel.grid.minor = element_blank())
age_plot_4

pdf( "E:/group/report/round19/age_temporal.pdf", width = 10, height =10)
age_plot_4
dev.off()



##############################################################################
#Age tertile spline plots
age_1<-readRDS("E:/Group/Ollies_workspace/pspline_func_runs/age54_nchild0.rds")
age_2<-readRDS("E:/Group/Ollies_workspace/pspline_func_runs/age54_nchild1.rds")
age_3<-readRDS("E:/Group/Ollies_workspace/pspline_func_runs/age55_nchild0.rds")
age_4<-readRDS("E:/Group/Ollies_workspace/pspline_func_runs/age55_nchild1.rds")

###############################################
age_1[[6]]$age_char <- "18-54: No children"
age_2[[6]]$age_char <- "18-54: Children"
age_3[[6]]$age_char <- "55+: No children"
age_4[[6]]$age_char <- "55+: Children"


age_plot <- rbind(age_1[[6]],
                  age_2[[6]],
                  age_3[[6]],
                  age_4[[6]])

age_1[[5]]$age_char <- "18-54: No children"
age_2[[5]]$age_char <- "18-54: Children"
age_3[[5]]$age_char <- "55+: No children"
age_4[[5]]$age_char <- "55+: Children"


age_data <- rbind(age_1[[5]],
                  age_2[[5]],
                  age_3[[5]],
                  age_4[[5]])

age_plot <- as.data.frame(age_plot)

age_plot$age_char <- factor(age_plot$age_char, levels= c("18-54: No children","18-54: Children",
                                                         "55+: No children", "55+: Children"))
age_data$age_char <- factor(age_data$age_char, levels= c("18-54: No children","18-54: Children",
                                                         "55+: No children", "55+: Children"))


mindate1<-as.Date("2020-05-01")
mindate8<-as.Date("2021-1-06")



maxdate7<-as.Date("2020-12-03")
maxdate13 <- as.Date("2021-07-12")

mindate14 <- as.Date("2021-09-09")

#age_plot3<-age_plot[age_plot$d_comb >=mindate1 &age_plot$d_comb<=maxdate7,]
#age_plot2<-age_plot[age_plot$d_comb >=mindate8 &age_plot$d_comb<=maxdate13,]
age_plot1<-age_plot[age_plot$d_comb >=mindate14 &age_plot$d_comb<=max(age_plot$d_comb),]

age_plot_4<-ggplot(data = age_plot1, aes(x = d_comb, y=inv_logit(p)*100, ymin=inv_logit(lb_2.5)*100, ymax=inv_logit(ub_97.5)*100))+
  geom_line(aes(color= age_char))+
  geom_ribbon(alpha=0.3, aes( fill = age_char))+
  geom_ribbon(aes(x = d_comb, y=inv_logit(p)*100, ymin=inv_logit(lb_25)*100, ymax=inv_logit(ub_75)*100, fill = age_char),alpha=0.3)+  
  coord_cartesian(ylim=c(0.1,3), xlim=c(mindate14, max(age_plot1$d_comb)))+
  theme_bw(base_size = 18)+ 
  xlab("Day of swab")+
  ylab("Prevalence (%)")+
  scale_x_date(date_breaks = "2 week", date_labels = "%d %b")+
  scale_y_log10()+
  scale_color_brewer(type = "q", palette = 6)+
  scale_fill_brewer(type = "q", palette = 6)+
  labs(color = "Age", fill = "Age")+
  geom_point(data = age_data, aes(x=d_comb, y=p*100, ymin=lb*100, ymax=ub*100, color=age_char), alpha=0.5)+
  geom_errorbar(data = age_data, aes(x=d_comb, y=p*100, ymin=lb*100, ymax=ub*100, color=age_char), alpha=0.5)+
  facet_wrap(vars(age_char), ncol=2)+
  theme(legend.position = "bottom",
        strip.text = element_blank(),
        panel.grid.minor = element_blank())
age_plot_4

pdf( "E:/group/report/round15/adult_child_temporal.pdf", width = 10, height =10)
age_plot_4
dev.off()



########################################################################################################
##############################################################################
#Age tertile spline plots
age_1<-readRDS("E:/Group/Ollies_workspace/pspline_func_runs/age_new1a.rds")
age_2<-readRDS("E:/Group/Ollies_workspace/pspline_func_runs/age_new1b.rds")
age_3<-readRDS("E:/Group/Ollies_workspace/pspline_func_runs/age_new1c.rds")
age_4<-readRDS("E:/Group/Ollies_workspace/pspline_func_runs/age_new1d.rds")
age_5<-readRDS("E:/Group/Ollies_workspace/pspline_func_runs/age_new1e.rds")
age_6<-readRDS("E:/Group/Ollies_workspace/pspline_func_runs/age_new1f.rds")
age_7<-readRDS("E:/Group/Ollies_workspace/pspline_func_runs/age_new1g.rds")
age_8<-readRDS("E:/Group/Ollies_workspace/pspline_func_runs/age_new1h.rds")
age_9<-readRDS("E:/Group/Ollies_workspace/pspline_func_runs/age_new1i.rds")


###############################################
age_1[[6]]$age_char <- "5 to 11"
age_2[[6]]$age_char <- "12 to 17"
age_3[[6]]$age_char <- "18 to 24"
age_4[[6]]$age_char <- "25 to 34"
age_5[[6]]$age_char <- "35 to 44"
age_6[[6]]$age_char <- "45 to 54"
age_7[[6]]$age_char <- "55 to 64"
age_8[[6]]$age_char <- "65 to 74"
age_9[[6]]$age_char <- "75 and over"


age_plot <- rbind(age_1[[6]],
                  age_2[[6]],
                  age_3[[6]],
                  age_4[[6]],
                  age_5[[6]],
                  age_6[[6]],
                  age_7[[6]],
                  age_8[[6]],
                  age_9[[6]])



age_1[[5]]$age_char <- "5 to 11"
age_2[[5]]$age_char <- "12 to 17"
age_3[[5]]$age_char <- "18 to 24"
age_4[[5]]$age_char <- "25 to 34"
age_5[[5]]$age_char <- "35 to 44"
age_6[[5]]$age_char <- "45 to 54"
age_7[[5]]$age_char <- "55 to 64"
age_8[[5]]$age_char <- "65 to 74"
age_9[[5]]$age_char <- "75 and over"


age_data <- rbind(age_1[[5]],
                  age_2[[5]],
                  age_3[[5]],
                  age_4[[5]],
                  age_5[[5]],
                  age_6[[5]],
                  age_7[[5]],
                  age_8[[5]],
                  age_9[[5]])


age_plot <- as.data.frame(age_plot)

age_plot$age_char <- factor(age_plot$age_char, levels= c("5 to 11","12 to 17","18 to 24","25 to 34","35 to 44","45 to 54","55 to 64","65 to 74","75 and over"))
age_data$age_char <- factor(age_data$age_char, levels= c("5 to 11","12 to 17","18 to 24","25 to 34","35 to 44","45 to 54","55 to 64","65 to 74","75 and over"))


mindate1<-as.Date("2020-05-01")
mindate8<-as.Date("2021-1-06")



maxdate7<-as.Date("2020-12-03")
maxdate13 <- as.Date("2021-07-12")

mindate14 <- as.Date("2021-09-09")

#age_plot3<-age_plot[age_plot$d_comb >=mindate1 &age_plot$d_comb<=maxdate7,]
#age_plot2<-age_plot[age_plot$d_comb >=mindate8 &age_plot$d_comb<=maxdate13,]
age_plot1<-age_plot[age_plot$d_comb >=mindate14 &age_plot$d_comb<=max(age_plot$d_comb),]

age_plot_4<-ggplot(data = age_plot1, aes(x = d_comb, y=inv_logit(p)*100, ymin=inv_logit(lb_2.5)*100, ymax=inv_logit(ub_97.5)*100))+
  geom_line()+
  geom_ribbon(alpha=0.3 )+
  geom_ribbon(aes(x = d_comb, y=inv_logit(p)*100, ymin=inv_logit(lb_25)*100, ymax=inv_logit(ub_75)*100),alpha=0.3)+  
  coord_cartesian(ylim=c(0.1,30), xlim=c(as.Date("2021-11-20"), max(age_plot1$d_comb)))+
  theme_bw(base_size = 18)+ 
  xlab("Day of swab")+
  ylab("Prevalence (%)")+
  scale_x_date(date_breaks = "1 month", date_labels = "%b")+
  scale_y_log10()+
  geom_point(data = age_data, aes(x=d_comb, y=p*100, ymin=lb*100, ymax=ub*100), alpha=0.5)+
  geom_errorbar(data = age_data, aes(x=d_comb, y=p*100, ymin=lb*100, ymax=ub*100), alpha=0.5)+
  facet_wrap(vars(age_char), ncol=3)+
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())
age_plot_4

pdf( "E:/group/report/round19/age_temporalRev.pdf", width = 10, height =12)
age_plot_4
dev.off()


