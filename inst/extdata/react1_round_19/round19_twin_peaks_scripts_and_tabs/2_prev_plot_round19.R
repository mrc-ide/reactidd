#' # Weighted & unweighted (when the amount of data is not enought to run weighted analysis)  

#' prevalence plots for REACT-1
#' Author: Haowei Wang
#' 
#' html generated with the knitr package by spin(<this file>)
#'
#' ## Preamble and script setup
Sys.time()
knitr.table.format = "markdown"

#' First clear the environment of variables
rm(list=ls(all=TRUE))
setwd("E:/Group")

#' Pull in packages needed
source("E:/Group/functions/load_packages.R")

pkgs <- c("prevalence", "mgcv", "mgcViz", "MASS", "dplyr",
          "tidyr", "forcats", "ggplot2", "qpcR","survey",
          "readr")
load_packages(pkgs)


#' Set output directory
setwd("E:/Group/report/round19")

#' ## Load data and add variables to working dataframe
output_tag <- "03_Apr"

col2 <-  c("#CCCCCC","#FD8D3C")

col4 <-  c("#FDD0A2", "#FDAE6B", "#F16913", "#D94801")

wt_out <- readRDS(paste0("wt_prev_outputs_for_plots_",output_tag,".rds"))
prev_tab_r19_interim <- wt_out[[1]]

# uw_out_r18 <- read_csv("E:/Group/report/round18/ut_prev_tab_raw_06_Mar.csv") %>%
#   filter(Variable %in% c("age_group_char3", "region"))
# uw_out_r19 <- read_csv(paste0("ut_prev_tab_raw_",output_tag,".csv")) %>%
#   filter(Variable %in% c("age_group_char3", "region"))
# 
# uw_out <- full_join(uw_out_r18, uw_out_r19, by = c("Variable", "Category"))

########## Weighted plots - r18r19 ###########
#' by age
df_wt_age <- prev_tab_r19_interim  %>%
  filter(level %in% c("5-11", "12-17", "18-24", "25-34",
                      "35-44", "45-54", "55-64", "65-74","75+")) %>%
  dplyr::select(level,
                prevalence_18, upper_18, lower_18,
                prevalence_19, upper_19, lower_19) %>%
  pivot_longer(., cols = 2:7, names_to = "stat", values_to = "value") %>%
  separate(stat, c("stat", "round"), sep = "_") %>%
  mutate(round = case_when(
    round == "18" ~ "Round 18",
    round == "19" ~ "Round 19"),
    level = factor(level, levels = c("5-11", "12-17", "18-24", "25-34",
                                     "35-44", "45-54", "55-64", "65-74","75+"))
  ) %>%
  pivot_wider(names_from = "stat", values_from = "value") %>%
  mutate(round = factor(round, levels = c("Round 18", "Round 19")),
         prevalence = prevalence*100,
         lower = lower*100,
         upper = upper*100) %>%
  rename(Round = round)


age_wt <- ggplot(data = df_wt_age,
                 aes(x = level,
                     y = prevalence,
                     fill=Round)) +
  geom_bar(stat= "identity",color="black", position = position_dodge2())+
  geom_errorbar(aes(x=level,
                    y=prevalence,
                    ymin = lower,
                    ymax = upper),
                width = .2, position = position_dodge(.9),color = "grey60")+
  xlab("Age Group") +
  ylab("Weighted Prevalence (%)") +
  theme_bw()+
  # scale_fill_manual(values= col2)+
  scale_fill_manual(values= col2, 
                    labels = c("Round 18 (08 Feb to 01 Mar 2022)", "Round 19 (08 Mar to 31 Mar 2022)"),
                    name = "Round")


age_wt + theme(legend.title = element_text(size = 14),
               legend.text = element_text(size = 14),
               text = element_text(size = 14),
               legend.position = c(0.8,0.85))

# save file (png and pdf)
ggsave(filename = paste0("wt_prev_age_bar_plot_",output_tag,".png"), width = 12, height = 8, dpi = 320,
       scale = 1)
ggsave(filename = paste0("wt_prev_age_bar_plot_",output_tag,".pdf"), width = 12, height = 8, dpi = 320,
       scale = 1)


#' by region
df_wt_region <- prev_tab_r19_interim %>%
  filter(level %in% c("North West", "Yorkshire and The Humber", "North East",
                      "West Midlands","East Midlands", "East of England",
                      "South West", "London", "South East")) %>%
  dplyr::select(level,                 
                prevalence_18, upper_18, lower_18,
                prevalence_19, upper_19, lower_19) %>%
  pivot_longer(., cols = 2:7, names_to = "stat", values_to = "value") %>%
  separate(stat, c("stat", "round"), sep = "_") %>%
  mutate(round = case_when(
    round == "18" ~ "Round 18",
    round == "19" ~ "Round 19"),
    level = factor(level, levels = c("North West", "North East", "Yorkshire and The Humber", 
                                     "West Midlands", "East Midlands", "East of England",
                                     "South West", "London", "South East"))
  ) %>%
  pivot_wider(names_from = "stat", values_from = "value") %>%
  mutate(round = factor(round, levels = c("Round 18", "Round 19")),
         prevalence = prevalence*100,
         lower = lower*100,
         upper = upper*100) %>%
  rename(Round = round)


region_wt <- ggplot(data = df_wt_region,
                    aes(x = level,
                        y = prevalence,
                        fill=Round)) +
  geom_bar(stat= "identity", position = position_dodge( ),color="black")+
  geom_errorbar(aes(x=level,
                    y=prevalence,
                    ymin = lower,
                    ymax = upper),
                width = .2, position = position_dodge(.9),color = "grey60")+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  coord_cartesian(ylim = c(0,10))+
  # scale_y_continuous(breaks = c(0.0, 0.5, 1.0, 1.5, 2.0, 2.5),
  #                    labels = c(0.0, 0.5, 1.0, 1.5, 2.0, 2.5))+
  xlab("Region") +
  ylab("Weighted Prevalence (%)") +
  theme_bw()+
  scale_fill_manual(values= col2, 
                    labels = c("Round 18 (08 Feb to 01 Mar 2022)", "Round 19 (08 Mar to 31 Mar 2022)"),
                    name = "Round")

region_wt + theme(legend.title = element_text(size = 14),
                  legend.text = element_text(size = 14),
                  text = element_text(size = 14),
                  legend.position = c(0.80,0.91))


# save file (png and pdf)
ggsave(filename = paste0("wt_prev_region_bar_plot_",output_tag,".png"), width = 13, height = 8, dpi = 320,
       scale = 1)
ggsave(filename = paste0("wt_prev_region_bar_plot_",output_tag,".pdf"), width = 13, height = 8, dpi = 320,
       scale = 1)

########### round 16 to 19 ########
wt_r1617_rds <- readRDS("E:/Group/report/round17/wt_prev_outputs_for_plots_25_Jan.rds")
wt_out_r1617 <-wt_r1617_rds[[1]] 

prev_tab_r16r19 <- left_join(wt_out_r1617, prev_tab_r19_interim, by = "level")

#' by age groups
df_wt_age_r16r19 <- prev_tab_r16r19 %>%
  filter(level %in% c("5-11", "12-17", "18-24", "25-34",
                      "35-44", "45-54", "55-64", "65-74","75+")) %>%
  dplyr::select(level,
                prevalence_16, upper_16, lower_16,
                prevalence_17, upper_17, lower_17,
                prevalence_18, upper_18, lower_18,
                prevalence_19, upper_19, lower_19) %>%
  pivot_longer(., cols = 2:13, names_to = "stat", values_to = "value") %>%
  separate(stat, c("stat", "round"), sep = "_") %>%
  mutate(round = case_when(
    round == "16" ~ "Round 16",
    round == "17" ~ "Round 17",
    round == "18" ~ "Round 18",
    round == "19" ~ "Round 19"),
    level = factor(level, levels = c("5-11", "12-17", "18-24", "25-34",
                                     "35-44", "45-54", "55-64", "65-74","75+")),
    round = factor(round, levels = c("Round 16","Round 17", "Round 18", "Round 19"))
  ) %>%
  pivot_wider(names_from = "stat", values_from = "value")%>%
  mutate(prevalence = prevalence*100,
         lower = lower*100,
         upper = upper*100) %>%
  rename(Round = round)


age_wt_r16r19 <- ggplot(data = df_wt_age_r16r19,
                        aes(x = level,
                            y = prevalence,
                            fill=Round)) +
  geom_bar(stat= "identity",color="black", position = position_dodge2())+
  geom_errorbar(aes(x=level,
                    y=prevalence,
                    ymin = lower,
                    ymax = upper),
                width = .2, position = position_dodge(.9),color = "grey60")+
  # coord_cartesian(ylim = c(0,10))+
  # scale_y_continuous(breaks = c(0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3, 3.5, 4.0),
  #                    labels = c(0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3, 3.5,4.0))+
  xlab("Age Group") +
  ylab("Weighted Prevalence (%)") +
  theme_bw()+
  scale_fill_manual(values= col4, 
                    labels = c("Round 16 (23 Nov to 14 Dec 2021)", "Round 17 (05 Jan to 20 Jan 2022)",
                               "Round 18 (08 Feb to 01 Mar 2022)", "Round 19 (08 Mar to 31 Mar 2022)"),
                    name = "Round")



age_wt_r16r19 + theme(legend.title = element_text(size = 14),
                      legend.text = element_text(size = 14),
                      text = element_text(size = 14),
                      legend.position = c(0.83,0.88))

# save file (png and pdf)
ggsave(filename = paste0("wt_prev_age_bar_plot_r16r19_",output_tag,".png"), width = 13, height = 9, dpi = 320,
       scale = 1)
ggsave(filename = paste0("wt_prev_age_bar_plot_r16r19_",output_tag,".pdf"), width = 13, height = 9, dpi = 320,
       scale = 1)

#' by region

df_wt_region_r16r19 <-  prev_tab_r16r19 %>%
  filter(level %in% c("North West", "Yorkshire and The Humber", "North East",
                      "West Midlands","East Midlands", "East of England",
                      "South West", "London", "South East")) %>%
  dplyr::select(level,
                prevalence_16, upper_16, lower_16,
                prevalence_17, upper_17, lower_17,
                prevalence_18, upper_18, lower_18,
                prevalence_19, upper_19, lower_19) %>%
  pivot_longer(., cols = 2:13, names_to = "stat", values_to = "value") %>%
  separate(stat, c("stat", "round"), sep = "_") %>%
  mutate(round = case_when(
    round == "16" ~ "Round 16",
    round == "17" ~ "Round 17",
    round == "18" ~ "Round 18",
    round == "19" ~ "Round 19"),
    level = factor(level, levels = c("North West", "North East", "Yorkshire and The Humber", 
                                     "West Midlands", "East Midlands", "East of England",
                                     "South West", "London", "South East")),
    round = factor(round, levels = c("Round 16","Round 17", "Round 18", "Round 19"))
  ) %>%
  pivot_wider(names_from = "stat", values_from = "value") %>%
  mutate(prevalence = prevalence*100,
         lower = lower*100,
         upper = upper*100) %>%
  rename(Round = round)

region_wt_r16r19 <- ggplot(data = df_wt_region_r16r19,
                           aes(x = level,
                               y = prevalence,
                               fill=Round)) +
  geom_bar(stat= "identity",color="black", position = position_dodge2())+
  #scale_fill_gradient(low = "blue",high = "red")+
  geom_errorbar(aes(x=level,
                    y=prevalence,
                    ymin = lower,
                    ymax = upper),
                width = .2, position = position_dodge(.9),color = "grey60")+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  coord_cartesian(ylim = c(0,10))+
  # scale_y_continuous(breaks = c(0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3, 3.5, 4.0),
  #                    labels = c(0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3, 3.5,4.0))+
  xlab("Region") +
  ylab("Weighted Prevalence (%)") +
  theme_bw()+
  scale_fill_manual(values= col4, 
                    labels = c("Round 16 (23 Nov to 14 Dec 2021)", "Round 17 (05 Jan to 20 Jan 2022)",
                               "Round 18 (08 Feb to 01 Mar 2022)", "Round 19 (08 Mar to 31 Mar 2022)"),
                    name = "Round")


region_wt_r16r19 + theme(legend.title = element_text(size = 14),
                         legend.text = element_text(size = 14),
                         text = element_text(size = 14),
                         legend.position = c(0.86,0.92))


# save file (png and pdf)
ggsave(filename = paste0("wt_prev_region_bar_plot_r16r19_",output_tag,".png"), width = 15, height = 10, dpi = 320,
       scale = 1)
ggsave(filename = paste0("wt_prev_region_bar_plot_r16r19_",output_tag,".pdf"), width = 15, height = 10, dpi = 320,
       scale = 1)



########## Unweighted plots - r18r19 ###########

#' by age groups
df_uw_age <- uw_out %>%
  dplyr::select(Variable, Category,
                Prevalence_r18, Upper_r18, Lower_r18,
                Prevalence_r19, Upper_r19, Lower_r19) %>%
  filter(Variable == "age_group_char3") %>%
  pivot_longer(., cols = 3:8, names_to = "stat", values_to = "value") %>%
  separate(stat, c("stat", "round"), sep = "_") %>%
  mutate(round = case_when(
    round == "r18" ~ "Round 18",
    round == "r19" ~ "Round 19 partial"),
    Category = factor(Category, levels = c("5-11", "12-17", "18-24", "25-34",
                                           "35-44", "45-54", "55-64", "65-74","75+")),
    round = factor(round, levels = c("Round 18","Round 19 partial"))
  ) %>%
  pivot_wider(names_from = "stat", values_from = "value")%>%
  mutate(Prevalence = Prevalence*100,
         Lower = Lower*100,
         Upper = Upper*100)%>%
  rename(Round = round)


age_uw <- ggplot(data = df_uw_age,
                 aes(x = Category,
                     y = Prevalence,
                     fill=Round)) +
  geom_bar(stat= "identity",color="black", position = position_dodge2())+
  geom_errorbar(aes(x=Category,
                    y=Prevalence,
                    ymin = Lower,
                    ymax = Upper),
                width = .2, position = position_dodge(.9),color = "grey60")+
  coord_cartesian(ylim = c(0,12))+
  scale_y_continuous(breaks = c(0, 3, 6, 9, 12),
                     labels = c(0, 3, 6, 9, 12))+
  xlab("Age Group") +
  ylab("Unweighted Prevalence (%)") +
  theme_bw()+
  scale_fill_manual(values= col2, 
                    labels = c("Round 18 (08 Feb to 01 Mar 2022)", "Round 19 partial (08 Mar to 15 Mar 2022)"),
                    name = "Round")


age_uw + theme(legend.title = element_text(size = 14),
               legend.text = element_text(size = 14),
               text = element_text(size = 14),
               legend.position = c(0.80,0.85))

# save file (png and pdf)
ggsave(filename = paste0("uw_prev_age_bar_plot_",output_tag,".png"), width = 12, height = 8, dpi = 320,
       scale = 1)
ggsave(filename = paste0("uw_prev_age_bar_plot_",output_tag,".pdf"), width = 12, height = 8, dpi = 320,
       scale = 1)

#' by region

df_uw_region <-  uw_out %>%
  filter(Variable == "region") %>%
  dplyr::select(Variable, Category,
                Prevalence_r18, Upper_r18, Lower_r18,
                Prevalence_r19, Upper_r19, Lower_r19) %>%
  pivot_longer(., cols = 3:8, names_to = "stat", values_to = "value") %>%
  separate(stat, c("stat", "round"), sep = "_") %>%
  mutate(round = case_when(
    round == "r18" ~ "Round 18",
    round == "r19" ~ "Round 19 partial"),
    Category = factor(Category, levels = c("North West", "North East", "Yorkshire and The Humber", 
                                           "West Midlands", "East Midlands", "East of England",
                                           "South West", "London", "South East")),
    round = factor(round, levels = c("Round 18","Round 19 partial"))
  ) %>%
  pivot_wider(names_from = "stat", values_from = "value") %>%
  mutate(Prevalence = Prevalence*100,
         Lower = Lower*100,
         Upper = Upper*100)%>%
  rename(Round = round)



region_uw <- ggplot(data = df_uw_region,
                    aes(x = Category,
                        y = Prevalence,
                        fill=Round)) +
  geom_bar(stat= "identity",color="black", position = position_dodge2())+
  #scale_fill_gradient(low = "blue",high = "red")+
  geom_errorbar(aes(x=Category,
                    y=Prevalence,
                    ymin = Lower,
                    ymax = Upper),
                width = .2, position = position_dodge(.9),color = "grey60")+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  coord_cartesian(ylim = c(0,8))+
  # scale_y_continuous(breaks = c(0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3, 3.5, 4.0),
  #                    labels = c(0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3, 3.5,4.0))+
  xlab("Region") +
  ylab("Unweighted Prevalence (%)") +
  theme_bw()+
  scale_fill_manual(values= col2, 
                    labels = c("Round 18 (08 Feb to 01 Mar 2022)", "Round 19 partial (08 Mar to 15 Mar 2022)"),
                    name = "Round")


region_uw + theme(legend.title = element_text(size = 14),
                  legend.text = element_text(size = 14),
                  text = element_text(size = 14),
                  legend.position = c(0.80,0.91))


# save file (png and pdf)
ggsave(filename = paste0("uw_prev_region_bar_plot_",output_tag,".png"), width = 12, height = 8, dpi = 320,
       scale = 1)
ggsave(filename = paste0("uw_prev_region_bar_plot_",output_tag,".pdf"), width = 12, height = 8, dpi = 320,
       scale = 1)

