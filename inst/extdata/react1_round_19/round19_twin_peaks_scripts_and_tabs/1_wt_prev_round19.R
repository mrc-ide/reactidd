#' # Weighted prevalence estimation for REACT-1 round 19
#' Author: Haowei Wang
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
          "readr", "scales")
load_packages(pkgs)

#' Source any functions from the local file
source("E:/Group/functions/add_conf_ints.R")
source("E:/Group/functions/make_tables.R")
source("E:/Group/functions/overall_prev.R")


#' ## Load data and add variables to working dataframe

#' Set output directory
setwd("E:/Group/report/round19/")

df_round19g <- readRDS("E:/Group/saved_objects/rep19.rds")


#' specify result variable
res_param <- "estbinres"

#' specify a tag for outputs
output_tag <- "03_Apr"


###########################################
#' ## Overall Prevalence - weighted
# remove NAs in estbinres

df_round19g_b <- df_round19g %>% filter(!is.na(estbinres)) %>% mutate(group = "Overall")


#' remove NAs in weights

df_round19g_b <- df_round19g_b %>% filter(!is.na(wt_antigen))
df_round19g_b$sympt_cat <- factor(df_round19g_b$sympt_cat, levels = c("Classic COVID symptoms", "Other symptoms",
                                                                      "No symptoms","NA"))


df_round19g_b  <- df_round19g_b  %>% 
  mutate(vax_status_cat = ifelse(is.na(vax_status_cat), "NA", vax_status_cat)) %>% 
  mutate(vax_status_cat = factor(vax_status_cat, levels = c("Not vaccinated", "One does", "Two does",
                                                            "Unknown does", "NA")),
         vax_status_noDate = factor(vax_status_noDate, levels = c("Not vaccinated", "One does", "Two does",
                                                                  "Unknown does", "NA")),
         vax_status_noDate_v2 = factor(vax_status_noDate_v2, levels = c("Not vaccinated", "One does", "Two does",
                                                                        "Three does", "Unknown does", "NA"))
         # vax_wane = factor(vax_wane, levels = c("Unvaccinated", "1 dose", "2 dose < 3 months",  "2 dose 3-6 months",
         #                                        "2 dose > 6 months",  "NA"))
  ) %>% 
  mutate(covidcon_char = ifelse(is.na(covidcon_char),"NA",as.character(covidcon_char))) %>% 
  mutate(covidcon_char = factor(covidcon_char,
                                levels = c("No", "Yes, contact with confirmed/tested COVID-19 case",
                                           "Yes, contact with suspected COVID-19 case","NA"))
  ) %>% 
  mutate(shield2_char = case_when(shield2 == 1 ~ "Yes",
                                  shield2 == 2 ~ "No",
                                  is.na(shield2) == TRUE ~ "sheiding - unknown"),
         indmask_char = case_when(indmask == 1 ~ "Always",
                                  indmask == 2 ~ "Sometimes",
                                  indmask == 3 ~ "Hardly ever",
                                  indmask == 4 ~ "Never",
                                  indmask %in% c(-92,-91,-77,5) == TRUE ~ "mask - unknown",
                                  is.na(indmask) == TRUE ~ "mask - unknown"))


dclus19g <- svydesign(id = ~id, strata = ~ lacode, weights = ~ wt_antigen, data = df_round19g_b, nest = TRUE)


##]####### overall prevalence ######
# wt_prev_o_r13 <-  svyby(~estbinres, by = ~group, design = dclus13, FUN = svyciprop, vartype = "ci") %>% rename(level = group)
wt_prev_o_r19 <-  svyby(~estbinres, by = ~group, design = dclus19g, FUN = svyciprop, vartype = "ci") %>% rename(level = group)

# wt_prev_o_r16Dec <-  svyby(~estbinres, by = ~group, design = dclus16g_Dec1To6, FUN = svyciprop, vartype = "ci") %>% rename(level = group)
# write.csv(wt_prev_o_r16Dec, paste0("wt_overall_prev_Dec1To6_",output_tag,".csv"), row.names = FALSE)

wt_prev_tab <- bind_rows(wt_prev_o_r19,
                         .id = "round") %>%
  rename(wt_prev = estbinres,
         lower = ci_l,
         upper = ci_u) %>%
  mutate(round = c("19"))

write.csv(wt_prev_tab, paste0("wt_overall_prev_tab_",output_tag,".csv"), row.names = FALSE)


############ core variables ##############

prev_tab_g_r19 <-  svyby(~estbinres, by = ~gender_char,     design = dclus19g, FUN = svyciprop, vartype = "ci") %>% rename(level = gender_char)
# prev_tab_a_r19 <-  svyby(~estbinres, by = ~age_group_char,  design = dclus19g, FUN = svyciprop, vartype = "ci") %>% rename(level = age_group_char)
prev_tab_a3_r19 <-  svyby(~estbinres, by = ~age_group_char3,  design = dclus19g, FUN = svyciprop, vartype = "ci") %>% rename(level = age_group_char3)
prev_tab_r_r19 <-  svyby(~estbinres, by = ~region,          design = dclus19g, FUN = svyciprop, vartype = "ci") %>% rename(level = region)
prev_tab_w_r19 <-  svyby(~estbinres, by = ~work_new_alt,    design = dclus19g, FUN = svyciprop, vartype = "ci") %>% rename(level = work_new_alt)
prev_tab_e_r19 <-  svyby(~estbinres, by = ~ethnic_new_char, design = dclus19g, FUN = svyciprop, vartype = "ci") %>% rename(level = ethnic_new_char)
prev_tab_hh_r19 <- svyby(~estbinres, by = ~hh_size_cat,     design = dclus19g, FUN = svyciprop, vartype = "ci") %>% rename(level = hh_size_cat)
prev_tab_c_r19 <-  svyby(~estbinres, by = ~covidcon_char,   design = dclus19g, FUN = svyciprop, vartype = "ci") %>% rename(level = covidcon_char)
prev_tab_h_r19 <-  svyby(~estbinres, by = ~hosp_char,       design = dclus19g, FUN = svyciprop, vartype = "ci") %>% rename(level = hosp_char)
prev_tab_s_r19 <-  svyby(~estbinres, by = ~sympt_cat,       design = dclus19g, FUN = svyciprop, vartype = "ci") %>% rename(level = sympt_cat)
prev_tab_d_r19 <-  svyby(~estbinres, by = ~imd_quintile,    design = dclus19g, FUN = svyciprop, vartype = "ci") %>% rename(level = imd_quintile)
prev_tab_nchild_r19 <-  svyby(~estbinres, by = ~nchild2,    design = dclus19g, FUN = svyciprop, vartype = "ci") %>% rename(level = nchild2)
prev_tab_urban_r19 <-  svyby(~estbinres, by = ~urban_char,    design = dclus19g, FUN = svyciprop, vartype = "ci") %>% rename(level = urban_char)
prev_tab_sheild2_r19 <-  svyby(~estbinres, by = ~shield2_char,    design = dclus19g, FUN = svyciprop, vartype = "ci") %>% rename(level = shield2_char)
prev_tab_indmask_r19 <-  svyby(~estbinres, by = ~indmask_char,    design = dclus19g, FUN = svyciprop, vartype = "ci") %>% rename(level = indmask_char)
prev_tab_vs_r19 <- svyby(~estbinres, by = ~vax_status_noDate_v2, design = dclus19g, FUN = svyciprop, vartype = "ci") %>% rename(level = vax_status_noDate_v2)

prev_tab_ar_r19 <- svyby(~estbinres, by = ~region + age_group_char3, design = dclus19g, FUN = svyciprop, vartype = "ci")

prev_tab_ar_r19 <- prev_tab_ar_r19 %>% 
  rename(age = age_group_char3,
         wt_prevalence = estbinres,
         wt_lower = ci_l,
         wt_upper = ci_u)

write.csv(prev_tab_ar_r19, "wt_prev_tab_ar_r19.csv", row.names = FALSE)

prev_tab_d_r19$level <- as.character(prev_tab_d_r19$level)

#' complete weighted table by covs
prev_tab_r19 <- bind_rows(
  prev_tab_g_r19,
  prev_tab_a3_r19,
  prev_tab_r_r19,
  prev_tab_w_r19,
  prev_tab_e_r19,
  prev_tab_hh_r19,
  prev_tab_nchild_r19,
  prev_tab_c_r19,
  prev_tab_h_r19,
  prev_tab_s_r19,
  prev_tab_d_r19,
  prev_tab_vs_r19,
  prev_tab_urban_r19,
  prev_tab_sheild2_r19,
  prev_tab_indmask_r19) %>%
  rename(prevalence_19 = estbinres,
         lower_19 = ci_l,
         upper_19 = ci_u)
rownames(prev_tab_r19) <- c(1:dim(prev_tab_r19)[1])

prev_tab_r19_ar <- prev_tab_r19 %>%
  filter(level %in% c("5-11", "12-17", "18-24", "25-34",
                      "35-44", "45-54", "55-64", "65-74", "75+",
                      "South East", "North East", "North West", "Yorkshire and The Humber",
                      "East Midlands", "West Midlands", "East of England", "London",
                      "South West"))


prev_tab_r19 <- bind_rows(prev_tab_r19[1:3,], #gender
                          # age
                          prev_tab_r19[4:12,,],
                          # region
                          prev_tab_r19[15,],prev_tab_r19[14,],prev_tab_r19[16,],prev_tab_r19[18,],
                          prev_tab_r19[17,],prev_tab_r19[19,],prev_tab_r19[21,],prev_tab_r19[20,],
                          prev_tab_r19[13,],
                          # urban
                          prev_tab_r19[c(69,67:68), ],
                          # work type
                          prev_tab_r19[23,], prev_tab_r19[24,], prev_tab_r19[22,],
                          prev_tab_r19[25,],prev_tab_r19[26,],
                          # ethnic
                          prev_tab_r19[27:32,],
                          # household size & # nchild & hospital contact
                          prev_tab_r19[33:51,],
                          # shield
                          prev_tab_r19[c(72,70:71),],
                          #mask
                          prev_tab_r19[c(73,77,74,76,75),],
                          # syptom
                          prev_tab_r19[52:55,],
                          # depreviation
                          prev_tab_r19[56:60,],
                          # vax
                          prev_tab_r19[66,], prev_tab_r19[61:65,])

# prev_tab_r19 <- bind_rows(prev_tab_r19[1:2,], #gender
#                           # age
#                           prev_tab_r19[3:11,],
#                           # region
#                           prev_tab_r19[12:20,],
#                           # work type
#                           prev_tab_r19[22,], prev_tab_r19[23,], prev_tab_r19[21,],
#                           prev_tab_r19[24,],prev_tab_r19[25,],
#                           # ethnic
#                           prev_tab_r19[26:31,],
#                           # household size & nchild & hospital contact
#                           prev_tab_r19[32:44,],
#                           # syptom
#                           prev_tab_r19[51:54,],
#                           # depreviation
#                           prev_tab_r19[55:59,],
#                           # vax
#                           prev_tab_r19[65,], prev_tab_r19[60:64,])

rownames(prev_tab_r19) <- c(1:dim(prev_tab_r19)[1])

prev_tab_r19_formatted <- prev_tab_r19 %>% 
  mutate(prevalence_19 = percent(prevalence_19, accuracy = 0.01), 
         lower_19 = percent(lower_19, accuracy = 0.01),
         upper_19 = percent(upper_19, accuracy = 0.01)) %>% 
  mutate(weighted_prevalence = paste0(prevalence_19, " (", lower_19, ", ", upper_19, ")"))

write.csv(prev_tab_r19_formatted, paste0("wt_prev_tab_by_covs_formatted_r19_",output_tag,".csv"),
          row.names = FALSE)

write.csv(prev_tab_r19, paste0("wt_prev_tab_by_covs_r19_",output_tag,".csv"),
          row.names = FALSE)



###### plot output

prev_tab_complete_r18 <- read.csv("E:/Group/report/round18/wt_prev_tab_by_covs_r18_06_Mar.csv")

prev_tab_r18_ar <- prev_tab_complete_r18  %>%
  filter(level %in% c("5-11", "12-17", "18-24", "25-34", 
                      "35-44", "45-54", "55-64", "65-74","75+",
                      "South East", "North East", "North West", "Yorkshire and The Humber", 
                      "East Midlands", "West Midlands", "East of England", "London",
                      "South West"))

prev_tab_r19_interim <- left_join(prev_tab_r18_ar, prev_tab_r19_ar, by = "level")


rtn <- list(prev_tab_r19_interim)

saveRDS(rtn, file = paste0("wt_prev_outputs_for_plots_",output_tag,".rds"))

