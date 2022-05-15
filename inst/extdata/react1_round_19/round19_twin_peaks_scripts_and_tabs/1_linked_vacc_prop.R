#' html generated with the knitr package by spin(<this file>)
#' Author: Haowei Wang
#' Unwieghted prevalence for round 19
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
          "scales" ,"reshape2")
load_packages(pkgs)

#' Source any functions from the local file
source("E:/Group/functions/add_conf_ints.R")
source("E:/Group/functions/make_tables.R")
source("E:/Group/functions/overall_prev.R")
# source("E:/Group/functions/model_table.R")
# source("E:/Group/functions/df_make.R")
# source("E:/Group/functions/get_gam_deets.R")
# source("E:/Group/functions/get_or_tables.R")


#' ## Load data and add variables to working dataframe

#' Set output directory
setwd("E:/Group/report/round19")

#' ## Read in the data and organise a little 
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

dfRes8 <- dfRes8 %>% 
  mutate(vacc_res = ifelse(vax_status_noDate_v2 %in% c("One does", "Two does", "Three does"), 1,
                           ifelse(vax_status_noDate_v2 %in% c("Not vaccinated", "Unknown does"), 0, NA))) %>% 
  filter(!is.na(vacc_res))%>% filter(!is.na(estbinres))

dfRes9 <- dfRes9 %>% 
  mutate(vacc_res = ifelse(vax_status_noDate_v2 %in% c("One does", "Two does", "Three does"), 1,
                           ifelse(vax_status_noDate_v2 %in% c("Not vaccinated", "Unknown does"), 0, NA))) %>% 
  filter(!is.na(vacc_res))%>% filter(!is.na(estbinres))

dfRes10 <- dfRes10 %>% 
  mutate(vacc_res = ifelse(vax_status_noDate_v2 %in% c("One does", "Two does", "Three does"), 1,
                           ifelse(vax_status_noDate_v2 %in% c("Not vaccinated", "Unknown does"), 0, NA))) %>% 
  filter(!is.na(vacc_res))%>% filter(!is.na(estbinres))

dfRes11 <- dfRes11 %>% 
  mutate(vacc_res = ifelse(vax_status_noDate_v2 %in% c("One does", "Two does", "Three does"), 1,
                           ifelse(vax_status_noDate_v2 %in% c("Not vaccinated", "Unknown does"), 0, NA))) %>% 
  filter(!is.na(vacc_res))%>% filter(!is.na(estbinres))

dfRes12 <- dfRes12 %>% 
  mutate(vacc_res = ifelse(vax_status_noDate_v2 %in% c("One does", "Two does", "Three does"), 1,
                           ifelse(vax_status_noDate_v2 %in% c("Not vaccinated", "Unknown does"), 0, NA))) %>% 
  filter(!is.na(vacc_res))%>% filter(!is.na(estbinres))

dfRes13 <- dfRes13 %>% 
  mutate(vacc_res = ifelse(vax_status_noDate_v2 %in% c("One does", "Two does", "Three does"), 1,
                           ifelse(vax_status_noDate_v2 %in% c("Not vaccinated", "Unknown does"), 0, NA))) %>% 
  filter(!is.na(vacc_res))%>% filter(!is.na(estbinres))

dfRes14 <- dfRes14 %>% 
  mutate(vacc_res = ifelse(vax_status_noDate_v2 %in% c("One does", "Two does", "Three does"), 1,
                           ifelse(vax_status_noDate_v2 %in% c("Not vaccinated", "Unknown does"), 0, NA))) %>% 
  filter(!is.na(vacc_res))%>% filter(!is.na(estbinres))

dfRes15 <- dfRes15 %>% 
  mutate(vacc_res = ifelse(vax_status_noDate_v2 %in% c("One does", "Two does","Three does"), 1,
                           ifelse(vax_status_noDate_v2 %in% c("Not vaccinated", "Unknown does"), 0, NA))) %>% 
  filter(!is.na(vacc_res))%>% filter(!is.na(estbinres))

dfRes16 <- dfRes16 %>% 
  mutate(vacc_res = ifelse(vax_status_noDate_v2 %in% c("One does", "Two does", "Three does"), 1,
                           ifelse(vax_status_noDate_v2 %in% c("Not vaccinated", "Unknown does"), 0, NA))) %>% 
  filter(!is.na(vacc_res))%>% filter(!is.na(estbinres))

dfRes17 <- dfRes17 %>% 
  mutate(vacc_res = ifelse(vax_status_noDate_v2 %in% c("One does", "Two does", "Three does"), 1,
                           ifelse(vax_status_noDate_v2 %in% c("Not vaccinated", "Unknown does"), 0, NA))) %>% 
  filter(!is.na(vacc_res))%>% filter(!is.na(estbinres))

dfRes18 <- dfRes18 %>% 
  mutate(vacc_res = ifelse(vax_status_noDate_v2 %in% c("One does", "Two does", "Three does"), 1,
                           ifelse(vax_status_noDate_v2 %in% c("Not vaccinated", "Unknown does"), 0, NA))) %>% 
  filter(!is.na(vacc_res))%>% filter(!is.na(estbinres))

dfRes19 <- dfRes19 %>% 
  mutate(vacc_res = ifelse(vax_status_noDate_v2 %in% c("One does", "Two does", "Three does"), 1,
                           ifelse(vax_status_noDate_v2 %in% c("Not vaccinated", "Unknown does"), 0, NA))) %>% 
  filter(!is.na(vacc_res))%>% filter(!is.na(estbinres))

####### Proportion vaccination ############

covs <- c("age_group_char3")

#' Make the prevalence tables for the above covariates (unweighted)
perc = FALSE
sig_figs = 6

prop_tables_r8 <- make_tables(dat = dfRes8, covariates = covs, sens = 1, spec = 1, method = "exact",
                              result_var = "vacc_res", suffix = "r8", sf = sig_figs, percent = perc) %>%
  bind_rows(.) %>%
  rename("Unweighted_proportion" = "Prevalence") %>% 
  mutate(Round = "Round 8") %>% 
  relocate(Round, .before = 1) 

prop_tables_r9 <- make_tables(dat = dfRes9, covariates = covs, sens = 1, spec = 1, method = "exact",
                              result_var = "vacc_res", suffix = "r9", sf = sig_figs, percent = perc) %>%
  bind_rows(.) %>%
  rename("Unweighted_proportion" = "Prevalence") %>% 
  mutate(Round = "Round 9") %>% 
  relocate(Round, .before = 1)

prop_tables_r10 <- make_tables(dat = dfRes10, covariates = covs, sens = 1, spec = 1, method = "exact",
                               result_var = "vacc_res", suffix = "r10", sf = sig_figs, percent = perc) %>%
  bind_rows(.) %>%
  rename("Unweighted_proportion" = "Prevalence") %>% 
  mutate(Round = "Round 10") %>% 
  relocate(Round, .before = 1)

prop_tables_r11 <- make_tables(dat = dfRes11, covariates = covs, sens = 1, spec = 1, method = "exact",
                               result_var = "vacc_res", suffix = "r11", sf = sig_figs, percent = perc) %>%
  bind_rows(.) %>%
  rename("Unweighted_proportion" = "Prevalence") %>% 
  mutate(Round = "Round 11") %>% 
  relocate(Round, .before = 1)

prop_tables_r12 <- make_tables(dat = dfRes12, covariates = covs, sens = 1, spec = 1, method = "exact",
                               result_var = "vacc_res", suffix = "r12", sf = sig_figs, percent = perc) %>%
  bind_rows(.) %>%
  rename("Unweighted_proportion" = "Prevalence") %>% 
  mutate(Round = "Round 12") %>% 
  relocate(Round, .before = 1)

prop_tables_r13 <- make_tables(dat = dfRes13, covariates = covs, sens = 1, spec = 1, method = "exact",
                               result_var = "vacc_res", suffix = "r13", sf = sig_figs, percent = perc) %>%
  bind_rows(.) %>%
  rename("Unweighted_proportion" = "Prevalence") %>% 
  mutate(Round = "Round 13") %>% 
  relocate(Round, .before = 1)

prop_tables_r14 <- make_tables(dat = dfRes14, covariates = covs, sens = 1, spec = 1, method = "exact",
                               result_var = "vacc_res", suffix = "r14", sf = sig_figs, percent = perc) %>%
  bind_rows(.) %>%
  rename("Unweighted_proportion" = "Prevalence") %>% 
  mutate(Round = "Round 14") %>% 
  relocate(Round, .before = 1)

prop_tables_r15 <- make_tables(dat = dfRes15, covariates = covs, sens = 1, spec = 1, method = "exact",
                               result_var = "vacc_res", suffix = "r15", sf = sig_figs, percent = perc) %>%
  bind_rows(.) %>%
  rename("Unweighted_proportion" = "Prevalence") %>% 
  mutate(Round = "Round 15") %>% 
  relocate(Round, .before = 1)

prop_tables_r16 <- make_tables(dat = dfRes16, covariates = covs, sens = 1, spec = 1, method = "exact",
                               result_var = "vacc_res", suffix = "r16", sf = sig_figs, percent = perc) %>%
  bind_rows(.) %>%
  rename("Unweighted_proportion" = "Prevalence") %>% 
  mutate(Round = "Round 16") %>% 
  relocate(Round, .before = 1)

prop_tables_r17 <- make_tables(dat = dfRes17, covariates = covs, sens = 1, spec = 1, method = "exact",
                               result_var = "vacc_res", suffix = "r17", sf = sig_figs, percent = perc) %>%
  bind_rows(.) %>%
  rename("Unweighted_proportion" = "Prevalence") %>% 
  mutate(Round = "Round 17") %>% 
  relocate(Round, .before = 1)

prop_tables_r18 <- make_tables(dat = dfRes18, covariates = covs, sens = 1, spec = 1, method = "exact",
                               result_var = "vacc_res", suffix = "r18", sf = sig_figs, percent = perc) %>%
  bind_rows(.) %>%
  rename("Unweighted_proportion" = "Prevalence") %>% 
  mutate(Round = "Round 18") %>% 
  relocate(Round, .before = 1)

prop_tables_r19 <- make_tables(dat = dfRes19, covariates = covs, sens = 1, spec = 1, method = "exact",
                               result_var = "vacc_res", suffix = "r19", sf = sig_figs, percent = perc) %>%
  bind_rows(.) %>%
  rename("Unweighted_proportion" = "Prevalence") %>% 
  mutate(Round = "Round 19") %>% 
  relocate(Round, .before = 1)

ut_prop_tables <- rbind(prop_tables_r8, prop_tables_r9, prop_tables_r10,
                        prop_tables_r11, prop_tables_r12, prop_tables_r13,
                        prop_tables_r14, prop_tables_r15, prop_tables_r16,
                        prop_tables_r17, prop_tables_r18, prop_tables_r19) %>% 
  mutate(unweighted_prop_full = paste0(percent(Unweighted_proportion, accuracy = 0.01), " (", 
                                       percent(Lower, accuracy = 0.01), ", ", 
                                       percent(Upper, accuracy = 0.01), ")"))

########## weighted proportion vaccination ##########
df_round8_v <- dfRes8 %>% filter(!is.na(estbinres)) %>% mutate(group = "Overall") %>% 
  filter(!is.na(wt_antigen)) %>% 
  filter(!is.na(vacc_res))

df_round9_v <- dfRes9 %>% filter(!is.na(estbinres)) %>% mutate(group = "Overall") %>% 
  filter(!is.na(wt_antigen)) %>% 
  filter(!is.na(vacc_res))

df_round10_v <- dfRes10 %>% filter(!is.na(estbinres)) %>% mutate(group = "Overall") %>% 
  filter(!is.na(wt_antigen)) %>% 
  filter(!is.na(vacc_res))

df_round11_v <- dfRes11 %>% filter(!is.na(estbinres)) %>% mutate(group = "Overall") %>% 
  filter(!is.na(wt_antigen)) %>% 
  filter(!is.na(vacc_res))

df_round12_v <- dfRes12 %>% filter(!is.na(estbinres)) %>% mutate(group = "Overall") %>% 
  filter(!is.na(wt_antigen)) %>% 
  filter(!is.na(vacc_res))

df_round13_v <- dfRes13 %>% filter(!is.na(estbinres)) %>% mutate(group = "Overall") %>% 
  filter(!is.na(wt_antigen)) %>% 
  filter(!is.na(vacc_res))

df_round14_v <- dfRes14 %>% filter(!is.na(estbinres)) %>% mutate(group = "Overall") %>% 
  filter(!is.na(wt_antigen)) %>% 
  filter(!is.na(vacc_res))

df_round15_v <- dfRes15 %>% filter(!is.na(estbinres)) %>% mutate(group = "Overall") %>% 
  filter(!is.na(wt_antigen)) %>% 
  filter(!is.na(vacc_res))

df_round16_v <- dfRes16 %>% filter(!is.na(estbinres)) %>% mutate(group = "Overall") %>% 
  filter(!is.na(wt_antigen)) %>% 
  filter(!is.na(vacc_res))

df_round17_v <- dfRes17 %>% filter(!is.na(estbinres)) %>% mutate(group = "Overall") %>% 
  filter(!is.na(wt_antigen)) %>% 
  filter(!is.na(vacc_res))

df_round18_v <- dfRes18 %>% filter(!is.na(estbinres)) %>% mutate(group = "Overall") %>% 
  filter(!is.na(wt_antigen)) %>% 
  filter(!is.na(vacc_res))

df_round19_v <- dfRes19 %>% filter(!is.na(estbinres)) %>% mutate(group = "Overall") %>% 
  filter(!is.na(wt_antigen)) %>% 
  filter(!is.na(vacc_res))

dclus8v <- svydesign(id = ~id, strata = ~ lacode, weights = ~ wt_antigen, data = df_round8_v, nest = TRUE)
dclus9v <- svydesign(id = ~id, strata = ~ lacode, weights = ~ wt_antigen, data = df_round9_v, nest = TRUE)
dclus10v <- svydesign(id = ~id, strata = ~ lacode, weights = ~ wt_antigen, data = df_round10_v, nest = TRUE)
dclus11v <- svydesign(id = ~id, strata = ~ lacode, weights = ~ wt_antigen, data = df_round11_v, nest = TRUE)
dclus12v <- svydesign(id = ~id, strata = ~ lacode, weights = ~ wt_antigen, data = df_round12_v, nest = TRUE)
dclus13v <- svydesign(id = ~id, strata = ~ lacode, weights = ~ wt_antigen, data = df_round13_v, nest = TRUE)
dclus14v <- svydesign(id = ~id, strata = ~ lacode, weights = ~ wt_antigen, data = df_round14_v, nest = TRUE)
dclus15v <- svydesign(id = ~id, strata = ~ lacode, weights = ~ wt_antigen, data = df_round15_v, nest = TRUE)
dclus16v <- svydesign(id = ~id, strata = ~ lacode, weights = ~ wt_antigen, data = df_round16_v, nest = TRUE)
dclus17v <- svydesign(id = ~id, strata = ~ lacode, weights = ~ wt_antigen, data = df_round17_v, nest = TRUE)
dclus18v <- svydesign(id = ~id, strata = ~ lacode, weights = ~ wt_antigen, data = df_round18_v, nest = TRUE)
dclus19v <- svydesign(id = ~id, strata = ~ lacode, weights = ~ wt_antigen, data = df_round19_v, nest = TRUE)


prev_tab_vacc_a3_r8 <-  svyby(~vacc_res, by = ~age_group_char3,  design = dclus8v, FUN = svyciprop, vartype = "ci") %>% rename(level = age_group_char3)
prev_tab_vacc_a3_r9 <-  svyby(~vacc_res, by = ~age_group_char3,  design = dclus9v, FUN = svyciprop, vartype = "ci") %>% rename(level = age_group_char3)
prev_tab_vacc_a3_r10 <-  svyby(~vacc_res, by = ~age_group_char3,  design = dclus10v, FUN = svyciprop, vartype = "ci") %>% rename(level = age_group_char3)
prev_tab_vacc_a3_r11 <-  svyby(~vacc_res, by = ~age_group_char3,  design = dclus11v, FUN = svyciprop, vartype = "ci") %>% rename(level = age_group_char3)
prev_tab_vacc_a3_r12 <-  svyby(~vacc_res, by = ~age_group_char3,  design = dclus12v, FUN = svyciprop, vartype = "ci") %>% rename(level = age_group_char3)
prev_tab_vacc_a3_r13 <-  svyby(~vacc_res, by = ~age_group_char3,  design = dclus13v, FUN = svyciprop, vartype = "ci") %>% rename(level = age_group_char3)
prev_tab_vacc_a3_r14 <-  svyby(~vacc_res, by = ~age_group_char3,  design = dclus14v, FUN = svyciprop, vartype = "ci") %>% rename(level = age_group_char3)
prev_tab_vacc_a3_r15 <-  svyby(~vacc_res, by = ~age_group_char3,  design = dclus15v, FUN = svyciprop, vartype = "ci") %>% rename(level = age_group_char3)
prev_tab_vacc_a3_r16 <-  svyby(~vacc_res, by = ~age_group_char3,  design = dclus16v, FUN = svyciprop, vartype = "ci") %>% rename(level = age_group_char3)
prev_tab_vacc_a3_r17 <-  svyby(~vacc_res, by = ~age_group_char3,  design = dclus17v, FUN = svyciprop, vartype = "ci") %>% rename(level = age_group_char3)
prev_tab_vacc_a3_r18 <-  svyby(~vacc_res, by = ~age_group_char3,  design = dclus18v, FUN = svyciprop, vartype = "ci") %>% rename(level = age_group_char3)
prev_tab_vacc_a3_r19 <-  svyby(~vacc_res, by = ~age_group_char3,  design = dclus19v, FUN = svyciprop, vartype = "ci") %>% rename(level = age_group_char3)


prev_tab_vacc_a3_r8 <- rbind(c("5-11", 0, 0, 0),prev_tab_vacc_a3_r8)
prev_tab_vacc_a3_r13 <- rbind(c("5-11", 0, 0, 0),prev_tab_vacc_a3_r13)
prev_tab_vacc_a3_r14 <- rbind(c("5-11", 0, 0, 0),prev_tab_vacc_a3_r14)
prev_tab_vacc_a3_r15 <- rbind(c("5-11", 0, 0, 0),prev_tab_vacc_a3_r15)
prev_tab_vacc_a3_r16 <- rbind(c("5-11", 0, 0, 0),prev_tab_vacc_a3_r16)
prev_tab_vacc_a3_r17 <- rbind(c("5-11", 0, 0, 0),prev_tab_vacc_a3_r17)
prev_tab_vacc_a3_r18 <- rbind(c("5-11", 0, 0, 0),prev_tab_vacc_a3_r18)
prev_tab_vacc_a3_r19 <- rbind(c("5-11", 0, 0, 0),prev_tab_vacc_a3_r19)

prev_tab_vacc_a3_r8 <- prev_tab_vacc_a3_r8 %>% 
  mutate(Round = "Round 8") %>% 
  relocate(Round, .before = 1) %>% 
  rename(Category = level)

prev_tab_vacc_a3_r9 <- prev_tab_vacc_a3_r9 %>% 
  mutate(Round = "Round 9") %>% 
  relocate(Round, .before = 1) %>% 
  rename(Category = level)

prev_tab_vacc_a3_r10 <- prev_tab_vacc_a3_r10 %>% 
  mutate(Round = "Round 10") %>% 
  relocate(Round, .before = 1) %>% 
  rename(Category = level)

prev_tab_vacc_a3_r11 <- prev_tab_vacc_a3_r11 %>% 
  mutate(Round = "Round 11") %>% 
  relocate(Round, .before = 1) %>% 
  rename(Category = level)

prev_tab_vacc_a3_r12 <- prev_tab_vacc_a3_r12 %>% 
  mutate(Round = "Round 12") %>% 
  relocate(Round, .before = 1) %>% 
  rename(Category = level)

prev_tab_vacc_a3_r13 <- prev_tab_vacc_a3_r13 %>% 
  mutate(Round = "Round 13") %>% 
  relocate(Round, .before = 1) %>% 
  rename(Category = level)

prev_tab_vacc_a3_r14 <- prev_tab_vacc_a3_r14 %>% 
  mutate(Round = "Round 14") %>% 
  relocate(Round, .before = 1) %>% 
  rename(Category = level)

prev_tab_vacc_a3_r15 <- prev_tab_vacc_a3_r15 %>% 
  mutate(Round = "Round 15") %>% 
  relocate(Round, .before = 1) %>% 
  rename(Category = level)

prev_tab_vacc_a3_r16 <- prev_tab_vacc_a3_r16 %>% 
  mutate(Round = "Round 16") %>% 
  relocate(Round, .before = 1) %>% 
  rename(Category = level)

prev_tab_vacc_a3_r17 <- prev_tab_vacc_a3_r17 %>% 
  mutate(Round = "Round 17") %>% 
  relocate(Round, .before = 1) %>% 
  rename(Category = level)

prev_tab_vacc_a3_r18 <- prev_tab_vacc_a3_r18 %>% 
  mutate(Round = "Round 18") %>% 
  relocate(Round, .before = 1) %>% 
  rename(Category = level)

prev_tab_vacc_a3_r19 <- prev_tab_vacc_a3_r19 %>% 
  mutate(Round = "Round 19") %>% 
  relocate(Round, .before = 1) %>% 
  rename(Category = level)


wt_prop_tables <- rbind(prev_tab_vacc_a3_r8, prev_tab_vacc_a3_r9, prev_tab_vacc_a3_r10,
                        prev_tab_vacc_a3_r11, prev_tab_vacc_a3_r12, prev_tab_vacc_a3_r13,
                        prev_tab_vacc_a3_r14, prev_tab_vacc_a3_r15, prev_tab_vacc_a3_r16,
                        prev_tab_vacc_a3_r17, prev_tab_vacc_a3_r18, prev_tab_vacc_a3_r19) %>% 
  mutate(vacc_res = as.numeric(vacc_res),
         ci_l = as.numeric(ci_l),
         ci_u = as.numeric(ci_u)) %>% 
  mutate(weighted_prop = paste0(percent(vacc_res, accuracy = 0.01), " (", 
                                percent(ci_l, accuracy = 0.01), ", ", 
                                percent(ci_u, accuracy = 0.01), ")"))

prop_table_all <- full_join(ut_prop_tables , wt_prop_tables, by = c("Round", "Category"))

prop_table_all_under18 <- prop_table_all %>% 
  filter(Category %in% c("5-11", "12-17"))

write.csv(prop_table_all, "./linked_vacc_prop_r8_r19.csv")
write.csv(prop_table_all_under18, "./linked_vacc_prop_under18_r8_r19.csv")

