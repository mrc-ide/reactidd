
#' First clear the environment of variables
rm(list=ls(all=TRUE))
# get root director of project
root.dir <- getwd()
setwd(dir = "/Group/react2_study5/report_phases_combined/projects/omicron_symptom_profiling/")
outpath <- paste0(root.dir,"/output/")
figpath <-  paste0(root.dir,"/plots/")



source("E:/Group/functions/load_packages.R", local = T)
source("E:/Group/functions/full_join_multiple_datasets.R", local = T)
source("E:/Group/functions/wrangle_cleaner_functions.R", local = T)
source("E:/Group/functions/cats_and_covs.R", local = T)
source("E:/Group/react2_study5/report_phases_combined/projects/function_scripts/create_subfolder.R", local = T)
source("E:/Group/react2_study5/report_phases_combined/projects/function_scripts/forest_plot.R", local = T)
source("E:/Group/react2_study5/report_phases_combined/projects/function_scripts/save_styled_table.R", local = T)
source("E:/Group/react2_study5/report_phases_combined/projects/function_scripts/stability_selection.R", local = T)

# 
# 
# source("E:/Group/react2_study5/report_phases_combined/projects/delta_symptom_prediction/code/00_bits_and_pieces.R", local = T)
# source("E:/Group/react2_study5/report_phases_combined/projects/delta_symptom_prediction/code/00_functions.R", local = T)
# source("E:/Group/react2_study5/report_phases_combined/projects/symptom_prediction_children/code/00_bits_and_pieces.R", local = T)


#' Pull in packages needed
package.list <- c("prevalence","mgcv","knitr","MASS","kableExtra","table1","dplyr",
                  "tidyr", "pheatmap","scales","OverReact","ggstance",
                  "ggplot2","ggsci", "RColorBrewer", "tidyverse", "lubridate", 
                  "readr","ggthemes", "questionr", "foreach", "doParallel","withr",
                  "patchwork","randomcoloR","focus"
)

load_packages(package.list)


# Import REACT-1 data -----------------------------------------------------
source("E:/Group/react2_study5/report_phases_combined/projects/omicron_symptom_profiling/code/00_bits_and_pieces.R", local = T)
source("E:/Group/react2_study5/report_phases_combined/projects/omicron_symptom_profiling//code/00_data_prep.R", local = T)

# create subfolder
createMySubfolder(subfolderName = "waning_analysis")


# Date edit
dfRes <- dfRes %>% mutate(vaccinethird=as.Date(vaccinethird, format = "%m/%d/%Y"),
                          days_since_boost = as.numeric(d_comb-vaccinethird),
                          impact_a_lot=case_when(covidability7==1 ~ 1,
                                                  T ~ 0),
                          impact_a_lot_or_a_litte=case_when(covidability7==1 ~ 1,
                                                            covidability7==2 ~ 1,
                                                 T ~ 0),
                          BA2 = case_when(variant_inferred_detail=="BA.2 (Omicron)" ~ 1,
                                          variant_inferred_detail=="BA.1 (Omicron)" ~ 0,
                                          T ~ NA_real_)
                          )

# create symptom count variable
dfRes$symptom_count_26 <- rowSums(dfRes[,covid_yesnos], na.rm=T)
dfRes$symptom_count_4 <- rowSums(dfRes[,sympnames_type_df$symptom_code[1:4]], na.rm=T)


table(dfRes$covidability7)
# subset data
dfRes_pos=dfRes %>% filter(variant_inferred_detail%in% c("BA.1 (Omicron)","BA.2 (Omicron)"), round %in%c(17:19))
dfRes_pos_boost=dfRes %>% filter(variant_inferred_detail%in% c("BA.1 (Omicron)","BA.2 (Omicron)"), 
                                 round %in%c(17:19), !covidability7%in% c(4:5),
                                 vaccdose==3,days_since_boost>=14)


# Waning analysis for infection in r17-19, BA1/BA2 ----------------------
unique(dfRes_pos_boost$age_group_named)
dfRes_pos_boost$age_group_named <- factor(dfRes_pos_boost$age_group_named,
                                          levels = c("18-24" ,"25-34","35-44","45-54","55-64","65-74", "74+"))

mod_glm_severe = glm(formula = as.formula("impact_a_lot ~ age_group_named + sex + BA2 + days_since_boost +I(days_since_boost^2)+ round"),
              data = dfRes_pos_boost,family = "binomial")
mod_glm_a_lot_or_a_little = glm(formula = 
                                  as.formula("impact_a_lot_or_a_litte ~ age_group_named + sex + BA2 + days_since_boost+
                                             I(days_since_boost^2)+ round"),
                     data = dfRes_pos_boost,family = "binomial")


jtools::summ(mod_glm_severe,exp=T)
jtools::summ(mod_glm_a_lot_or_a_little,exp=T)

# output severe table
tab_severity=jtools::summ(mod_glm_severe,exp=T)$coeftable %>% as.data.frame() %>% 
  select(-`z val.`)
tab_severity$OR_with_CI=paste0(round(tab_severity$`exp(Est.)`,2), ", [", round(tab_severity$`2.5%`,2),",",
                               round(tab_severity$`97.5%`,2),"]")
tab_severity$variable = rownames(tab_severity)
tab_severity <- tab_severity %>% select(variable,OR_with_CI,p)
colnames(tab_severity)= c("Independent variable","Odds ratio with CI","P-value")
# tab_severity <- tab_severity %>% select(variable, everything())
tab_severity$`Independent variable` <- c("Intercept","Age group: 25-34",
                                         "Age group: 35-44",
                                         "Age group: 45-54",
                                         "Age group: 55-64",
                                         "Age group: 65-74",
                                         "Age group: 74+",
                                         "Sex: male",
                                         "BA.2",
                                         "Days since booster vaccine",
                                         "[Days since booster vaccine]^2",
                                         "Round 18",
                                         "Round 19"
                                         )
# output severe table
tab_severity_2=jtools::summ(mod_glm_a_lot_or_a_little,exp=T)$coeftable %>% as.data.frame() %>% 
  select(-`z val.`)
colnames(tab_severity_2)= c("Odds ratio","Lower (95% CI)","Upper (95% CI)","P-value")
tab_severity_2$variable = rownames(tab_severity_2)
tab_severity_2 <- tab_severity_2 %>% select(variable, everything())

# Loop over all symptoms --------------------------------------------------

sympFunc <- function(symp){
  f=as.formula(paste0(symp," ~ age_group_named + sex + BA2 + days_since_boost + round"))
  mod = glm(formula = f,data = dfRes_pos_boost,family = "binomial")
  tab=jtools::summ(mod,exp=T)$coeftable %>% as.data.frame()
  tab$symptom = sympnames_type_df$symptom[sympnames_type_df$symptom_code==symp]
  return(tab[4,])
}

# run over all symptoms
allmods=lapply(covid_yesnos,sympFunc)
names(allmods)=sympnames_type_df$symptom[1:26]
allmods=bind_rows(allmods) %>% select(-`z val.`)
colnames(allmods)= c("Odds ratio","Lower (95% CI)","Upper (95% CI)","P-value","Symptom")
allmods <- allmods %>% select(Symptom, everything())


savePrettyExcelWorkbook(listOfTables = list(tab_a_lot=tab_severity, tab_a_little_or_a_lot=tab_severity_2,
                                            tab_all_symptoms=allmods),
                        workbookName = "severity_modelling",outpath = outpath)
