#' First clear the environment of variables
rm(list=ls(all=TRUE))
setwd("E:/Group/react2_study5/report_phases_combined/projects/react_1_antigen_symptom_prediction/")
outpath <- "E:/Group/react2_study5/report_phases_combined/projects/react_1_antigen_symptom_prediction/output/"
figpath <- "E:/Group/react2_study5/report_phases_combined/projects/react_1_antigen_symptom_prediction/plots/"

#' Source any functions from the local file
source("E:/Group/functions/load_packages.R")
source("E:/Group/functions/hamming_distance.R")
source("E:/Group/functions/cats_and_covs.R")
source("E:/Group/functions/crosstab.R")
source("E:/Group/functions/relevel_all_categorical_vars.R")
source("E:/Group/functions/wrangle_cleaner_functions.R", local = T)
source("E:/Group/functions/full_join_multiple_datasets.R", local = T)
source("E:/Group/react2_study5/report_phases_combined/projects/clustering/code/00_functions.R")
source("E:/Group/react2_study5/report_phases_combined/projects/functions/00_xgboost_funcs.R")

source("E:/Group/react2_study5/report_phases_combined/projects/functions/00_longcovid_calculator.R")
source("E:/Group/react2_study5/report_phases_combined/projects/functions/00_LCA_optimiser.R")
source("E:/Group/react2_study5/report_phases_combined/projects/functions/00_lasso_stability_funcs.R")
source("E:/Group/react2_study5/report_phases_combined/projects/functions/00_cif_funcs.R")
source("E:/Group/react2_study5/report_phases_combined/projects/functions/univariate_models.R")

#' Pull in packages needed
package.list <- c("prevalence","mgcv","knitr","MASS","kableExtra","table1","dplyr","factoextra","tableone","networkD3",
                  "tidyr","forcats", "cluster", "fpc", "mclust", "pheatmap","FactoMineR", "NbClust","clValid","plotly",
                  "ggplot2","ggsci", "RColorBrewer", "tidyverse", "lubridate", "egg", "poLCA", "Rcpp","xml2","splitstackshape",
                  "fs", "later", "promises","proxy","dendextend", "ComplexHeatmap","circlize","doSNOW","ClusterR","htmlwidgets",
                  "readr","ggthemes", "questionr", "gridExtra", "foreach", "doParallel", "patchwork",
                  "purrr", "httr", "htmltools","ggalluvial","datapasta","xgboost","SHAPforxgboost"
)

load_packages(package.list)



############################ Run data import script ########################
r8_include_switch = F
source("E:/Group/react2_study5/report_phases_combined/projects/react_1_antigen_symptom_prediction/code/00_data_prep.R")
############################ Run data import script ########################

# Univariate analysis -----------------------------------------------------
mod.dat <- dfRes%>% filter(!is.na(estbinres), !is.na(age), !is.na(gender)) %>% 
  select(estbinres,all_of(covid_yesnos), age, gender, round, interim_label, age_group_pred)

results_all <- univariateModels(mod.dat, adjustments = c("age", "gender"), modelname = "all")
results_u18 <- univariateModels(mod.dat %>% filter(age_group_pred == 1), adjustments = c("age", "gender"), 
                                modelname = "5-17")
results_18_55 <- univariateModels(mod.dat %>% filter(age_group_pred == 2), adjustments = c("age", "gender"),
                                  modelname = "18_54")
results_55_plus <- univariateModels(mod.dat %>% filter(age_group_pred == 3), adjustments = c("age", "gender"),
                                    modelname = "55_plus")

results_u18$crude$modname <- "5-17"
results_18_55$crude$modname <- "18-54"
results_55_plus$crude$modname <- "55 and over"
results_comb <- rbind(results_all$crude,results_u18$crude, results_18_55$crude, results_55_plus$crude)

plotUnivResults <- function(results_comb){
  dodge=position_dodge(width=0.6)
  #### unadjusted results
  p1 <- results_comb %>% 
    filter(modname!="all") %>% 
    filter(var%in%c("Loss or change of sense of smell","Loss or change of sense of taste","Fever",
                    "New persistent cough")) %>%
    mutate(modname=factor(modname, levels = unique(modname)),
           var=factor(var, levels = unique(var)),
           classic=ifelse(var%in%c("Loss or change of sense of smell","Loss or change of sense of taste","Fever",
                                   "New persistent cough"), "Classic symptoms", "Other symptoms")) %>% 
    ggplot(aes(x=reorder(var, coef), y= coef, col =(modname))) +
    geom_point(position = dodge) +
    geom_errorbar(aes(ymin = lower, ymax = upper),position = dodge, 
                  size=0.3, width=0.5) +
    coord_flip() +
    facet_wrap(.~classic, scales = "free", ncol = 1) +
    theme_bw() +
    geom_hline(yintercept = 1, linetype = "dashed")+
    # ylim(c(0,100))+
    scale_color_manual(values = myCols[c(1:3)])+
    labs(x="Symptom", y = "Odds ratio", col="Age \rstratum") +
    scale_y_log10(breaks = c(1,10,100),limits = c(0.3,100)) +
    guides(col = guide_legend(reverse=T)) +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          strip.background = element_rect(fill = "white"),
          strip.text = element_text(face = "bold"),
          axis.text.y = element_text(size = 10))
  p1
  
  #### unadjusted results
  p2 <- results_comb %>% 
    filter(modname!="all") %>% 
    filter(!var%in%c("Loss or change of sense of smell","Loss or change of sense of taste","Fever",
                     "New persistent cough")) %>%
    mutate(modname=factor(modname, levels = unique(modname)),
           var=factor(var, levels = unique(var)),
           classic=ifelse(var%in%c("Loss or change of sense of smell","Loss or change of sense of taste","Fever",
                                   "New persistent cough"), "Classic symptoms", "Other symptoms")) %>% 
    ggplot(aes(x=reorder(var, coef), y= coef, col =(modname))) +
    geom_point(position = dodge) +
    geom_errorbar(aes(ymin = lower, ymax = upper),position = dodge, 
                  size=0.3, width=0.5) +
    coord_flip() +
    facet_wrap(.~classic, scales = "free", ncol = 1) +
    theme_bw() +
    geom_hline(yintercept = 1, linetype = "dashed")+
    # ylim(c(0,65))+
    scale_color_manual(values = myCols[c(1:3)])+
    labs(x="Symptom", y = "Odds ratio", col="Age \rstratum") +
    scale_y_log10(breaks = c(1,10,100),limits = c(0.3,100))  +
    guides(col = guide_legend(reverse=T)) +
    theme(
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(face = "bold"),
      axis.text.y = element_text(size = 10)
        )
  
  p2
  
  pcomb <- p1/p2 + plot_layout(guides = "collect", heights = c(4,22))
  return(pcomb)
  
}

pcomb <- plotUnivResults(results_comb)
pcomb

### Save plot
ggsave(filename = paste0(figpath, "uni_ORs_crude_include_asymptomatics.png"), plot = pcomb, 
       width = 7, height=8, 
       dpi = 300, units = "in")
### Save plot
ggsave(filename = paste0(figpath, "uni_ORs_crude_include_asymptomatics.pdf"), plot = pcomb, 
       width = 7, height=8, 
       dpi = 300, units = "in")


### Save univariate results
write_csv(results_comb, paste0(outpath,"univariate_results_crude_include_asymptomatics.csv"))

### Save univariate results
saveRDS(results_comb, paste0(outpath,"univariate_results_crude_include_asymptomatics.rds"))


# As above, with adjusted results -----------------------------------------

results_comb_2 <- rbind(results_all$adj,results_u18$adj, results_18_55$adj, results_55_plus$adj)

pcomb_adj <- plotUnivResults(results_comb_2)
pcomb_adj






# Excluding asymptomatics --------------------------------------------------


mydata <- dfRes %>% filter(!is.na(estbinres), !is.na(age), !is.na(gender), symptomatic == 1) %>% 
  select(estbinres,all_of(covid_yesnos), age, gender, round, interim_label, age_group_pred)

results_all <- univariateModels(mydata, adjustments = c("age", "gender"), modelname = "all")
results_u18 <- univariateModels(mydata %>% filter(age_group_pred == 1), adjustments = c("age", "gender"), 
                                modelname = "5-17")
results_18_55 <- univariateModels(mydata %>% filter(age_group_pred == 2), adjustments = c("age", "gender"),
                                  modelname = "18_54")
results_55_plus <- univariateModels(mydata %>% filter(age_group_pred == 3), adjustments = c("age", "gender"),
                                    modelname = "55_plus")

results_u18$crude$modname <- "5-17"
results_18_55$crude$modname <- "18-54"
results_55_plus$crude$modname <- "55 and over"
results_comb_symp_only <- rbind(results_all$crude,results_u18$crude, results_18_55$crude, results_55_plus$crude)

results_comb_symp_only_adj <- rbind(results_all$adj,results_u18$adj, results_18_55$adj, results_55_plus$adj)

pcomb_symp_only <- plotUnivResults(results_comb_symp_only)
pcomb_symp_only
ggsave(filename = paste0(figpath, "uni_ORs_crude_exclude_asymptomatics.png"), plot = pcomb_symp_only, 
       width = 7, height=8, 
       dpi = 300, units = "in")


### Save univariate results
write_csv(results_comb_symp_only, paste0(outpath,"univariate_results_crude_exclude_asymptomatics.csv"))



