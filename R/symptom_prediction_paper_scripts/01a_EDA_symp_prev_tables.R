#' First clear the environment of variables
rm(list=ls(all=TRUE))
setwd("E:/Group/react2_study5/report_phases_combined/projects/react_1_antigen_symptom_prediction/")
outpath <- "E:/Group/react2_study5/report_phases_combined/projects/react_1_antigen_symptom_prediction/round_compare_analysis/output/"
figpath <- "E:/Group/react2_study5/report_phases_combined/projects/react_1_antigen_symptom_prediction/round_compare_analysis/plots/"

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
source("E:/Group/functions/prev_difference_test.R", local = TRUE)
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
                  "purrr", "httr", "htmltools","ggalluvial","datapasta","xgboost","SHAPforxgboost",
                  "scales", "stringr"
)

load_packages(package.list)



############################ Run data import script ########################
r8_include_switch=T
source("E:/Group/react2_study5/report_phases_combined/projects/react_1_antigen_symptom_prediction/code/00_data_prep.R")
############################ Run data import script ########################
outpath <- "E:/Group/react2_study5/report_phases_combined/projects/react_1_antigen_symptom_prediction/round_compare_analysis/output/"
figpath <- "E:/Group/react2_study5/report_phases_combined/projects/react_1_antigen_symptom_prediction/round_compare_analysis/plots/"
source("E:/Group/react2_study5/report_phases_combined/projects/functions/univariate_models.R")

### new subfolder
dir.create(paste0(outpath,"symp_prevs_2_8_vs_8"), showWarnings = F)
mydata_symp$one_of_22 <- as.numeric(rowSums(mydata_symp[,varnames[c(5:26)]]) > 0 & mydata_symp$one_of_four == 0)

#### save data for r7&8 only
tabdat <- mydata_symp %>% filter(round %in% c(2:8))

## drop unused levels
tabdat$r2_7_vs_r8 <-  as.factor(ifelse(tabdat$round_split_8_all == 1, "Round 8", "Round 2-7"))
### select vars for checking
myVars <- c("symptomatic","one_of_four","one_of_22",varnames)
### Convert relevant variables to factor
tabdat[myVars] <- lapply(tabdat[myVars], factor)

### Save positive and negative dfs
tabdat_pos <- tabdat %>% filter(estbinres == 1)
tabdat_neg <- tabdat %>% filter(estbinres == 0)

dat <- tabdat_pos
table(tabdat$round_split_8_all,tabdat$symptomatic)



### quick function to run table ones
createMyTableOne <- function(dat, stratum = "r2_7_vs_r8"){
  
  ### Create table 1
  tab1 <- CreateTableOne(vars = myVars,
                         data = dat,
                         strata = stratum)
  ### 'print' table 1
  tab1mat <- print(tab1, quote=FALSE, noSpaces = TRUE, printToggle = FALSE, showAllLevels = TRUE) %>% 
    as.data.frame()
  
  ### Some wrangling
  pvals <- tab1mat$p
  pvals <- pvals[!is.na(pvals) & pvals!= ""] 
  
  ### filter negatives
  tab1mat <- tab1mat %>% filter(level != 0)
  
  ### Add sympnames
  tab1mat$level <- c("Full cohort", "Symptomatic",
                     "Classic symtoms",
                     "Non-classic symptoms only",varnames)
  ### add test result
  tab1mat$p <- c(NA_real_,pvals)
  
  tab1mat$symptom_class <- c(NA_character_,
                             "Symptomatic",
                             "Classic symtoms",
                             "Non-classic symptoms only",
                             rep("Classic",4),
                             rep("Cold/flu like",6),
                             rep("Headache /dizziness",2),
                             rep("Gastrointestinal",4),
                             rep("Respiratory",3),
                             "Cold/flu like",
                             rep("Fatigue",3),
                             rep("Somatic",3)
                             
  )
  
  
  ### remove redundant colum
  tab1mat <- tab1mat %>% select(-test)
  return(tab1mat)
}


### Run table ones 7 vs 8
tab1_all <- createMyTableOne(tabdat)
tab1_pos <- createMyTableOne(tabdat_pos)
tab1_neg <- createMyTableOne(tabdat_neg)

write.csv(tab1_all, paste0(outpath, "symp_prevs/symp_and_asymp_prev_all_r2_7_vs_r8.csv"))
write.csv(tab1_pos, paste0(outpath, "symp_prevs/symp_and_asymp_prev_pos_r2_7_vs_r8.csv"))
write.csv(tab1_neg, paste0(outpath, "symp_prevs/symp_and_asymp_prev_neg_r2_7_vs_r8.csv"))

