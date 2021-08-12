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
source("E:/Group/react2_study5/report_phases_combined/projects/functions/barbara_functions/penalisation_functions.R")


#' Pull in packages needed
package.list <- c("prevalence","mgcv","knitr","MASS","kableExtra","table1","dplyr","factoextra","tableone","networkD3",
                  "tidyr","forcats", "cluster", "fpc", "mclust", "pheatmap","FactoMineR", "NbClust","clValid","plotly",
                  "ggplot2","ggsci", "RColorBrewer", "tidyverse", "lubridate", "egg", "poLCA", "Rcpp","xml2","splitstackshape",
                  "fs", "later", "promises","proxy","dendextend", "ComplexHeatmap","circlize","doSNOW","ClusterR","htmlwidgets",
                  "readr","ggthemes", "questionr", "gridExtra", "foreach", "doParallel","glmnet","caTools","pROC","patchwork",
                  "purrr", "httr", "htmltools","ggalluvial","datapasta","xgboost",
                  "SHAPforxgboost", "plotROC","ggupset", "UpSetR", "ComplexUpset"
)

load_packages(package.list)


############################ Run data import script ########################
r8_include_switch = F
source("E:/Group/react2_study5/report_phases_combined/projects/react_1_antigen_symptom_prediction/code/00_data_prep.R")
############################ Run data import script ########################



### 
ages <- list(under_18 = c(0:17),from_18_to_54 = c(18:54), over_55 = c(55:200))
agenames <- c("under_18","from_18_to_54","over_55")
agelabels <- c("Age 0-17","Age 18-54","Age 55+")
varnames_classic <- c("loss_of_sense_of_smell","loss_of_sense_of_taste","new_persistent_cough","fever")






# Complex Upsetr ----------------------------------------------------------

posdat <- mydata %>% filter(estbinres==1, !is.na(age)) %>% 
  as.data.frame() %>% select(varnames_classic,age) %>% 
  rename_all(funs(gsub("[[:punct:]]", " ",make.names(c(varnames_classic, "age")))))
prop.table(table(rowSums(posdat[,1:2])== 2))

### Plot big 4 symps (top 10 sets)
up_big_4_pos <- ComplexUpset::upset(posdat,
                                    intersect = list("fever", "new persistent cough",
                                                     "loss or change of sense of taste","loss or change of sense of smell"),
                                    max_degree = 4, sort_sets = F, sort_intersections = "descending",
                                    min_size = 0,
                                    wrap = T,
                                    themes=upset_modify_themes(
                                      list('intersections_matrix' = theme(text = element_text(size=20)))
                                    ),
                                    base_annotations = list('Size'=(intersection_size(counts=F,
                                                                                       mode='inclusive_union')),
                                                            # with manual aes specification:
                                                            'Intersection size'=intersection_size(text_mapping=aes(label=paste0(round(
                                                              !!get_size_mode('exclusive_intersection')/nrow(posdat) * 100
                                                            ), '%')))),
                                    queries=list(
                                      upset_query(set =c("fever", "new persistent cough",
                                                        "loss or change of sense of taste","loss or change of sense of smell"),
                                                  color = 'blue',
                                                  fill = 'blue',
                                                  only_components = c('intersections_matrix',
                                                                      'Intersection size')),
                                      # upset_query(group="fever", color='blue'),
                                      # upset_query(group="new persistent cough", color='blue'),
                                      # upset_query(group="loss or change of sense of taste", color='blue'),
                                      # upset_query(group="loss or change of sense of smell", color='blue'),
                                      upset_query(set="fever", fill='blue'),
                                      upset_query(set="new persistent cough", fill='blue'),
                                      upset_query(set="loss or change of sense of taste", fill='blue'),
                                      upset_query(set="loss or change of sense of smell", fill='blue')
                                    )
                                    )
                                    
up_big_4_pos     


### Negatives
negdat <- mydata %>% filter(estbinres==0, !is.na(age)) %>% 
  as.data.frame() %>% select(varnames_classic,age) %>% 
  rename_all(funs(gsub("[[:punct:]]", " ",make.names(c(varnames_classic, "age")))))
prop.table(table(rowSums(posdat[,1:2])== 2))

### Plot big 4 symps (top 10 sets)
up_big_4_neg <- ComplexUpset::upset(negdat,
                                    intersect = list("fever", "new persistent cough",
                                                     "loss or change of sense of taste","loss or change of sense of smell"),
                                    max_degree = 4, sort_sets = F, sort_intersections = "descending",
                                    min_size = 0,
                                    wrap = T,
                                    themes=upset_modify_themes(
                                      list('intersections_matrix' = theme(text = element_text(size=20)))
                                    ),
                                    base_annotations = list('Size'=(intersection_size(counts=F,
                                                                                      mode='inclusive_union')),
                                                            # with manual aes specification:
                                                            'Intersection size'=intersection_size(text_mapping=aes(label=paste0(round(
                                                              !!get_size_mode('exclusive_intersection')/nrow(negdat) * 100
                                                            ), '%')))),
                                    queries=list(
                                      upset_query(set =c("fever", "new persistent cough",
                                                         "loss or change of sense of taste","loss or change of sense of smell"),
                                                  color = 'blue',
                                                  fill = 'blue',
                                                  only_components = c('intersections_matrix',
                                                                      'Intersection size')),
                                      # upset_query(group="fever", color='blue'),
                                      # upset_query(group="new persistent cough", color='blue'),
                                      # upset_query(group="loss or change of sense of taste", color='blue'),
                                      # upset_query(group="loss or change of sense of smell", color='blue'),
                                      upset_query(set="fever", fill='blue'),
                                      upset_query(set="new persistent cough", fill='blue'),
                                      upset_query(set="loss or change of sense of taste", fill='blue'),
                                      upset_query(set="loss or change of sense of smell", fill='blue')
                                    )
)



png(paste0(figpath, "upset_big4_pos.png"), width = 17, height=10, res = 300, units = "in")
up_big_4_pos
dev.off()
pdf(paste0(figpath, "upset_big4_pos.pdf"), width = 17, height=10)
up_big_4_pos
dev.off()



png(paste0(figpath, "upset_big4_neg.png"), width = 17, height=10, res = 300, units = "in")
up_big_4_neg
dev.off()
pdf(paste0(figpath, "upset_big4_neg.pdf"), width = 17, height=10)
up_big_4_neg
dev.off()




# Now a version with flippin loads of variables ---------------------------

varnames_spaced <- gsub("_", " ", varnames)
varnames_classic_spaced <- gsub("_", " ", varnames_classic)

posdat <- mydata %>% filter(estbinres==1, !is.na(age)) %>% 
  as.data.frame() %>% select(varnames,age) %>% 
  rename_all(funs(gsub("[[:punct:]]", " ",make.names(c(varnames, "age")))))
prop.table(table(rowSums(posdat[,1:2])== 2))






### Plot big 4 symps (top 10 sets)
up_all_pos <- ComplexUpset::upset(posdat,
                                    intersect = list(varnames_spaced),
                                    max_degree = 10, sort_sets = F, sort_intersections = "descending",
                                    min_size = 3,
                                    wrap = T,
                                    themes=upset_modify_themes(
                                      list('intersections_matrix' = theme(text = element_text(size=20)))
                                    ),
                                    base_annotations = list('Size'=(intersection_size(counts=F,
                                                                                      mode='inclusive_union')),
                                                            # with manual aes specification:
                                                            'Intersection size'=intersection_size(text_mapping=aes(label=paste0(round(
                                                              !!get_size_mode('exclusive_intersection')/nrow(posdat) * 100
                                                            ), '%')))),
                                    queries=list(
                                      upset_query(set =c("fever", "new persistent cough",
                                                         "loss or change of sense of taste","loss or change of sense of smell"),
                                                  color = 'blue',
                                                  fill = 'blue',
                                                  only_components = c('intersections_matrix',
                                                                      'Intersection size')),
                                      # upset_query(group="fever", color='blue'),
                                      # upset_query(group="new persistent cough", color='blue'),
                                      # upset_query(group="loss or change of sense of taste", color='blue'),
                                      # upset_query(group="loss or change of sense of smell", color='blue'),
                                      upset_query(set="fever", fill='blue'),
                                      upset_query(set="new persistent cough", fill='blue'),
                                      upset_query(set="loss or change of sense of taste", fill='blue'),
                                      upset_query(set="loss or change of sense of smell", fill='blue')
                                    )
)

up_all_pos     


### Negatives
negdat <- mydata %>% filter(estbinres==0, !is.na(age)) %>% 
  as.data.frame() %>% select(varnames,age) %>% 
  rename_all(funs(gsub("[[:punct:]]", " ",make.names(c(varnames, "age")))))
prop.table(table(rowSums(posdat[,1:2])== 2))

### Plot big 4 symps (top 10 sets)
up_all_neg <- ComplexUpset::upset(negdat,
                                    intersect = list(varnames_spaced),
                                    max_degree = 10, sort_sets = F, sort_intersections = "descending",
                                    min_size = 200,
                                    wrap = T,
                                    themes=upset_modify_themes(
                                      list('intersections_matrix' = theme(text = element_text(size=20)))
                                    ),
                                    base_annotations = list('Size'=(intersection_size(counts=F,
                                                                                      mode='inclusive_union')),
                                                            # with manual aes specification:
                                                            'Intersection size'=intersection_size(text_mapping=aes(label=paste0(round(
                                                              !!get_size_mode('exclusive_intersection')/nrow(negdat) * 100
                                                            ), '%')))),
                                    queries=list(
                                      upset_query(set =c("fever", "new persistent cough",
                                                         "loss or change of sense of taste","loss or change of sense of smell"),
                                                  color = 'blue',
                                                  fill = 'blue',
                                                  only_components = c('intersections_matrix',
                                                                      'Intersection size')),
                                      # upset_query(group="fever", color='blue'),
                                      # upset_query(group="new persistent cough", color='blue'),
                                      # upset_query(group="loss or change of sense of taste", color='blue'),
                                      # upset_query(group="loss or change of sense of smell", color='blue'),
                                      upset_query(set="fever", fill='blue'),
                                      upset_query(set="new persistent cough", fill='blue'),
                                      upset_query(set="loss or change of sense of taste", fill='blue'),
                                      upset_query(set="loss or change of sense of smell", fill='blue')
                                    )
)

up_all_neg    


png(paste0(figpath, "upset_all_pos.png"), width = 20, height=17, res = 300, units = "in")
up_all_pos
dev.off()
pdf(paste0(figpath, "upset_all_pos.pdf"), width = 20, height=17)
up_all_pos
dev.off()



png(paste0(figpath, "upset_all_neg.png"), width = 20, height=17, res = 300, units = "in")
up_all_neg
dev.off()
pdf(paste0(figpath, "upset_all_neg.pdf"), width = 20, height=17)
up_all_neg
dev.off()

