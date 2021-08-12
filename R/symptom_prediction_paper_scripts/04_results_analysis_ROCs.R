## Import model results and analyse ###


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
                  "purrr", "httr", "htmltools","ggalluvial","datapasta","xgboost","SHAPforxgboost", "plotROC"
)

load_packages(package.list)


############################ Run data import script ########################
r8_include_switch=F
source("E:/Group/react2_study5/report_phases_combined/projects/react_1_antigen_symptom_prediction/code/00_data_prep.R")
############################ Run data import script ########################



load(paste0(outpath, "m_b_stab_results.rdata"))

agelabels <- c("Age 0-17","Age 18-54","Age 55+")




# Run spec / tp analyses --------------------------------------------------


## function to calculate proportion of positives identified
calculateStats <- function(data, data_subset = NULL, age_subset){
  if(!is.null(data_subset)){
    dat <- data %>% filter(Data %in% data_subset, age_group %in% age_subset)
  }else{
    dat <- data %>% filter(age %in% age_subset, preds_new <= thresh)
  }
  n = nrow(dat)
  res = data.frame(thresh=seq(0,1,0.001), npredpos = NA_real_, 
                   npos= NA_real_, tp= NA_real_, fp= NA_real_,
                   tn= NA_real_,fn= NA_real_,sens= NA_real_, spec= NA_real_,
                   prop_pos_identified = NA_real_)
  
  t <- 0.05
  i <- 1
  for (t in seq(0,1,0.001)){
    print(t)
    npredpos = sum(dat$preds_new >= t)
    npredneg = sum(dat$preds_new < t)
    npos = sum(dat$estbinres)
    nneg = n-sum(dat$estbinres)
    prev = npos/(npos+nneg)
    tp = sum(dat$preds_new >= t &  dat$estbinres == 1)
    fp = sum(dat$preds_new >= t &  dat$estbinres == 0)
    tn = sum(dat$preds_new < t &  dat$estbinres == 0)
    fn = sum(dat$preds_new < t &  dat$estbinres == 1)
    prop_pos_identified = tp/npos
    sens = tp/(tp+fn)
    spec = tn/(tn+fp)
    res[i,] <- c(t,npredpos, npos, tp,fp,tn,fn, sens, spec,prop_pos_identified)
    i <- i+1
  }
  
  return(res)
}




### Age group 1
res_1 <- calculateStats(data=reclass_df_comb,
                        data_subset="Train",
                        age_subset=agenames[[1]])

### Age group 2
res_2 <- calculateStats(data=reclass_df_comb,
                        data_subset="Train",
                        age_subset=agenames[[2]])

### Age group 3
res_3 <- calculateStats(data=reclass_df_comb,
                        data_subset="Train",
                        age_subset=agenames[[3]])

###### TEST data

### Age group 1
res_4 <- calculateStats(data=reclass_df_comb,
                        data_subset="Test",
                        age_subset=agenames[[1]])

### Age group 2
res_5 <- calculateStats(data=reclass_df_comb,
                        data_subset="Test",
                        age_subset=agenames[[2]])

### Age group 3
res_6 <- calculateStats(data=reclass_df_comb,
                        data_subset="Test",
                        age_subset=agenames[[3]])






# ROC curves --------------------------------------------------------------



### Function to plot nice roc curves
plotMyROC <- function(preds_select, preds_classic, actual, age_stratum, preds_one_of_four){
  
  ### Create ROC objects
  roc.obj <- roc(actual, preds_select)
  roc.obj.2 <- roc(actual,preds_classic)
  roc.obj.3 <- roc(actual,preds_one_of_four)
  
  ### create CI objects
  ci.obj <- ci.sp(roc.obj, sensitivities = seq(0,1,0.01), boot.n = 500, progress ="none")
  ci.obj.2 <- ci.sp(roc.obj.2, sensitivities = seq(0,1,0.01), boot.n = 500, progress ="none")
  ci.obj.3 <- ci.sp(roc.obj.3, sensitivities = seq(0,1,0.01), boot.n = 500, progress ="none")
  
  
  ### Create plot
  plot(roc(actual, preds_select), print.auc = F, col = rgb(1,0.1,0.1, alpha = 1),
       xlim = c(1,0), ylim = c(0,1))
  plot(ci.obj, type = "shape", col = rgb(1,0,0, alpha = 0.1),
       bg = rgb(1,0,0, alpha = 0.3),
       lwd = 0.01,
       lty = 1)
  
  plot(ci.obj.2, type = "shape", col = rgb(0.1,0.1,1, alpha = 0.2),
       lwd = 0.01,
       lty = 1)
  
  plot(ci.obj.3, type = "shape", col = rgb(0.1,0.1,0.1, alpha = 0.2),
       lwd = 0.01,
       lty = 1)
  
  lines(roc(actual,preds_classic),lty=1, col = rgb(0.1,0.1,1, alpha = 1))
  lines(roc(actual,preds_one_of_four),lty=2, col = rgb(0.1,0.1,0.1, alpha = 1))
  lines(roc(actual, preds_select), print.auc = F, col = rgb(1,0.1,0.1, alpha = 1),
        xlim = c(1,0), ylim = c(0,1))
  
  
  title(str_to_sentence(gsub("_", " ", age_stratum)), line =3)
  auc.1 <- pROC::ci.auc(roc.obj, conf.level = 0.95, method="bootstrap",
                        boot.n =500,boot.stratified = T, reuse.auc=T, progress ="none")
  auc.2 <- pROC::ci.auc(roc.obj.2, conf.level = 0.95, method="bootstrap",
                        boot.n =500,boot.stratified = T, reuse.auc=T, progress ="none")
  auc.3 <- pROC::ci.auc(roc.obj.3, conf.level = 0.95, method="bootstrap",
                        boot.n =500,boot.stratified = T, reuse.auc=T, progress ="none")
  
  legend("bottomright", inset = c(0,0.05),lty=c(1,1,1), bty="n",
         col = c(col = rgb(1,0.1,0.1, alpha = 1), rgb(0.1,0.1,1, alpha = 1), 
                 rgb(0.1,0.1,0.1, alpha = 1)),
         legend=c(paste0("Stability model: AUC = ",
                         round(auc.1[2],3)," [",
                         round(auc.1[1],3),"-",
                         round(auc.1[3],3),"]"),
                  paste0("Four symptom model: AUC = ",
                         round(auc.2[2],3)," [",
                         round(auc.2[1],3),"-",
                         round(auc.2[3],3),"]"),
                  paste0("Any of four: AUC = ",
                         round(auc.3[2],3)," [",
                         round(auc.3[1],3),"-",
                         round(auc.3[3],3),"]")))
  
}




### Create combined plots
png(paste0(figpath, "mein_buhl_stability_increment_ROC_optimal_train_combined.png"),
    width = 15,height = 5, units = "in", res = 300)
layout.matrix <- matrix(c(1,2,3), nrow=1, ncol=3)
graphics::layout(mat = layout.matrix, heights = c(2), widths = c(2,2,2), respect = layout.matrix)
# layout.show(6)
plotMyROC(preds_select = reclass_df_comb %>% filter(Data == "Train", age_group == agenames[[1]]) %>% select(preds_new) %>% unlist(),
          preds_classic = reclass_df_comb %>% filter(Data == "Train", age_group == agenames[[1]]) %>% select(preds) %>% unlist(),
          preds_one_of_four =reclass_df_comb %>% filter(Data == "Train", age_group == agenames[[1]]) %>% select(one_of_four) %>% unlist(),
          actual = reclass_df_comb %>% filter(Data == "Train", age_group == agenames[[1]]) %>% select(estbinres) %>% unlist(),
          age_stratum = agelabels[[1]]
)


plotMyROC(preds_select = reclass_df_comb %>% filter(Data == "Train", age_group == agenames[[2]]) %>% select(preds_new) %>% unlist(),
          preds_classic = reclass_df_comb %>% filter(Data == "Train", age_group == agenames[[2]]) %>% select(preds) %>% unlist(),
          preds_one_of_four =reclass_df_comb %>% filter(Data == "Train", age_group == agenames[[2]]) %>% select(one_of_four) %>% unlist(),
          actual = reclass_df_comb %>% filter(Data == "Train", age_group == agenames[[2]]) %>% select(estbinres) %>% unlist(),
          age_stratum = agelabels[[2]]
)


plotMyROC(preds_select = reclass_df_comb %>% filter(Data == "Train", age_group == agenames[[3]]) %>% select(preds_new) %>% unlist(),
          preds_classic = reclass_df_comb %>% filter(Data == "Train", age_group == agenames[[3]]) %>% select(preds) %>% unlist(),
          preds_one_of_four =reclass_df_comb %>% filter(Data == "Train", age_group == agenames[[3]]) %>% select(one_of_four) %>% unlist(),
          actual = reclass_df_comb %>% filter(Data == "Train", age_group == agenames[[3]]) %>% select(estbinres) %>% unlist(),
          age_stratum = agelabels[[2]]
)
# 
# ### Plot second row of TP pickups
# 
# plot(res_1$spec, res_1$prop_pos_identified, type = "s",col = rgb(1,0,0, alpha = 1),
#      xlab = "Specificity", ylab = "Proportion of positives identified",
#      xlim = rev(range(res_1$spec)))
# 
# 
# plot(res_2$spec, res_2$prop_pos_identified, type = "s",col = rgb(1,0,0, alpha = 1),
#      xlab = "Specificity", ylab = "Proportion of positives identified",
#      xlim = rev(range(res_2$spec)))
# 
# 
# plot(res_3$spec, res_3$prop_pos_identified, type = "s",col = rgb(1,0,0, alpha = 1),
#      xlab = "Specificity", ylab = "Proportion of positives identified",
#      xlim = rev(range(res_3$spec)))

dev.off()


### Create combined plots
pdf(paste0(figpath, "mein_buhl_stability_increment_ROC_optimal_train_combined.pdf"),
    width = 15,height = 5, units = "in", res = 300)
layout.matrix <- matrix(c(1,2,3), nrow=1, ncol=3)
graphics::layout(mat = layout.matrix, heights = c(2), widths = c(2,2,2), respect = layout.matrix)
# layout.show(6)
plotMyROC(preds_select = reclass_df_comb %>% filter(Data == "Train", age_group == agenames[[1]]) %>% select(preds_new) %>% unlist(),
          preds_classic = reclass_df_comb %>% filter(Data == "Train", age_group == agenames[[1]]) %>% select(preds) %>% unlist(),
          preds_one_of_four =reclass_df_comb %>% filter(Data == "Train", age_group == agenames[[1]]) %>% select(one_of_four) %>% unlist(),
          actual = reclass_df_comb %>% filter(Data == "Train", age_group == agenames[[1]]) %>% select(estbinres) %>% unlist(),
          age_stratum = agelabels[[1]]
)


plotMyROC(preds_select = reclass_df_comb %>% filter(Data == "Train", age_group == agenames[[2]]) %>% select(preds_new) %>% unlist(),
          preds_classic = reclass_df_comb %>% filter(Data == "Train", age_group == agenames[[2]]) %>% select(preds) %>% unlist(),
          preds_one_of_four =reclass_df_comb %>% filter(Data == "Train", age_group == agenames[[2]]) %>% select(one_of_four) %>% unlist(),
          actual = reclass_df_comb %>% filter(Data == "Train", age_group == agenames[[2]]) %>% select(estbinres) %>% unlist(),
          age_stratum = agelabels[[2]]
)


plotMyROC(preds_select = reclass_df_comb %>% filter(Data == "Train", age_group == agenames[[3]]) %>% select(preds_new) %>% unlist(),
          preds_classic = reclass_df_comb %>% filter(Data == "Train", age_group == agenames[[3]]) %>% select(preds) %>% unlist(),
          preds_one_of_four =reclass_df_comb %>% filter(Data == "Train", age_group == agenames[[3]]) %>% select(one_of_four) %>% unlist(),
          actual = reclass_df_comb %>% filter(Data == "Train", age_group == agenames[[3]]) %>% select(estbinres) %>% unlist(),
          age_stratum = agelabels[[2]]
)
# 
# ### Plot second row of TP pickups
# 
# plot(res_1$spec, res_1$prop_pos_identified, type = "s",col = rgb(1,0,0, alpha = 1),
#      xlab = "Specificity", ylab = "Proportion of positives identified",
#      xlim = rev(range(res_1$spec)))
# 
# 
# plot(res_2$spec, res_2$prop_pos_identified, type = "s",col = rgb(1,0,0, alpha = 1),
#      xlab = "Specificity", ylab = "Proportion of positives identified",
#      xlim = rev(range(res_2$spec)))
# 
# 
# plot(res_3$spec, res_3$prop_pos_identified, type = "s",col = rgb(1,0,0, alpha = 1),
#      xlab = "Specificity", ylab = "Proportion of positives identified",
#      xlim = rev(range(res_3$spec)))

dev.off()







# Test set ----------------------------------------------------------------






### Create combined plots
png(paste0(figpath, "mein_buhl_stability_increment_ROC_optimal_test_combined.png"),
    width = 15,height = 5, units = "in", res = 300)
layout.matrix <- matrix(c(1,2,3), nrow=1, ncol=3)
graphics::layout(mat = layout.matrix, heights = c(2), widths = c(2,2,2), respect = layout.matrix)
# layout.show(6)
plotMyROC(preds_select = reclass_df_comb %>% filter(Data == "Test", age_group == agenames[[1]]) %>% select(preds_new) %>% unlist(),
          preds_classic = reclass_df_comb %>% filter(Data == "Test", age_group == agenames[[1]]) %>% select(preds) %>% unlist(),
          preds_one_of_four =reclass_df_comb %>% filter(Data == "Test", age_group == agenames[[1]]) %>% select(one_of_four) %>% unlist(),
          actual = reclass_df_comb %>% filter(Data == "Test", age_group == agenames[[1]]) %>% select(estbinres) %>% unlist(),
          age_stratum = agelabels[[1]]
)


plotMyROC(preds_select = reclass_df_comb %>% filter(Data == "Test", age_group == agenames[[2]]) %>% select(preds_new) %>% unlist(),
          preds_classic = reclass_df_comb %>% filter(Data == "Test", age_group == agenames[[2]]) %>% select(preds) %>% unlist(),
          preds_one_of_four =reclass_df_comb %>% filter(Data == "Test", age_group == agenames[[2]]) %>% select(one_of_four) %>% unlist(),
          actual = reclass_df_comb %>% filter(Data == "Test", age_group == agenames[[2]]) %>% select(estbinres) %>% unlist(),
          age_stratum = agelabels[[2]]
)


plotMyROC(preds_select = reclass_df_comb %>% filter(Data == "Test", age_group == agenames[[3]]) %>% select(preds_new) %>% unlist(),
          preds_classic = reclass_df_comb %>% filter(Data == "Test", age_group == agenames[[3]]) %>% select(preds) %>% unlist(),
          preds_one_of_four =reclass_df_comb %>% filter(Data == "Test", age_group == agenames[[3]]) %>% select(one_of_four) %>% unlist(),
          actual = reclass_df_comb %>% filter(Data == "Test", age_group == agenames[[3]]) %>% select(estbinres) %>% unlist(),
          age_stratum = agelabels[[2]]
)
# 
# ### Plot second row of TP pickups
# 
# plot(res_4$spec, res_4$prop_pos_identified, type = "s",col = rgb(1,0,0, alpha = 1),
#      xlab = "Specificity", ylab = "Proportion of positives identified",
#      xlim = rev(range(res_4$spec)))
# 
# 
# plot(res_5$spec, res_5$prop_pos_identified, type = "s",col = rgb(1,0,0, alpha = 1),
#      xlab = "Specificity", ylab = "Proportion of positives identified",
#      xlim = rev(range(res_5$spec)))
# 
# 
# plot(res_6$spec, res_6$prop_pos_identified, type = "s",col = rgb(1,0,0, alpha = 1),
#      xlab = "Specificity", ylab = "Proportion of positives identified",
#      xlim = rev(range(res_6$spec)))

dev.off()



##### Now pdf version 3####

### Create combined plots
pdf(paste0(figpath, "mein_buhl_stability_increment_ROC_optimal_test_combined.pdf"),
    width = 15,height = 5, units = "in", res = 300)
layout.matrix <- matrix(c(1,2,3), nrow=1, ncol=3)
graphics::layout(mat = layout.matrix, heights = c(2), widths = c(2,2,2), respect = layout.matrix)
# layout.show(6)
plotMyROC(preds_select = reclass_df_comb %>% filter(Data == "Test", age_group == agenames[[1]]) %>% select(preds_new) %>% unlist(),
          preds_classic = reclass_df_comb %>% filter(Data == "Test", age_group == agenames[[1]]) %>% select(preds) %>% unlist(),
          preds_one_of_four =reclass_df_comb %>% filter(Data == "Test", age_group == agenames[[1]]) %>% select(one_of_four) %>% unlist(),
          actual = reclass_df_comb %>% filter(Data == "Test", age_group == agenames[[1]]) %>% select(estbinres) %>% unlist(),
          age_stratum = agelabels[[1]]
)


plotMyROC(preds_select = reclass_df_comb %>% filter(Data == "Test", age_group == agenames[[2]]) %>% select(preds_new) %>% unlist(),
          preds_classic = reclass_df_comb %>% filter(Data == "Test", age_group == agenames[[2]]) %>% select(preds) %>% unlist(),
          preds_one_of_four =reclass_df_comb %>% filter(Data == "Test", age_group == agenames[[2]]) %>% select(one_of_four) %>% unlist(),
          actual = reclass_df_comb %>% filter(Data == "Test", age_group == agenames[[2]]) %>% select(estbinres) %>% unlist(),
          age_stratum = agelabels[[2]]
)


plotMyROC(preds_select = reclass_df_comb %>% filter(Data == "Test", age_group == agenames[[3]]) %>% select(preds_new) %>% unlist(),
          preds_classic = reclass_df_comb %>% filter(Data == "Test", age_group == agenames[[3]]) %>% select(preds) %>% unlist(),
          preds_one_of_four =reclass_df_comb %>% filter(Data == "Test", age_group == agenames[[3]]) %>% select(one_of_four) %>% unlist(),
          actual = reclass_df_comb %>% filter(Data == "Test", age_group == agenames[[3]]) %>% select(estbinres) %>% unlist(),
          age_stratum = agelabels[[2]]
)
# 
# ### Plot second row of TP pickups
# 
# plot(res_4$spec, res_4$prop_pos_identified, type = "s",col = rgb(1,0,0, alpha = 1),
#      xlab = "Specificity", ylab = "Proportion of positives identified",
#      xlim = rev(range(res_4$spec)))
# 
# 
# plot(res_5$spec, res_5$prop_pos_identified, type = "s",col = rgb(1,0,0, alpha = 1),
#      xlab = "Specificity", ylab = "Proportion of positives identified",
#      xlim = rev(range(res_5$spec)))
# 
# 
# plot(res_6$spec, res_6$prop_pos_identified, type = "s",col = rgb(1,0,0, alpha = 1),
#      xlab = "Specificity", ylab = "Proportion of positives identified",
#      xlim = rev(range(res_6$spec)))

dev.off()

