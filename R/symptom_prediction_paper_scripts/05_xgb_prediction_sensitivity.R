### Data prep for React 1 ###



#' First clear the environment of variables
rm(list=ls(all=TRUE))
setwd("E:/Group/react2_study5/report_phases_combined/projects/react_1_antigen_symptom_prediction/")
outpath <- "E:/Group/react2_study5/report_phases_combined/projects/react_1_antigen_symptom_prediction/output/"
figpath <- "E:/Group/react2_study5/report_phases_combined/projects/react_1_antigen_symptom_prediction/plots/"

#' Source any functions from the local file
source("E:/Group/functions/load_packages.R")
source("E:/Group/functions/hamming_distance.R")
source("E:/Group/functions/cats_and_covs.R")
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
                  "ggplot2","ggsci", "RColorBrewer", "tidyverse", "lubridate", 
                  "egg", "poLCA", "Rcpp","xml2","splitstackshape",
                  "fs", "later", "promises","proxy","dendextend", 
                  "ComplexHeatmap","circlize","doSNOW","ClusterR","htmlwidgets",
                  "readr","ggthemes", "questionr", "gridExtra", "foreach", 
                  "doParallel","glmnet","pROC","DiagrammeR",
                  "purrr", "httr", "htmltools","ggalluvial","rsvg",
                  "datapasta","xgboost","SHAPforxgboost", "patchwork"
)

load_packages(package.list)


############################ Run data import script ########################
r8_include_switch = F
source("E:/Group/react2_study5/report_phases_combined/projects/react_1_antigen_symptom_prediction/code/00_data_prep.R")
############################ Run data import script ########################


### 
ages <- list(under_18 = c(0:17),from_18_to_54 = c(18:54), over_55 = c(55:200))
agenames <- c("under_18","from_18_to_54","over_55")

pred.results.list  <-list()
i <- 2

for (i in 1:3){
  
  
  ### create training data
  Xdata <- mydata_train %>% filter(age %in% ages[[i]]) %>% select(varnames) %>% as.matrix()
  Ydata <- mydata_train %>% filter(age %in% ages[[i]]) %>% select(estbinres)  %>% as.matrix()
  Xdata_test <- mydata_test %>% filter( age %in% ages[[i]]) %>% select(varnames) %>% as.matrix()
  Ydata_test <- mydata_test %>% filter( age %in% ages[[i]]) %>% select(estbinres)  %>% as.matrix()
  
  ### train signle tree XGB model
  obj="binary:logistic"
  eval_metric = "logloss"
  
  one_tree_mod <- xgboost(
    data = Xdata,
    label =  Ydata,
    nrounds =1,
    params = list(max_depth =3),
    print_every_n = 10,
    eval_metric = eval_metric
    # eval_metric = "error",
  )
  
 
    plt <- xgb.plot.tree(model = one_tree_mod, trees = 0,show_node_id = T, render = F)
    xgb.plot.tree(model = one_tree_mod, trees = 0,show_node_id = T, render = T)
    export_graph(plt, paste0(figpath,"decision_tree_", agenames[[i]], ".png"))
    # 

    
  ### Run full XGB model
  t0 <- Sys.time()
  dtrain <- xgboost::xgb.DMatrix(data = Xdata, label = Ydata)
  estop=10
  nfold = 3
  max_depth_list <- c(2:7)
  best_logloss <- 10
  best_max_depth <- 3
  for(maxd in max_depth_list){
    params_phil <- list(booster ="gbtree", 
                        eta = 0.05,
                        gamma = 0,
                        max_depth = maxd, obj  = "binary:logistic"
                        # min_child_weight = 2, subsample = 0.64,
                        # colsample_bytree = 0.71
    )
    results.xgboost.cv <- xgboost::xgb.cv(data = dtrain,
                                          nrounds = 100,
                                          eval_metric = "logloss",
                                          nfold = nfold,
                                          params = params_phil,
                                          print_every_n = 10,
                                          early_stopping_rounds = estop,
                                          maximize = F)
  
    
    if(min(results.xgboost.cv$evaluation_log$test_logloss_mean) <best_logloss){
      best_logloss <- min(results.xgboost.cv$evaluation_log$test_logloss_mean)
      best_max_depth <- maxd
      best_nrounds <- results.xgboost.cv$best_iteration
    }
  }
  best_max_depth
  best_nrounds
  
  params_phil <- list(booster ="gbtree", 
                      eta = 0.05, 
                      gamma = 0,
                      max_depth = best_max_depth, 
                      obj  = "binary:logistic"
                      # min_child_weight = 2, subsample = 0.64,
                      # colsample_bytree = 0.71
  )
  # min(results.xgboost.cv$evaluation_log$test_logloss_mean)
  
  
  results.xgboost <- xgboost::xgb.train(data = dtrain,
                                        nrounds = best_nrounds,
                                 eval_metric = "logloss",
                                 params = params_phil,
                                 maximize = F)
  
  
  print(Sys.time() - t0)
  
  
  
  ### Lasso model
  lasso.mod <- cv.glmnet(x = Xdata, y = Ydata, family = "binomial", alpha = 1)
  lambdaopt <- lasso.mod$lambda.1se
  
  ### create results df
  df.results <- data.frame(model = c("Single tree", "Boosted trees", "Lasso"), 
                           age_group = NA_character_,
                           auc_train_lower = NA_real_, 
                           auc_train= NA_real_, 
                           auc_train_upper= NA_real_,
                           auc_test_lower = NA_real_, 
                           auc_test= NA_real_, 
                           auc_test_upper= NA_real_)
  
  mods <-list(one_tree_mod,results.xgboost, lasso.mod)
x <- 3
    for (x in 1:3){
    ### get predictions
      if(x==3){
        preds_train <- predict(object = mods[[x]],newx   = Xdata, type = "link", s = lambdaopt)
        preds_test <- predict(object = mods[[x]],newx = Xdata_test, type = "link", s = lambdaopt)
      }else{
        preds_train <- predict(object = mods[[x]],newdata = Xdata)
        preds_test <- predict(object = mods[[x]],newdata = Xdata_test)
      }
    
    
    ### Get metrics (train)
    auc_train <- pROC::auc(response=as.vector(Ydata), predictor=as.vector(preds_train))
    roc_train <- pROC::roc(response=as.vector(Ydata), predictor=as.vector(preds_train))
    auc_ci_train <- pROC::ci.auc(roc_train, conf.level = 0.95, method="bootstrap",
                                 boot.n =1000,boot.stratified = T, reuse.auc=T)
    ### Get metrics (test)
    auc_test <- pROC::auc(response=as.vector(Ydata_test), predictor=as.vector(preds_test))
    roc_test<- pROC::roc(response=as.vector(Ydata_test), predictor=as.vector(preds_test))
    auc_ci_test <- pROC::ci.auc(roc_test, conf.level = 0.95, method="bootstrap",
                                boot.n =1000,boot.stratified = T, reuse.auc=T)
    
    
    res <- c(agenames[[i]],round(c(auc_ci_train,auc_ci_test),4))
    df.results[x,2:ncol(df.results)] <- res
    
    
  }
  pred.results.list[[i]] <- df.results
}

df.results.all <- do.call(rbind,pred.results.list)

df.results.all.train <- df.results.all[,c("model","age_group","auc_train_lower", "auc_train","auc_train_upper" )]
df.results.all.test <- df.results.all[,c("model","age_group","auc_test_lower", "auc_test","auc_test_upper" )]

names(df.results.all.train) <- names(df.results.all.test)  <- c("model","age_group","auc_lower", "auc","auc_upper" )
df.results.all.train$Data = "Train"
df.results.all.test$Data = "Test"

### rebind
df.results.all <- rbind(df.results.all.train,df.results.all.test)
df.results.all <- df.results.all %>% mutate(auc = as.numeric(auc),
                          auc_lower = as.numeric(auc_lower),
                          auc_upper = as.numeric(auc_upper))


### Load up Lasso AUC data
lasso_under_18 <- readRDS(paste0(outpath,"mein_buhl_stability_increment_under_18_results_auc.rds"))
lasso_18_54 <- readRDS(paste0(outpath,"mein_buhl_stability_increment_from_18_to_54_results_auc.rds"))
lasso_over_55 <- readRDS(paste0(outpath,"mein_buhl_stability_increment_over_55_results_auc.rds"))

### wrangle
lasso_under_18 <- lasso_under_18[lasso_under_18$symptom == lasso_under_18[nrow(lasso_under_18),]$symptom,c("auc_lower","auc", "auc_upper",  "Data")]
lasso_under_18$age_group <- agenames[[1]]
lasso_under_18$model <- "Lasso stability"
lasso_18_54 <- lasso_18_54[lasso_18_54$symptom == lasso_18_54[nrow(lasso_18_54),]$symptom,c("auc_lower","auc", "auc_upper",  "Data")]
lasso_18_54$age_group <- agenames[[2]]
lasso_18_54$model <- "Lasso stability"
lasso_over_55 <- lasso_over_55[lasso_over_55$symptom == lasso_over_55[nrow(lasso_over_55),]$symptom,c("auc_lower","auc", "auc_upper",  "Data")]
lasso_over_55$age_group <- agenames[[3]]
lasso_over_55$model <- "Lasso stability"


### Bind
df.results.all <- rbind(df.results.all,lasso_under_18,lasso_18_54,lasso_over_55)


### Plot results
dodge=position_dodge(width=0.5)
p1 <- df.results.all %>% 
  filter(!model %in% c("Single tree","Lasso"), Data != "Round_8") %>% 
  mutate(age_group_new = case_when(age_group == "under_18" ~ "5-17",
                                   age_group == "from_18_to_54"~"18-54",
                                   TRUE ~ "55 plus"),
         age_group_new = factor(age_group_new, levels = unique(age_group_new))) %>%
  mutate(model = factor(model, levels = unique(model))) %>%
  ggplot(aes(x=model, y= auc, col=Data)) +
  geom_point(position = dodge) +
  ylim(0.45, 0.9) +
  facet_wrap(.~ age_group_new) +
  geom_errorbar(aes(ymin = auc_lower, ymax =auc_upper), 
                width=0.5, size =0.5,position = dodge) +
  scale_color_manual(values = myCols[c(8,2)]) +
  labs(x="Model", y = "AUC") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust = 1),
          strip.background = element_rect(fill = "white"),
          strip.text = element_text(face = "bold"),
          axis.text.y = element_text(size = 10)
        )


p1
ggsave(filename = paste0(figpath, "xgboost_tree_forest_results_auc.png"), 
       plot = p1,width = 7, height=7, dpi = 300, units = "in")
ggsave(filename = paste0(figpath, "xgboost_tree_forest_results_auc.pdf"), 
       plot = p1,width = 7, height=7, dpi = 300, units = "in")


