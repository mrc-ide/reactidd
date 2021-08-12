### Run models ###



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
r8_include_switch = F
source("E:/Group/react2_study5/report_phases_combined/projects/react_1_antigen_symptom_prediction/code/00_data_prep.R")
############################ Run data import script ########################


### 
ages <- list(under_18 = c(0:17),from_18_to_54 = c(18:54), over_55 = c(55:200))
agenames <- c("under_18","from_18_to_54","over_55")
agelabels <- c("Age 5-17","Age 18-54","Age 55+")
varnames_classic <- c("loss_or_change_of_sense_of_smell","loss_or_change_of_sense_of_taste","new_persistent_cough","fever")
pred.results.list  <-list()
test.preds.list <- list()
r8.preds.list <- list()
train.preds.list <- list()
ids <- list(train = list(),
            test= list(),
            r8=list())
plots.list <- list()
i <- 1
auc_one_of_four <- list(train = c(0.628, 0.710,0.691),
                        test = c(0.641, 0.730, 0.661))

for (i in 1:3){
  
  ### create training data
  Xdata <- mydata_train %>% filter(age %in% ages[[i]]) %>% select(varnames) %>% as.matrix()
  Ydata <- mydata_train %>% filter(age %in% ages[[i]]) %>% select(estbinres)  %>% as.matrix()
  Xdata_test <- mydata_test %>% filter( age %in% ages[[i]]) %>% select(varnames) %>% as.matrix()
  Ydata_test <- mydata_test %>% filter( age %in% ages[[i]]) %>% select(estbinres)  %>% as.matrix()
  Xdata_r8 <- mydata_r8 %>% filter(age %in% ages[[i]]) %>% select(varnames) %>% as.matrix()
  Ydata_r8 <- mydata_r8 %>% filter( age %in% ages[[i]]) %>% select(estbinres)  %>% as.matrix()
  
  ### Save IDs
  ids$train[[i]] <- mydata_train %>% filter(age %in% ages[[i]]) %>% 
    select(new_id) %>% unlist() %>% 
    as.numeric()
  ids$test[[i]] <- mydata_test %>% filter(age %in% ages[[i]]) %>% 
    select(new_id) %>% unlist() %>% 
    as.numeric()
  ids$r8[[i]] <- mydata_r8 %>% filter(age %in% ages[[i]]) %>% 
    select(new_id) %>% unlist() %>% 
    as.numeric()
  
  
  ### Running stability selection
  
  t0=Sys.time()
  out=CalibrateRegression(xdata=Xdata, ydata=Ydata, family="binomial", K=1000, pi_list=seq(0.51,0.99,by=0.01))
  t1=Sys.time()
  print(as.numeric(difftime(t1,t0, units="secs")))
  
  
  ### Calibration (Meinshausen and Buhlmann)
  
  ## Parameters to set
  
  PFER_thr=5 # allowing 5 falsely selected variables
  hat_pi=0.51 # arbitrary threshold in selection proportion 
  hat_pi_true=0.50
  
  ## Calibration 
  
  N_block=ncol(out$selprop)
  K=out$params$K
  PFER=rep(NA, nrow(out$selprop))
  for (k in 1:nrow(out$selprop)){
    tmpselprop=out$selprop[k,]
    q=round(sum(tmpselprop))
    if (any(tmpselprop>=hat_pi)){
      pi=min(tmpselprop[tmpselprop>=hat_pi])
      PFER[k]=ComputePFER(q=q, pi=pi, N=N_block, K=K, method="MB")
    }
  }
  
  hat_lambda=out$Lambda[max(which(PFER<PFER_thr))]
  theta_lambda=which.min(abs(out$Lambda[,1]-hat_lambda))
  
  selprop=out$selprop[theta_lambda,]
  
  # theta_pi=which(out$params$pi_list==hat_pi)
  # hat_lambda=out$Lambda[max(which(out$PFER_2d[,theta_pi]<PFER_thr))]
  # theta_lambda=which.min(abs(out$Lambda[,1]-hat_lambda))
  
  # selprop=out$selprop[theta_lambda,]
  average_beta=apply(out$Beta[theta_lambda,,],1,FUN=function(x){mean(x[x!=0])})
  se_beta=apply(out$Beta[theta_lambda,,],1,FUN=function(x){sd(x[x!=0])}) / sqrt(100)
  
  
  ### Just to check
  
  a=apply(out$Beta[theta_lambda,,],1,FUN=function(x){sum(x>0)})
  b=apply(out$Beta[theta_lambda,,],1,FUN=function(x){sum(x<0)})
  # Just to make sure that good consistency in sign of the beta coefficients for the variables with high selprop
  plot(a/(a+b), selprop, las=1,
       xlab="Proportion of positive beta among non-zero betas", ylab="Selection Proportion")
  
  
  ### Running the model with average betas
  
  selprop_nonzero=selprop[selprop>0]
  myorder=names(selprop_nonzero)[sort.list(selprop_nonzero, decreasing = TRUE)]
  myselprops=selprop_nonzero[sort.list(selprop_nonzero, decreasing = TRUE)]
  num_include = length(myselprops[myselprops>=hat_pi_true])
  myaverage_beta=average_beta[selprop>0][sort.list(selprop_nonzero, decreasing = TRUE)]
  myse_beta=se_beta[selprop>0][sort.list(selprop_nonzero, decreasing = TRUE)]
  
  # Just one iteration to see 
  k=1
  mymodel=glm(Ydata~offset(Xdata[,myorder[1:k],drop=FALSE]%*%matrix(average_beta[myorder[1:k]], ncol=1)), family="binomial")
  coef(mymodel) # estimated intercept (does not matter for AUC)
  mymodel$fitted.values # predicted probabilities
  myroc=roc(response=Ydata, predictor=mymodel$fitted.values)
  myroc$auc
  ci.auc(myroc)
  
  # Note: the linear predictors are the sum of the intercept and the offset
  a=coef(mymodel)+offset(Xdata[,myorder[1:k],drop=FALSE]%*%matrix(average_beta[myorder[1:k]], ncol=1))
  mymodel$linear.predictors
  plot(a, mymodel$linear.predictors)
  
  
  df.results <- data.frame(symptom = rep(myorder,3), 
                           selection_proportion = rep(myselprops,3),
                           mean_beta = rep(myaverage_beta,3),
                           se_beta = rep(myse_beta,3),
                           auc_lower = NA_real_, 
                           auc= NA_real_, 
                           auc_upper= NA_real_,
                           Data = NA_character_)
  
  ## For loop incrementally adding the predictors
  myaucs=NULL
  for (k in 1:length(myorder)){
    mymodel=glm(Ydata~offset(Xdata[,myorder[1:k],drop=FALSE]%*%matrix(average_beta[myorder[1:k]], ncol=1)), family="binomial")
    # myroc=roc(response=Ydata, predictor=mymodel$fitted.values)
    auc_train <- pROC::auc(response=as.vector(Ydata), predictor=as.vector(mymodel$fitted.values))
    roc_train <- pROC::roc(response=as.vector(Ydata), predictor=as.vector(mymodel$fitted.values))
    auc_ci_train <- pROC::ci.auc(roc_train, conf.level = 0.95, method="bootstrap",
                                 boot.n =1000,boot.stratified = T, reuse.auc=T, progress = "none")
    df.results[k,5:7] <- auc_ci_train
    df.results[k,8] <- "Train"
    if(k==num_include){
      train.preds.list[[i]] <- data.frame(preds=as.vector(mymodel$fitted.values), test_result =as.vector(Ydata))
    }
  }
  
  # plot(myaucs, las=1, pch=19, xlab="Order", ylab="AUC")
  
  
  ### Test set
  
  ## For loop incrementally adding the predictors
  myaucs=NULL
  for (k in 1:length(myorder)){
    mymodel=glm(Ydata_test~offset(Xdata_test[,myorder[1:k],drop=FALSE]%*%matrix(average_beta[myorder[1:k]], ncol=1)), family="binomial")
    # myroc=roc(response=Ydata_test, predictor=mymodel$fitted.values)
    auc_test <- pROC::auc(response=as.vector(Ydata_test), predictor=as.vector(mymodel$fitted.values))
    roc_test<- pROC::roc(response=as.vector(Ydata_test), predictor=as.vector(mymodel$fitted.values))
    auc_ci_test <- pROC::ci.auc(roc_test, conf.level = 0.95, method="bootstrap",
                                boot.n =1000,boot.stratified = T, reuse.auc=T, progress = "none")
    df.results[(k+length(myorder)),5:7] <- auc_ci_test
    df.results[(k+length(myorder)),8] <- "Test"
    
    if(k==num_include){
      test.preds.list[[i]] <- data.frame(preds=as.vector(mymodel$fitted.values), test_result =as.vector(Ydata_test))
    }
  }
  
  
  ### Round 8!!!

  ## For loop incrementally adding the predictors
  myaucs=NULL
  for (k in 1:length(myorder)){
    mymodel=glm(Ydata_r8~offset(Xdata_r8[,myorder[1:k],drop=FALSE]%*%matrix(average_beta[myorder[1:k]], ncol=1)), family="binomial")
    # myroc=roc(response=Ydata_test, predictor=mymodel$fitted.values)
    auc_test <- pROC::auc(response=as.vector(Ydata_r8), predictor=as.vector(mymodel$fitted.values))
    roc_test<- pROC::roc(response=as.vector(Ydata_r8), predictor=as.vector(mymodel$fitted.values))
    auc_ci_test <- pROC::ci.auc(roc_test, conf.level = 0.95, method="bootstrap",
                                boot.n =1000,boot.stratified = T, reuse.auc=T, progress = "none")
    df.results[(k+(2*length(myorder))),5:7] <- auc_ci_test
    df.results[(k+(2*length(myorder))),8] <- "Round_8"
    
    if(k==num_include){
      r8.preds.list[[i]] <- data.frame(preds=as.vector(mymodel$fitted.values), test_result =as.vector(Ydata_r8))
    }
  }
  
  
  saveRDS(df.results, paste0(outpath, "mein_buhl_stability_increment_",agenames[[i]],"_results_auc.rds"))
  pred.results.list[[i]] <- df.results
  
}


names(mydata_train)

# As above but for classic symptoms only logistic model ----------------------------------

varnames_classic <- c("loss_or_change_of_sense_of_smell","loss_or_change_of_sense_of_taste","new_persistent_cough","fever")
pred.results.list.classic  <-list()
test.preds.list.classic <- list()
train.preds.list.classic <- list()
r8.preds.list.classic <- list()

intercept.check <- list()
i <- 1
for (i in 1:3){
  
  
  ### create training data
  Xdata <- mydata_train %>% filter(age %in% ages[[i]]) %>% select(varnames_classic) %>% as.matrix()
  Ydata <- mydata_train %>% filter(age %in% ages[[i]]) %>% select(estbinres)  %>% as.matrix()
  traindata <- mydata_train %>% filter(age %in% ages[[i]])
  Xdata_test <- mydata_test %>% filter( age %in% ages[[i]]) %>% select(varnames_classic) %>% as.matrix()
  Ydata_test <- mydata_test %>% filter( age %in% ages[[i]]) %>% select(estbinres)  %>% as.matrix()
  testdata <- mydata_test %>% filter( age %in% ages[[i]])
  Xdata_r8 <- mydata_r8 %>% filter(age %in% ages[[i]]) %>% select(varnames_classic) %>% as.matrix()
  Ydata_r8 <- mydata_r8 %>% filter( age %in% ages[[i]]) %>% select(estbinres)  %>% as.matrix()
  testdata_r8 <- mydata_r8 %>% filter( age %in% ages[[i]])
  
  
  ### Run unpenalised model
  
  f <- as.formula("estbinres ~ loss_or_change_of_sense_of_smell+loss_or_change_of_sense_of_taste+new_persistent_cough+fever")
  mod.fit <- glm(formula = f, family="binomial", data =  traindata)
  
  ### create resultss df
  df.results <- data.frame(symptom = rep(varnames_classic,3), 
                           selection_proportion = NA_real_,
                           mean_beta = rep(coef(mod.fit)[2:5],3),
                           se_beta = NA_real_,
                           auc_lower = NA_real_, 
                           auc= NA_real_, 
                           auc_upper= NA_real_,
                           Data = NA_character_)
  
  ## For loop incrementally adding the predictors
  myaucs=NULL
  for (k in 1:length(varnames_classic)){
    # mymodel=glm(Ydata~offset(Xdata[,myorder[1:k],drop=FALSE]%*%matrix(average_beta[myorder[1:k]], ncol=1)), family="binomial")
    # myroc=roc(response=Ydata, predictor=mymodel$fitted.values)
    auc_train <- pROC::auc(response=as.vector(Ydata), predictor=predict(object = mod.fit, newdata = traindata))
    roc_train <- pROC::roc(response=as.vector(Ydata), predictor=predict(object = mod.fit, newdata = traindata))
    auc_ci_train <- pROC::ci.auc(roc_train, conf.level = 0.95, method="bootstrap",
                                 boot.n =1000,boot.stratified = T, reuse.auc=T, progress = "none")
    df.results[k,5:7] <- auc_ci_train
    df.results[k,8] <- "Train"
    train.preds.list.classic[[i]] <- data.frame(preds=predict(object = mod.fit, newdata = traindata, type = "response"), test_result =as.vector(Ydata))
    
  }
  
  # plot(myaucs, las=1, pch=19, xlab="Order", ylab="AUC")
  
  
  ### Test set
  
  ## For loop incrementally adding the predictors
  myaucs=NULL
  for (k in 1:length(varnames_classic)){
    auc_test <- pROC::auc(response=as.vector(Ydata_test), predictor=predict(object = mod.fit, newdata = testdata))
    roc_test <- pROC::roc(response=as.vector(Ydata_test), predictor=predict(object = mod.fit, newdata = testdata))
    auc_ci_test <- pROC::ci.auc(roc_test, conf.level = 0.95, method="bootstrap",
                                boot.n =1000,boot.stratified = T, reuse.auc=T, progress = "none")
    df.results[k+length(varnames_classic),5:7] <- auc_ci_test
    df.results[k+length(varnames_classic),8] <- "Test"
    
    test.preds.list.classic[[i]] <- data.frame(preds=predict(object = mod.fit, newdata = testdata, type = "response"), test_result =as.vector(Ydata_test))
  }
  
  
  ### Round 8!!!
  
  ## For loop incrementally adding the predictors
  myaucs=NULL
  for (k in 1:length(varnames_classic)){
    auc_test <- pROC::auc(response=as.vector(Ydata_r8), predictor=predict(object = mod.fit, newdata = testdata_r8))
    roc_test <- pROC::roc(response=as.vector(Ydata_r8), predictor=predict(object = mod.fit, newdata = testdata_r8))
    auc_ci_test <- pROC::ci.auc(roc_test, conf.level = 0.95, method="bootstrap",
                                boot.n =1000,boot.stratified = T, reuse.auc=T, progress = "none")
    df.results[k+(2*length(varnames_classic)),5:7] <- auc_ci_test
    df.results[k+(2*length(varnames_classic)),8] <- "Round_8"
    
    r8.preds.list.classic[[i]] <- data.frame(preds=predict(object = mod.fit, newdata = testdata_r8, type = "response"), test_result =as.vector(Ydata_r8))
  }
  
  
  saveRDS(df.results, paste0(outpath, "mein_buhl_stability_increment_",agenames[[i]],"_results_auc_CLASSIC.rds"))
  pred.results.list.classic[[i]] <- df.results
}



# Reclassification analysis --------------------------------------------------------



### Train data
reclass_dfs <- list()
for (i in 1:3){
  
  preds_classic <- train.preds.list.classic[[i]]
  preds_new <- train.preds.list[[i]]
  preds_comb <- cbind(preds_classic, preds_new = preds_new$preds)
  
  df <- preds_comb %>% 
    mutate(test_result_1 = ifelse(test_result == 1, "Positive", "Negative"),
           Prediction_improve = case_when((preds < preds_new) & test_result == 1 ~ "Improved",
                                          (preds>preds_new) & test_result == 1 ~ "Not improved",
                                          (preds > preds_new) & test_result == 0 ~ "Improved",
                                          (preds < preds_new) & test_result == 0 ~ "Not improved"),
           Prediction_change =case_when((preds < preds_new) ~ "Higher",
                                        (preds>preds_new)  ~ "Lower")
    ) 
  
  df$age <- agenames[[i]]
  df$new_id <- ids$train[[i]]
  reclass_dfs[[i]] <- df
  
}

reclass_df_comb_train <- do.call(rbind, reclass_dfs)


### Test data
reclass_dfs <- list()
for (i in 1:3){
  preds_classic <- test.preds.list.classic[[i]]
  preds_new <- test.preds.list[[i]]
  preds_comb <- cbind(preds_classic, preds_new = preds_new$preds)
  
  df <- preds_comb %>% 
    mutate(test_result_1 = ifelse(test_result == 1, "Positive", "Negative"),
           Prediction_improve = case_when((preds < preds_new) & test_result == 1 ~ "Improved",
                                          (preds>preds_new) & test_result == 1 ~ "Not improved",
                                          (preds > preds_new) & test_result == 0 ~ "Improved",
                                          (preds < preds_new) & test_result == 0 ~ "Not improved"),
           Prediction_change =case_when((preds < preds_new) ~ "Higher",
                                        (preds>preds_new)  ~ "Lower")
    ) 
  
  df$age <- agenames[[i]]
  df$new_id <- ids$test[[i]]
  reclass_dfs[[i]] <- df
  
}

### Bind final file
reclass_df_comb_test <- do.call(rbind, reclass_dfs)



### Test data
reclass_dfs <- list()
for (i in 1:3){
  preds_classic <- r8.preds.list.classic[[i]]
  preds_new <- r8.preds.list[[i]]
  preds_comb <- cbind(preds_classic, preds_new = preds_new$preds)
  
  df <- preds_comb %>% 
    mutate(test_result_1 = ifelse(test_result == 1, "Positive", "Negative"),
           Prediction_improve = case_when((preds < preds_new) & test_result == 1 ~ "Improved",
                                          (preds>preds_new) & test_result == 1 ~ "Not improved",
                                          (preds > preds_new) & test_result == 0 ~ "Improved",
                                          (preds < preds_new) & test_result == 0 ~ "Not improved"),
           Prediction_change =case_when((preds < preds_new) ~ "Higher",
                                        (preds>preds_new)  ~ "Lower")
    ) 
  
  df$age <- agenames[[i]]
  df$new_id <- ids$r8[[i]]
  reclass_dfs[[i]] <- df
  
}

reclass_df_comb_r8 <- do.call(rbind, reclass_dfs)

### join and rename a couple of things
reclass_df_comb_test$Data = "Test"
reclass_df_comb_train$Data = "Train"
reclass_df_comb_r8$Data = "Round_8"
reclass_df_comb_test <- reclass_df_comb_test %>% rename(age_group = age)
reclass_df_comb_train <- reclass_df_comb_train %>% rename(age_group = age)
reclass_df_comb_r8 <- reclass_df_comb_r8 %>% rename(age_group = age)


### Combine all
reclass_df_comb <- do.call(rbind,list(right_join(mydata, reclass_df_comb_train, by = "new_id"),
                         right_join(mydata, reclass_df_comb_test, by = "new_id"),
                         right_join(mydata_r8, reclass_df_comb_r8, by = "new_id")))


reclass_df_comb <- reclass_df_comb %>% mutate(
  tfpf = factor(case_when(test_result == 1 & Prediction_change == "Higher" ~ "TP",
                          test_result == 1 & Prediction_change == "Lower" ~ "FN",
                          test_result == 0 & Prediction_change == "Higher" ~ "FP",
                          test_result == 0 & Prediction_change == "Lower" ~ "TN"), 
                levels = c("TP", "FN", "TN", "FP"))
)


### Add variable for one of four classical symptoms
reclass_df_comb$one_of_four <- as.numeric(rowSums(reclass_df_comb[, varnames_classic]) > 0)




### Save results
save(varnames_classic ,
     reclass_df_comb,
     test.preds.list.classic, 
     train.preds.list.classic, 
     intercept.check,
     ages ,
     agenames ,
     pred.results.list,
     test.preds.list,
     train.preds.list,
     ids, file = paste0(outpath, "m_b_stab_results.rdata"))
load(paste0(outpath, "m_b_stab_results.rdata"))


# Export results lists ----------------------------------------------------



### Add new variable - any of 'our' symptoms

symptoms_1=varnames[c(1,2,11,20,4)]
symptoms_2=varnames[c(1,3,4,26,13,2,20)]
symptoms_3=varnames[c(1,13,3,4,20,2)]

### Count first
reclass_df_comb[reclass_df_comb$age_group == agenames[[1]],"count_of_our_vars"] <- rowSums(reclass_df_comb[reclass_df_comb$age_group == agenames[[1]],symptoms_1])
reclass_df_comb[reclass_df_comb$age_group == agenames[[2]],"count_of_our_vars"] <- rowSums(reclass_df_comb[reclass_df_comb$age_group == agenames[[2]],symptoms_2])
reclass_df_comb[reclass_df_comb$age_group == agenames[[3]],"count_of_our_vars"] <- rowSums(reclass_df_comb[reclass_df_comb$age_group == agenames[[3]],symptoms_3])


### Now binary
reclass_df_comb$any_of_our_vars <- as.numeric(reclass_df_comb$count_of_our_vars > 0)

forBarbara <- reclass_df_comb %>% select(preds,preds_new,test_result_1,Prediction_improve,
                                         Prediction_change,age_group,Data,tfpf,
                                         one_of_four) %>% 
  rename(StabilityModelPreds = preds_new,
         fourSymptomModelPreds = preds 
  )
forBarbara_ext <- reclass_df_comb %>% 
  rename(StabilityModelPreds = preds_new,
         fourSymptomModelPreds = preds 
  )

saveRDS(forBarbara, paste0(outpath, "predictions.rds"))
saveRDS(forBarbara_ext, paste0(outpath, "predictions_detailed.rds"))


# join plots --------------------------------------------------------------


### PLot loop 
i <- 1
abc <- c("A: ","B: ", "C: ")
for (i in 1:3){
  
  ### Plot
  df.results <- pred.results.list[[i]]
  faces <- ifelse(df.results$symptom %in% varnames_classic, "bold", "plain")
  
  dodge=position_dodge(width=0.5)
  p1 <- df.results %>% 
    filter(Data =="Test") %>% 
    mutate(symptom = stringr::str_to_sentence(gsub(pattern =  "_", " ", symptom))) %>% 
    mutate(symptom = factor(symptom, levels = unique(symptom))) %>% 
    ggplot(aes(x=symptom, y= auc, col = selection_proportion > hat_pi_true)) +
    geom_point(position = dodge) +
    ylim(0.5, 0.85) +
    geom_errorbar(aes(ymin = auc_lower, ymax =auc_upper), 
                  width=0.5, size =0.5,position = dodge) +
    scale_color_manual(values = c(myCols[c(5)], "black")) +
    labs(x="Symptoms (ordered by selection proportion)", y = "AUC in test data") +
    theme_bw() +
    geom_hline(yintercept = auc_one_of_four$test[[i]], linetype = "dashed",
               col = myCols[c(7)]) +
    theme(axis.text.x = element_text(angle=70, 
                                     hjust = 1,
                                     face = faces,
                                     color = ifelse(df.results$selection_proportion > hat_pi_true,
                                                    "black", "grey50")),
          legend.position = "none")+
    annotate("text", x = 16, col = myCols[c(7)], hjust =0,
             y = auc_one_of_four$test[[i]] + 0.01,
             label =paste0("AUC:",auc_one_of_four$test[[i]]),
             size = 3)
  p1
  
  p2 <- df.results %>% 
    filter(Data=="Train" ) %>% 
    mutate(symptom = factor(symptom, levels = unique(symptom))) %>% 
    ggplot(aes(x=symptom, y= selection_proportion, col = selection_proportion > hat_pi_true)) +
    geom_col( width=0.01) +
    geom_point() +
    geom_hline(yintercept = hat_pi_true, linetype = "dashed", col = myCols[1]) +
    theme_minimal() +
    scale_color_manual(values = c(myCols[c(5)], "black"))  +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "none"
    ) +
    labs(y = "Selection prop.") 
  p2
  
  p3 <- df.results %>% 
    filter(Data=="Train" ) %>% 
    mutate(symptom = factor(symptom, levels = unique(symptom)),
           mean_OR=exp(mean_beta),
           Sign = ifelse(mean_beta>0, ">1", "<1")) %>% 
    ggplot(aes(x=symptom, y= mean_beta, col = Sign)) +
    geom_point() +
    geom_col( width=0.01) +
    # geom_errorbar(aes(ymin = mean_beta - 1.96*se_beta, ymax =mean_beta + 1.96*se_beta), 
    #               width=0.2, size =0.5) +
    scale_colour_manual(values = myCols[c(7,4)]) +
    theme_minimal() +
    ylim(-0.3,2.5) +
    geom_hline(yintercept = 0, linetype = "dashed")+
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "none"
    ) +
    labs(y = "Mean Log OR", fill = "OR") + ggtitle(paste0(abc[[i]], agelabels[[i]]))
  
  p3
  
  
  plot.comb <- p3/p2/p1 + plot_layout(heights = c(1,1,3))
  plot.comb
  plots.list[[i]] <- plot.comb
  ggsave(filename = paste0(figpath, "mein_buhl_stability_increment_",agenames[[i]],".png"), 
         plot = plot.comb,width = 9, height=10, dpi = 300, units = "in")
  
  
  
}

### Use patchwork to create one huge plot
megaplot <- wrap_plots(plots.list[[1]] , 
           plots.list[[2]],
           plots.list[[3]],
           ncol =3, 
           widths = c(26,23,24)) + 
  plot_layout(guides = "collect")

megaplot

ggsave(filename = paste0(figpath, "mein_buhl_stability_combined_plot.png"), 
       plot = megaplot,width = 17, height=10, dpi = 300, units = "in")

ggsave(filename = paste0(figpath, "mein_buhl_stability_combined_plot.pdf"), 
       plot = megaplot,width = 17, height=10, dpi = 300, units = "in")




# Output mean betas -------------------------------------------------------

df_comb <- bind_rows(pred.results.list, .id = "id")
df_comb$beta_lower <- df_comb$mean_beta - 1.96*df_comb$se_beta
df_comb$beta_upper <- df_comb$mean_beta + 1.96*df_comb$se_beta

write.csv(df_comb, paste0(outpath, "stab_sel_results.csv"))



