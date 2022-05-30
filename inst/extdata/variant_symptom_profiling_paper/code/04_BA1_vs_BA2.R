
#' First clear the environment of variables
rm(list=ls(all=TRUE))

# get root director of project
root.dir <-"E:/Group/react2_study5/report_phases_combined/projects/omicron_symptom_profiling/"


outpath <- paste0(root.dir,"/output/")
figpath <-  paste0(root.dir,"/plots/")



source("E:/Group/functions/load_packages.R", local = T)
source("E:/Group/functions/full_join_multiple_datasets.R", local = T)
source("E:/Group/functions/wrangle_cleaner_functions.R", local = T)
source("E:/Group/functions/cats_and_covs.R", local = T)
source("E:/Group/react2_study5/report_phases_combined/projects/function_scripts/create_subfolder.R", local = T)
source("E:/Group/react2_study5/report_phases_combined/projects/function_scripts/forest_plot.R", local = T)
source("E:/Group/react2_study5/report_phases_combined/projects/function_scripts/save_styled_table.R", local = T)
source("E:/Group/react2_study5/report_phases_combined/projects/function_scripts/catboost_and_shap.R", local = T)
source("E:/Group/react2_study5/report_phases_combined/projects/function_scripts/stability_selection.R", local = T)

# 
# 
# source("E:/Group/react2_study5/report_phases_combined/projects/delta_symptom_prediction/code/00_bits_and_pieces.R", local = T)
# source("E:/Group/react2_study5/report_phases_combined/projects/delta_symptom_prediction/code/00_functions.R", local = T)
# source("E:/Group/react2_study5/report_phases_combined/projects/symptom_prediction_children/code/00_bits_and_pieces.R", local = T)


#' Pull in packages needed
package.list <- c("prevalence","mgcv","knitr","MASS","kableExtra","table1","dplyr",
                  "tidyr", "pheatmap","catboost","datapasta",
                  "ggplot2","ggsci", "RColorBrewer", "tidyverse", "lubridate", 
                  "readr","ggthemes", "questionr", "foreach", "doParallel","withr",
                  "patchwork","randomcoloR","focus","ComplexHeatmap"
)

load_packages(package.list)


# Import REACT-1 data -----------------------------------------------------
source("E:/Group/react2_study5/report_phases_combined/projects/omicron_symptom_profiling/code/00_bits_and_pieces.R", local = T)
source("E:/Group/react2_study5/report_phases_combined/projects/omicron_symptom_profiling/code/00_data_prep.R", local = T)

# create subfolder
createMySubfolder(subfolderName = "bA1_bA2")

# check lineage categories
table(dfRes$variant_inferred_detail)

# clean up key symptom variables
dfRes <- dfRes %>% 
  mutate_at(covid_yesnos_firstsymp,binaryCleaner_1_0) %>% 
  mutate_at(covid_yesnos_month,binaryCleaner_1_0)  %>% 
  mutate(doublevaxxed_or_boosted=case_when(vax_status_number>=2 ~1,
                                           TRUE ~0),
         sex=factor(sex)) %>% 
  mutate_at(c(sympnames_type_df$symptom_code,covid_yesnos_firstsymp,covid_yesnos_month,
              "sex","doublevaxxed_or_boosted"),as.factor) 


# Filter data to typed BA1/BA2
dfRes_18 <- dfRes %>% filter((variant_inferred_detail%in% c("BA.1 (Omicron)","BA.2 (Omicron)") & estbinres==1) | 
                               (estbinres==0 & round %in% c(17:19)))

# Check
dfRes_18$variant_inferred_detail %>% table(dfRes_18$estbinres,exclude="none")  #all good



# Univariate --------------------------------------------------------------

# Define function
runUnivariate <- function(mydat, title="BA.2 vs BA.1",rounds =c(17:19), adj_level =5,
                          outcome = "estbinres",
                          subtitle = "Odds ratios from logistic regression, adjusted for age, sex, round and vaccination status"){
  
  # OMICRON #
  joint_adjustment_vars = c("age_group_named", "sex", "vax_status_number","round")
  uiv_omicron =OverReact::ModelMakerMulti(dat = mydat %>% filter(round %in% rounds),
                                          list_of_variables_of_interest = sympnames_type_df$symptom_code,
                                          outcome = outcome,
                                          cov_name_list = cov_name_list,
                                          joint_adjustment_vars = joint_adjustment_vars)
  # extract plottable dat
  plot_df <- uiv_omicron$plot_output
  plot_df$Category <- plot_df$Variable
  plot_df$Variable <- "Symptom in past week"
  
  
  
  # Get summary df with average vals by group
  mean_ors <- plot_df %>% group_by(Category) %>% 
    summarise(mean_or=mean(OR, na.rm=T)) %>% 
    arrange(-mean_or) %>% 
    ungroup()

  
  # join with original data
  plot_df_withmean <- plot_df %>% left_join(mean_ors) %>% filter(adjustment==5) %>% 
    arrange(-OR) %>% filter(!is.na(OR))
  
  
  bold_vect=rev(case_when(plot_df_withmean$Category %in% sympnames_type_df$symptom[c(27)] ~ "bold",
                          T ~ "plain"))
  
  # plot_df_withmean$variant <- variant
  p_compare <- plotReactForest(univ_df_plot = plot_df_withmean,adjustment_numbers = adj_level,
                               adjustment_descriptions = c("Age, sex, vaccination status adjusted"),
                               insignificant_results_greyed_out = T,
                               legend.position = "bottom", 
                               palette = "cool")
  p_compare <- p_compare + labs(title = title, 
                                subtitle = subtitle) +
    theme(legend.position = "none",
          axis.text.y = element_text(face = bold_vect))
  
  return(list(plot=p_compare,
              results=plot_df_withmean))
  
}


age_levs=c("35-44","18-24","25-34" ,"45-54" , "55-64","65-74","74+" )


### BA2 vs BA1
univ_ba2_ba1 <- runUnivariate(mydat = dfRes_18 %>% mutate(ba2_ba1=case_when(variant_inferred_detail== "BA.2 (Omicron)" ~ 1,
                                                                            T ~ 0),
                                                          age_group_named=factor(age_group_named,levels = age_levs)) %>% 
                                filter(estbinres==1), outcome = "ba2_ba1",
                              subtitle = NULL)

# save plot
OverReact::saveREACTplot(p = univ_ba2_ba1$plot,figpath = figpath,filename = "univ_ba2_ba1",width = 4,height = 6,savePDF = F)
OverReact::saveREACTtable(tab = univ_ba2_ba1$results,outpath = outpath,
                          filename = "univ_ba2_ba1", save_rds = T)



### Delta vs BA1
univ_ba1_delta <- runUnivariate(mydat = dfRes %>% mutate(ba1_delta=case_when(variant_inferred_detail== "BA.1 (Omicron)" ~ 1,
                                                                           grepl("B.1.617.2",react_lineage) ~ 0,
                                                                           grepl("AY",react_lineage) ~ 0,
                                                                            T ~ NA_real_)) %>% 
                                filter(estbinres==1, round %in% c(15:16), !is.na(ba1_delta)), 
                              outcome = "ba1_delta",title = "BA.1 vs Delta",rounds =c(15:16),
                              subtitle = NULL)
univ_ba1_delta$plot

# save plot
OverReact::saveREACTplot(p = univ_ba1_delta$plot,figpath = figpath,filename = "univ_ba1_delta",width = 4,height = 6,savePDF = F)


### Save workbook with output
savePrettyExcelWorkbook(listOfTables = list(univ_ba1_delta=univ_ba1_delta$results,
                                            univ_ba2_ba1=univ_ba2_ba1$results), workbookName = "univariate_delta_ba1_ba2",
                        outpath = outpath)



















# Stability selection -----------------------------------------------------



# Function that does everything
doEverything = function(mydat_for_analysis,
                        yname="alpha_delta",
                        analysis_name="alpha_vs_delta",
                        myK=1000,
                        myseed=123,
                        subfolder=NULL,
                        adjustments=NULL, 
                        positive_coefs_only = F,
                        pi_list=seq(0.6,0.95,0.01),
                        ...
){
  
  
  
  
  
  # create folder
  if(!is.null(subfolder)){
    dir.create(paste0(outpath,subfolder,"/"),showWarnings = F)
    dir.create(paste0(figpath,subfolder,"/"),showWarnings = F)
  }
  
  if(is.null(pi_list)){
    pi_list=seq(0.6,0.9,0.01)
  }
  
  
  # Wrangle data ------------------------------------------------------------
  
  set.seed(myseed)
  # filter to complete cases
  predictors=c(adjustments,covid_yesnos)
  mydat_for_analysis <- mydat_for_analysis[complete.cases(mydat_for_analysis[,predictors]),]
  
  mydat_train = mydat_for_analysis %>% filter(test_train_split =="train")
  mydat_test = mydat_for_analysis %>% filter(test_train_split == "test")
  penalty.factor=c(rep(0,length(adjustments)),rep(1,length(covid_yesnos)))
  if(positive_coefs_only){
    lower.limits =c(rep(-Inf,length(adjustments)),rep(0,length(covid_yesnos)))
  }else{
    lower.limits=c(rep(-Inf,length(adjustments)),rep(-Inf,length(covid_yesnos)))
  }
  
  X = mydat_train[,predictors]
  y = pull(mydat_train,yname)
  X_test=mydat_test[,predictors]
  colnames(X) <- colnames(X_test)<- c(paste0(adjustments," [unpenalised]"),react1_sympnames[1:26])
  y_test = pull(mydat_test,yname)
  
  # Run analysis ------------------------------------------------------------
  
  stab_out <- focus::VariableSelection(xdata = X, 
                                       ydata = y,
                                       penalty.factor = penalty.factor,
                                       pi_list = pi_list,
                                       K = myK,
                                       seed = myseed,
                                       family = "binomial",
                                       PFER_method = "MB",
                                       lower.limits = lower.limits,
                                       # n_cores = 30,  
                                       verbose = T,
                                       resampling = "bootstrap",
                                       output_data = T
  )
  
  # Save results #
  OverReact::saveREACT(stab_out,outpath = paste0(outpath,subfolder,"/"),filename = paste0(analysis_name, "_stab_results"))
  
  # Calibrationplot #
  png(paste0(figpath,subfolder,"/",analysis_name,"_calbration_plot.png"),width = 7,height = 7, units="in", res = 300)
  par(mar=c(7,5,7,6))
  focus::CalibrationPlot(stab_out)
  dev.off()
  
  
  
  results_df <- getIncrementalSummary(stab_out = stab_out,
                                      ydata = y_test,
                                      xdata = X_test,
                                      K = 100,
                                      n_thr = min(10,ncol(X_test)),
                                      family = "binomial",
                                      reorder_by_normalised_beta = F)
  opt_iter=which.max(stab_out$S)
  
  
  
  
  
  
  
  myplot=plotStabResultsStripedFlipped(results_df, opt_thresh=stab_out$P[opt_iter], 
                                       stab_out = stab_out,
                                       plotOnlyTopN = 20,
                                       reorder_by_normalised_beta = F)
  myplot
  
  
  ### PLOT WITHOUT AUCs
  myplot_2panel=plotStabResultsStripedFlipped(results_df, opt_thresh=stab_out$P[opt_iter], 
                                              stab_out = stab_out,plotOnlyTopN = 15,
                                              reorder_by_normalised_beta = F,plot_aucs = F)
  myplot_2panel
  
  OverReact::saveREACTplot(p = myplot,figpath = paste0(figpath,subfolder,"/"),
                           filename = paste0(analysis_name, "_stability_selection_plot"),
                           width = 10,height = 3+nrow(results_df)*0.08,savePDF = F)
  
  OverReact::saveREACTplot(p = myplot_2panel,figpath = paste0(figpath,subfolder,"/"),
                           filename = paste0(analysis_name, "_stability_selection_plot_two_Panel"),
                           width = 7,height = 3+nrow(results_df)*0.08,savePDF = F)
  
  
  return(list(results=results_df,
              raw_results = stab_out,
              plot=myplot))
  
}

# numeric sex
dfRes_18$sex = dfRes_18$gender - 1 
# scaled age
dfRes_18$age_scaled = (dfRes_18$age -mean(dfRes_18$age,na.rm=T))/sd(dfRes_18$age,na.rm=T)
# Numeric round
dfRes_18$round_num = as.numeric(as.character(dfRes_18$round))
dfRes$round_num = as.numeric(as.character(dfRes$round))


# all positives in r17 (implicitly all omicron) - adjusted on vaccine status
res_round_ba1_ba2=doEverything(mydat_for_analysis = 
                               dfRes_18 %>% mutate(ba2_ba1=case_when(variant_inferred_detail== "BA.2 (Omicron)" ~ 1,
                                                                     T ~ 0)) %>% 
                               filter(estbinres==1), 
                             yname = "ba2_ba1",
                             analysis_name = "r17_19_BA1_age_sex_vax_round", 
                             myK = 1000,
                             myseed = 122, 
                             subfolder = "omicron",
                             positive_coefs_only = T,
                             adjustments = c("age_scaled","sex", "vax_status_number","round_num"),
                             pi_list=seq(0.6,0.99,0.01)
)

res_round_ba1_ba2$plot
res_round_ba1_ba2$results


# all positives in r17 (implicitly all omicron) - adjusted on vaccine status
res_ba1_ba2=doEverything(mydat_for_analysis = 
                           dfRes_18 %>% mutate(ba2_ba1=case_when(variant_inferred_detail== "BA.2 (Omicron)" ~ 1,
                                                                 T ~ 0)) %>% 
                           filter(estbinres==1), 
                         yname = "ba2_ba1",
                         analysis_name = "r17_19_BA1_age_sex_vax", 
                         myK = 1000,
                         myseed = 122, 
                         subfolder = "omicron",
                         positive_coefs_only = T,
                         adjustments = c("age_scaled","sex", "vax_status_number"),
                         pi_list=seq(0.6,0.99,0.01)
)

res_ba1_ba2$plot
res_ba1_ba2$results
