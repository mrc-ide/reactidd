
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
source("E:/Group/react2_study5/report_phases_combined/projects/function_scripts/stability_selection.R", local = T)

# 
# 
# source("E:/Group/react2_study5/report_phases_combined/projects/delta_symptom_prediction/code/00_bits_and_pieces.R", local = T)
# source("E:/Group/react2_study5/report_phases_combined/projects/delta_symptom_prediction/code/00_functions.R", local = T)
# source("E:/Group/react2_study5/report_phases_combined/projects/symptom_prediction_children/code/00_bits_and_pieces.R", local = T)


#' Pull in packages needed
package.list <- c("prevalence","mgcv","knitr","MASS","kableExtra","table1","dplyr",
                  "tidyr", "pheatmap",
                  "ggplot2","ggsci", "RColorBrewer", "tidyverse", "lubridate", 
                  "readr","ggthemes", "questionr", "foreach", "doParallel","withr",
                  "patchwork","randomcoloR","focus","ComplexHeatmap"
)

load_packages(package.list)


# Import REACT-1 data -----------------------------------------------------

source("E:/Group/react2_study5/report_phases_combined/projects/omicron_symptom_profiling//code/00_bits_and_pieces.R", local = T)
source("E:/Group/react2_study5/report_phases_combined/projects/omicron_symptom_profiling//code/00_data_prep.R", local = T)

# create subfolder
createMySubfolder(subfolderName = "stability_selection")



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
    lower.limits=NULL
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



# Some tinkering ----------------------------------------------------------

# numeric sex
dfRes$sex = dfRes$gender - 1 
# scaled age
dfRes$age_scaled = (dfRes$age -mean(dfRes$age,na.rm=T))/sd(dfRes$age,na.rm=T)



# Omicron BA1 only  -----------------------------------------------------------




# all positives in r17 (implicitly all omicron) - adjusted on vaccine status
res_omicron_BA1=doEverything(mydat_for_analysis = dfRes %>% 
                               filter(round %in% c(17:19) & 
                                        (variant_inferred_detail ==  "BA.1 (Omicron)" | estbinres==0) # filter for BA1/negative
                                      ),
                             yname = "estbinres",
                             analysis_name = "BA1_omicron_r17_19_pos_v_neg_age_sex_vax_adjusted", 
                             myK = 1000,
                             myseed = 123, 
                             subfolder = "BA1_omicron",
                             positive_coefs_only = T,
                             adjustments = c("age_scaled","sex", "vax_status_number"),
                             pi_list=seq(0.92,0.998,0.002)
)

res_omicron_BA1$plot
res_omicron_BA1$results



# Omicron BA2 only  -----------------------------------------------------------

# BA2
res_omicron_BA2=doEverything(mydat_for_analysis = dfRes %>% 
                               filter(round %in% c(17:19) & 
                                        (variant_inferred_detail ==  "BA.2 (Omicron)" | estbinres==0) # filter for BA1/negative
                               ),
                             yname = "estbinres",
                             analysis_name = "BA2_omicron_r17_19_pos_v_neg_age_sex_vax_adjusted", 
                             myK = 1000,
                             myseed = 123, 
                             subfolder = "BA2_omicron",
                             positive_coefs_only = T,
                             adjustments = c("age_scaled","sex", "vax_status_number"),
                             pi_list=seq(0.92,0.998,0.002)
)

res_omicron_BA2$plot
res_omicron_BA2$results









# Delta only  -----------------------------------------------------------


# all positives in r17 (implicitly all delta) - adjusted on vaccine status
res_delta_r13_15=doEverything(mydat_for_analysis = dfRes %>% 
                               mutate(doublevaxxed=case_when(vax_status_number==2 ~1,
                                                             TRUE ~0),
                                      boosted=case_when(vax_status_number==3 ~1,
                                                        TRUE ~0)
                               ) %>% 
                               filter(round %in% c(13:15)),
                             yname = "estbinres",
                             analysis_name = "r13_15_pos_v_neg_age_sex_vax_adjusted", 
                             myK = 1000,myseed = 122, 
                             subfolder = "delta",
                             positive_coefs_only = T,
                             adjustments = c("age_scaled","sex", "vax_status_number"),
                             pi_list=seq(0.92,0.998,0.002)
)


res_delta_r13_15$plot
res_delta_r13_15$results




# Alpha only  -----------------------------------------------------------


# all positives in r17 (implicitly all delta) - adjusted on vaccine status
res_alpha_r8_10=doEverything(mydat_for_analysis = dfRes %>% 
                                mutate(doublevaxxed_or_boosted=case_when(vax_status_number>=2 ~1,
                                                                         TRUE ~0)
                                ) %>% 
                                filter(round %in% c(8:10)),
                              yname = "estbinres",
                              analysis_name = "r8_10_pos_v_neg_age_sex_vax_adjusted", 
                              myK = 1000,
                             myseed = 123, 
                              subfolder = "alpha",
                             positive_coefs_only = T,
                              adjustments = c("age_scaled","sex", "vax_status_number"),
                             pi_list=seq(0.92,0.998,0.002)
)

# 
# # all positives in r17 (implicitly all delta) - adjusted on vaccine status
# res_alpha_r8_10_novax=doEverything(mydat_for_analysis = dfRes %>% 
#                                filter(round %in% c(8:10)),
#                              yname = "estbinres",
#                              analysis_name = "r8_10_pos_v_neg_age_sex_adjusted", 
#                              myK = 1000,
#                              myseed = 123, 
#                              subfolder = "alpha",
#                              positive_coefs_only = T,
#                              adjustments = c("age_scaled","sex"))
# 
# res_alpha_r8_10_novax$results
# res_alpha_r8_10$results



# Wildtype only  -----------------------------------------------------------


# all positives in r17 (implicitly all delta) - adjusted on vaccine status
res_alpha_r2_7=doEverything(mydat_for_analysis = dfRes %>% 
                               mutate(doublevaxxed_or_boosted=case_when(vax_status_number>=2 ~1,
                                                                        TRUE ~0)
                               ) %>% 
                               filter(round %in% c(2:7)),
                             yname = "estbinres",
                             analysis_name = "r2_7_pos_v_neg_age_sex_adjusted", 
                             myK = 1000,
                             myseed = 123, 
                            positive_coefs_only = T,
                             subfolder = "wildtype",
                             adjustments = c("age_scaled","sex"),
                            pi_list=seq(0.92,0.998,0.002)
)


res_alpha_r2_7$plot
res_alpha_r2_7$results







# Save all results --------------------------------------------------------
wild_sx=res_alpha_r2_7$results$variable[res_alpha_r2_7$results$selected]
wild_sx=wild_sx[!grepl(" \\[",wild_sx,ignore.case = T)]
alpha_sx=res_alpha_r8_10$results$variable[res_alpha_r8_10$results$selected]
alpha_sx=alpha_sx[!grepl(" \\[",alpha_sx,ignore.case = T)]
delta_sx=res_delta_r13_15$results$variable[res_delta_r13_15$results$selected]
delta_sx=delta_sx[!grepl(" \\[",delta_sx,ignore.case = T)]
ba1_sx = res_omicron_BA1$results$variable[res_omicron_BA1$results$selected]
ba1_sx=ba1_sx[!grepl(" \\[",ba1_sx,ignore.case = T)]
ba2_sx = res_omicron_BA2$results$variable[res_omicron_BA2$results$selected]
ba2_sx=ba2_sx[!grepl(" \\[",ba2_sx,ignore.case = T)]




### save all names and codes
selected_sx = list(codes=list(wildtype=sympnames_type_df$symptom_code[sympnames_type_df$symptom %in% wild_sx],
                              alpha=sympnames_type_df$symptom_code[sympnames_type_df$symptom %in% alpha_sx],
                              delta=sympnames_type_df$symptom_code[sympnames_type_df$symptom %in% delta_sx],
                              BA1_omicron=sympnames_type_df$symptom_code[sympnames_type_df$symptom %in% ba1_sx],
                              BA2_omicron=sympnames_type_df$symptom_code[sympnames_type_df$symptom %in% ba2_sx]),
                   sx=list(wildtype=wild_sx,alpha=alpha_sx,delta=delta_sx,ba1=ba1_sx,ba2=ba2_sx))

## output
OverReact::saveREACT(file = selected_sx,filename = "selected_sx",outpath = outpath)

# save full results tables
savePrettyExcelWorkbook(listOfTables = list(wildtype=res_alpha_r2_7$results,
                                            alpha=res_alpha_r8_10$results,
                                            delta=res_delta_r13_15$results,
                                            ba1=res_omicron_BA1$results,
                                            ba2=res_omicron_BA2$results),
                        workbookName = "stability_selection_summary_data",
                        outpath = outpath)



# Simple plot -------------------------------------------------------------

# Add variant types
res_alpha_r2_7$results$Variant="Wild-type"
res_alpha_r8_10$results$Variant="Alpha"
res_delta_r13_15$results$Variant="Delta"
res_omicron_BA1$results$Variant="BA.1"
res_omicron_BA2$results$Variant="BA.2"

# bind all
res_all=rbind(res_alpha_r2_7$results,res_alpha_r8_10$results,res_delta_r13_15$results,
              res_omicron_BA1$results, res_omicron_BA2$results) %>% 
  filter(!grepl("unpenalised",variable)) %>% 
  mutate(Variant=factor(Variant, levels = c("Wild-type","Alpha","Delta","BA.1","BA.2")))

#  plot
pal="bnt"
p_selprop=res_all %>% ggplot(aes(x=100*selprop, y=reorder(variable,selprop), col =selected,fill =selected)) +
  geom_point(size=1.5) +
  geom_col(width=0.03)+
  facet_wrap(.~Variant, ncol=5) +
  OverReact::theme_react() +
  scale_fill_manual(values=c("#9D9D9D","#D24000")) +
  scale_color_manual(values=c("#9D9D9D","#D24000"))+
  theme(legend.position = "none") +
  labs(x="Selection proportion (% of 1000 models)", y="")
p_selprop
OverReact::imperial_palettes
# save
OverReact::saveREACTplot(p = p_selprop,figpath = figpath,
                         filename = "selection_proportions_panel",
                         width = 9,height = 5,savePDF = F)










# Create Venn plot --------------------------------------------------------
adj <- c("age_scaled [unpenalised]","sex [unpenalised]","vax_status_number [unpenalised]" )
# get selected variables from each wave
selvars_omicron=setdiff(res_omicron_r17_19$results$variable[res_omicron_r17_19$results$selected],adj)
selvars_delta=setdiff(res_delta_r13_15$results$variable[res_delta_r13_15$results$selected],adj)
selvars_alpha=setdiff(res_alpha_r8_10$results$variable[res_alpha_r8_10$results$selected],adj)
selvars_wildtype=setdiff(res_alpha_r2_7$results$variable[res_alpha_r2_7$results$selected],adj)

# create Venn diagram
venn_plot<- venn::venn(list("Wild type" = selvars_wildtype,"Alpha" = selvars_alpha,
                            "Delta" = selvars_delta,"Omicron" = selvars_omicron
                                      ),
                   ggplot=T, zcolor = myCols[c(1,3,6,7)])





# Cross evaluation --------------------------------------------------------


# sub function to evaluate performance on other variants
evalPerformance <- function(stab_out,
                            n_thr = 10,
                            myK=1000,
                            myseed=123,
                            yname = "estbinres",
                            positive_coefs_only = T,
                            pi_list=seq(0.8,0.99,0.01),
                            reorder_by_normalised_beta = F,
                            omirounds=c(17)){
  
  set.seed(myseed)
  
  
  # create list of all rounds
  roundlist=list(wildtype=c(2:7),
                 alpha=c(8:10),
                 delta=c(13:15),
                 omicron=omirounds)
  adjustlist=list(wildtype=c("age_scaled","sex"),
                  alpha=c("age_scaled","sex", "vax_status_number"),
                  delta=c("age_scaled","sex", "vax_status_number"),
                  omicron=c("age_scaled","sex", "vax_status_number"))
  roundlist_res=list()
  
  
  # selection proportions from trained model
  selprop=focus::SelectionProportions(stab_out)  
  opt_iter=which.max(stab_out$S)
  opt_thresh=stab_out$P[opt_iter]
  selvars=names(selprop)[selprop>=opt_thresh]
  
  
  i=1
  for(i in 1:4){
    print(paste0("Evaluating on ",names(roundlist)[[i]]))
    # create data sets
    adjustments=adjustlist[[i]]
    mydat_for_analysis=dfRes %>% filter(round %in% roundlist[[i]]) %>% 
      filter(test_train_split == "test")%>% 
      select(all_of(yname),adjustments,covid_yesnos)
    predictors=c(adjustments,covid_yesnos)
    mydat_for_analysis <- mydat_for_analysis[complete.cases(mydat_for_analysis[,predictors]),]
    # mydat_test = mydat_for_analysis %>% filter(test_train_split =="test")
    ydata=pull(mydat_for_analysis,yname)
    X_test=mydat_for_analysis %>% select(all_of(predictors))
    colnames(X_test)<- c(adjustments,react1_sympnames[1:26])
    xdata=X_test[,c(adjustments,selvars)]
    lower.limits =c(rep(-Inf,length(adjustments)),rep(0,length(selvars)))
    
    perf <- ExplanatoryPerformanceMW(xdata = xdata[,,drop=FALSE], ydata = ydata, stability = NULL, 
                                     family = "binomial", 
                                     K = 100, tau = 0.8, seed = 1, 
                                     n_thr = NULL, ij_method = FALSE, 
                                     time = 1000, lower.limits = lower.limits)
    
    
    roundlist_res[[i]] <- perf$AUC
  }
  return(roundlist_res)
}






# Evaluate on other variants
perf_omi=evalPerformance(stab_out = res_omicron_r17_19$raw_results,
                         myK = 100,myseed = 123,
                         positive_coefs_only = T,
                         n_thr = 10,yname = "estbinres",
                         reorder_by_normalised_beta = F)



# Evaluate on other variants
perf_delta=evalPerformance(stab_out = res_delta_r13_15$raw_results,
                           myK = 1000,myseed = 123,
                           positive_coefs_only = T,
                           n_thr = 10,yname = "estbinres",
                           reorder_by_normalised_beta = F)


# Evaluate on other variants
perf_alpha=evalPerformance(stab_out = res_alpha_r8_10$raw_results,
                           myK = 1000,myseed = 123,
                           positive_coefs_only = T,
                           n_thr = 10,yname = "estbinres",
                           reorder_by_normalised_beta = F)
perf_alpha

# Evaluate on other variants
perf_wild=evalPerformance(stab_out = res_alpha_r2_7$raw_results,
                          myK = 1000,myseed = 123,
                          positive_coefs_only = T,
                          n_thr = 10,yname = "estbinres",
                          reorder_by_normalised_beta = F)

# get matrix of final x-performance
perf_mat = matrix(data = NA,nrow = 4,ncol = 4,
                  dimnames = list(c("Wildtype","Alpha","Delta","Omicron"),
                                  c("Wildtype","Alpha","Delta","Omicron")))


# create list of all results
pers_all=list(perf_wild,perf_alpha,perf_delta,perf_omi)
for(i in 1:4){
  perf=pers_all[[i]]
  for(j in 1:4){
    perf_mat[i,j] <- round(mean(perf[[j]]),3)
  }
}

# rows = trained on
# columns = evaluated on

### PLot heatmap
col_fun=circlize::colorRamp2(breaks = c(0.65,0.75,0.85), colors = c("white",myCols[[2]],myCols[[1]]))
auc_heatmap=ComplexHeatmap::Heatmap(name = "AUC",
                        perf_mat,
                        cluster_rows = F,cluster_columns = F,
                        cell_fun = function(j,i,x,y,width, height, fill){
                          grid.text(sprintf("%.3f", perf_mat[i,j]), 
                                    x, y, gp=gpar(fontsize=10, col = 
                                                    ifelse(perf_mat[i,j]>0.7,"white","black")))
                        },
                        row_title = "Variables selected on",
                        row_title_gp = gpar(fontsize=16, fontface="bold"),
                        column_title_gp = gpar(fontsize=16, fontface="bold"),
                        row_names_side = "left",
                        column_names_side = "top",
                        
                        column_title = "Models evaluated on",
                        column_gap = unit(5,"mm"),
                        rect_gp = gpar(col = "white", lwd = 1),
                        col = col_fun)

auc_heatmap

png(filename = paste0(figpath,"cross_evaluation_auc_heatmap_omi_17_18.png"),
    width = 5,height = 4.3,units = "in",res = 300)
auc_heatmap
dev.off()



# With rounds 17/18 for Omicron -------------------------------------------






# Evaluate on other variants
perf_2_omi=evalPerformance(stab_out = res_omicron_r17_18$raw_results,
                         myK = 100,myseed = 123,
                         positive_coefs_only = T,
                         n_thr = 10,yname = "estbinres",
                         reorder_by_normalised_beta = F,omirounds=c(17,18))



# Evaluate on other variants
perf_2_delta=evalPerformance(stab_out = res_delta_r13_15$raw_results,
                           myK = 1000,myseed = 123,
                           positive_coefs_only = T,
                           n_thr = 10,yname = "estbinres",
                           reorder_by_normalised_beta = F,omirounds=c(17,18))


# Evaluate on other variants
perf_2_alpha=evalPerformance(stab_out = res_alpha_r8_10$raw_results,
                           myK = 1000,myseed = 123,
                           positive_coefs_only = T,
                           n_thr = 10,yname = "estbinres",
                           reorder_by_normalised_beta = F,omirounds=c(17,18))
perf_2_alpha

# Evaluate on other variants
perf_2_wild=evalPerformance(stab_out = res_alpha_r2_7$raw_results,
                          myK = 1000,myseed = 123,
                          positive_coefs_only = T,
                          n_thr = 10,yname = "estbinres",
                          reorder_by_normalised_beta = F,omirounds=c(17,18))

# get matrix of final x-performance
perf_2_mat = matrix(data = NA,nrow = 4,ncol = 4,
                  dimnames = list(c("Wildtype","Alpha","Delta","Omicron"),
                                  c("Wildtype","Alpha","Delta","Omicron")))


# create list of all results
pers_all=list(perf_2_wild,perf_2_alpha,perf_2_delta,perf_2_omi)
for(i in 1:4){
  perf=pers_all[[i]]
  for(j in 1:4){
    perf_2_mat[i,j] <- round(mean(perf[[j]]),3)
  }
}

# rows = trained on
# columns = evaluated on

### PLot heatmap
col_fun=circlize::colorRamp2(breaks = c(0.65,0.75,0.85), colors = c("white",myCols[[2]],myCols[[1]]))
auc_heatmap=ComplexHeatmap::Heatmap(name = "AUC",
                                    perf_2_mat,
                                    cluster_rows = F,cluster_columns = F,
                                    cell_fun = function(j,i,x,y,width, height, fill){
                                      grid.text(sprintf("%.3f", perf_2_mat[i,j]), 
                                                x, y, gp=gpar(fontsize=10, col = 
                                                                ifelse(perf_2_mat[i,j]>0.7,"white","black")))
                                    },
                                    row_title = "Variables selected on",
                                    row_title_gp = gpar(fontsize=16, fontface="bold"),
                                    column_title_gp = gpar(fontsize=16, fontface="bold"),
                                    row_names_side = "left",
                                    column_names_side = "top",
                                    
                                    column_title = "Models evaluated on",
                                    column_gap = unit(5,"mm"),
                                    rect_gp = gpar(col = "white", lwd = 1),
                                    col = col_fun)

auc_heatmap

png(filename = paste0(figpath,"cross_evaluation_auc_heatmap_omi_17_18.png"),
    width = 5,height = 4.3,units = "in",res = 300)
auc_heatmap
dev.off()






