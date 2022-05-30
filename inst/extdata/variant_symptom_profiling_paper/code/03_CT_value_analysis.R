
#' First clear the environment of variables
rm(list=ls(all=TRUE))
# get root director of project
root.dir <- getwd()
# setwd(dir = "/Group/react2_study5/report_phases_combined/projects/omicron_symptom_profiling/")
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
                  "tidyr", "pheatmap","scales","survival","Epi",
                  "ggplot2","ggsci", "RColorBrewer", "tidyverse", "lubridate", 
                  "readr","ggthemes", "questionr", "foreach", "doParallel","withr",
                  "patchwork","randomcoloR","focus"
)

load_packages(package.list)


# Import REACT-1 data -----------------------------------------------------

source("E:/Group/react2_study5/report_phases_combined/projects/omicron_symptom_profiling/code/00_bits_and_pieces.R", local = T)
source("E:/Group/react2_study5/report_phases_combined/projects/omicron_symptom_profiling//code/00_data_prep.R", local = T)

# create subfolder
createMySubfolder(subfolderName = "ct_value_analysis")


# Distribution plot -------------------------------------------------------


### Set up comparisons for ggstats
waves= dfRes$variant_inferred_detail%>% unique() %>% as.character()
waves <- waves[!is.na(waves)][4:5]
mycomparisons <- list(c(waves[[1]],waves[[2]])
                      )

# Boxplot
p_boxplot=dfRes %>% 
  dplyr::filter(ct1 >=0, !is.na(ct1), estbinres==1, !is.na(variant_inferred_detail) & 
                  variant_inferred_detail %in% waves) %>% 
  ggplot(aes( y=ct1, x=variant_inferred_detail)) +
  geom_boxplot() +
  ggbeeswarm::geom_quasirandom(alpha=0.05) +
  OverReact::theme_react() +
  scale_x_discrete(labels =function(x) str_wrap(x, width=15)) +
  ggpubr::stat_compare_means(comparisons = mycomparisons,size=3) +
  ggpubr::stat_compare_means(label.y = 90,size=3) +
  labs(y="N-gene CT value", x="", col ="")

p_boxplot

# save
OverReact::saveREACTplot(p = p_boxplot,figpath = figpath,
                         filename = "ct_values_boxplot_n_gene",
                         width = 3,height = 5,savePDF = F)


# Boxplot
p_boxplot_e=dfRes %>% filter(ct2 >=0, !is.na(ct2), estbinres==1, variant_inferred_detail %in% waves &
                               !is.na(variant_inferred_detail)) %>% 
  ggplot(aes( y=ct2, x=variant_inferred_detail)) +
  geom_boxplot() +
  ggbeeswarm::geom_quasirandom(alpha=0.05) +
  OverReact::theme_react() +
  scale_x_discrete(labels =function(x) str_wrap(x, width=15)) +
  ggpubr::stat_compare_means(comparisons = mycomparisons,size=3) +
  ggpubr::stat_compare_means(label.y = 90,size=3) +
  labs(y="E-gene CT value", x="", col ="")

p_boxplot_e

# save
OverReact::saveREACTplot(p = p_boxplot_e,figpath = figpath,
                         filename = "ct_values_boxplot_e_gene",
                         width = 3,height = 5,savePDF = F)



# Add symp/asymp plot: N-gene -----------------------------------------------------


# medians symp
meds_symp=dfRes %>% filter(ct1 >=0, !is.na(ct1), estbinres==1,
                           symptomatic == 1,
                           !is.na(variant_inferred_detail) & variant_inferred_detail %in% waves) %>% 
  group_by(variant_inferred_detail) %>% 
  summarise(med=round(median(ct1, na.rm=T),2))
meds_symp$col = c("white","white")

# medians asymp
meds_asymp=dfRes %>% filter(ct1 >=0, !is.na(ct1), estbinres==1,
                            symptomatic == 0,
                            !is.na(variant_inferred_detail)&variant_inferred_detail %in% waves) %>% 
  group_by(variant_inferred_detail) %>% 
  summarise(med=round(median(ct1, na.rm=T),2))
meds_asymp$col = c("white",rep("black",1))

# Boxplot
p_boxplot_symp=dfRes %>% filter(ct1 >=0, !is.na(ct1), estbinres==1,
                                symptomatic == 1,
                                !is.na(variant_inferred_detail) & 
                                  variant_inferred_detail %in% waves) %>% 
  ggplot(aes( y=ct1, x=variant_inferred_detail)) +
  geom_boxplot() +
  ggbeeswarm::geom_quasirandom(alpha=0.05) +
  OverReact::theme_react() +
  scale_x_discrete(labels =function(x) str_wrap(x, width=12)) +
  ggpubr::stat_compare_means(comparisons = mycomparisons,size=3) +
  ggpubr::stat_compare_means(label.y = 55,size=2.5) +
  geom_text(data=meds_symp,aes(x=variant_inferred_detail, y=med, label = med, col=col), 
            size=2.5, vjust=-0.4)+
  ylim(c(0,60))+
  scale_colour_manual(values = c("white"))+  
  labs(y="N-gene Ct value",   x="", col ="", title = "Symptomatic")

p_boxplot_symp

# Boxplot
p_boxplot_asymp=dfRes %>% filter(ct1 >=0, !is.na(ct1), estbinres==1,
                                 symptomatic == 0,
                                 !is.na(variant_inferred_detail) & 
                                   variant_inferred_detail %in% waves) %>% 
  ggplot(aes( y=ct1, x=variant_inferred_detail)) +
  geom_boxplot() +
  ggbeeswarm::geom_quasirandom(alpha=0.05) +
  OverReact::theme_react() +
  scale_x_discrete(labels =function(x) str_wrap(x, width=12)) +
  ggpubr::stat_compare_means(comparisons = mycomparisons,size=3) +
  ggpubr::stat_compare_means(label.y = 55,size=2.5) +
  geom_text(data=meds_asymp,aes(x=variant_inferred_detail, y=med, label = med, col=col), 
            size=2.5, vjust=-0.4)+
  ylim(c(0,60))+
  scale_colour_manual(values = c("black","black"))+
  labs(y="N-gene Ct value", x="", col ="", title = "Asymptomatic")

p_boxplot_asymp




# combine into one master plot
p_boxplot_all <- p_boxplot_symp+p_boxplot_asymp&theme(legend.position = "none")
p_boxplot_all


# save
OverReact::saveREACTplot(p = p_boxplot_all,figpath = figpath,
                         filename = "ct_values_boxplot_symptom_compare",
                         width = 4.85,height = 5,savePDF = F)




# Add symp/asymp plot: E-gene -----------------------------------------------------

# medians symp
meds_symp_e=dfRes %>% filter(ct2 >=0, !is.na(ct2), estbinres==1,
                             symptomatic == 1,
                             !is.na(variant_inferred_detail)) %>% 
  group_by(variant_inferred_detail) %>% 
  summarise(med=round(median(ct2, na.rm=T),2),
            mean=round(mean(ct2, na.rm=T),2))
meds_symp_e$col = c(rep("black",3),"white","white")

# medians asymp
meds_asymp_e=dfRes %>% filter(ct2 >=0, !is.na(ct2), estbinres==1,
                              symptomatic == 0,
                              !is.na(variant_inferred_detail)) %>% 
  group_by(variant_inferred_detail) %>% 
  summarise(med=round(median(ct2, na.rm=T),2))
meds_asymp_e$col = c(rep("black",5))

# Boxplot
p_boxplot_symp_e=dfRes %>% filter(ct2 >=0, !is.na(ct2), estbinres==1,
                                  symptomatic == 1,
                                  !is.na(variant_inferred_detail)) %>% 
  ggplot(aes( y=ct2, x=variant_inferred_detail)) +
  geom_boxplot() +
  ggbeeswarm::geom_quasirandom(alpha=0.05) +
  OverReact::theme_react() +
  scale_x_discrete(labels =function(x) str_wrap(x, width=15)) +
  ggpubr::stat_compare_means(comparisons = mycomparisons,size=3) +
  ggpubr::stat_compare_means(label.y = 85,size=2.5) +
  geom_text(data=meds_symp_e,aes(x=variant_inferred_detail, y=med, label = med, col=col), 
            size=2.5, vjust=-0.4)+
  scale_colour_manual(values = c("black","white"))+  labs(y="E-gene Ct value", x="", 
                                                          col ="", title = "Symptomatic")

p_boxplot_symp_e

# Boxplot
p_boxplot_asymp_e=dfRes %>% filter(ct2 >=0, !is.na(ct2), estbinres==1,
                                   symptomatic == 0,
                                   !is.na(variant_inferred_detail)) %>% 
  ggplot(aes( y=ct2, x=variant_inferred_detail)) +
  geom_boxplot() +
  ggbeeswarm::geom_quasirandom(alpha=0.05) +
  OverReact::theme_react() +
  scale_x_discrete(labels =function(x) str_wrap(x, width=15)) +
  ggpubr::stat_compare_means(comparisons = mycomparisons,size=3) +
  ggpubr::stat_compare_means(label.y = 85,size=2.5) +
  geom_text(data=meds_asymp_e,aes(x=variant_inferred_detail, y=med, label = med, col=col), 
            size=2.5, vjust=-0.4)+
  scale_colour_manual(values = c("black","white"))+
  labs(y="E-gene Ct value", x="", col ="", title = "Asymptomatic")

p_boxplot_asymp_e




# combine into one master plot
p_boxplot_all_e <- p_boxplot_symp_e+p_boxplot_asymp_e&theme(legend.position = "none")
p_boxplot_all_e
# save
OverReact::saveREACTplot(p = p_boxplot_all_e,figpath = figpath,
                         filename = "ct_values_boxplot_symptom_compare_e_gene",
                         width = 8,height = 5,savePDF = F)






# Stability analysis ------------------------------------------------------


# Function that does everything
doEverything = function(mydat_for_analysis,
                        yname="alpha_delta",
                        analysis_name="alpha_vs_delta",
                        myK=1000,
                        myseed=123,
                        subfolder=NULL,
                        adjustments=NULL, 
                        positive_coefs_only = F,
                        symps=covid_yesnos,
                        pi_list,
                        
                        ...
){
  
  
  # create folder
  if(!is.null(subfolder)){
    dir.create(paste0(outpath,subfolder,"/"),showWarnings = F)
    dir.create(paste0(figpath,subfolder,"/"),showWarnings = F)
  }
  
  
  # Wrangle data ------------------------------------------------------------
  
  set.seed(myseed)
  # filter to complete cases
  predictors=c(adjustments,symps)
  mydat_for_analysis <- mydat_for_analysis[complete.cases(mydat_for_analysis[,predictors]),]
  mydat_train = mydat_for_analysis %>% filter(test_train_split =="train")
  mydat_test = mydat_for_analysis %>% filter(test_train_split == "test")
  penalty.factor=c(rep(0,length(adjustments)),rep(1,length(symps)))
  if(positive_coefs_only){
    lower.limits =c(rep(-Inf,length(adjustments)),rep(0,length(symps)))
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
                                       family = "gaussian",
                                       PFER_method = "MB",
                                       # lower.limits = lower.limits,
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
                                      K = 20,
                                      n_thr = min(10,ncol(X_test)),
                                      family = "gaussian",
                                      reorder_by_normalised_beta = F)
  opt_iter=which.max(stab_out$S)
  
  myplot=plotStabResultsStripedFlipped(results_df, opt_thresh=stab_out$P[opt_iter], 
                                       stab_out = stab_out,
                                       family =  "gaussian",
                                       loss = "Q-squared",
                                       plotOnlyTopN = 15,
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






# Omicron only  -----------------------------------------------------------

dfRes$vax_status_number %>% table(dfRes$round, exclude="none")
dfRes$vax_status_number[dfRes$vax_status_number<0] <- NA_real_


# all positives in r17 (implicitly all omicron) - adjusted on vaccine status
res_omicron_ba1=doEverything(mydat_for_analysis = dfRes %>% 
                               mutate(doublevaxxed=case_when(vax_status_number==2 ~1,
                                                             TRUE ~0),
                                      boosted=case_when(vax_status_number==3 ~1,
                                                        TRUE ~0)
                               ) %>% 
                               filter(round %in% c(17:19), estbinres==1,variant_inferred_detail == "BA.1 (Omicron)"),
                             yname = "ct1",
                             analysis_name = "BA1_r17_19_pos_v_neg_age_sex_vax_adjusted", 
                             myK = 1000,
                             myseed = 122, 
                             subfolder = "ba1_omicron",
                             positive_coefs_only = F,
                             pi_list=seq(0.90,0.998,0.002),
                             adjustments = c("age_scaled","sex", "vax_status_number"))

res_omicron_ba1$plot


# BA2 ---------------------------------------------------------------------


# all positives in r17 (implicitly all omicron) - adjusted on vaccine status
res_omicron_ba2=doEverything(mydat_for_analysis = dfRes %>% 
                               mutate(doublevaxxed=case_when(vax_status_number==2 ~1,
                                                             TRUE ~0),
                                      boosted=case_when(vax_status_number==3 ~1,
                                                        TRUE ~0)
                               ) %>% 
                               filter(round %in% c(17:19), estbinres==1,variant_inferred_detail == "BA.2 (Omicron)"),
                             yname = "ct1",
                             analysis_name = "BA2_r17_19_pos_v_neg_age_sex_vax_adjusted", 
                             myK = 1000,
                             myseed = 122, 
                             subfolder = "ba2_omicron",
                             positive_coefs_only = F,
                             pi_list=seq(0.90,0.998,0.002),
                             adjustments = c("age_scaled","sex", "vax_status_number"))

res_omicron_ba2$plot




# Delta only  -----------------------------------------------------------
table(dfRes$variant_inferred_detail)


# all positives in r17 (implicitly all delta) - adjusted on vaccine status
res_delta_r13_15=doEverything(mydat_for_analysis = dfRes %>% 
                                mutate(doublevaxxed=case_when(vax_status_number==2 ~1,
                                                              TRUE ~0),
                                       boosted=case_when(vax_status_number==3 ~1,
                                                         TRUE ~0)
                                ) %>% 
                                filter(variant_inferred_detail == "Rounds 12-15 (Delta)", estbinres==1),
                              yname = "ct1",
                              analysis_name = "r13_15_pos_v_neg_age_sex_vax_adjusted", 
                              myK = 1000,myseed = 122, 
                              subfolder = "delta",
                              positive_coefs_only = F,
                              adjustments = c("age_scaled","sex", "vax_status_number"))


res_delta_r13_15$plot




# Alpha only  -----------------------------------------------------------

# all positives in r17 (implicitly all delta) - adjusted on vaccine status
res_alpha_r8_10=doEverything(mydat_for_analysis = dfRes %>% 
                               mutate(doublevaxxed_or_boosted=case_when(vax_status_number>=2 ~1,
                                                                        TRUE ~0)
                               ) %>% 
                               filter(variant_inferred_detail == "Rounds 8-10 (Alpha)", estbinres==1),
                             yname = "ct1",
                             analysis_name = "r8_10_pos_v_neg_age_sex_vax_adjusted", 
                             myK = 1000,
                             myseed = 123, 
                             subfolder = "alpha",
                             positive_coefs_only = F,
                             pi_list=seq(0.89,0.998,0.002),
                             adjustments = c("age_scaled","sex", "vax_status_number"))




# Wildtype only  -----------------------------------------------------------


# all positives in r17 (implicitly all delta) - adjusted on vaccine status
res_alpha_r2_7=doEverything(mydat_for_analysis = dfRes %>% 
                              mutate(doublevaxxed_or_boosted=case_when(vax_status_number>=2 ~1,
                                                                       TRUE ~0)
                              ) %>% 
                              filter(variant_inferred_detail == "Rounds 2-7 (Wild type)",round %in% c(2:7), estbinres==1),
                            yname = "ct1",
                            analysis_name = "r2_7_pos_v_neg_age_sex_adjusted", 
                            myK = 1000,
                            myseed = 123, 
                            positive_coefs_only = F,
                            subfolder = "wildtype",
                            adjustments = c("age_scaled","sex"))


res_alpha_r2_7$plot
res_alpha_r2_7$results




# Univariate --------------------------------------------------------------

dfRes$age_group_named %>% table()

### Function to run all univariate
runUniv <- function(mydata, adjustments = c("age_group_named","sex","vax_status_number"),
                    symps=covid_yesnos){
  mydata$symptomatic_26 = as.numeric(rowSums(mydata[,symps],na.rm = T)>0)
  univ_df <- data.frame(Variable ="Symptoms",Category = sympnames_type_df$symptom,Beta = NA_real_, 
                        Lower=NA_real_,Upper=NA_real_, 
                        Symp_split = c(rep("Individual symptoms",26),rep("Any", 1)))
  
  for (i in 1:27){
    mydata$predictor = as.numeric(unlist(mydata[,sympnames_type_df$symptom_code[[i]]]))
    f=as.formula(paste("ct1 ~ predictor +", paste(adjustments,collapse = "+")))
    mod_lm <- lm(f, data = mydata)
    modsum=jtools::summ(mod_lm) 
    univ_df$Beta[[i]] <- modsum$coeftable[2,1]
    univ_df$Lower[[i]] <- modsum$coeftable[2,1] - 1.96*modsum$coeftable[2,2]
    univ_df$Upper[[i]] <- modsum$coeftable[2,1] + 1.96*modsum$coeftable[2,2]
  }
  return(univ_df)
}

## Run over all waves
## Wildtype
univ_w <- runUniv(mydata = dfRes  %>% 
                    filter(round %in% c(2:7), estbinres==1),
                  adjustments = c("age_group_named","sex"))
univ_w$model =  univ_w$adjustment ="Wildtype"

## Alpha
univ_a <- runUniv(mydata = dfRes  %>% 
                    filter(round %in% c(8:10), estbinres==1))
univ_a$model = univ_a$adjustment = "Alpha"

## Delta
univ_d <- runUniv(mydata = dfRes  %>% 
                    filter(round %in% c(12:15), estbinres==1))
univ_d$model =  univ_d$adjustment ="Delta"

## BA1
univ_ba1 <- runUniv(mydata = dfRes  %>% 
                    filter(variant_inferred_detail == "BA.1 (Omicron)", estbinres==1))
univ_ba1$model =  univ_ba1$adjustment ="BA.1"


## BA2
univ_ba2 <- runUniv(mydata = dfRes  %>% 
                    filter(variant_inferred_detail == "BA.2 (Omicron)", estbinres==1))
univ_ba2$model =  univ_ba2$adjustment ="BA.2"




## Omicron
univ_omi <- runUniv(mydata = dfRes  %>% 
                      filter(round %in%c(17:19), estbinres==1))
univ_omi$model =  univ_omi$adjustment ="Omicron"





### univ_all
univ_all <- rbind(univ_w,univ_a,univ_d,univ_ba1,univ_ba2)
# univ_all$adjustment = "1"
adjustment_numbers=adjustment_descriptions=c("Wildtype","Alpha","Delta","BA.1","BA.2")
insignificant_results_greyed_out=strip_borders=T

univ_all$P_value = case_when(univ_all$Lower*univ_all$Upper <0 ~ 1,
                             T ~ 0)
univ_all <- univ_all
stripes_df <- univ_all %>%  
  group_by(Category) %>% 
  summarise(mean_or=mean(Beta, na.rm=T)) %>% 
  arrange(-mean_or) %>% 
  ungroup() %>% 
  mutate(stripe_col = case_when(row_number() %%2 ==0 ~ "NA",
                                TRUE ~ "grey90"))

stripes_df$Category <- factor(stripes_df$Category, levels=unique(stripes_df$Category))

# create dummy df of adjustment levels for simple joining rather than complex case_When
adj_df <- data.frame(adjustment = as.character(adjustment_numbers), 
                     adjust_desc = factor(adjustment_descriptions,levels = adjustment_descriptions))
dodger = position_dodge2(width=0.9,reverse = T)
alpha_val=0
strip_text_size = 8

# create bold label vector
boldlabel = (c(rep("plain",19),"bold","plain","plain",rep("plain",6)))
max(univ_all$Upper)
# make plot
p_univ <- univ_all %>% 
  dplyr::filter(adjustment %in% adjustment_numbers) %>% 
  dplyr::mutate(Category=factor(Category, levels = unique(Category)),
                Variable=factor(Variable, levels = unique(Variable)),
                shape=case_when(!insignificant_results_greyed_out ~ "15",
                                P_value >0.05 ~ "0",
                                TRUE ~ "15")) %>%
  left_join(adj_df) %>% 
  left_join(stripes_df) %>% 
  ggplot(aes(x=(Category), y = Beta, col = adjust_desc, group = shape)) +
  geom_tile(data = stripes_df, aes(x=Category, y =1, height = Inf,
                                   fill = stripe_col),inherit.aes = F,
            alpha=0.8,
            col = if(strip_borders) "grey70" else NA,
            linetype = "dashed",
            show.legend = F) +
  scale_fill_manual(values = c("white","grey96")) +
  ggnewscale::new_scale_fill() +
  geom_hline(yintercept = 0, linetype = "dashed", col = "grey50", size=0.1 )+
  geom_point(position = dodger, size=0.9,aes(shape = shape)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), 
                position = dodger, width = 0.9,
                show.legend = F, size=0.2) +
  coord_flip() +
  scale_shape_manual("group",values = c("0" = 0, "15" = 15),guide = "none") +
  # scale_y_continuous(trans = "log10",breaks = x_seq) +
  # scale_y_continuous(trans = "log2", #breaks = x_seq, 
  #                    labels =function(x) MASS::fractions(x))+ 
  # scale_color_manual(values = myCols) +
  OverReact::scale_color_imperial() +
  OverReact::theme_react(strip_text_size = strip_text_size) +
  # ggforce::facet_col(.~Symp_split,  scales = "free", space = "free") +
  labs(x="", y="Change in Ct value (adjusted)", col ="") +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        panel.grid.major.x = element_line(size = rel(0.1),linetype = "dashed")
  )+ 
  theme(axis.text.y =element_text(face = boldlabel)) +
  guides(col=guide_legend(nrow=2,byrow=T))
p_univ


# save
OverReact::saveREACTplot(p = p_univ,figpath = figpath,
                         filename = "univariate_ct_vals",
                         width = 5,height = 6,savePDF = F)



# Show BA1/BA2 only -------------------------------------------------------
univ_all <- univ_all
stripes_df <- univ_all  %>% 
  dplyr::filter(adjustment %in% c("BA.1","BA.2"))%>%  
  group_by(Category) %>% 
  summarise(mean_or=mean(Beta, na.rm=T)) %>% 
  arrange(-mean_or) %>% 
  ungroup() %>% 
  mutate(stripe_col = case_when(row_number() %%2 ==0 ~ "NA",
                                TRUE ~ "grey90"))

stripes_df$Category <- factor(stripes_df$Category, levels=unique(stripes_df$Category))

# create dummy df of adjustment levels for simple joining rather than complex case_When
adj_df <- data.frame(adjustment = as.character(adjustment_numbers), 
                     adjust_desc = factor(adjustment_descriptions,levels = adjustment_descriptions))
dodger = position_dodge2(width=0.9,reverse = T)
alpha_val=0
strip_text_size = 8

# create bold label vector
boldlabel = (c(rep("plain",19),"bold","plain","plain",rep("plain",6)))
max(univ_all$Upper)

# create bold label vector
boldlabel = c(rep("plain",26),"bold")
max(univ_all$Upper)

# make plot
p_univ_omi <- univ_all %>% 
  dplyr::filter(adjustment %in% c("BA.1","BA.2")) %>% 
  dplyr::mutate(Category=factor(Category, levels = unique(Category)),
                Variable=factor(Variable, levels = unique(Variable)),
                shape=case_when(!insignificant_results_greyed_out ~ "15",
                                P_value >0.05 ~ "0",
                                TRUE ~ "15")) %>%
  left_join(adj_df) %>% 
  left_join(stripes_df) %>% 
  ggplot(aes(x=(Category), y = Beta, col = adjust_desc, group = shape)) +
  geom_tile(data = stripes_df, aes(x=Category, y =1, height = Inf,
                                   fill = stripe_col),inherit.aes = F,
            alpha=0.8,
            col = if(strip_borders) "grey70" else NA,
            linetype = "dashed",
            show.legend = F) +
  scale_fill_manual(values = c("white","grey96")) +
  ggnewscale::new_scale_fill() +
  geom_hline(yintercept = 0, linetype = "dashed", col = "grey50", size=0.1 )+
  geom_point(position = dodger, size=0.9,aes(shape = shape)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), 
                position = dodger, width = 0.9,
                show.legend = F, size=0.2) +
  coord_flip() +
  scale_shape_manual("group",values = c("0" = 0, "15" = 15),guide = "none") +
  # scale_y_continuous(trans = "log10",breaks = x_seq) +
  # scale_y_continuous(trans = "log2", #breaks = x_seq, 
  #                    labels =function(x) MASS::fractions(x))+ 
  # scale_color_manual(values = myCols) +
  OverReact::scale_color_imperial() +
  OverReact::theme_react(strip_text_size = strip_text_size) +
  # ggforce::facet_col(.~Symp_split,  scales = "free", space = "free") +
  labs(x="", y="Change in Ct value (adjusted)", col ="") +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        panel.grid.major.x = element_line(size = rel(0.1),linetype = "dashed")
  )+ 
  theme(axis.text.y =element_text(face = boldlabel))+
  scale_colour_manual(values = c( "#EC7300", "#960078" ))
p_univ_omi


# save
OverReact::saveREACTplot(p = p_univ_omi,figpath = figpath,
                         filename = "univariate_ct_vals_BA1_BA2",
                         width = 5,height = 5,savePDF = F)




# Omicron phase only ------------------------------------------------------



# univ_omi$adjustment = "1"
adjustment_numbers=adjustment_descriptions=c("Wildtype","Alpha","Delta","BA.1","BA.2")
insignificant_results_greyed_out=strip_borders=T

univ_omi$P_value = case_when(univ_omi$Lower*univ_omi$Upper <0 ~ 1,
                             T ~ 0)
univ_omi <- univ_omi
stripes_df <- univ_omi %>%  
  group_by(Category) %>% 
  summarise(mean_or=mean(Beta, na.rm=T)) %>% 
  arrange(-mean_or) %>% 
  ungroup() %>% 
  mutate(stripe_col = case_when(row_number() %%2 ==0 ~ "NA",
                                TRUE ~ "grey90"))

stripes_df$Category <- factor(stripes_df$Category, levels=unique(stripes_df$Category))

# create dummy df of adjustment levels for simple joining rather than complex case_When
adj_df <- data.frame(adjustment = as.character(adjustment_numbers), 
                     adjust_desc = factor(adjustment_descriptions,levels = adjustment_descriptions))
dodger = position_dodge2(width=0.9,reverse = T)
alpha_val=0
strip_text_size = 8

# create bold label vector
boldlabel = (c(rep("plain",23),"bold","plain","plain",rep("plain",4)))
col ="#00ACD7"

# make plot
p_univ <- univ_omi %>% 
  dplyr::mutate(Category=factor(Category, levels = unique(Category)),
                Variable=factor(Variable, levels = unique(Variable)),
                shape=case_when(!insignificant_results_greyed_out ~ "15",
                                P_value >0.05 ~ "0",
                                TRUE ~ "15")) %>%
  left_join(adj_df) %>% 
  left_join(stripes_df) %>% 
  ggplot(aes(x=(Category), y = Beta, col = adjustment, group = shape)) +
  geom_tile(data = stripes_df, aes(x=Category, y =1, height = Inf,
                                   fill = stripe_col),inherit.aes = F,
            alpha=0.8,
            col = if(strip_borders) "grey70" else NA,
            linetype = "dashed",
            show.legend = F) +
  scale_fill_manual(values = c("white","grey96")) +
  ggnewscale::new_scale_fill() +
  geom_hline(yintercept = 0, linetype = "dashed", col = "grey50", size=0.1 )+
  geom_point(position = dodger, size=0.9,aes(shape = shape), col =col ) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), 
                position = dodger, width = 0.3,
                show.legend = F, size=0.2, col =col ) +
  coord_flip() +
  scale_shape_manual("group",values = c("0" = 0, "15" = 15),guide = "none") +
  OverReact::scale_color_imperial() +
  OverReact::theme_react(strip_text_size = strip_text_size) +
  labs(x="", y="Change in Ct value (adjusted)", col ="") +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        panel.grid.major.x = element_line(size = rel(0.1),linetype = "dashed")
  )+ 
  theme(axis.text.y =element_text(face = boldlabel)) +
  guides(col=guide_legend(nrow=2,byrow=T))
p_univ

OverReact::imperial_palettes
# save
OverReact::saveREACTplot(p = p_univ,figpath = figpath,
                         filename = "univariate_ct_vals",
                         width = 5,height = 5,savePDF = F)

