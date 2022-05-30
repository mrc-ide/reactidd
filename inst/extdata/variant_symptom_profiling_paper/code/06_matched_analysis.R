
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
                  "tidyr", "pheatmap","catboost","datapasta","ggtext","Rcpp",
                  "ggplot2","ggsci", "RColorBrewer", "tidyverse", "lubridate", 
                  "readr","ggthemes", "questionr", "foreach", "doParallel","withr",
                  "patchwork","randomcoloR","focus","ComplexHeatmap"
)
load_packages(package.list)


# Import REACT-1 data -----------------------------------------------------
source("E:/Group/react2_study5/report_phases_combined/projects/omicron_symptom_profiling/code/00_bits_and_pieces.R", local = T)
source("E:/Group/react2_study5/report_phases_combined/projects/omicron_symptom_profiling/code/00_data_prep.R", local = T)

# create subfolder
createMySubfolder(subfolderName = "matched_analysis")



# Functions ---------------------------------------------------------------


## Manual matching ##
generateMatchedData <- function(data,
                                case_passcodes,
                                control_passcodes,
                                age_offset=5,
                                match_ratio=2,
                                myseed = 123,
                                strata_offset=0){
  set.seed(myseed)
  data <- data[!duplicated(data$u_passcode),]
  
  # clean up vax status
  data$vax_status_number[data$vax_status_number<0] <- NA_real_
  data=data %>% filter(u_passcode %in% c(case_passcodes,control_passcodes))
  data_case=data %>% filter(u_passcode %in% c(case_passcodes))
  data_control=data %>% filter(u_passcode %in% c(control_passcodes))
  internal_case_passcodes=case_passcodes
  internal_control_passcodes=control_passcodes
  match_pcs <- c()
  # matched_data=data.frame()
  
  
  # get summary data frame to speed up loop
  dat_summary=data_case %>% filter(u_passcode %in% c(case_passcodes)) %>% 
    group_by(round,age,vax_status_number) %>% 
    summarise(n=n())
  
  pb=progress::progress_bar$new(format = " Running matching [:bar] :percent eta: :eta",
                                width = 100,clear = F,
                                total = nrow(dat_summary))
  
  pb$tick(0)
  
  strata_count=1
  res_list=list()
  # run loop
  for(i in 1:nrow(dat_summary)){
    pb$tick()
    
    pc <- internal_case_passcodes[[i]]
    age_o <- dat_summary$age[[i]]
    vax_o <- dat_summary$vax_status_number[[i]]
    round_o <- dat_summary$round[[i]]
    
    
    dat_temp_control <- data_control %>% filter(round == round_o,
                                                age>=(age_o-age_offset), age<=(age_o+age_offset))
    
    
    dat_temp_case <- data_case %>% filter(round == round_o,
                                          # is.na(vax_status_number),
                                          # vax_status_number==vax_o,
                                          age==age_o)
    
    
    
    
    if(!is.na(vax_o)){
      dat_temp_control <- dat_temp_control %>% filter(vax_status_number ==vax_o)
      dat_temp_case <- dat_temp_case %>% filter(vax_status_number ==vax_o)
      
    }else{
      dat_temp_control <- dat_temp_control %>% filter(is.na(vax_status_number))
      dat_temp_case <- dat_temp_case %>% filter(is.na(vax_status_number))
      
    }
    
    mpc=sample(x = 1:nrow(dat_temp_control),size = dat_summary$n[[i]]*match_ratio,replace = T)
    # internal_control_passcodes=setdiff(internal_control_passcodes,mpc)
    # match_pcs <- c(match_pcs,mpc)
    contdat=dat_temp_control[mpc,]
    if(nrow(contdat)>0){
      dat_temp_case$strata=contdat$strata=strata_count:(strata_count+nrow(contdat)-1)
      strata_count=strata_count+nrow(contdat)
      res_list[[i]] <- rbind(dat_temp_case,contdat)
    }
    
  }
  
  # bind all rows
  matched_data=bind_rows(res_list)
  
  
  # output new data set
  return(matched_data)
  
}



# Generate matched data ---------------------------------------------------



# Round 17
ba2_17_pcs <- dfRes$u_passcode[dfRes$variant_inferred_detail=="BA.2 (Omicron)" & !is.na(dfRes$variant_inferred_detail) & dfRes$round ==17]
ba1_17_pcs <- dfRes$u_passcode[dfRes$variant_inferred_detail=="BA.1 (Omicron)" & !is.na(dfRes$variant_inferred_detail) & dfRes$round ==17]

### Generate matched data
r17_matched <- generateMatchedData(data = dfRes %>% filter(round == 17, !is.na(vax_status_number)) ,
                                     case_passcodes = ba2_17_pcs,
                                     control_passcodes = ba1_17_pcs,
                                     age_offset = 5,
                                     match_ratio = 1,
                                     myseed = 123)


# Round 18
ba2_18_pcs <- dfRes$u_passcode[dfRes$variant_inferred_detail=="BA.2 (Omicron)" & !is.na(dfRes$variant_inferred_detail) & dfRes$round ==18]
ba1_18_pcs <- dfRes$u_passcode[dfRes$variant_inferred_detail=="BA.1 (Omicron)" & !is.na(dfRes$variant_inferred_detail) & dfRes$round ==18]

### Generate matched data
r18_matched <- generateMatchedData(data = dfRes %>% filter(round == 18, !is.na(vax_status_number)) ,
                                   case_passcodes = ba2_18_pcs,
                                   control_passcodes = ba1_18_pcs,
                                   age_offset = 5,
                                   match_ratio = 1,
                                   myseed = 123,
                                   strata_offset = max(r17_matched$strata))


# Round 19
ba2_19_pcs <- dfRes$u_passcode[dfRes$variant_inferred_detail=="BA.2 (Omicron)" & !is.na(dfRes$variant_inferred_detail) & dfRes$round ==19]
ba1_19_pcs <- dfRes$u_passcode[dfRes$variant_inferred_detail=="BA.1 (Omicron)" & !is.na(dfRes$variant_inferred_detail) & dfRes$round ==19]

### Generate matched data
r19_matched <- generateMatchedData(data = dfRes %>% filter(round == 19,!is.na(vax_status_number)) ,
                                   case_passcodes = ba1_19_pcs,
                                   control_passcodes = ba2_19_pcs,
                                   age_offset = 5,
                                   match_ratio = 1,
                                   myseed = 123,
                                   strata_offset = max(r18_matched$strata))

# bind all matched data
matched_data = rbind(r17_matched,r18_matched,r19_matched)
matched_data$symptomatic_26 = as.numeric(rowSums(matched_data[,covid_yesnos],na.rm = T)>0)
matched_data$ba2=case_when(matched_data$variant_inferred_detail=="BA.2 (Omicron)" ~ 1,
                           T~0)



# Run analysis ------------------------------------------------------------


univ_df <- data.frame(Variable ="Symptoms",Category = sympnames_type_df$symptom,OR = NA_real_, 
                      Lower=NA_real_,Upper=NA_real_, Symp_split = c(rep("Individual symptoms",26),rep("Any", 1)))

for ( i in 1:27){
  matched_data$predictor = as.numeric(unlist(matched_data[,sympnames_type_df$symptom_code[[i]]]))
  mod_clogistic <- Epi::clogistic(ba2 ~ predictor, strata = strata, data = matched_data)
  coef=mod_clogistic$coefficients
  se=mod_clogistic$var^0.5
  univ_df$OR[[i]] <- exp(coef)
  univ_df$Lower[[i]] <-exp(coef-1.96*se)
  univ_df$Upper[[i]] <-exp(coef+1.96*se)
  
}

univ_df$sig = case_when(univ_df$Lower <1 & univ_df$Upper <1 ~ 1,
                        univ_df$Lower >1 & univ_df$Upper >1 ~ 1,
                        T~0)


univ_df$P_value = case_when(univ_df$Lower <1 & univ_df$Upper <1 ~ 0.01,
                            univ_df$Lower >1 & univ_df$Upper >1 ~ 0.01,
                            T~1)
univ_df$adjustment ="1"
univ_df$model =1


# create bold label vector
univ_df <- univ_df %>% arrange(-OR)
boldlabel=rev(case_when(univ_df$Category=="Any of 26 symptoms" ~ "bold",
                        T~"plain"))

p_univ <- plotReactForest(univ_df_plot = univ_df,
                          adjustment_numbers = 1,
                          adjustment_descriptions = "",
                          legend.position = "none",
                          alpha_val = 0.5) + 
  theme(axis.text.y =element_text(face = boldlabel))

p_univ

# save
OverReact::saveREACTplot(p = p_univ,figpath = figpath,
                         filename = "univariate_matched_conditional_logistic_ba1_ba2",
                         width = 5,height = 6,savePDF = F)








# Combined plot with standard univariate ----------------------------------

# load 'normal' regression results
univ_results <- readRDS("E:/Group/react2_study5/report_phases_combined/projects/omicron_symptom_profiling//output/bA1_bA2/univ_ba2_ba1.rds")
univ_results$model="Logistic regression adjusted on age, sex, vaccine, and round"
univ_results$suffix="unmatched"

univ_df$model="Conditional logistic regression on matched data"
univ_df$suffix="matched"

univ_comb_df=bind_rows(univ_results,univ_df)
setdiff(names(univ_df),names(univ_results))
univ_comb_df$Variable="Symptom in past week"

# univ_comb_df$adjustment=univ_comb_df$
p_forest=plotReactForest(univ_df_plot = univ_comb_df,adjustment_numbers = c(1,5),
                adjustment_descriptions = c("Conditional logistic regression on matched data",
                                            "Logistic regression adjusted on age, sex, vaccine, and round"
                                            ),
                palette = "two_col_grey_pool",arrangeByOR = T,
                )+
  guides(col=guide_legend(nrow=2)) +
  scale_colour_manual(values = c(myCols[[3]],myCols[[1]]))

p_forest


# Build df for scatter plot
df_for_scatter=univ_comb_df %>% pivot_wider(id_cols = Category,
                             names_from = suffix,
                             values_from = c(OR,Lower,Upper)) %>% 
  mutate(col_sig=case_when(Lower_unmatched >1 & Lower_matched >1 ~"p<0.05 in both bodels",
                           Lower_unmatched >1  ~"p<0.05 in unmatched model",
                           T ~"p>0.05 in both models"
  ),
  text_sig = Lower_unmatched >1,
  label=ifelse(text_sig, Category,"")) 

# get colour labels for forest plot

labs_df=data.frame(Category=ggplot_build(p_forest)$layout$panel_params[[1]]$y$get_labels(), order=1:27)
labs_df <- labs_df %>% left_join(df_for_scatter) %>% 
  mutate(cols=case_when(grepl("unmatched",col_sig ) ~ "orange",
                            grepl("p<0.05",col_sig ) ~ "firebrick2",
                            T ~"grey40"))
cols_sig=rev(labs_df$cols)
bold_sig=rev(ifelse(cols_sig=="grey40","plain","bold"))


# add colour to forest
p_forest <- p_forest +theme(axis.text.y = element_text(face = bold_sig,
                                                       colour = rev(cols_sig )))
p_forest

df_for_scatter$Upper_matched %>% max()
p_scatter= df_for_scatter%>% 
  ggplot(aes(x=OR_unmatched,OR_matched,
             size=factor(ifelse(Category=="Any of 26 symptoms",2,1)), 
             col=col_sig,
             label=label)) +
  geom_abline(intercept = 0,slope = 1, linetype="dashed", size=0.2)+
  geom_vline(xintercept = 1,linetype="dashed", size=0.2, col="grey60")+
  geom_hline(yintercept = 1,linetype="dashed", size=0.2, col="grey60")+
  geom_errorbar(aes(xmin=Lower_unmatched,xmax=Upper_unmatched), size=0.1) +
  geom_errorbar(aes(ymin=Lower_matched,ymax=Upper_matched), size=0.1) +
  geom_point() +
  scale_size_manual(values = c(2,2),guide=F)+
  scale_color_manual(values = rev(c("grey40","orange","firebrick2"))) +
  ggrepel::geom_label_repel(force = 4,label.padding = 0.15,min.segment.length = 0.01,
                            show.legend =F) +
  
  scale_y_continuous(trans = "log2", #breaks = x_seq, 
                     labels =function(x) MASS::fractions(x),limits =c(0.55,4.3) )+ 
  scale_x_continuous(trans = "log2", #breaks = x_seq, 
                     labels =function(x) MASS::fractions(x),limits =c(0.55,4.3))+ 
  OverReact::theme_react() +
  labs(x="Odds ratio (logistic regression)", 
       y="Odds ratio (matched conditional logistic regression)",
       col="") +
  theme(legend.position = "bottom") +
  guides(col=guide_legend(nrow=3))
p_scatter

p_combined=p_forest+p_scatter +plot_layout(widths = c(2,5))


OverReact::saveREACTplot(p = p_scatter,figpath = figpath,
                         filename = "univ_matched_unmatched_compare_scatter",width = 5,
                         height = 6
                         )


OverReact::saveREACTplot(p = p_combined,figpath = figpath,
                         filename = "univ_matched_unmatched_compare_panel",width = 9,
                         height = 6
)

