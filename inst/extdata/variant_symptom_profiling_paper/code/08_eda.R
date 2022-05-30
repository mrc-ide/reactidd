
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




#' Pull in packages needed
package.list <- c("prevalence","mgcv","knitr","MASS","kableExtra","table1","dplyr",
                  "tidyr", "pheatmap","OverReact",
                  "ggplot2","ggsci", "RColorBrewer", "tidyverse", "lubridate", 
                  "readr","ggthemes", "questionr", "foreach", "doParallel","withr",
                  "patchwork","randomcoloR","focus"
)

load_packages(package.list)


# Import REACT-1 data -----------------------------------------------------
source("E:/Group/react2_study5/report_phases_combined/projects/omicron_symptom_profiling/code/00_bits_and_pieces.R",
       local = T)
source("E:/Group/react2_study5/report_phases_combined/projects/omicron_symptom_profiling/code/00_data_prep.R", 
       local = T)


# create subfolder
createMySubfolder(subfolderName = "eda")

# Table one ---------------------------------------------------------------
dfRes$symptomatic_26

dfRes_tab <- dfRes %>% 
  filter(round %in% c(2:7,8:10,13:15,16,17,18,19)) %>% 
  mutate(round = as.numeric(round),
         wave = factor(case_when(round <8 ~ "Rounds 2-7 (Wild type)",
                          round <11 ~ "Rounds 8-10 (Alpha)",
                          round %in% c(13:15) ~ "Rounds 13-15 (Delta)",
                          round%in% c(17) ~ "Round 17 (BA.1 Omicron)",
                          round%in% c(19) ~ "Rounds 19 (BA.2 Omicron)"
                          ),
                       levels = c("Rounds 2-7 (Wild type)",
                                  "Rounds 8-10 (Alpha)",
                                  "Rounds 13-15 (Delta)",
                                  "Round 17 (BA.1 Omicron)",
                                  "Rounds 19 (BA.2 Omicron)")),
         wave_omi = factor(case_when(round <8 ~ "Rounds 2-7 (Wild type)",
                                 round <11 ~ "Rounds 8-10 (Alpha)",
                                 round %in% c(13:15) ~ "Rounds 13-15 (Delta)",
                                 round%in% c(17:19) ~ "Rounds 17-19 (Omicron)",
         ),
         levels = c("Rounds 2-7 (Wild type)",
                    "Rounds 8-10 (Alpha)",
                    "Rounds 13-15 (Delta)",
                    "Rounds 17-19 (Omicron)")),
         
         brediff_cat = case_when(brediff==2 ~ "Yes, Even when sitting/lying down",
                                 brediff==1 ~ "Yes, affecting normal activities",
                                 brediff==3 ~ "No, affecting normal activities",
                                 T ~ NA_character_
                                 ),
         covidability_cat = factor(case_when(symptomatic_26==0 ~ "No symptoms reported",
                                               covidability7 == 1 ~ "A lot",
                                      covidability7 == 2 ~ "A little",
                                      covidability7 == 3 ~ "Not at all",
                                      covidability7 %in% c(4,5) ~ "Don't know / PNA / Non-reponse",
                                      T ~ "Don't know / PNA / Non-reponse",
                                      ), 
                                   levels = c("A lot","A little","Not at all",
                                                    "Don't know / PNA / Non-reponse","No symptoms reported")))

# create symptom count variable
dfRes_tab$symptom_count_26 <- rowSums(dfRes_tab[,covid_yesnos], na.rm=T)
dfRes_tab$symptom_count_4 <- rowSums(dfRes_tab[,sympnames_type_df$symptom_code[1:4]], na.rm=T)

dfRes_tab$seekmed_cat %>% table(dfRes_tab$wave, exclude="none") %>% prop.table() *100
dfRes_tab$covidability7[dfRes_tab$estbinres==1] %>% 
  table(dfRes_tab$symptomatic_26[dfRes_tab$estbinres==1], exclude="none") %>% 
  prop.table(2) *100



dfRes_tab$one_of_four <- as.factor(dfRes_tab$one_of_four)
dfRes_tab$symptomatic <- as.factor(dfRes_tab$symptomatic)
dfRes_tab$estbinres_char <- ifelse(dfRes_tab$estbinres==1,"Yes","No")
dfRes_tab$round <- as.character(dfRes_tab$round)
dfRes_tab$ethnic_new %>% table(dfRes_tab$round)
table(dfRes$round,dfRes$covida_cat)
cov_name_list$prior_inf = "Prior COVID-19 infection"
cov_name_list$symptom_count_26 = "Number of reported symptoms"
cov_name_list$symptom_count_4 = "Number of reported symptoms (classic symptoms only)"
cov_name_list$vax_status_noDate_v2 = "Vaccination status"
cov_name_list$covidability_cat="Symptoms affecting day-to-day activities"
cov_name_list$seekmed_cat="Sought medical attention for symptoms"
cov_name_list$brediff_cat="Breathing difficulties"

# define list of rowvars
rowvar_list=c("all_participants", "sex","age_group_named","ethnic_new","estbinres_char","prior_inf",
              "vax_status_noDate_v2",
              "symptomatic", "symptom_count_26","seekmed_cat","covidability_cat",
              "days_since_symptom_onset",covid_yesnos,sympnames_type_df$first_symptom_code[1:26])

# Run on all participants
table_one <- OverReact::crossTabMulti(dat = dfRes_tab,rowvar_list =  rowvar_list,
                                      colvar = "wave",cov_names = cov_name_list,
                                      confint = T,include_percentages = T,comma_thousands = T,
                                      rowwise_precentages = F,statistical_test = F)





# Tinker
table_one <- table_one %>% filter(Category!=0)
switchindex <- table_one$Category==1
table_one$Category[switchindex] <- table_one$Variable[switchindex]
table_one$Variable[switchindex] <- "Symptoms"

### Replace NaNs
table_one[table_one=="0 (NaN%)"] <- "0"
table_one$Variable[(nrow(table_one)-25):nrow(table_one)] <- "First reported symptoms"



# Table one just omicron --------------------------------------------------

# Run on all participants
table_one_omicron <- OverReact::crossTabMulti(dat = dfRes_tab,rowvar_list =  rowvar_list,
                                      colvar = "wave_omi",cov_names = cov_name_list,
                                      confint = T,include_percentages = T,comma_thousands = T,
                                      rowwise_precentages = F,statistical_test = F)
# Tinker
table_one_omicron <- table_one_omicron %>% filter(Category!=0)
switchindex <- table_one_omicron$Category==1
table_one_omicron$Category[switchindex] <- table_one_omicron$Variable[switchindex]
table_one_omicron$Variable[switchindex] <- "Symptoms"

### Replace NaNs
table_one_omicron[table_one_omicron=="0 (NaN%)"] <- "0"
table_one_omicron$Variable[(nrow(table_one_omicron)-25):nrow(table_one_omicron)] <- "First reported symptoms"




# Run on positives --------------------------------------------------------

dfRes %>% filter(estbinres==1) %>% group_by(variant_inferred_detail) %>% 
  summarise(mean_days=mean(days_since_symptom_onset,na.rm=T))

# Run on only positives
table_one_pos <- OverReact::crossTabMulti(dat = dfRes_tab %>% filter( estbinres==1),rowvar_list =  rowvar_list,
                                      colvar = "variant_inferred_detail",cov_names = cov_name_list,
                                      confint = T,include_percentages = T,comma_thousands = T,
                                      rowwise_precentages = F,statistical_test = F)
# Tinker
table_one_pos <- table_one_pos %>% filter(Category!=0)
switchindex <- table_one_pos$Category==1
table_one_pos$Category[switchindex] <- table_one_pos$Variable[switchindex]
table_one_pos$Variable[switchindex] <- "Symptoms"
### Replace NaNs
table_one_pos[table_one_pos=="0 (NaN%)"] <- "0"
# Add first symptom name
table_one_pos$Variable[(nrow(table_one_pos)-25):nrow(table_one_pos)] <- "First reported symptoms"

# Run on only positives
table_one_pos_60_plus <- OverReact::crossTabMulti(dat = dfRes_tab %>% filter( estbinres==1, age>=60),
                                                  rowvar_list =  rowvar_list,
                                          colvar = "variant_inferred_detail",cov_names = cov_name_list,
                                          confint = F,include_percentages = T,
                                          rowwise_precentages = F,statistical_test = F)
# Tinker
table_one_pos_60_plus <- table_one_pos_60_plus %>% filter(Category!=0)
switchindex <- table_one_pos_60_plus$Category==1
table_one_pos_60_plus$Category[switchindex] <- table_one_pos_60_plus$Variable[switchindex]
table_one_pos_60_plus$Variable[switchindex] <- "Symptoms"
### Replace NaNs
table_one_pos_60_plus[table_one_pos_60_plus=="0 (NaN%)"] <- "0"



# FIrst symptoms ----------------------------------------------------------

# add first symptoms to covnamelist
cov_name_list[sympnames_type_df$first_symptom_code[1:26]] <- sympnames_type_df$symptom[1:26]

# Sort out symptoms
dfRes_tab <- dfRes_tab %>% 
  mutate_at(all_of(sympnames_type_df$first_symptom_code[1:26]),binaryCleaner_1_0) 


### Replace NA with 0 in symptoms
dfRes_tab[sympnames_type_df$first_symptom_code[1:26]][is.na(dfRes_tab[sympnames_type_df$first_symptom_code[1:26]])] <- 0

# Run on only positives
table_one_pos_first_symp <- OverReact::crossTabMulti(dat = dfRes_tab %>% 
                                                       filter( estbinres==1, symptomatic==1),
                                                     rowvar_list =  
                                                       sympnames_type_df$first_symptom_code[1:26],
                                          colvar = "variant_inferred_detail",
                                          cov_names = cov_name_list,
                                          confint = F,
                                          include_percentages = T,
                                          rowwise_precentages = F,statistical_test = F) 

table_one_pos_first_symp <- table_one_pos_first_symp %>% filter(Category!=0)


### save workbooks
savePrettyExcelWorkbook(listOfTables = list(tab1_all=table_one, table_one_omicron=table_one_omicron,
                                            tab1_pos=table_one_pos,
                                            tab1_pos_60_plus=table_one_pos_60_plus,
                                            tab1_first_symptom=table_one_pos_first_symp),
                        workbookName = "table_one",outpath = outpath)




# Plot symptom prevalence per wave ----------------------------------------

dfRes_tab$symptomatic_26 = as.numeric(rowSums(dfRes_tab[,covid_yesnos],na.rm = T)>0)

# Run on only positives
table_one_pos <- OverReact::crossTabMulti(dat = dfRes_tab %>% filter( estbinres==1),rowvar_list =  
                                                 c("symptomatic_26", covid_yesnos),
                                          colvar = "variant_inferred_detail",
                                          cov_names = cov_name_list,confint = T,include_percentages = T,
                                          rowwise_precentages = F,statistical_test = F) 

# table_one_pos_plot <- OverReact::makeXtabPlottable(myxtab = table_one_pos_plot)

table_one_pos


table_one_pos <- table_one_pos%>% 
  filter(Category!=0) %>% 
  select(-Category, -Total) %>% 
  pivot_longer(cols = -Variable)

table_one_pos_plot <- table_one_pos

# convert % to numeric
table_one_pos_plot$lower <- as.numeric(OverReact::xtabPercentageExtractor(mystring = 
                                                                                 table_one_pos_plot$value,
                                                                               lookbehind = "\\[",
                                                                          lookahead = "[-]"))
table_one_pos_plot$upper <- as.numeric(OverReact::xtabPercentageExtractor(mystring = 
                                                                            table_one_pos_plot$value,
                                                                          lookbehind = "[-]",
                                                                          lookahead = "\\]"))

table_one_pos_plot$percentage <- as.numeric(OverReact::xtabPercentageExtractor(mystring = 
                                                                                 table_one_pos_plot$value,
                                                                               lookbehind = "\\(",
                                                                               lookahead = "\\%"))



dodge_mult=0.8
dodger <- position_dodge2(width = dodge_mult,reverse = F)

stripes_df <- table_one_pos_plot %>% arrange(-percentage) %>% 
  filter(name=="Rounds 2-7 (Wild type)") %>% 
  left_join(sympnames_type_df, by = c("Variable"= "symptom"))%>% 
  group_by(symptom_type) %>% 
  mutate(stripe_col = case_when(row_number() %%2 ==0 ~ "NA", TRUE ~ "grey90")
  )


### Get summary df to fix column widths
summ <- table_one_pos_plot %>% 
  left_join(sympnames_type_df, by = c("Variable"= "symptom")) %>% 
  group_by(symptom_type) %>% 
  summarise(n=n())

# join
sympnames_type_df <- sympnames_type_df %>% left_join(summ)
# set col width
sympnames_type_df$col_width=dodge_mult*sympnames_type_df$n/max(sympnames_type_df$n)

# PLot
p_prevs=table_one_pos_plot %>% 
  left_join(sympnames_type_df, by = c("Variable"= "symptom")) %>% 
  mutate(name = factor(name, levels = unique(name)),
         symptom_type = factor(symptom_type, levels = c("Overall",  "Smell/taste","Respiratory", 
                                                "Coryzal", "Gastrointestinal","Fatigue", "Other"))) %>% 
  ggplot() +
  # geom_tile(data = stripes_df , 
  #           aes(x=reorder(Variable,percentage), y =1, height = Inf,
  #                                  fill = stripe_col),
  #           alpha=0.8,
  #           col = "grey70",
  #           linetype = "dashed",
  #           show.legend = F) +
  # scale_fill_manual(values = c("white","grey96")) +
  # ggnewscale::new_scale_fill() +
  geom_bar(aes(x=reorder(Variable,percentage), y= percentage, fill = name, width=col_width),
           stat = "identity",
           position = dodger, col = "black", size = 0.01) +
  geom_errorbar(position = dodger,
                aes(x=reorder(Variable,percentage),ymin= lower, ymax=upper, width=col_width),
                size=0.3) +
  theme_react(strip_text_size = 10) +
  scale_x_discrete(labels = function(x) str_wrap(x, width =17)) +
  scale_y_continuous(breaks = scales::breaks_width(10)) +
  # scale_fill_brewer(palette = "Reds") +
  OverReact::scale_fill_imperial(palette = "default") +
  # coord_flip() +
  # facet_grid(scales = "free",rows = "symptom_type",space = "fixed",shrink = T) +
  ggforce::facet_col(facets = "symptom_type",scales = "free",space = "free") +
  labs(x="", y="% of PCR positive respondents with symptom in past week", fill = "") +
  theme(legend.position = "bottom",
        panel.spacing = unit(1, "cm"))
p_prevs


# save
OverReact::saveREACTplot(p = p_prevs,figpath = figpath,filename = "symptom_prevalence_by_variant",
                         width = 7.5,height = 10)




# As above for negatives --------------------------------------------------


# Run on only negitives
table_one_neg_plot <- OverReact::crossTabMulti(dat = dfRes_tab %>% filter( estbinres==0),rowvar_list =  
                                                 c("symptomatic_26", covid_yesnos),
                                               colvar = "wave",
                                               cov_names = cov_name_list,confint = T,include_percentages = T,
                                               rowwise_precentages = F,statistical_test = F) 

# table_one_neg_plot <- OverReact::makeXtabPlottable(myxtab = table_one_neg_plot)

table_one_neg_plot


table_one_neg_plot <- table_one_neg_plot%>% 
  filter(Category!=0) %>% 
  select(-Category, -Total) %>% 
  pivot_longer(cols = -Variable)

# convert % to numeric
table_one_neg_plot$lower <- as.numeric(OverReact::xtabPercentageExtractor(mystring = 
                                                                            table_one_neg_plot$value,
                                                                          lookbehind = "\\[",
                                                                          lookahead = "[-]"))
table_one_neg_plot$upper <- as.numeric(OverReact::xtabPercentageExtractor(mystring = 
                                                                            table_one_neg_plot$value,
                                                                          lookbehind = "[-]",
                                                                          lookahead = "\\]"))

table_one_neg_plot$percentage <- as.numeric(OverReact::xtabPercentageExtractor(mystring = 
                                                                                 table_one_neg_plot$value,
                                                                               lookbehind = "\\(",
                                                                               lookahead = "\\%"))



dodge_mult=0.8
dodger <- position_dodge2(width = dodge_mult,reverse = F)

stripes_df <- table_one_neg_plot %>% arrange(-percentage) %>% 
  filter(name=="Rounds 2-7 (Wild type)") %>% 
  left_join(sympnames_type_df, by = c("Variable"= "symptom"))%>% 
  group_by(symptom_type) %>% 
  mutate(stripe_col = case_when(row_number() %%2 ==0 ~ "NA", TRUE ~ "grey90")
  )


### Get summary df to fix column widths
summ <- table_one_neg_plot %>% 
  left_join(sympnames_type_df, by = c("Variable"= "symptom")) %>% 
  group_by(symptom_type) %>% 
  summarise(n=n())

# join
sympnames_type_df <- sympnames_type_df %>% left_join(summ)
# set col width
sympnames_type_df$col_width=dodge_mult*sympnames_type_df$n/max(sympnames_type_df$n)

# PLot
p_prevs_neg=table_one_neg_plot %>% 
  left_join(sympnames_type_df, by = c("Variable"= "symptom")) %>% 
  mutate(name = factor(name, levels = unique(name)),
         symptom_type = factor(symptom_type, levels = c("Overall",  "Smell/taste","Respiratory", 
                                                        "Coryzal", "Gastrointestinal","Fatigue", "Other"))) %>% 
  ggplot() +
  # geom_tile(data = stripes_df , 
  #           aes(x=reorder(Variable,percentage), y =1, height = Inf,
  #                                  fill = stripe_col),
  #           alpha=0.8,
  #           col = "grey70",
  #           linetype = "dashed",
  #           show.legend = F) +
  # scale_fill_manual(values = c("white","grey96")) +
  # ggnewscale::new_scale_fill() +
  geom_bar(aes(x=reorder(Variable,percentage), y= percentage, fill = name, width=col_width),
           stat = "identity",
           position = dodger, col = "black", size = 0.01) +
  geom_errorbar(position = dodger,
                aes(x=reorder(Variable,percentage),ymin= lower, ymax=upper, width=col_width),
                size=0.3) +
  theme_react(strip_text_size = 10) +
  scale_x_discrete(labels = function(x) str_wrap(x, width =17)) +
  scale_y_continuous(breaks = scales::breaks_width(10)) +
  # scale_fill_brewer(palette = "Reds") +
  OverReact::scale_fill_imperial(palette = "default") +
  # coord_flip() +
  # facet_grid(scales = "free",rows = "symptom_type",space = "fixed",shrink = T) +
  ggforce::facet_col(facets = "symptom_type",scales = "free",space = "free") +
  labs(x="", y="% of PCR negative respondents with symptom in past week", fill = "") +
  theme(legend.position = "bottom",
        panel.spacing = unit(1, "cm"))
p_prevs_neg


# save
OverReact::saveREACTplot(p = p_prevs_neg,figpath = figpath,
                         filename = "symptom_prevalence_by_variant_pcr_neg",
                         width = 7.5,height = 10)




# Plot comparing severity -------------------------------------------------

# get plot dat to save compute
plotdat_severity=dfRes_tab %>% filter(!is.na(covidability_cat),
                                      !is.na(variant_inferred_detail),
                                      estbinres==1)

p_sev <-  plotdat_severity %>% 
  group_by(variant_inferred_detail,covidability_cat) %>% 
  summarise(n=n()) %>% 
  group_by(variant_inferred_detail) %>% 
  mutate(percent=100*n/sum(n)) %>% 
    ggplot(aes(x=variant_inferred_detail, y=percent, fill=covidability_cat)) +
  geom_col(position="dodge") +
  OverReact::theme_react() +
  scale_fill_brewer(palette = "Reds",direction = -1) +
  labs(x="", "%")
p_sev






# First symptom coocurrence -----------------------------------------------

