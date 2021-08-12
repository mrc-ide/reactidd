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

#' Pull in packages needed
package.list <- c("prevalence","mgcv","knitr","MASS","kableExtra","table1","dplyr","factoextra","tableone","networkD3",
                  "tidyr","forcats", "cluster", "fpc", "mclust", "pheatmap","FactoMineR", "NbClust","clValid","plotly",
                  "ggplot2","ggsci", "RColorBrewer", "tidyverse", "lubridate", "egg", "poLCA", "Rcpp","xml2","splitstackshape",
                  "fs", "later", "promises","proxy","dendextend", "ComplexHeatmap","circlize","doSNOW","ClusterR","htmlwidgets",
                  "readr","ggthemes", "questionr", "gridExtra", "foreach", "doParallel", 
                  "purrr", "httr", "htmltools","ggalluvial","datapasta","xgboost","SHAPforxgboost"
)

load_packages(package.list)



############################ Run data import script ########################
r8_include_switch=F
source("E:/Group/react2_study5/report_phases_combined/projects/react_1_antigen_symptom_prediction/code/00_data_prep.R")
############################ Run data import script ########################




# Responses by date  ---------------------------------------------------------------------



### Initial plot of response dates by round
p.date <- mydata_symp %>% 
  filter(d_comb >  "2020-05-01") %>% 
  ggplot(aes(x=d_comb, fill =factor(round))) + 
  geom_histogram(col = "black") + 
  theme_bw() +
  labs(x = "Date of response", y = "Count", fill = "Round") +
  scale_fill_brewer() +
  scale_x_date(date_labels = c("%d %b"), date_breaks = "2 weeks") +
  theme(axis.text.x = element_text(angle =45, hjust = 0.5, vjust = 0.5))
ggsave(filename = paste0(figpath, "date_of_swab_by_round.png"), 
       plot = p.date,width = 10, height=7, dpi = 300, units = "in")



# Symptomaticness by age --------------------------------------------------

### Proportions among all respondents
age_symp_tab <- prop.table(table(mydata_symp$symptomatic, mydata_symp$age_group_pred), 2) %>% as.data.frame.matrix()
names(age_symp_tab) <- c("age_5_17", "age_18_54", "age_55_plus")
write_csv(age_symp_tab,paste0(outpath,"age_symptomatic_xtab_prop.csv"))

### Counts among all respondents
age_symp_tab <- (table(mydata_symp$symptomatic, mydata_symp$age_group_pred)) %>% as.data.frame.matrix()
names(age_symp_tab) <- c("age_5_17", "age_18_54", "age_55_plus")
write_csv(age_symp_tab,paste0(outpath,"age_symptomatic_xtab.csv"))

### Proportions among positive
mydata_symp_pos <- mydata_symp %>% filter(estbinres==1)
age_symp_tab <- prop.table(table(mydata_symp_pos$symptomatic, mydata_symp_pos$age_group_pred), 2) %>% as.data.frame.matrix()
names(age_symp_tab) <- c("age_5_17", "age_18_54", "age_55_plus")
write_csv(age_symp_tab,paste0(outpath,"age_symptomatic_xtab_positive.csv"))


# Master symptom prevalence table -----------------------------------------

###' Master table one: symptom prevalence stratified by age and positivity


### Create function to create table so we can do it on R8 after
createMasterTable <- function(dat){
  agenames <- c("5-17","from_18_to_54","over_55")
  age_strat_results_list <- list()
  for (i in 1:3){
    print(i)
    dat_temp <- dat %>% filter(age_bin==i)
    df_all_age <- crossTabMulti(symps = c("one_of_four",covid_yesnos),
                                sympnames = c("one_of_four", react1_sympnames[1:26]), 
                                dat = dat_temp,
                                var1 = "estbinres",  
                                catscovs = NULL,
                                var1_levels = c("Negative", "Positive"),
                                var2_levels = c("No", "Yes"))
    names(df_all_age) <- c("Symptom", paste0(c("Negative_","Positive_"),agenames[[i]] ))
    age_strat_results_list[[i]] <- df_all_age
  }
  
  ### join DFs
  age_strat_results_df <- age_strat_results_list %>% purrr::reduce(left_join)
  
  #### Get top level number
  df <- table(dat$estbinres, dat$age_bin) %>% as.data.frame.matrix()
  df_new <- cbind(df[1,1],df[2,1],df[1,2],df[2,2],df[1,3],df[2,3])
  
  ### Join
  age_strat_results_df <- rbind(c("Full cohort",df_new),age_strat_results_df)
  age_strat_results_df
}

#### subset data frames
dfRes <- dfRes %>% filter(!is.na(estbinres), !is.na(age), !is.na(gender)) %>% 
  select(estbinres,all_of(covid_yesnos), age, gender, round, interim_label, age_group_pred, one_of_four) %>% 
  mutate(age_bin = case_when(age <18 ~ 1,
                             age <55 ~ 2,
                             age >=55 ~ 3,
                             TRUE ~ NA_real_))
dfRes_r8 <- dfRes_r8 %>% filter(!is.na(estbinres), !is.na(age), !is.na(gender)) %>% 
  select(estbinres,all_of(covid_yesnos), age, gender, round, interim_label, age_group_pred, one_of_four) %>% 
  mutate(age_bin = case_when(age <18 ~ 1,
                             age <55 ~ 2,
                             age >=55 ~ 3,
                             TRUE ~ NA_real_))

### Run tables
age_strat_results_df_2_7 <- createMasterTable(dfRes)
age_strat_results_df_8 <- createMasterTable(dfRes_r8)

write.csv(age_strat_results_df_2_7, paste0(outpath, "symptom_prev_age_positivity_2_7.csv"))
write.csv(age_strat_results_df_8, paste0(outpath, "symptom_prev_age_positivity_8.csv"))





# Distribution of symptom counts ------------------------------------------

varnames_classic <- varnames[c(1:4)]
mydata$classic_count <- rowSums(mydata[, varnames_classic])
mydata$allsymp_count <- rowSums(mydata[, varnames])


pdist <- mydata %>% pivot_longer(c(classic_count,allsymp_count)) %>% 
  mutate(Result=ifelse(estbinres==1, "Positive", "Negative"),
         symptoms = ifelse(name=="classic_count", "Classic four symptoms", "All 26 symptoms")) %>% 
  ggplot(aes(x=value,  fill = factor(Result))) +
  geom_bar(stat = "count") +
  theme_bw() +
  scale_colour_manual(values = myCols[c(2,8)]) +
  scale_fill_manual(values = myCols[c(2,8)]) +
  facet_wrap(.~symptoms*Result, scales = "free") +
  labs(x="Number of symptoms", "Count of participants") +
  theme(legend.position = "none",
        # axis.text.x = element_blank(),
              # axis.title.x = element_blank(),
              strip.background = element_rect(fill="white"),
              strip.text = element_text(face = "bold")
        )
pdist


ggsave(filename = paste0(figpath, "symptom_count_distributions.png"), 
       plot = pdist,width = 8, height=10, dpi = 300, units = "in")

### Save csvs with equivalent data
df_dist_classic <- mydata %>% group_by(estbinres) %>% count(classic_count)
df_dist_all <- mydata %>% group_by(estbinres) %>% count(allsymp_count)


write.csv(df_dist_classic,paste0(outpath, "symp_count_classic.csv"))
write.csv(df_dist_all,paste0(outpath, "symp_count_all.csv"))







# Cross tabs --------------------------------------------------------------


# cross tabs among all respondents

### Sex
df_all_sex <- crossTabMulti(symps = covid_yesnos, 
                            sympnames = react1_sympnames[1:26], 
                            dat = dfRes,
                            var1 = "gender",  
                            catscovs = NULL,
                            var1_levels = c("Male", "Female"),
                            var2_levels = c("No", "Yes"))

### Age
df_all_age <- crossTabMulti(symps = covid_yesnos, 
                            sympnames = react1_sympnames[1:26], 
                            dat = dfRes,
                            var1 = "age_group_r1",  
                            catscovs = NULL,
                            var1_levels = c("<18",cat_cov_list$age_group$Levels),
                            var2_levels = c("No", "Yes"))

### Age
df_all_round <- crossTabMulti(symps = covid_yesnos, 
                              sympnames = react1_sympnames[1:26], 
                              dat = dfRes,
                              var1 = "round",  
                              catscovs = NULL,
                              var1_levels = c(2:7),
                              var2_levels = c("No", "Yes"))

### All
# create_dummy var
dfRes$dummy <- 1
df_all <- crossTabMulti(symps = covid_yesnos, 
                        sympnames = react1_sympnames[1:26], 
                        dat = dfRes,
                        var1 = "dummy",  
                        catscovs = NULL,
                        var1_levels = "Full cohort",
                        var2_levels = c("No", "Yes"))


write.csv(df_all_sex, paste0(outpath, "symptom_prev_sex.csv"))
write.csv(df_all_age, paste0(outpath, "symptom_prev_age.csv"))
write.csv(df_all_round, paste0(outpath, "symptom_prev_round.csv"))
write.csv(df_all, paste0(outpath, "symptom_prev.csv"))








# Cross tabs among pcr positive only --------------------------------------

dfRes_pos <- dfRes %>% filter(estbinres==1)

### Sex
df_pos_sex <- crossTabMulti(symps = covid_yesnos, 
                            sympnames = react1_sympnames[1:26], 
                            dat = dfRes_pos,
                            var1 = "gender",  
                            catscovs = NULL,
                            var1_levels = c("Male", "Female"),
                            var2_levels = c("No", "Yes"))

### Age
df_pos_age <- crossTabMulti(symps = covid_yesnos, 
                            sympnames = react1_sympnames[1:26], 
                            dat = dfRes_pos,
                            var1 = "age_group_r1",  
                            catscovs = NULL,
                            var1_levels = c("<18",cat_cov_list$age_group$Levels),
                            var2_levels = c("No", "Yes"))


### Age
df_pos_round <- crossTabMulti(symps = covid_yesnos, 
                              sympnames = react1_sympnames[1:26], 
                              dat = dfRes_pos,
                              var1 = "round",  
                              catscovs = NULL,
                              var1_levels = c(2:7),
                              var2_levels = c("No", "Yes"))

### All
# create_dummy var
dfRes_pos$dummy <- 1
df_pos_all <- crossTabMulti(symps = covid_yesnos, 
                        sympnames = react1_sympnames[1:26], 
                        dat = dfRes_pos,
                        var1 = "dummy",  
                        catscovs = NULL,
                        var1_levels = "Full cohort",
                        var2_levels = c("No", "Yes"))



write.csv(df_pos_sex, paste0(outpath, "symptom_prev_pos_sex.csv"))
write.csv(df_pos_age, paste0(outpath, "symptom_prev_pos_age.csv"))
write.csv(df_pos_round, paste0(outpath, "symptom_prev_pos_round.csv"))
write.csv(df_pos_all, paste0(outpath, "symptom_prev_pos.csv"))



# Prev among neg ----------------------------------------------------------


dfRes_neg <- dfRes %>% filter(estbinres==0)

### Sex
df_neg_sex <- crossTabMulti(symps = covid_yesnos, 
                            sympnames = react1_sympnames[1:26], 
                            dat = dfRes_neg,
                            var1 = "gender",  
                            catscovs = NULL,
                            var1_levels = c("Male", "Female"),
                            var2_levels = c("No", "Yes"))

### Age
df_neg_age <- crossTabMulti(symps = covid_yesnos, 
                            sympnames = react1_sympnames[1:26], 
                            dat = dfRes_neg,
                            var1 = "age_group_r1",  
                            catscovs = NULL,
                            var1_levels = c("<18",cat_cov_list$age_group$Levels),
                            var2_levels = c("No", "Yes"))


### Age
df_neg_round <- crossTabMulti(symps = covid_yesnos, 
                              sympnames = react1_sympnames[1:26], 
                              dat = dfRes_neg,
                              var1 = "round",  
                              catscovs = NULL,
                              var1_levels = c(2:7),
                              var2_levels = c("No", "Yes"))

### All
# create_dummy var
dfRes_neg$dummy <- 1
df_neg_all <- crossTabMulti(symps = covid_yesnos, 
                            sympnames = react1_sympnames[1:26], 
                            dat = dfRes_neg,
                            var1 = "dummy",  
                            catscovs = NULL,
                            var1_levels = "Full cohort",
                            var2_levels = c("No", "Yes"))


write.csv(df_neg_sex, paste0(outpath, "symptom_prev_neg_sex.csv"))
write.csv(df_neg_age, paste0(outpath, "symptom_prev_neg_age.csv"))
write.csv(df_neg_round, paste0(outpath, "symptom_prev_neg_round.csv"))
write.csv(df_neg_all, paste0(outpath, "symptom_prev_neg.csv"))


# Plot distribution of responses ------------------------------------------

table(df.comb$interim_label)
df.comb %>% 
  mutate(date = case_when(!is.na(d_swab) ~ d_swab,
                          TRUE ~ d_col)) %>% 
  filter(date > "2020-06-01") %>% 
  ggplot(aes(x=date, fill = factor(round))) + geom_histogram() +
  theme_bw() +
  scale_x_date(breaks = "1 month", date_minor_breaks = "1 week", date_labels = "%B") +
  labs(fill = "Round", x="Date of swab")


# Prevalence by region ----------------------------------------------------



### add symptomatic variable
df.comb$symptomatic <- as.numeric(rowSums(df.comb[,covid_yesnos], na.rm = T) > 0)
prop.table(table(df.comb$symptomatic ))


preg <- dfRes %>% 
  filter(!is.na(estbinres)) %>%
  group_by(region, round, estbinres) %>% 
  summarise(n = n(),
            perc_symptomatic = (sum(symptomatic, na.rm = T)/ n),
            se_prop = sqrt((perc_symptomatic*(1-perc_symptomatic))/n),
            lower = max(0,perc_symptomatic-(1.96*se_prop)),
            upper = min(1, perc_symptomatic+(1.96*se_prop))) %>% 
  # mutate(perc_symptomatic = 100*(symptomatic)/nrow(dfRes[dfRes$estbinres == 1,])) %>% 
  ggplot(aes(x=round, y= perc_symptomatic, col=factor(estbinres))) + 
  geom_point() +
  geom_errorbar(aes(ymin= lower, ymax =upper)) +
  ylim(c(0,1)) +
  facet_wrap(.~region) +
  theme_bw() + labs(y = "Proportion with symptom", x="Round",
                    col = "Test \rresult")

ggsave(filename = paste0(figpath, "symptom_prevalence_by_region.png"), 
       plot = preg,width = 10, height=10, dpi = 300, units = "in")


colnames(dfRes)[13:38] <- tolower(react1_sympnames[1:26])

p <- dfRes%>% 
  filter(!is.na(estbinres)) %>%
  pivot_longer(cols= tolower(react1_sympnames[1:26])) %>% 
  group_by(name,round,estbinres) %>% 
  summarise(n = n(),
            perc_symptomatic = (sum(value, na.rm = T)/ n),
            se_prop = sqrt((perc_symptomatic*(1-perc_symptomatic))/n),
            lower = max(0,perc_symptomatic-(1.96*se_prop)),
            upper = min(1, perc_symptomatic+(1.96*se_prop))) %>% 
  # mutate(perc_symptomatic = 100*(symptomatic)/nrow(dfRes[dfRes$estbinres == 1,])) %>% 
  ggplot(aes(x=round, y= perc_symptomatic, col =factor(estbinres))) + 
  geom_point() +
  geom_errorbar(aes(ymin= lower, ymax =upper)) +
  ylim(c(0,0.3)) +
  facet_wrap(.~name) +
  theme_bw() + labs(y = "Proportion with symptom", x="Round",
                    col = "Test \rresult")
p
ggsave(filename = paste0(figpath, "symptom_prevalence_over_time.png"), 
       plot = p,width = 15, height=12, dpi = 300, units = "in")



# Pairwise PCR prevalence -------------------------------------------------

df.pairwise.prevs <- data.frame(matrix(nrow=26,ncol=26, dimnames = list(react1_sympnames[c(1:26)],
                                                                        react1_sympnames[c(1:26)])))

df.pairwise.ORs <- data.frame(matrix(nrow=26,ncol=26, dimnames = list(react1_sympnames[c(1:26)],
                                                                      react1_sympnames[c(1:26)])))
 
for(i in 1:26){
  print(i)
  for(j in 1:26){
    var1=covid_yesnos[[i]]
    var2=covid_yesnos[[j]]
    prev=dfRes %>% summarise(n=n(),
                             npos=sum(dfRes[,var1]*dfRes[,var2]),
                             prev=npos/n)
    df.pairwise.prevs[i,j] <- prev$prev*100
    
    ### Model
    dfRes$newvar <- as.numeric(unlist(dfRes[,var1] * dfRes[,var2]))
    # f <- as.formula(paste("estbinres ~ (", paste0(c(var1,var2), collapse = " + "),")^2"))
    f <- as.formula("estbinres ~ newvar")
    mod <- glm(formula = f, data = dfRes, family = "binomial")
    modsumm <- jtools::summ(mod)
    coefs <- modsumm$coeftable
    or <- exp(coefs[2,1])
    df.pairwise.ORs[i,j] <- or
  }
}


### plot heatmap
p <- ComplexHeatmap::Heatmap(df.pairwise.prevs)
p
p2 <- ComplexHeatmap::Heatmap(df.pairwise.ORs,name = "Odds \nratio")
p2
table(dfRes$newvar)
bu <- df.pairwise.ORs

