setwd("E:/Group/react2_study5/report_phases_combined/projects/react_1_antigen_symptom_prediction/")

source("E:/Group/functions/full_join_multiple_datasets.R", local = T)
source("E:/Group/functions/load_packages.R", local = T)

source("E:/Group/functions/relevel_all_categorical_vars.R")
source("E:/Group/functions/wrangle_cleaner_functions.R", local = T)
source("E:/Group/functions/cats_and_covs.R", local = T)

source("E:/Group/functions/full_join_multiple_datasets.R", local = T)
source("E:/Group/react2_study5/report_phases_combined/projects/clustering/code/00_functions.R")


# Pull in packages needed
package.list <- c("prevalence","mgcv","knitr","MASS","kableExtra","table1","dplyr","factoextra","tableone","networkD3",
                  "tidyr","forcats", "cluster", "fpc", "mclust", "pheatmap","FactoMineR", "NbClust","clValid","plotly",
                  "ggplot2","ggsci", "RColorBrewer", "tidyverse", "lubridate", "egg", "poLCA", "Rcpp","xml2","splitstackshape",
                  "fs", "later", "promises","proxy","dendextend", "ComplexHeatmap","circlize","doSNOW","ClusterR","htmlwidgets",
                  "readr","ggthemes", "questionr", "gridExtra", "foreach", "doParallel",
                  "purrr", "httr", "htmltools","ggalluvial","datapasta","xgboost","SHAPforxgboost"
)

load_packages(package.list)



# Import REACT-1 data -----------------------------------------------------

r1 <- readRDS("E:/Group/react2_study5/saved_objects/react_1_saved_objects/rep1.rds")
r2 <- readRDS("E:/Group/react2_study5/saved_objects/react_1_saved_objects/rep2.rds")
r3 <- readRDS("E:/Group/react2_study5/saved_objects/react_1_saved_objects/rep3.rds")
r4 <- readRDS("E:/Group/react2_study5/saved_objects/react_1_saved_objects/rep4.rds")
r5 <- readRDS("E:/Group/react2_study5/saved_objects/react_1_saved_objects/rep5.rds")
r6 <- readRDS("E:/Group/react2_study5/saved_objects/react_1_saved_objects/rep6.rds")
r7 <- readRDS("E:/Group/react2_study5/saved_objects/react_1_saved_objects/rep7.rds")
r8 <- readRDS("E:/Group/react2_study5/saved_objects/react_1_saved_objects/rep8.rds")
#' We will exclude R1 because it doesn't include the 'do you have symptoms now?' questions
#' 
#
#' #' bind together to create data set for both rounds
dfRes   <- bind_rows(r1,r2,r3,r4,r5,r6,r7,r8, .id = "round") %>%
  mutate(round = as.factor(round))
### filter r1
dfRes <- dfRes %>% filter(round != 1)

# df_all <- joinMyDatasets(list(r2,r3,r4,r5,r6,r7, r8))
rm(r1,r2,r3,r4,r5,r6,r7,r8)

prop.table(table(dfRes$estbinres, exclude = "none"))
table(dfRes$round)

### Exclude the missing esbinres data as we will never use this
dfRes <- dfRes %>% filter(!is.na(estbinres))

# Wrangle -----------------------------------------------------------------

### This segment of the code should be insertable into all the REACT-1 wrangling scripts ###

### define lists of variable names
covid_yesnos <- dfRes %>% select(symptnowaw_1:symptnowaw_26) %>% colnames()
covid_yesnos_1month <- dfRes %>% select(sympt_any1_1:sympt_any3_9) %>% colnames()

### Replace NA with 0 in symptoms
dfRes[covid_yesnos][is.na(dfRes[covid_yesnos])] <- 0

### clean up the data
dfRes <- dfRes %>% 
  mutate_at(all_of(covid_yesnos),binaryCleaner_1_0)  %>% 
  mutate_at(c("gender","age", "ethnic","smokenow", "smokecig", 
              "empl", "educ","covidcon"),
            continuousCleaner)  %>% 
  mutate(age_group_pred = case_when(age <18 ~ 1,
                                    age <=54 ~ 2,
                                    age >54  ~ 3,
                                    TRUE ~NA_real_),
         round_split = case_when(round %in% 2:6 ~ 0,
                                 round %in% 7:8 ~ 1),
         round_split_8_all = case_when(round %in% 2:7 ~ 0,
                                       round  == 8 ~ 1))

### Add a variable for any one of the 26 reported symptoms
dfRes$symptomatic <- as.numeric(rowSums(dfRes[,covid_yesnos]) > 0)

### Add variable for one of four classical symptoms
dfRes$one_of_four <- as.numeric(rowSums(dfRes[, covid_yesnos[c(1:4)]]) > 0)


### Wrangling segment ends ###


# Modelling symptoms only -------------------------------------------------------

### create predictor and output data
# take r8_switch from parent file
if(!r8_include_switch){
  dfRes_r8 <- dfRes %>% filter(round == 8)
  dfRes <- dfRes %>% filter(round != 8)
}

exclusions_gender <- sum(is.na(dfRes$gender))
exclusions_age <- sum(is.na(dfRes$age))

print(paste("There were",exclusions_gender, "exclusions for missing data in gender"))
print(paste("There were",exclusions_age, "exclusions for missing data in age"))

mydata_symp <- dfRes %>% filter(!is.na(age), !is.na(gender)) %>% 
  select(estbinres,all_of(covid_yesnos), age, gender, round, interim_label, ethnic_new,
         smokenow, u_passcode, symptomatic,age_group_pred,d_comb,
         one_of_four, round_split, round_split_8_all) 

### add variable names
colnames(mydata_symp)[2:27] <- react1_sympnames[1:26]

### clean with Janitor
mydata_symp <- janitor::clean_names(mydata_symp)
varnames <- names(mydata_symp)[2:27]

### create predictor and output data
mydata <- mydata_symp %>% filter(symptomatic == 1) 
mydata$new_id <- 1:nrow(mydata)

### Create test-train split
set.seed(123)
ttsplit <- caret::createDataPartition( y = mydata$estbinres,times = 1,p = 0.7,list = F)
mydata_train <- mydata[ttsplit,]
mydata_test <- mydata[-ttsplit,]


### Wrangle R8 separately
if(!r8_include_switch){

mydata_symp_r8 <- dfRes_r8 %>% filter(!is.na(age), !is.na(gender)) %>% 
  select(estbinres,all_of(covid_yesnos), age, gender, round, interim_label, ethnic_new,
         smokenow, u_passcode, symptomatic,age_group_pred,d_comb,
         one_of_four, round_split, round_split_8_all) 

### Add indices
ind <- 1+nrow(mydata)
mydata_symp_r8$new_id <- ind:(ind+nrow(mydata_symp_r8)-1)

### add variable names
colnames(mydata_symp_r8)[2:27] <- react1_sympnames[1:26]

### clean with Janitor
mydata_symp_r8 <- janitor::clean_names(mydata_symp_r8)

### create predictor and output data
mydata_r8 <- mydata_symp_r8 %>% filter(symptomatic == 1) 

}

