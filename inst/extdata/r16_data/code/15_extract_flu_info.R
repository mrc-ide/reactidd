rm(list = ls(all = TRUE))
setwd("E:/Group/report/round16/Scripts/")

# Loading required packages
source("functions/load_packages.R")
pkgs <- c(
  "prevalence", "mgcv", "mgcViz", "MASS", "dplyr",
  "tidyr", "forcats", "ggplot2", "qpcR", "survey", "reshape2",
  "openxlsx", "ggvenn"
)
load_packages(pkgs)

# Source any functions from the local file
source("functions/add_conf_ints.R")
source("functions/make_tables.R")
source("functions/overall_prev.R")
source("functions/formatting_functions.R")


## Parametrisation

# Paths to files
round_id <- 16
data_file <- paste0("E:/Group/saved_objects/rep", round_id, "_flu.rds")
output_file <- "E:/Group/report/round16/Tables/Prevalence_influenza"
output_tag <- Sys.Date()

annot_file <- "E:/Group/report/round16/Parameters/Variable_names.xlsx"
recoding_file <- "E:/Group/report/round16/Parameters/Recoding.xlsx"
recoding_from_cont_file <- "E:/Group/report/round16/Parameters/Recoding_from_continuous.xlsx"

# Modelling options
weighted <- FALSE # whether weighted prevalences should be included or not

# Copying output files directly to transfer folder
direct_export <- TRUE

# Variable for test results
res_param <- "estbinres"

# Variable for weights
if (weighted) {
  weight_params <- c("id", "lacode", "wt_antigen")
  names(weight_params) <- c("id", "strata", "weights")
} else {
  weight_params <- NULL
}

# Updating output file name
if (weighted) {
  output_file <- paste0(output_file, "_weighted")
} else {
  output_file <- paste0(output_file, "_unweighted")
}

# Variables for stratification
covs <- c(
  "u_passcode",
  "gender_char", "age", "region",
  "work_new_alt", "ethnic_new_char",
  "hh_size_cat", "covidcon_char", "sympt_cat",
  "nchild2",
  "imd_quintile", "vax_status_noDate_v2",
  "influenzaa", "influenzab",
  "fluvacc"
)


## Loading and preparing the data

# Loading the data
df_round <- data.frame(readRDS(data_file))
rownames(df_round)=df_round$u_passcode

table(!is.na(df_round$estbinres), !is.na(df_round$influenzaa))
table(!is.na(df_round$estbinres), !is.na(df_round$influenzab))
table(!is.na(df_round$estbinres), !is.na(df_round$influenzaacpvalue))
table(!is.na(df_round$estbinres), !is.na(df_round$influenzabcpvalue))
table(!is.na(df_round$influenzaa), !is.na(df_round$influenzaacpvalue))

table(df_round$influenzaa, !is.na(df_round$influenzaacpvalue), useNA="always")
table(df_round$influenzab, !is.na(df_round$influenzabcpvalue), useNA="always")

table(df_round$influenzaa, (df_round$influenzaacpvalue==0), useNA="always")
table(df_round$influenzab, (df_round$influenzabcpvalue==0), useNA="always")

table(df_round$influenzaa, df_round$influenzaacpvalue, useNA="always")
table(df_round$influenzab, df_round$influenzabcpvalue, useNA="always")

sheet=1
for (sheet in 1:2){
  mydata=read.xlsx("../Data/FLU Confirmation Results 2021 01 03.xlsx", sheet = sheet)
  myu_passcode=mydata$Sample.Name
  myu_passcode=gsub("UK", "", gsub(" \\(.*", "", myu_passcode))
  extracted=df_round[myu_passcode, c("influenzaa", "influenzab", "influenzaacpvalue", "influenzabcpvalue")]
  mydata=cbind(mydata, extracted)
  mydata=mydata[,3:ncol(mydata)]
  write.xlsx(mydata, file=paste0("../Data/Flu_sheet",sheet,".xlsx"))
}


