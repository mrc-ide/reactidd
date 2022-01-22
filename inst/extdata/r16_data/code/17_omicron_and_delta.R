rm(list = ls(all = TRUE))
setwd("E:/Group/report/round16/Scripts/")

# Loading required packages
source("functions/load_packages.R")
pkgs <- c(
  "prevalence", "mgcv", "mgcViz", "MASS", "dplyr",
  "tidyr", "forcats", "ggplot2", "qpcR", "survey", "reshape2",
  "openxlsx"
)
load_packages(pkgs)

# Source any functions from the local file
source("functions/add_conf_ints.R")
source("functions/make_tables.R")
source("functions/overall_prev.R")
source("functions/formatting_functions.R")


## Parametrisation

# Modelling options
weighted <- TRUE # whether weighted prevalences should be included or not

# Copying output files directly to transfer folder
direct_export <- TRUE

# Choice of rounds
round_id <- c(16)

# Paths to files
overall_file <- "E:/Group/report/round16/Tables/Overall_prevalence"
output_file <- "E:/Group/report/round16/Tables/Prevalence"

output_tag <- Sys.Date()
annot_file <- "E:/Group/report/round16/Parameters/Variable_names.xlsx"
template_file <- "E:/Group/report/round16/Parameters/Table2a_urban"
template_sheet1 <- "Table2a"
template_sheet2 <- "Table2b"

recoding_file <- "E:/Group/report/round16/Parameters/Recoding.xlsx"
recoding_from_cont_file <- "E:/Group/report/round16/Parameters/Recoding_from_continuous.xlsx"

# Fetching corresponding template
if (weighted) {
  template_file <- paste0(template_file, "_weighted.xlsx")
  output_file <- paste0(output_file, "_weighted")
} else {
  template_file <- paste0(template_file, "_unweighted.xlsx")
  output_file <- paste0(output_file, "_unweighted")
}

# Variable for test results
res_param <- "estbinres"

# Variable for weights
if (weighted) {
  weight_params <- c("id", "lacode", "wt_antigen")
  names(weight_params) <- c("id", "strata", "weights")
} else {
  weight_params <- NULL
}

# Variables for stratification
covs <- c(
  "gender_char", "age", "region",
  "work_new_alt", "ethnic_new_char",
  "hh_size_cat", "covidcon_char", "sympt_cat",
  "nchild2",
  "urban",
  "imd_quintile", "vax_status_noDate_v2"
)

# Defining the column widths / row heights
column_widths <- c(5.5, 5.5, 19)
if (weighted) {
  column_widths <- c(column_widths, 19)
}
row_height <- 15


## Computing the prevalences

# Path to the data
data_file <- paste0("E:/Group/saved_objects/rep", round_id, ".rds")

## Checking parameters

# Checking required sheets are in template
sheet_names <- getSheetNames(template_file)
tocheck <- c(template_sheet1, template_sheet2) %in% sheet_names
if (!all(tocheck)) {
  stop(paste0(
    "Sheets are not found in the template: ",
    c(template_sheet1, template_sheet2)[!tocheck]
  ))
}


## Loading and preparing the data

# Loading the data
df_round <- data.frame(readRDS(data_file))
rownames(df_round)=df_round$u_passcode

df16 <- readRDS("E:/Group/saved_objects/rep16_lineage.rds")
rownames(df16)=df16$u_passcode

ids=df16$u_passcode
df_round=df_round[ids,]

ids_omicron=df16$u_passcode[which(df16$react_lineage=="BA.1")]
table(df16[ids_omicron,"react_lineage"])

ids_delta=df16$u_passcode[which(df16$react_lineage!="BA.1")]
table(df16[ids_delta,"react_lineage"])

ct1_omicron=df_round[ids_omicron, "ct2"]
ct1_omicron=ct1_omicron[which(ct1_omicron!=0)]
mean(ct1_omicron)
sd(ct1_omicron)

ct1_delta=df_round[ids_delta, "ct2"]
ct1_delta=ct1_delta[which(ct1_delta!=0)]
mean(ct1_delta)
sd(ct1_delta)

t.test(ct1_delta, ct1_omicron)

write.xlsx(cbind(ids_omicron), 
           overwrite=TRUE,
           "E:/Group/report/round16/Tables/u_passcode_omicron.xlsx")

write.xlsx(cbind(ids_delta), 
           overwrite=TRUE,
           "E:/Group/report/round16/Tables/u_passcode_delta.xlsx")



