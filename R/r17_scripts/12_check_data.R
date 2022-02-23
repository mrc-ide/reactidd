rm(list = ls(all = TRUE))

# Choice of rounds
round_ids <- c(17)
round_id=round_ids[1]

# Setting working directory
setwd(paste0("E:/Group/report/round",round_id))

# Loading required packages
source("Scripts/functions/load_packages.R")
pkgs <- c(
  "prevalence", "mgcv", "mgcViz", "MASS", "dplyr",
  "tidyr", "forcats", "ggplot2", "qpcR", "survey", "reshape2",
  "openxlsx"
)
load_packages(pkgs)

# Source any functions from the local file
source("Scripts/functions/add_conf_ints.R")
source("Scripts/functions/make_tables.R")
source("Scripts/functions/overall_prev.R")
source("Scripts/functions/formatting_functions.R")


## Parametrisation

# Modelling options
weighted <- TRUE # whether weighted prevalences should be included or not

# Copying output files directly to transfer folder
direct_export <- TRUE

# Paths to files
overall_file <- "Tables/Overall_prevalence"
output_file <- "Tables/Prevalence"

output_tag <- Sys.Date()
annot_file <- "Parameters/Variable_names.xlsx"
template_file <- "Parameters/Table2a_urban"
template_sheet1 <- "Table2a"
template_sheet2 <- "Table2b"

recoding_file <- "Parameters/Recoding.xlsx"
recoding_from_cont_file <- "Parameters/Recoding_from_continuous.xlsx"

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

table(df_round$influenzab, df_round$influenzaa, df_round$estbinres, useNA = "always")
table(is.na(df_round$influenzab), is.na(df_round$influenzaa))
table(df_round$influenzab, df_round$influenzaa, df_round$res, useNA = "always")

ids=which(is.na(df_round$estbinres)&(!is.na(df_round$influenzaa)))
tmp=df_round[ids,c("estbinres","res","ct1","ct2","influenzaa", "influenzaacpvalue", "influenzab", "influenzabcpvalue")]
write.xlsx(tmp, "Data/Missing_covid_test_result.xlsx", row.names=TRUE, col.names=TRUE)





positives=read.xlsx("Data/20220113_IPSOS W17 Positive Samples Manifest_Shipment_3.xlsx")
positives$Sample_ID=gsub("UK", "", positives$Sample_ID)

positives=cbind(positives, 
                df_round[positives$Sample_ID, c("u_passcode", "estbinres", "res", "ct2", "ct1", 
                                                "influenzaa", "influenzaacpvalue", "influenzab", "influenzabcpvalue")])

write.xlsx(positives, "Data/20220113_IPSOS W17 Positive Samples Manifest_Shipment_3_comparison.xlsx")

par(mfrow=c(1,2))
plot(positives$`Cp_E-Gene`, positives$ct2, 
     xlim=c(0,50), ylim=c(0,50),
     panel.first=abline(0,1,col="red"))
plot(positives$`Cp_N-Gene`, positives$ct1, 
     xlim=c(0,50), ylim=c(0,50),
     panel.first=abline(0,1,col="red"))

table(positives$estbinres)

ids=which(as.character(df_round$d_comb)%in%names(table(df_round[positives$Sample_ID, "d_comb"])))
table(df_round[ids, "estbinres"])

tmp=positives[which(is.na(positives$estbinres)),]
tmp=tmp[!is.na(tmp$Sample_ID),]

write.xlsx(tmp, "Data/Extracted_missing.xlsx")

ids_missing=read.table("Data/id_missing.csv", sep=",")

table(tmp$Sample_ID%in%ids_missing$V2)

positives[positives$Sample_ID%in%ids_missing$V2,]

write.xlsx(positives[positives$Sample_ID%in%ids_missing$V2,], "Data/Extracted_missing_40.xlsx")
