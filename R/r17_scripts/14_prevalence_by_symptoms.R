# Script to produce the table of prevalences
# Outputs: 1 xlsx file with overall prevalence in the round
# under "Overall_prevalence" with round ID and date;
# and 1 xlsx file with prevalence by category
# under "Prevalence_weighted" with round ID and date
# Files are stored in "E:/Group/report/round17/Tables"
# and automatically copied to the transfer folder

# To change age categories:
# Step 1: change the thresholds in recoding_from_cont_file (see path below)
# Step 2: manually update the template in template_file (see path below)


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

# Age bounds
age_lower_bound <- 17
age_upper_bound <- 200

# Modelling options
weighted <- FALSE # whether weighted prevalences should be included or not

# Copying output files directly to transfer folder
direct_export <- TRUE

# Paths to files
overall_file <- "Tables/Overall_prevalence"
output_file <- "Tables/Prevalence_by_symptoms"
output_file=paste0(output_file, "_", age_lower_bound, "_", age_upper_bound)

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
  "sympt_cat",
  paste0("symptnowaw_", 1:26)
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

# Adding variable introduced in round 15
if (!"vax_status_noDate_v2" %in% colnames(df_round)) {
  df_round$vax_status_noDate_v2 <- df_round$vax_status_noDate
}

# Removing missing in estbinres
df_round <- df_round %>%
  filter(!is.na(estbinres)) %>%
  mutate(group = "Overall")

if (weighted) {
  # Removing missing in weights
  df_round <- df_round %>% filter(!is.na(wt_antigen))
}

df_round <- df_round %>%
  mutate(vax_status_cat = ifelse(is.na(vax_status_cat), "NA", vax_status_cat)) %>%
  # mutate(
  #   vax_wane = ifelse(is.na(vax_wane), "NA", vax_wane),
  #   rm_dip = ifelse(rm_dip==-1, "NA", rm_dip),
  #   rm_dip2 = case_when(rm_dip == 1 ~ "1",
  #                       rm_dip == 2 ~ "2",
  #                       rm_dip %in% c(3:12)  ~ "3+",
  #                       rm_dip == "NA" ~ "NA")) %>%
  mutate(
    # vax_status_cat = factor(vax_status_cat, levels = c("Not vaccinated", "One does", "Two does",
    #                                                    "Unknown does", "NA")),
    # vax_wane = factor(vax_wane, levels = c("Unvaccinated", "1 dose", "2 dose < 3 months",  "2 dose 3-6 months",
    #                                        "2 dose > 6 months",  "NA")),
    vax_status_noDate = factor(vax_status_noDate, levels = c(
      "Not vaccinated", "One does", "Two does",
      "Unknown does", "NA"
    )),
    vax_status_noDate_v2 = factor(vax_status_noDate_v2, levels = c(
      "Not vaccinated", "One does", "Two does",
      "Three does", "Unknown does", "NA"
    ))
  ) %>%
  mutate(covidcon_char = ifelse(is.na(covidcon_char), "NA", as.character(covidcon_char))) %>%
  mutate(covidcon_char = factor(covidcon_char,
                                levels = c(
                                  "Yes, contact with confirmed/tested COVID-19 case",
                                  "Yes, contact with suspected COVID-19 case",
                                  "No", "NA"
                                )
  ))


## Prevalence stratified by covariates

# Extracting covariate names
tmp <- read.xlsx(annot_file)
covs_names <- tmp[, 2]
names(covs_names) <- tmp[, 1]
covs_names <- covs_names[covs]

# Stratifying by age
df_round=filter(df_round, 
                age > as.numeric(age_lower_bound) & 
                  age < as.numeric(age_upper_bound))

# Removing unused variables
df_round <- df_round[, c(res_param, covs, weight_params)]

# Recoding categorical variables
covs_to_recode <- getSheetNames(recoding_file)
covs_to_recode <- intersect(names(covs_names), covs_to_recode)
if (length(covs_to_recode) > 0) {
  for (i in 1:length(covs_to_recode)) {
    recoding <- read.xlsx(recoding_file, sheet = covs_to_recode[i])
    recoding[which(is.na(recoding[, 1])), 1] <- "NA"
    renaming <- recoding[, 2]
    names(renaming) <- recoding[, 1]
    x <- as.character(df_round[, covs_to_recode[i]])
    x[is.na(x)] <- "NA"
    x <- factor(x, levels = names(renaming), labels = renaming)
    df_round[, covs_to_recode[i]] <- x
  }
}

# Recoding symptom status
column_ids=grep("symptnowaw", colnames(df_round))
for (k in column_ids){
  if(all(sort(unique(df_round[,k]))==c(0,1))){
    df_round[,k]=as.character(df_round[,k])
    ids=which(is.na(df_round[,k]))
    df_round[ids,k]="Unknown"
    df_round[,k]=factor(df_round[,k], levels=c(1,0,"Unknown"), labels=c("Yes","No","Unknown"))
  } else {
    stop("Not binary.")
  }
}

# Make the prevalence tables for the above covariates using Vivi's code (unweighted)
system.time({
  mytable <- ExtractPrevalence(
    df_round = df_round,
    covs = covs, covs_names = covs_names,
    res_param = res_param, weighted = FALSE
  )
})

if (weighted) {
  system.time({
    mytable_weighted <- ExtractPrevalence(
      df_round = df_round,
      covs = covs, covs_names = covs_names,
      res_param = res_param,
      weight_params = weight_params, weighted = TRUE
    )
  })
  tmp <- mytable
  tmp[rownames(mytable_weighted), 5] <- mytable_weighted[, 3]
  mytable <- cbind(mytable, tmp[, 5])
}

# Saving output table
colnames(mytable)=c("Variable","","Positive","Total","Unweighted Prevalence")
write.xlsx(mytable, paste0(output_file, "_r", round_id, "_", output_tag, ".xlsx"),
           overwrite=TRUE, row.names=FALSE, col.names=TRUE)

# Copying output to transfer folder
if (direct_export) {
  file.copy(
    from = paste0(output_file, "_r", round_id, "_", output_tag, ".xlsx"),
    to = "T:/", overwrite = TRUE
  )
}


