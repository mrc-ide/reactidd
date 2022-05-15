# Script to produce the table of odds ratios
# Two models: adjusted on age and gender and mutually-adjusted
# Output: 1 xlsx file "Logistic_models" with round ID and date
# stored in "E:/Group/report/round19/Tables/"
# and automatically copied into the transfer folder

# FOR NIC::
# To add new categories:
# Step 1: to add new row in annot_file (see path) - with variable name and pretty name for report
# Step 2: manually update the template in template_file (see path below)
    # - to add in the new variable category and levels (be careful as this script excludes observations with n<5)
# Step 3: add a new sheet with variable category and levels to recoding_file (see path below)
# Step 4: add in new variable name to "covs" (line 60) / or remove those not adjusting
# Step 5: adjust the tidyverse mutates from line 118 and age filter on line 171


rm(list = ls(all = TRUE))

# Choice of round
round_id <- 19

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

# Paths to files
data_file <- paste0("E:/vaccination_data_temp/linkedR", round_id, "datANG.rds")
output_file <- "Tables/Logistic_models_vaxtype_linked"

output_tag <- Sys.Date()
annot_file <- "Parameters/Variable_names.xlsx"
template_file <- "Parameters/OR_table_boosted.xlsx"

recoding_file <- "Parameters/Recoding_routine.xlsx"
recoding_from_cont_file <- "Parameters/Recoding_from_continuous.xlsx"

# Copying output files directly to transfer folder
direct_export <- TRUE

# Variable for test results
res_param <- "estbinres"

# Variables for stratification
covs <- c(
  "gender_char", "age", "region",
  "work_new_alt", "ethnic_new_char",
  "urban",
  "hh_size_cat", "nchild2", "imd_quintile",
  "boosters_time_since_recent_jab",
  "boosters35_time_since_jab"
)

# Adjustment in base model
confounders <- c("gender_char", "age")

# Formatting
CI <- c(" (", ",", ")")


## Loading and preparing the data

df_round <- data.frame(readRDS(data_file))

# Adding variable introduced in round 15
if (!"vax_status_noDate_v2" %in% colnames(df_round)) {
  df_round$vax_status_noDate_v2 <- df_round$vax_status_noDate
}

# Removing missing in estbinres
df_round <- df_round %>%
  filter(!is.na(estbinres)) %>%
  mutate(group = "Overall")

df_round <- df_round %>%
  mutate(vax_status_cat = ifelse(is.na(vax_status_cat), "NA", vax_status_cat)) %>%
  mutate(
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
  )) %>%
  mutate(new_vax_brands = case_when(
    link_vaxdose == 3 & !is.na(link_vax1type) & !is.na(link_vax2type) & !is.na(link_boostertype) ~ paste0(link_vax1type, "/", link_vax2type, "/", link_boostertype),
    TRUE ~ NA_character_
    )) %>%
  mutate(new_vax_brands = as.factor(new_vax_brands))

# overall / 3 and 2 doses and 3 doses aged 35+
df_round <- df_round %>%
  mutate(time_since_recent_jab = case_when(
    link_vaxdose == 3 ~ as.numeric(d_comb) - link_vax3date,
    link_vaxdose == 2 ~ as.numeric(d_comb) - link_vax2date,
    link_vaxdose == 1 ~ as.numeric(d_comb) - link_vax1date,
    link_vaxdose == 0 ~ NA_real_)) %>%
  mutate(time_since_recent_jab = ifelse(time_since_recent_jab <= 30 , "0-30",
                                 ifelse(time_since_recent_jab  > 30 & time_since_recent_jab <= 60, "31-60",
                                 ifelse(time_since_recent_jab  > 60 & time_since_recent_jab <= 90, "61-90",
                                 ifelse(time_since_recent_jab  > 90 & time_since_recent_jab <= 120, "91-120",
                                 ifelse(time_since_recent_jab  > 120 & time_since_recent_jab <= 150, "121-150",
                                 ifelse(time_since_recent_jab  > 150 & time_since_recent_jab <= 180, "151-180",
                                 "180+"))))))) %>%
  mutate(time_since_recent_jab = factor(time_since_recent_jab,
         levels=c("0-30", "31-60", "61-90", "91-120", "121-150", "151-180", "180+"),ordered=FALSE)) %>%
  mutate(boosters_time_since_recent_jab = case_when(
    link_vaxdose == 3 ~ as.numeric(d_comb) - link_vax3date,
    TRUE ~ NA_real_)) %>%
  mutate(boosters_time_since_recent_jab = ifelse(boosters_time_since_recent_jab <= 30 , "0-30",
                                          ifelse(boosters_time_since_recent_jab  > 30 & boosters_time_since_recent_jab <= 60, "31-60",
                                          ifelse(boosters_time_since_recent_jab  > 60 & boosters_time_since_recent_jab <= 90, "61-90",
                                          ifelse(boosters_time_since_recent_jab  > 90 & boosters_time_since_recent_jab <= 120, "91-120",
                                          ifelse(boosters_time_since_recent_jab  > 120 & boosters_time_since_recent_jab <= 150, "121-150",
                                          ifelse(boosters_time_since_recent_jab  > 150 & boosters_time_since_recent_jab <= 180, "151-180",
                                          "180+"))))))) %>%
  mutate(boosters_time_since_recent_jab = factor(boosters_time_since_recent_jab,
         levels=c("0-30", "31-60", "61-90", "91-120", "121-150", "151-180", "180+"),ordered=FALSE)) %>%
  mutate(double_time_since_recent_jab = case_when(
    link_vaxdose == 2 ~ as.numeric(d_comb) - link_vax2date,
    TRUE ~ NA_real_)) %>%
  mutate(double_time_since_recent_jab = ifelse(double_time_since_recent_jab <= 30 , "0-30",
                                        ifelse(double_time_since_recent_jab  > 30 & double_time_since_recent_jab <= 60, "31-60",
                                        ifelse(double_time_since_recent_jab  > 60 & double_time_since_recent_jab <= 90, "61-90",
                                        ifelse(double_time_since_recent_jab  > 90 & double_time_since_recent_jab <= 120, "91-120",
                                        ifelse(double_time_since_recent_jab  > 120 & double_time_since_recent_jab <= 150, "121-150",
                                        ifelse(double_time_since_recent_jab  > 150 & double_time_since_recent_jab <= 180, "151-180",
                                        "180+"))))))) %>%
  mutate(double_time_since_recent_jab = factor(double_time_since_recent_jab,
         levels=c("0-30", "31-60", "61-90", "91-120", "121-150", "151-180", "180+"),ordered=FALSE)) %>%
  mutate(boosters35_time_since_jab = case_when(
    link_vaxdose == 3 & age >= 35 ~ as.numeric(d_comb) - link_vax3date,
    TRUE ~ NA_real_)) %>%
  mutate(boosters35_time_since_jab = ifelse(boosters35_time_since_jab <= 30 , "0-30",
                                     ifelse(boosters35_time_since_jab  > 30 & boosters35_time_since_jab <= 60, "31-60",
                                     ifelse(boosters35_time_since_jab  > 60 & boosters35_time_since_jab <= 90, "61-90",
                                     ifelse(boosters35_time_since_jab  > 90 & boosters35_time_since_jab <= 120, "91-120",
                                     ifelse(boosters35_time_since_jab  > 120 & boosters35_time_since_jab <= 150, "121-150",
                                     ifelse(boosters35_time_since_jab  > 150 & boosters35_time_since_jab <= 180, "151-180",
                                     "180+"))))))) %>%
  mutate(boosters35_time_since_jab = factor(boosters35_time_since_jab,
         levels=c("0-30", "31-60", "61-90", "91-120", "121-150", "151-180", "180+"),ordered=FALSE))

# limit analysis to adults
df_round <- df_round %>%
  filter(age > 17)

# Extracting covariate names
tmp <- read.xlsx(annot_file)
covs_names <- tmp[, 2]
names(covs_names) <- tmp[, 1]
covs_names <- covs_names[covs]

# Removing unused variables
df_round <- df_round[, c(res_param, covs)]

# Recoding variables
covs_to_recode <- getSheetNames(recoding_file)
covs_to_recode <- intersect(names(covs_names), covs_to_recode)
for (i in 1:length(covs_to_recode)) {
  recoding <- read.xlsx(recoding_file, sheet = covs_to_recode[i])
  recoding[which(is.na(recoding[, 1])), 1] <- "NA"
  renaming <- recoding[, 2]
  names(renaming) <- recoding[, 1]
  x <- as.character(df_round[, covs_to_recode[i]])
  print(table(x))
  x[is.na(x)] <- "NA"

  # Removing categories with less than 10 observations
  if (any(table(x) < 10)) {
    toremove <- names(table(x))[which(table(x) < 10)]
    print(paste0("Excluding category ", renaming[toremove]))
    x[which(x == toremove)] <- NA
    renaming <- renaming[!names(renaming) %in% toremove]
  }

  x <- factor(x, levels = names(renaming), labels = renaming)
  print(table(x, useNA = "always"))

  # Defining the reference level
  x <- relevel(x, ref = recoding[which(recoding[, 3] == "*"), 2])

  df_round[, covs_to_recode[i]] <- x
}

# Recoding continuous to categorical
covs_to_recode <- getSheetNames(recoding_from_cont_file)
covs_to_recode <- intersect(names(covs_names), covs_to_recode)
if (length(covs_to_recode) > 0) {
  for (i in 1:length(covs_to_recode)) {
    recoding <- read.xlsx(recoding_from_cont_file, sheet = covs_to_recode[i])
    x <- as.numeric(df_round[,covs_to_recode[i]])
    x <- cut(x, breaks = c(min(x, na.rm=TRUE) - 10, recoding[, 1]),
             labels = recoding[, 2])
    if (sum(is.na(x))>0){
      x=as.character(x)
      x[is.na(x)]="Unknown"
      x=as.factor(x)
    }

    # Defining the reference level
    x <- relevel(x, ref = recoding[which(recoding[, 3] == "*"), 2])
    print(table(x))

    df_round[, covs_to_recode[i]] <- x
  }
}

# df_round <- na.exclude(df_round)

# Specific recoding for r14/r15
if ("7+" %in% df_round$hh_size_cat) {
  df_round$hh_size_cat <- factor(df_round$hh_size_cat,
    levels = c(1:6, "7+"),
    labels = c(
      "1-2", "1-2",
      "3-5", "3-5", "3-5",
      "6+", "6+"
    )
  )
} else {
  df_round$hh_size_cat <- factor(df_round$hh_size_cat,
    levels = c(1:5, "6+"),
    labels = c(
      "1-2", "1-2",
      "3-5", "3-5", "3-5",
      "6+"
    )
  )
}

# Logistic regressions (base model)
i <- 1
mytable <- NULL
for (i in 1:length(covs)) {
  covariate <- covs[i]
  print(covariate)


  # Base model (adjusted on age and gender)
  myformula <- paste0(
    res_param, "~",
    paste(unique(c(covariate, confounders)), collapse = "+")
  )
  mymodel <- glm(as.formula(myformula), data = df_round, family = binomial(link = "logit"))
  coefs <- exp(mymodel$coefficients)
  coefs_ci <- exp(confint.default(mymodel))
  tmp <- cbind(
    coefs[grep(paste0("^", covariate), names(coefs))],
    coefs_ci[grep(paste0("^", covariate), rownames(coefs_ci)), , drop = FALSE]
  )
  base_model <- FormatCI(FormatOR(tmp))
  tmpnames <- names(coefs)[grep(paste0("^", covariate), names(coefs))]
  tmpnames <- gsub(covariate, "", tmpnames)

  myref <- levels(df_round[, covariate])[!levels(df_round[, covariate]) %in% tmpnames]
  base_model <- rbind("Ref", base_model)
  tmpnames <- c(myref[1], tmpnames)

  # Special case for unknown (e.g. for gender)
  if (length(myref)>1){
    base_model=c(base_model, rep(NA, length(myref)-1))
    tmpnames=c(tmpnames, myref[-1])
  }

  tmp_output <- cbind(
    c(covariate, rep("", length(tmpnames) - 1)),
    tmpnames, base_model
  )
  mytable <- rbind(mytable, tmp_output)
}
mytable_base <- mytable

# Logistic regression (full model)
myformula <- paste0(
  res_param, "~",
  paste(unique(covs)[1:9], collapse = "+")
)
print(myformula)
mymodel <- glm(as.formula(myformula), data = df_round, family = binomial(link = "logit"))

mytable <- NULL
for (i in 1:length(covs)) {
  covariate <- covs[i]
  print(covariate)
  coefs <- exp(mymodel$coefficients)
  coefs_ci <- exp(confint.default(mymodel))
  tmp <- cbind(
    coefs[grep(paste0("^", covariate), names(coefs))],
    coefs_ci[grep(paste0("^", covariate), rownames(coefs_ci)), , drop = FALSE]
  )
  base_model <- FormatCI(FormatOR(tmp))
  tmpnames <- names(coefs)[grep(paste0("^", covariate), names(coefs))]
  tmpnames <- gsub(covariate, "", tmpnames)

  myref <- levels(df_round[, covariate])[!levels(df_round[, covariate]) %in% tmpnames]
  base_model <- rbind("Ref", base_model)
  tmpnames <- c(myref[1], tmpnames)

  # Special case for unknown (e.g. for gender)
  if (length(myref)>1){
    base_model=c(base_model, rep(NA, length(myref)-1))
    tmpnames=c(tmpnames, myref[-1])
  }

  tmp_output <- cbind(
    c(covariate, rep("", length(tmpnames) - 1)),
    tmpnames, base_model
  )
  mytable <- rbind(mytable, tmp_output)
}
mytable_full <- mytable

mytable <- cbind(mytable_base, mytable_full[, 3])
mytable[, 1] <- covs_names[mytable[, 1]]
rownames(mytable) <- paste0(ExtendList(mytable[, 1]), "_", mytable[, 2])
mytable[which(is.na(mytable[, 1])), 1] <- ""

mysheet <- readWorkbook(template_file, sheet = 1)
mysheet <- mysheet[-nrow(mysheet), ]
rownames(mysheet) <- paste0(ExtendList(mysheet[, 1]), "_", mysheet[, 2])
rownames(mysheet) %in% rownames(mytable)
filled_sheet <- mytable[rownames(mysheet), ]
wb <- loadWorkbook(template_file)
writeData(wb,
  sheet = 1, x = filled_sheet,
  startCol = 1, startRow = 2, colNames = FALSE
)
saveWorkbook(wb,
  file = paste0(output_file, "_", round_id, "_", output_tag, ".xlsx"),
  overwrite = TRUE
)

if (direct_export) {
  file.copy(
    from = paste0(output_file, "_", round_id, "_", output_tag, ".xlsx"),
    to = "T:/", overwrite = TRUE
  )
}

