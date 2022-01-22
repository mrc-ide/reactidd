# Script to produce the table of odds ratios
# Two models: adjusted on age and gender and mutually-adjusted
# Output: 1 xlsx file "Logistic_models" with round ID and date
# stored in "E:/Group/report/round16/Tables/"
# and automatically copied into the transfer folder

# To change age categories:
# Step 1: change the thresholds in recoding_from_cont_file (see path below)
# Step 2: manually update the template in template_file (see path below)


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

# Loading functions
source("functions/add_conf_ints.R")
source("functions/make_tables.R")
source("functions/overall_prev.R")
source("functions/formatting_functions.R")


## Parametrisation

# Paths to files
round_id <- 16

data_file <- paste0("E:/Group/saved_objects/rep", round_id, ".rds")
output_file <- "E:/Group/report/round16/Tables/Logistic_models"

output_tag <- Sys.Date()
annot_file <- "E:/Group/report/round16/Parameters/Variable_names.xlsx"
template_file <- "E:/Group/report/round16/Parameters/OR_table_urban.xlsx"

recoding_file <- "E:/Group/report/round16/Parameters/Recoding.xlsx"
recoding_from_cont_file <- "E:/Group/report/round16/Parameters/Recoding_from_continuous.xlsx"

# Copying output files directly to transfer folder
direct_export <- TRUE

# Variable for test results
res_param <- "estbinres"

# Variables for stratification
covs <- c(
  "gender_char", "age", "region",
  "work_new_alt", "ethnic_new_char",
  "urban",
  "hh_size_cat", "nchild2", "imd_quintile"
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

# # Recoding continuous to categorical
covs_to_recode <- getSheetNames(recoding_from_cont_file)
covs_to_recode <- intersect(names(covs_names), covs_to_recode)
if (length(covs_to_recode) > 0) {
  for (i in 1:length(covs_to_recode)) {
    recoding <- read.xlsx(recoding_from_cont_file, sheet = covs_to_recode[i])
    x <- as.numeric(df_round$age)
    x <- cut(x, breaks = c(min(x) - 10, recoding[, 1]), labels = recoding[, 2])
    
    # Defining the reference level
    x <- relevel(x, ref = recoding[which(recoding[, 3] == "*"), 2])
    print(table(x))
    
    df_round[, covs_to_recode[i]] <- x
  }
}

df_round <- na.exclude(df_round)

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
  tmpnames <- c(myref, tmpnames)

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
  paste(unique(covs), collapse = "+")
)
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
  tmpnames <- c(myref, tmpnames)

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

