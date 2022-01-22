rm(list = ls())
setwd("E:/Group/report/round16/Scripts/")


library(magrittr)
library(dplyr)
library(data.table)
library(openxlsx)

## Parameters

# Path to output file
output_file <- "E:/Group/report/round16/Tables/Vaccine_efficiency"

# Choice of data
round_ids <- c(14,15,16)
linkage <- TRUE
direct_export=TRUE

# Outcome and predictor variables
outcome <- "estbinres"
predictor <- "vax_status"
effect_of_interest <- "vax_statusVaccinated - 2 doses"

# List of models: variables to adjust on
models <- list(
  m0 = NULL, # unadjusted
  m1 = c("gender_char", "age"), # adjusted on age and gender
  m2 = c("gender_char", "age", "imd_quintile", "region_id", "ethnic_new")
) # further adjusted on deprivation, region and ethnicity

# Models are further adjusted on round if looking at > 1 round
if (length(round_ids) > 1) {
  models <- lapply(models, FUN = function(x) {
    c("round", x)
  })
}

# Adding an interaction term with round if looking at > 1 round
if (length(round_ids) > 1) {
  interaction_term <- "vax_status*round"
} else {
  interaction_term <- NULL
}

# Stratification in the table
# strata_var_list: vector of variables to use for classification
# strata_list: list with categories to use for each stratum
strata_var_list <- rep("vax_type", 5)
strata_list <- list(
  "all",
  c("unvaccinated", "AZ"),
  c("unvaccinated", "Pfizer"),
  c("unvaccinated", "Moderna"),
  c("unvaccinated", "unknown")
) # for linked
# strata_list=list("all",
#                  "Yes, confirmed by a positive test",
#                  "No")

# Age bounds
age_lower_bound <- 17
age_upper_bound <- 65
output_file <- paste0(output_file, "_", age_lower_bound, "_", age_upper_bound)

# Analyses restricted to symptomatics
symptomatics <- FALSE
if (symptomatics) {
  output_file <- paste0(output_file, "_sympt")
}

# Formatting
CI <- c(" (", ", ", ")")

# Path to file for recoding of the categorical variables
recoding_file <- "E:/Group/report/round15/TmpBB/Recoding.xlsx"

for (id in 1:length(round_ids)) {
  round_id <- round_ids[id]

  if (linkage) {
    df_round <- readRDS(paste0("E:/dt20/linkedR", round_id, "datANG.rds"))

    vaxFrom1 <- 14 # Consider vaccinated from this many days after first dose
    vaxFrom2 <- 14
    if(round_id==16){
      df_round <- df_round %>%
        mutate(
          vax_status = link_vax_status_booster14days,
          vax_type = ifelse(link_vax1type == "Oxford", "AZ", link_vax1type),
          vax_type = ifelse(!(vax_type %in% c("AZ", "Pfizer", "Moderna", "unvaccinated")) & !vax_status == "unvaccinated", "unknown", vax_type),
          vax_type = ifelse(!(vax_type %in% c("AZ", "Pfizer", "Moderna", "unvaccinated", "unknown")), "unvaccinated", vax_type),
          # Overwrite if not vaccinated by swab:
          vax_type = ifelse(vax_status == "unvaccinated", "unvaccinated", vax_type)
        ) %>%
        filter(vax_status %in% c("unvaccinated", "2 doses"))
    }
    else if(round_id==15){
      df_round <- df_round %>%
          mutate(
          vax_status = link_vax_status_14days,
          vax_type = ifelse(link_vaxtype == "Oxford", "AZ", link_vaxtype),
          vax_type = ifelse(!(vax_type %in% c("AZ", "Pfizer", "Moderna", "unvaccinated")) & !vax_status == "unvaccinated", "unknown", vax_type),
          vax_type = ifelse(!(vax_type %in% c("AZ", "Pfizer", "Moderna", "unvaccinated", "unknown")), "unvaccinated", vax_type),
          # Overwrite if not vaccinated by swab:
          vax_type = ifelse(vax_status == "unvaccinated", "unvaccinated", vax_type)
          ) %>%
          filter(vax_status %in% c("unvaccinated", "2 doses"))
    }else{
      df_round <- df_round %>%
        mutate(
          vax_status = link_vax_status_2dose14days,
          vax_type = ifelse(link_vaxtype == "Oxford", "AZ", link_vaxtype),
          vax_type = ifelse(!(vax_type %in% c("AZ", "Pfizer", "Moderna", "unvaccinated")) & !vax_status == "unvaccinated", "unknown", vax_type),
          vax_type = ifelse(!(vax_type %in% c("AZ", "Pfizer", "Moderna", "unvaccinated", "unknown")), "unvaccinated", vax_type),
          # Overwrite if not vaccinated by swab:
          vax_type = ifelse(vax_status == "unvaccinated", "unvaccinated", vax_type)
        ) %>%
      filter(vax_status %in% c("unvaccinated", "2 doses"))
      
      
    }
    df_round$vax_status <- factor(df_round$vax_status,
      levels = c("unvaccinated", "2 doses"),
      labels = c("Unvaccinated", "Vaccinated - 2 doses")
      
    # df_round$vax_status <- factor(df_round$vax_status,
    #                               levels = c("unvaccinated", "1 dose", "2 doses", "3 doses"),
    #                               labels = c("Unvaccinated", "Vaccinated - 1 dose", "Vaccinated - 2 doses", "Vaccinated - 3 doses")
    #                                 
      
    )
    df_round$vax_type <- relevel(as.factor(df_round$vax_type), ref = "unvaccinated")
    if(round_id %in% c(13,14)){
      df_round$gender_char= factor(df_round$gender_char,levels=c('1','2'),labels=c('Male','Female'))
    }
    if(round_id==15){
      df_round$gender_char= factor(df_round$gender_char,levels=c('Male','Female'),labels=c('Male','Female')) 
    }
  } else {
    df_round <- data.frame(readRDS(paste0("E:/Group/saved_objects/rep", round_id, ".rds")))

    # Adding variable introduced in round 15
    if (!"vax_status_noDate_v2" %in% colnames(df_round)) {
      df_round$vax_status_noDate_v2 <- df_round$vax_status_noDate
    }

    # ids=intersect(colnames(dfl1), colnames(dfl2))
    # dfl=rbind(dfl1[,ids], dfl2[,ids])
    # dfl$round=as.factor(c(rep(paste0('Round', round_ids[1]), nrow(dfl1)),
    # rep(paste0('Round', round_ids[2]), nrow(dfl2))))

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

    # Excluding unknowns and 1 doses
    df_round$vax_status <- df_round$vax_status_noDate_v2
    df_round$vax_status <- factor(df_round$vax_status,
      levels = c("Not vaccinated", "Two does", "Three does"),
      labels = c("Unvaccinated", "Vaccinated - 2 doses", "Vaccinated - 3 doses")
    )
  }

  assign(paste0("dfl", id), df_round)
}

if (length(round_ids) ==1) {
  full_dataset <- dfl1
  full_dataset$round <- as.factor(rep(paste0("Round", round_ids[1]), nrow(dfl1)))
}

if (length(round_ids) ==2) {
  ids <- intersect(colnames(dfl1), colnames(dfl2))
  full_dataset <- rbind(dfl1[, ids], dfl2[, ids])
  full_dataset$round <- as.factor(c(
    rep(paste0("Round", round_ids[1]), nrow(dfl1)),
    rep(paste0("Round", round_ids[2]), nrow(dfl2))
  ))
}
if (length(round_ids)==3) {
  ids <- intersect(intersect(colnames(dfl1), colnames(dfl2)), colnames(dfl3))
  full_dataset <- rbind(dfl1[, ids], dfl2[, ids], dfl3[, ids])
  full_dataset$round <- as.factor(c(
    rep(paste0("Round", round_ids[1]), nrow(dfl1)),
    rep(paste0("Round", round_ids[2]), nrow(dfl2)),
    rep(paste0("Round", round_ids[3]), nrow(dfl3))
  ))
} 

full_dataset$sympt_cat=factor(full_dataset[,"sympt_cat"], 
                     levels=c(1,2,3,
                              "No symptoms",
                              "Classic COVID symptoms",
                              "Other symptoms",
                              "NA"),
                     labels=c("No symptoms",
                              "Classic COVID symptoms",
                              "Other symptoms",
                              "No symptoms",
                              "Classic COVID symptoms",
                              "Other symptoms",
                              "Unknown"))





# Filtering: looking only at 18-64
full_dataset <- filter(full_dataset, age > as.numeric(age_lower_bound) & age < as.numeric(age_upper_bound))
print("Age range:")
print(range(full_dataset$age))

if (symptomatics) {
  full_dataset <- full_dataset[which(full_dataset$sympt_cat %in% c("Classic COVID symptoms", "Other symptoms")), ]
}

annot_file <- "E:/Group/report/round15/TmpBB/Variable_names.xlsx"

# Extracting covariate names
covs <- unique(c(
  outcome, predictor,
  unique(unlist(models)),
  "round", strata_var_list)
)

tmp <- read.xlsx(annot_file)
covs_names <- tmp[, 2]
names(covs_names) <- tmp[, 1]
covs_names <- covs_names[covs]

# Recoding categorical variables
covs_to_recode <- getSheetNames(recoding_file)
covs_to_recode <- intersect(names(covs_names), covs_to_recode)
if (length(covs_to_recode) > 0) {
  for (i in 1:length(covs_to_recode)) {
    recoding <- read.xlsx(recoding_file, sheet = covs_to_recode[i])
    recoding[which(is.na(recoding[, 1])), 1] <- "NA"
    renaming <- recoding[, 2]
    names(renaming) <- recoding[, 1]
    x <- as.character(full_dataset[, covs_to_recode[i]])
    print(table(x))
    x[is.na(x)] <- "NA"
    x <- factor(x, levels = names(renaming), labels = renaming)
    print(table(x))
    full_dataset[, covs_to_recode[i]] <- x
  }
}

# List of variables to keep
full_dataset <- full_dataset[, covs]

full_dataset$vax_status_Round= (paste0(as.character(full_dataset$vax_status),"_",as.character(full_dataset$round)))


full_dataset$vax_status_Round=factor(full_dataset$vax_status_Round,
                                     levels=c("Unvaccinated_Round14","Vaccinated - 2 doses_Round14",
                                              "Unvaccinated_Round15","Vaccinated - 2 doses_Round15",
                                              "Unvaccinated_Round16","Vaccinated - 2 doses_Round16"),
                                     labels=c("Unvaccinated","Vaccinated - 2 doses_Round14",
                                              "Unvaccinated","Vaccinated - 2 doses_Round15",
                                              "Unvaccinated","Vaccinated - 2 doses_Round16"))



full_dataset$vax_status_Type= (paste0(as.character(full_dataset$vax_status),"_",as.character(full_dataset$vax_type)))


full_dataset$vax_status_Type=factor(full_dataset$vax_status_Round,
                                     levels=c("Unvaccinated_Round14","Vaccinated - 2 doses_Round14",
                                              "Unvaccinated_Round15","Vaccinated - 2 doses_Round15",
                                              "Unvaccinated_Round16","Vaccinated - 2 doses_Round16"),
                                     labels=c("Unvaccinated","Vaccinated - 2 doses_Round14",
                                              "Unvaccinated","Vaccinated - 2 doses_Round15",
                                              "Unvaccinated","Vaccinated - 2 doses_Round16"))



myfulltable <- NULL
for (stratum_id in 1:length(strata_list)) {
  mytable <- NULL
  for (model_id in paste0("m", 0:2)) {
    # Applying stratification
    mydata <- full_dataset
    if (strata_list[stratum_id] != "all") {
      mydata <- mydata[which(mydata[, strata_var_list[stratum_id]] %in% strata_list[[stratum_id]]), ]
      print(nrow(mydata))
    }

    # Computing the counts
    counts <- table(mydata[which(mydata[,predictor]=="Vaccinated - 2 doses"), outcome])

    # Running the model (no interaction)
    f <- paste0(outcome, " ~ ", paste0(c(predictor, models[[model_id]]), collapse = " + "))
    print(f)
    mymodel <- glm(as.formula(f), data = mydata, family = binomial(link = "logit"))

    # Extracting relevant coefficients
    VE <- cbind(
      1 - exp(mymodel$coefficients)[effect_of_interest],
      1 - exp(confint.default(mymodel))[effect_of_interest, , drop = FALSE]
    )
    VE_output <- paste0(
      paste0(formatC(VE[1, 1] * 100, format = "f", digits = 2), "%"),
      CI[1],
      paste0(formatC(VE[1, 3] * 100, format = "f", digits = 2), "%"),
      CI[2],
      paste0(formatC(VE[1, 2] * 100, format = "f", digits = 2), "%"),
      CI[3]
    )

    if (!is.null(interaction_term)) {
      # Running the model (with interaction)
      f <- paste0(f, " + ", interaction_term)
      print(f)
      mymodel_interaction <- glm(as.formula(f), data = mydata, family = binomial(link = "logit"))

      # Extracting relevant coefficients
      interaction <- summary(mymodel_interaction)$coefficients

      # Combining outputs from the two models
      output <- c(
        counts,
        VE_output,
        formatC(interaction[nrow(interaction), 4], format = "e", digits = 2)
      )
    } else {
      output <- c(counts, VE_output)
    }
    mytable <- rbind(mytable, output)
  }
  mytable <- cbind(
    c(paste(strata_list[[stratum_id]], collapse = " vs "), rep("", length(models) - 1)),
    paste0("Model ", seq(0, length(models) - 1)),
    mytable
  )
  myfulltable <- rbind(myfulltable, mytable)
}
colnames(myfulltable) <- c("Category", "Adjustment", "Test negatives", "Test positives", "Vaccine Effectiveness", "p-Interaction")

myfulltable

write.xlsx(myfulltable,
           paste0(output_file, "_", paste(round_ids, collapse=""), "_", Sys.Date(), ".xlsx"),
  rowNames = FALSE
)

if (direct_export) {
  file.copy(
    from = paste0(output_file, "_", paste(round_ids, collapse=""), "_", Sys.Date(), ".xlsx"),
    to = "T:/", overwrite = TRUE
  )
}

