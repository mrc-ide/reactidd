rm(list = ls())
setwd("E:/Group/report/round16/Scripts/")

library(magrittr)
library(dplyr)
library(data.table)
library(openxlsx)

## Parameters
lag=14
strict=TRUE
symptomatics <- "ANY" #FALSE, COVID, or ANY

# Path to output file
if(strict==TRUE){
  output_file <- paste0("E:/Group/report/round16/Tables/Strict_Children_efficiency")
}else{
  output_file <- paste0("E:/Group/report/round16/Tables/NonStrict_Children_efficiency")
  
}
# Choice of data
round_ids <- c(14,15,16)
linkage <- TRUE
direct_export=TRUE

# Outcome and predictor variables
outcome <- "estbinres2"
predictor <- "vax_status2"
if(strict==FALSE){
  effect_of_interest <- "vax_status2Vaccinated - 1/2 doses"
}else{
  effect_of_interest <- "vax_status2Vaccinated - 1 dose"
}
# List of models: variables to adjust on
models <- list(
  m0 = "age", # unadjusted
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


# Age bounds
age_lower_bound <- 11
age_upper_bound <- 18
output_file <- paste0(output_file, "_", age_lower_bound, "_", age_upper_bound)

# Analyses restricted to symptomatics
if (symptomatics!=FALSE) {
  output_file <- paste0(output_file, "_",symptomatics)
}

# Formatting
CI <- c(" (", ", ", ")")

# Path to file for recoding of the categorical variables
recoding_file <- "E:/Group/report/round15/TmpBB/Recoding.xlsx"




for (id in 1:length(round_ids)) {
  round_id <- round_ids[id]
  print(round_id)
  
  if (linkage) {
    if(round_id==15){
      df_round <- readRDS(paste0("E:/dt20/linkedR", round_id, "datANG.rds"))
    }
    if(round_id==16){
      df_round <- readRDS(paste0("E:/dt20/linkedR", round_id, "datANG.rds"))
    }else{
      df_round <- readRDS(paste0("E:/dt20/linkedR", round_id, "datANG.rds"))
    }
    vaxFrom1 <- 14 # Consider vaccinated from this many days after first dose
    vaxFrom2 <- 14
    if(round_id==16){
      df_round$link_vaxtype=df_round$link_vax1type
      df_round <- df_round %>%
        mutate(
          vax_status = link_vax_status_booster14days,
          vax_type = ifelse(link_vaxtype == "Oxford", "AZ", link_vaxtype),
          vax_type = ifelse(!(vax_type %in% c("AZ", "Pfizer", "Moderna", "unvaccinated")) & !vax_status == "unvaccinated", "unknown", vax_type),
          vax_type = ifelse(!(vax_type %in% c("AZ", "Pfizer", "Moderna", "unvaccinated", "unknown")), "unvaccinated", vax_type),
          # Overwrite if not vaccinated by swab:
          vax_type = ifelse(vax_status == "unvaccinated", "unvaccinated", vax_type)
        ) %>%
        # filter(vax_status %in% c("2 doses","3 doses"))
        filter(vax_status %in% c("unvaccinated","1 dose","2 doses")) # original
    } else if(round_id==15){
      df_round$link_vaxtype=df_round$link_boostertype
      df_round <- df_round %>%
        mutate(
          vax_status = link_vax_status_booster14days,
          vax_type = ifelse(link_vaxtype == "Oxford", "AZ", link_vaxtype),
          vax_type = ifelse(!(vax_type %in% c("AZ", "Pfizer", "Moderna", "unvaccinated")) & !vax_status == "unvaccinated", "unknown", vax_type),
          vax_type = ifelse(!(vax_type %in% c("AZ", "Pfizer", "Moderna", "unvaccinated", "unknown")), "unvaccinated", vax_type),
          # Overwrite if not vaccinated by swab:
          vax_type = ifelse(vax_status == "unvaccinated", "unvaccinated", vax_type)
        ) %>%
        # filter(vax_status %in% c("2 doses","3 doses"))
        filter(vax_status %in% c("unvaccinated","1 dose","2 doses")) # original
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
        filter(vax_status %in% c("unvaccinated","1 dose","2 doses"))
      
      
    }
    if(strict==FALSE){
      df_round$vax_status2 <- factor(df_round$vax_status,
                                     levels = c("unvaccinated","1 dose", "2 doses"),
                                     labels = c("unvaccinated","Vaccinated - 1/2 doses","Vaccinated - 1/2 doses"))
    }else{
    df_round$vax_status2 <- factor(df_round$vax_status,
                                   levels = c("unvaccinated","1 dose"),
                                   labels = c("unvaccinated","Vaccinated - 1 dose"))
      # df_round$vax_status2 <- factor(df_round$vax_status,
      #                                levels = c("unvaccinated","1 dose", "2 doses"),
      #                                labels = c("unvaccinated","Vaccinated - 1 dose", "Vaccinated - 2 doses"))
      
      
    }  
    
    # change definition of vaccination status 
    # below original definition
    #
    # df_round$vax_status <- factor(df_round$vax_status,
    #                               levels = c("unvaccinated", "1 dose", "2 doses", "3 doses"),
    #                               labels = c("Unvaccinated", "Vaccinated - 1 dose", "Vaccinated - 2 doses", "Vaccinated - 3 doses")
    #                                 
    
    
    df_round$vax_type <- relevel(as.factor(df_round$vax_type), ref = "AZ")
    if(round_id %in% c(13,14,15)){
      df_round$gender_char= factor(df_round$gender_char,levels=c('1','2'),labels=c('Male','Female'))
    }
#    if(round_id==15){
 #     df_round$gender_char= factor(df_round$gender_char,levels=c('Male','Female'),labels=c('Male','Female')) 
  #  }
    
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

full_dataset$MyCovida=factor(full_dataset$covida, levels=c(1,2,3,4),
                             labels=c("COVID19 case Yes, confirmed",
                                      "COVID19 case Yes, suspected",
                                      "COVID19 case Yes, suspected",
                                      "No"))


# Filtering: looking only at 35+
full_dataset <- filter(full_dataset, age > as.numeric(age_lower_bound) & age < as.numeric(age_upper_bound))
print("Age range:")
print(range(full_dataset$age))


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


full_dataset$estbinres2=full_dataset$estbinres

if (symptomatics=="COVID") {
  Id1=which(full_dataset$sympt_cat %in% c("No symptoms","Other symptoms"))
  full_dataset$estbinres2[Id1]=0
}

if (symptomatics=="ANY") {
  Id1=which(full_dataset$sympt_cat %in% c("No symptoms"))
  full_dataset$estbinres2[Id1]=0
}


annot_file <- "E:/Group/report/round15/TmpBB/Variable_names.xlsx"

# Extracting covariate names
covs <- unique(c(
  outcome, predictor,
  unique(unlist(models)),
  "round","vax_status")
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

mytable <- NULL
for (model_id in paste0("m", 0:2)) {
  # Applying stratification
  mydata <- full_dataset
  print(nrow(mydata))
  if(model_id=="m0"){
    if(is.null(interaction_term)){
      mytable <- c(table(full_dataset$vax_status2,full_dataset$estbinres)[1,],"-")
    }else{
      mytable <- c(table(full_dataset$vax_status2,full_dataset$estbinres)[1,],"-","-")
    }
  }
  # Computing the counts
  if(strict==FALSE){
    counts <- table(mydata[which(mydata[,predictor]=="Vaccinated - 1/2 doses"), outcome])
  }else{
    counts <- table(mydata[which(mydata[,predictor]=="Vaccinated - 1 dose"), outcome])
    
  }
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
if(is.null(interaction_term)) {
  colnames(mytable) <- c("Test negatives", "Test positives", "Vaccine Effectiveness")
  rownames(mytable) <- c("Unvaccinated",paste0("model",c(0:2)))
}else{
  colnames(mytable) <- c("Test negatives", "Test positives", "Vaccine Effectiveness","p-interaction")
  rownames(mytable) <- c("Unvaccinated",paste0("model",c(0:2)))
}
write.xlsx(mytable,
           paste0(output_file, "_", paste(round_ids, collapse=""), "_", Sys.Date(), ".xlsx"),
           rowNames = TRUE
)

if (direct_export) {
  file.copy(
    from = paste0(output_file, "_", paste(round_ids, collapse=""), "_", Sys.Date(), ".xlsx"),
    to = "T:/", overwrite = TRUE
  )
}

print(output_file)
