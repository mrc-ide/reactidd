rm(list = ls())
setwd("E:/Group/report/round16/Scripts/")

library(magrittr)
library(dplyr)
library(data.table)
library(openxlsx)

## Parameters
lag=14
# Path to output file
output_file <- paste0("E:/Group/report/round16/Tables/variant_3v2Dose_14_",lag)

# Choice of data
round_ids <- c(16)
linkage <- TRUE
direct_export=TRUE

# Outcome and predictor variables
outcome <- "variant_type"
predictor <- "vax_status"
effect_of_interest <- "vax_statusVaccinated - 3 doses"

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




# Age bounds
age_lower_bound <- 17
age_upper_bound <- 180
output_file <- paste0(output_file, "_", age_lower_bound, "_", age_upper_bound)

# Analyses restricted to symptomatics
symptomatics <- FALSE
if (symptomatics) {
  output_file <- paste0(output_file, "_sympt")
}

# Read RDS for omicron cases
df_lineage=readRDS(paste0("E:/Group/saved_objects/rep",round_ids,"_lineage.rds"))
ids_omicron=df_lineage$u_passcode[which(df_lineage$react_lineage=="BA.1")]
ids_delta=df_lineage$u_passcode[which(df_lineage$react_lineage!="BA.1")]


# Formatting
CI <- c(" (", ", ", ")")

# Path to file for recoding of the categorical variables
recoding_file <- "E:/Group/report/round15/TmpBB/Recoding.xlsx"




for (id in 1:length(round_ids)) {
  round_id <- round_ids[id]
  
  df_round <- readRDS(paste0("E:/dt20/linkedR", round_id, "datANG.rds"))
  df_round$lineage=rep(NA,nrow(df_round))
  
  df_round <- df_round %>% mutate(lineage = ifelse(
    u_passcode %in% ids_omicron, "omicron", ifelse(
    u_passcode %in% ids_delta, "delta", NA  
    )))
  
  # df_round[ids_omicron, "lineage"]="omicron"
  # df_round[ids_delta, "lineage"]="delta"

  vaxFrom1 <- 14 # Consider vaccinated from this many days after first dose
  vaxFrom2 <- 14
  if(round_id==16){
    df_round$link_vaxtype=df_round$link_vax1type
    if(lag==21){
      df_round <- df_round %>%
        mutate(
          vax_status = link_vax_status_booster21days,
          vax_type = ifelse(link_vaxtype == "Oxford", "AZ", link_vaxtype),
          vax_type = ifelse(!(vax_type %in% c("AZ", "Pfizer", "Moderna", "unvaccinated")) & !vax_status == "unvaccinated", "unknown", vax_type),
          vax_type = ifelse(!(vax_type %in% c("AZ", "Pfizer", "Moderna", "unvaccinated", "unknown")), "unvaccinated", vax_type),
          # Overwrite if not vaccinated by swab:
          vax_type = ifelse(vax_status == "unvaccinated", "unvaccinated", vax_type)
        ) #%>%
      #filter(vax_status %in% c("2 doses", "3 doses"))
    }
    if(lag==14){
      df_round <- df_round %>%
        mutate(
          vax_status = link_vax_status_booster14days,
          vax_type = ifelse(link_vaxtype == "Oxford", "AZ", link_vaxtype),
          vax_type = ifelse(!(vax_type %in% c("AZ", "Pfizer", "Moderna", "unvaccinated")) & !vax_status == "unvaccinated", "unknown", vax_type),
          vax_type = ifelse(!(vax_type %in% c("AZ", "Pfizer", "Moderna", "unvaccinated", "unknown")), "unvaccinated", vax_type),
          # Overwrite if not vaccinated by swab:
          vax_type = ifelse(vax_status == "unvaccinated", "unvaccinated", vax_type)
        ) %>%
      mutate(variant_type = ifelse(lineage == "delta", 0, ifelse(lineage == "omicron", 1, NA))) %>%
      filter(vax_status %in% c("2 doses", "3 doses"))
    }
    if(lag==7){
      df_round <- df_round %>%
        mutate(
          vax_status = link_vax_status_7days,
          vax_type = ifelse(link_vaxtype == "Oxford", "AZ", link_vaxtype),
          vax_type = ifelse(!(vax_type %in% c("AZ", "Pfizer", "Moderna", "unvaccinated")) & !vax_status == "unvaccinated", "unknown", vax_type),
          vax_type = ifelse(!(vax_type %in% c("AZ", "Pfizer", "Moderna", "unvaccinated", "unknown")), "unvaccinated", vax_type),
          # Overwrite if not vaccinated by swab:
          vax_type = ifelse(vax_status == "unvaccinated", "unvaccinated", vax_type)
        ) %>%
        filter(vax_status %in% c("2 doses", "3 doses"))
      
    }
    if(lag==0){
      df_round <- df_round %>%
        mutate(
          vax_status = link_vax_status_booster0days,
          vax_type = ifelse(link_vaxtype == "Oxford", "AZ", link_vaxtype),
          vax_type = ifelse(!(vax_type %in% c("AZ", "Pfizer", "Moderna", "unvaccinated")) & !vax_status == "unvaccinated", "unknown", vax_type),
          vax_type = ifelse(!(vax_type %in% c("AZ", "Pfizer", "Moderna", "unvaccinated", "unknown")), "unvaccinated", vax_type),
          # Overwrite if not vaccinated by swab:
          vax_type = ifelse(vax_status == "unvaccinated", "unvaccinated", vax_type)
        ) %>%
        filter(vax_status %in% c("2 doses", "3 doses"))
      
    }
  }else{
    df_round <- df_round %>%
      mutate(
        vax_status = link_vax_status_2dose14days,
        vax_type = ifelse(link_vaxtype == "Oxford", "AZ", link_vaxtype),
        vax_type = ifelse(!(vax_type %in% c("AZ", "Pfizer", "Moderna", "unvaccinated")) & !vax_status == "unvaccinated", "unknown", vax_type),
        vax_type = ifelse(!(vax_type %in% c("AZ", "Pfizer", "Moderna", "unvaccinated", "unknown")), "unvaccinated", vax_type),
        # Overwrite if not vaccinated by swab:
        vax_type = ifelse(vax_status == "unvaccinated", "unvaccinated", vax_type)
      )# %>%
    #filter(vax_status %in% c("2 doses", "3 doses"))
    
    
  }
  df_round$vax_status <- factor(df_round$vax_status,
                                levels = c("2 doses", "3 doses"),
                                labels = c("Vaccinated - 2 doses", "Vaccinated - 3 doses")
                                
                                #df_round$vax_status <- factor(df_round$vax_status,
                                #                             levels = c("unvaccinated", "1 dose", "2 doses", "3 doses"),
                                #                            labels = c("Unvaccinated", "Vaccinated - 1 dose", "Vaccinated - 2 doses", "Vaccinated - 3 doses")
                                
                                
  )
  df_round$vax_type <- relevel(as.factor(df_round$vax_type), ref = "AZ")
  if(round_id %in% c(13,14,15)){
    df_round$gender_char= factor(df_round$gender_char,levels=c('1','2'),labels=c('Male','Female'))
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

# Filtering: looking only at 35+
full_dataset <- filter(full_dataset, age > as.numeric(age_lower_bound) & age < as.numeric(age_upper_bound))
print("Age range:")
print(range(full_dataset$age))

if (symptomatics) {
  full_dataset <- full_dataset[which(full_dataset$sympt_cat %in% c("Classic COVID symptoms", "Other symptoms")), ]
}


annot_file <- "E:/Group/report/round15/TmpBB/Variable_names.xlsx"

#Only Consider double doses >6 months
# Id1=which(full_dataset$link_vaxrecenttoswab>180)
# Id1=which(full_dataset$link_vaxrecenttoswab>181 & full_dataset$vax_status=="Vaccinated - 2 doses")
# Id2=which(full_dataset$vax_status=="Vaccinated - 3 doses")
# Id=sort(c(Id1,Id2))
# 
# table(full_dataset[Id,]$vax_status)

# Extracting covariate names
covs <- unique(c(
  outcome, predictor,
  unique(unlist(models)),
  "round")
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
# full_dataset <- full_dataset[Id, covs]
full_dataset <- full_dataset[, covs] # for those with 2 dose only definition


mytable <- NULL
for (model_id in paste0("m", 0:2)) {
  # Applying stratification
  mydata <- full_dataset
  print(nrow(mydata))
  if(model_id=="m0"){
    mytable <- c(table(full_dataset$vax_status,full_dataset$variant_type)[1,],"-")
  }
  # Computing the counts
  counts <- table(mydata[which(mydata[,predictor]=="Vaccinated - 3 doses"), outcome])
  # Running the model (no interaction)
  f <- paste0(outcome, " ~ ", paste0(c(predictor, models[[model_id]]), collapse = " + "))
  print(f)
  
  mymodel <- glm(as.formula(f), data = mydata, family = binomial(link = "logit"))
  
  # Extracting relevant coefficients
  VE <- cbind(
    exp(mymodel$coefficients)[effect_of_interest],
    exp(confint.default(mymodel))[effect_of_interest, , drop = FALSE]
  )
  VE_output <- paste0(
    paste0(formatC(VE[1, 1] , format = "f", digits = 4)),
    CI[1],
    paste0(formatC(VE[1, 2] , format = "f", digits = 4)),
    CI[2],
    paste0(formatC(VE[1, 3] , format = "f", digits = 4)),
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
  colnames(mytable) <- c("Test negatives", "Test positives", "OR (95% CI)")
  rownames(mytable) <- c("Vaccinated- 2Doses",paste0("model",c(0:2)))
}else{
  colnames(mytable) <- c("Category", "Adjustment", "Test negatives", "Test positives", "OR (95% CI)", "p-Interaction")
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

