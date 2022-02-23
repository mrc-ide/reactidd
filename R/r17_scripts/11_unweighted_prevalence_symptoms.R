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
setwd(paste0("E:/Group/report/round", 17))

# Loading required packages
source("Scripts/functions/load_packages.R")
pkgs <- c(
  "prevalence", "mgcv", "mgcViz", "MASS", "dplyr",
  "tidyr", "forcats", "ggplot2", "qpcR", "survey", "reshape2",
  "openxlsx", "data.table"
)
load_packages(pkgs)

# Source any functions from the local file
source("Scripts/functions/add_conf_ints.R")
source("Scripts/functions/make_tables.R")
source("Scripts/functions/overall_prev.R")
source("Scripts/functions/formatting_functions.R")


## Parametrisation

# Modelling options
weighted <- FALSE # whether weighted prevalences should be included or not

# Copying output files directly to transfer folder
direct_export <- FALSE

# Paths to files
overall_file <- "Tables/overall_covid_prevalence_symptoms"
positive_file <- "Tables/covid_positive_symptoms"
omicron_file <- "Tables/covid_positive_symptoms_omicron_only"
delta_file <- "Tables/covid_positive_symptoms_delta_only"

output_tag <- Sys.Date()
annot_file <- "Parameters/Variable_names.xlsx"
template_file <- "Parameters/Table2a_urban"
template_sheet1 <- "Table2a"
template_sheet2 <- "Table2b"

recoding_file <- "Parameters/Recoding.xlsx"
recoding_from_cont_file <- "Parameters/Recoding_from_continuous.xlsx"


# Variable for test results
res_param <- "estbinres"

# Variable for weights
weight_params <- NULL


# Variables for stratification
# covs <- c(
#   "gender_char", "age", "region",
#   "work_new_alt", "ethnic_new_char",
#   "hh_size_cat", "covidcon_char", "sympt_cat",
#   "nchild2",
#   "urban",
#   "imd_quintile", "vax_status_noDate_v2"
# )

# Defining the column widths / row heights
# column_widths <- c(5.5, 5.5, 19)
# if (weighted) {
#   column_widths <- c(column_widths, 19)
# }
# row_height <- 15


## Computing the prevalences

data_file <- paste0("E:/Group/saved_objects/rep", round_id, ".rds")
df_round <- setDT(readRDS(data_file))

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

## Task
fake_symptoms_list <- c("covidsym2_1" , "covidsym2_2" , "covidsym2_3" ,
                   "covidsym2_4" , "covidsym2_5" , "covidsym2_6" ,
                   "covidsym2_7" , "covidsym2_8" , "covidsym2_9" ,
                   "covidsym2_10", "covidsym2_11", "covidsym2_12",
                   "covidsym2_13", "covidsym2_14",
                   "covidsym3_1" , "covidsym3_2" , "covidsym3_3" ,
                   "covidsym3_4" , "covidsym3_5" , "covidsym3_6" ,
                   "covidsym3_7" , "covidsym3_8" , "covidsym3_9" ,
                   "covidsym3_9_other", "covidsym3_10", "covidsym3_11" )

symptoms_list <- c("sympt_any1_1"   ,           "sympt_any1_2"    ,
"sympt_any1_3"   ,           "sympt_any1_4"    ,
"sympt_any2_1"   ,           "sympt_any2_2"    ,          "sympt_any2_3"     ,
"sympt_any2_4"   ,           "sympt_any2_5"    ,          "sympt_any2_6"     ,
"sympt_any2_7"   ,           "sympt_any2_8"    ,          "sympt_any2_9"     ,
"sympt_any2_10"  ,           "sympt_any2_11"   ,          "sympt_any2_12"    ,
"sympt_any2_13"  ,
"sympt_any3_1"   ,
"sympt_any3_2"   ,           "sympt_any3_3"    ,          "sympt_any3_4"     ,
"sympt_any3_5"   ,           "sympt_any3_6"    ,          "sympt_any3_7"     ,
"sympt_any3_8"   ,           "sympt_any3_9"    ,          "sympt_any3_10"    ,
"sympt_any3_11")

sympnames <- c("Loss or change of sense of smell", "Loss or change of sense of taste",
               "New persistent cough", "Fever",
               "Runny nose", "Sneezing", "Blocked nose",
               "Sore eyes", "Sore throat", "Hoarse voice",
               "Headache", "Dizziness", "Appetite loss",
               "Nausea/vomiting", "Diarrhoea", "Abdominal pain / belly ache",
               "Shortness of breath",
               "Tight chest",
               "Chest pain", "Chills", "Difficulty sleeping",
               "Tiredness", "Severe fatigue", "Numbness/tingling",
               "Heavy arms/legs", "Muscle ache", "No symptoms",
               "Leg swelling",
               "Don't know")



############################################
#### prevalence for round 17   (ALL)   #####
############################################

output_row <- NULL
temp_df <- NULL
mytable <- NULL
classic_df <- NULL

overall_prev_tab_exact <- NULL
classic_table_exact <- NULL

for (n in 1:length(symptoms_list)) {

  print(paste(n, symptoms_list[n], sympnames[n]))

  temp_df <- df_round

  temp_df <- as.data.frame(temp_df) %>%
    mutate(across(starts_with("sympt_any"),
                  ~recode(., `-92` = 0L,
                             `-91` = 0L,
                             `-77` = 0L))) %>%
    mutate(across(starts_with("sympt_any"),
                  ~replace_na(., 0)))

  temp_df <- as.data.frame(
    temp_df[
      with(temp_df,
           eval(parse(text=symptoms_list[n])) == 1 |
             eval(parse(text=symptoms_list[n])) == 0),])

  print(nrow(temp_df))

  # Unweighted
  overall_prev_tab_exact <- bind_rows(
    overall_prev(temp_df, method = "exact", outcome=symptoms_list[n]),
    .id = "Round") %>%
    mutate(Round = c(as.character(round_id)))

  overall_prev_tab_exact <- as.matrix(overall_prev_tab_exact)


  mytable <- cbind(
    sympnames[n],
    overall_prev_tab_exact[, 1, drop = FALSE],
    FormatCount(overall_prev_tab_exact[, 2:3, drop = FALSE]),
    FormatCI(FormatPrevalence(overall_prev_tab_exact[, 4:6, drop = FALSE]))
  )

  output_row <- rbind(output_row, mytable)

}

classic_df <- as.data.frame(df_round)
classic_df <- classic_df[, c("u_passcode", "sympt_cat")]
classic_df$sympt_cat <- classic_df$sympt_cat %>%
  recode("Classic COVID symptoms" = 1,
         "No symptoms" = 0,
         "Other symptoms" = 0,
         "NA" = 0 )

classic_table_exact <- as.matrix(
  overall_prev(classic_df, method = "exact", outcome="sympt_cat"))

classic_table <- cbind(
  "Classic symptoms",
  as.character(round_id),
  FormatCount(classic_table_exact[, 1:2, drop = FALSE]),
  FormatCI(FormatPrevalence(classic_table_exact[, 3:5, drop = FALSE]))
)

output_row <- rbind(classic_table, output_row)
colnames(output_row) <- c("Symptom", "Round", "Present", "Total COVID+ve", "Unweighted prevalence")


write.xlsx(output_row,
           colNames = TRUE,
           paste0(overall_file, "_r", round_id, "_", output_tag, ".xlsx")
)

if (direct_export) {
  file.copy(
    from = paste0(overall_file, "_r", round_id, "_", output_tag, ".xlsx"),
    to = "T:/", overwrite = TRUE
  )
}


############################################
#### prevalence for round 17 covid +ve #####
############################################

output_row <- NULL
temp_df <- NULL
mytable <- NULL
classic_df <- NULL

overall_prev_tab_exact <- NULL
classic_table_exact <- NULL

for (n in 1:length(symptoms_list)) {

  print(paste(n, symptoms_list[n], sympnames[n]))

  temp_df <- df_round
  temp_df <- temp_df[estbinres == 1, ]

  temp_df <- as.data.frame(temp_df) %>%
    mutate(across(starts_with("sympt_any"),
                  ~recode(., `-92` = 0L,
                             `-91` = 0L,
                             `-77` = 0L))) %>%
    mutate(across(starts_with("sympt_any"),
                  ~replace_na(., 0)))

  temp_df <- as.data.frame(
    temp_df[
      with(temp_df,
           eval(parse(text=symptoms_list[n])) == 1 |
           eval(parse(text=symptoms_list[n])) == 0),])


  print(nrow(temp_df))

  # Unweighted
  overall_prev_tab_exact <- bind_rows(
    overall_prev(temp_df, method = "exact", outcome=symptoms_list[n]),
    .id = "Round") %>%
      mutate(Round = c(as.character(round_id)))

  overall_prev_tab_exact <- as.matrix(overall_prev_tab_exact)


  mytable <- cbind(
    sympnames[n],
    overall_prev_tab_exact[, 1, drop = FALSE],
    FormatCount(overall_prev_tab_exact[, 2:3, drop = FALSE]),
    FormatCI(FormatPrevalence(overall_prev_tab_exact[, 4:6, drop = FALSE]))
  )

  output_row <- rbind(output_row, mytable)

}

classic_df <- df_round
classic_df <- classic_df[estbinres == 1, ]
classic_df <- as.data.frame(classic_df[, c("u_passcode", "sympt_cat")])
classic_df$sympt_cat <- classic_df$sympt_cat %>%
                recode("Classic COVID symptoms" = 1,
                       "No symptoms" = 0,
                       "Other symptoms" = 0,
                       "NA" = 0 )

classic_table_exact <- as.matrix(
  overall_prev(classic_df, method = "exact", outcome="sympt_cat"))

classic_table <- cbind(
  "Classic symptoms",
  as.character(round_id),
  FormatCount(classic_table_exact[, 1:2, drop = FALSE]),
  FormatCI(FormatPrevalence(classic_table_exact[, 3:5, drop = FALSE]))
  )

output_row <- rbind(classic_table, output_row)
colnames(output_row) <- c("Symptom", "Round", "Present", "Total COVID+ve", "Unweighted prevalence")


write.xlsx(output_row,
  colNames = TRUE,
  paste0(positive_file, "_r", round_id, "_", output_tag, ".xlsx")
)

if (direct_export) {
  file.copy(
    from = paste0(positive_file, "_r", round_id, "_", output_tag, ".xlsx"),
    to = "T:/", overwrite = TRUE
  )
}

######################################################
#### prevalence for round 17 covid +ve (OMICRON) #####
######################################################

output_row <- NULL
temp_df <- NULL
mytable <- NULL
classic_df <- NULL

overall_prev_tab_exact <- NULL
classic_table_exact <- NULL

df_lineage=readRDS(paste0("E:/Group/saved_objects/rep",round_id,"_lineage.rds"))

ids_omicron=df_lineage$u_passcode[which(df_lineage$react_lineage=="BA.1")]
ids_delta=df_lineage$u_passcode[which(df_lineage$react_lineage!="BA.1")]

df_round$lineage=rep(as.character(NA),nrow(df_round))
df_round[u_passcode %in% ids_omicron, lineage := "omicron"]
df_round[u_passcode %in% ids_delta  , lineage := "delta"  ]

df_omicron <- df_round[u_passcode %in% ids_omicron, ]
df_delta   <- df_round[u_passcode %in% ids_delta, ]

for (n in 1:length(symptoms_list)) {

  print(paste(n, symptoms_list[n], sympnames[n]))

  temp_df <- df_omicron
  temp_df <- temp_df[estbinres == 1, ]

  temp_df <- as.data.frame(temp_df) %>%
    mutate(across(starts_with("sympt_any"),
                  ~recode(., `-92` = 0L,
                             `-91` = 0L,
                             `-77` = 0L))) %>%
    mutate(across(starts_with("sympt_any"),
                  ~replace_na(., 0)))

  temp_df <- as.data.frame(
    temp_df[
      with(temp_df,
           eval(parse(text=symptoms_list[n])) == 1 |
             eval(parse(text=symptoms_list[n])) == 0),])


  print(nrow(temp_df))

  # Unweighted
  overall_prev_tab_exact <- bind_rows(
    overall_prev(temp_df, method = "exact", outcome=symptoms_list[n]),
    .id = "Round") %>%
    mutate(Round = c(as.character(round_id)))

  overall_prev_tab_exact <- as.matrix(overall_prev_tab_exact)


  mytable <- cbind(
    sympnames[n],
    overall_prev_tab_exact[, 1, drop = FALSE],
    FormatCount(overall_prev_tab_exact[, 2:3, drop = FALSE]),
    FormatCI(FormatPrevalence(overall_prev_tab_exact[, 4:6, drop = FALSE]))
  )

  output_row <- rbind(output_row, mytable)

}

classic_df <- df_omicron
classic_df <- classic_df[estbinres == 1, ]
classic_df <- classic_df[, c("u_passcode", "sympt_cat")]
classic_df$sympt_cat <- classic_df$sympt_cat %>%
  recode("Classic COVID symptoms" = 1,
         "No symptoms" = 0,
         "Other symptoms" = 0,
         "NA" = 0 )

classic_table_exact <- as.matrix(
  overall_prev(as.data.frame(classic_df), method = "exact", outcome="sympt_cat"))

classic_table <- cbind(
  "Classic symptoms",
  as.character(round_id),
  FormatCount(classic_table_exact[, 1:2, drop = FALSE]),
  FormatCI(FormatPrevalence(classic_table_exact[, 3:5, drop = FALSE]))
)

output_row <- rbind(classic_table, output_row)
colnames(output_row) <- c("Symptom", "Round", "Present", "Total COVID+ve", "Unweighted prevalence")


write.xlsx(output_row,
           colNames = TRUE,
           paste0(omicron_file, "_r", round_id, "_", output_tag, ".xlsx")
)

if (direct_export) {
  file.copy(
    from = paste0(omicron_file, "_r", round_id, "_", output_tag, ".xlsx"),
    to = "T:/", overwrite = TRUE
  )
}

######################################################
#### prevalence for round 17 covid +ve (DELTA)   #####
######################################################

output_row <- NULL
temp_df <- NULL
mytable <- NULL
classic_df <- NULL

overall_prev_tab_exact <- NULL
classic_table_exact <- NULL

for (n in 1:length(symptoms_list)) {

  print(paste(n, symptoms_list[n], sympnames[n]))

  temp_df <- df_delta
  temp_df <- temp_df[estbinres == 1, ]

  temp_df <- as.data.frame(temp_df) %>%
    mutate(across(starts_with("sympt_any"),
                  ~recode(., `-92` = 0L,
                          `-91` = 0L,
                          `-77` = 0L))) %>%
    mutate(across(starts_with("sympt_any"),
                  ~replace_na(., 0)))

  temp_df <- as.data.frame(
    temp_df[
      with(temp_df,
           eval(parse(text=symptoms_list[n])) == 1 |
             eval(parse(text=symptoms_list[n])) == 0),])


  print(nrow(temp_df))

  # Unweighted
  overall_prev_tab_exact <- bind_rows(
    overall_prev(temp_df, method = "exact", outcome=symptoms_list[n]),
    .id = "Round") %>%
    mutate(Round = c(as.character(round_id)))

  overall_prev_tab_exact <- as.matrix(overall_prev_tab_exact)


  mytable <- cbind(
    sympnames[n],
    overall_prev_tab_exact[, 1, drop = FALSE],
    FormatCount(overall_prev_tab_exact[, 2:3, drop = FALSE]),
    FormatCI(FormatPrevalence(overall_prev_tab_exact[, 4:6, drop = FALSE]))
  )

  output_row <- rbind(output_row, mytable)

}

classic_df <- df_delta
classic_df <- classic_df[estbinres == 1, ]
classic_df <- classic_df[, c("u_passcode", "sympt_cat")]
classic_df$sympt_cat <- classic_df$sympt_cat %>%
  recode("Classic COVID symptoms" = 1,
         "No symptoms" = 0,
         "Other symptoms" = 0,
         "NA" = 0 )

classic_table_exact <- as.matrix(
  overall_prev(as.data.frame(classic_df), method = "exact", outcome="sympt_cat"))

classic_table <- cbind(
  "Classic symptoms",
  as.character(round_id),
  FormatCount(classic_table_exact[, 1:2, drop = FALSE]),
  FormatCI(FormatPrevalence(classic_table_exact[, 3:5, drop = FALSE]))
)

output_row <- rbind(classic_table, output_row)
colnames(output_row) <- c("Symptom", "Round", "Present", "Total COVID+ve", "Unweighted prevalence")


write.xlsx(output_row,
           colNames = TRUE,
           paste0(delta_file, "_r", round_id, "_", output_tag, ".xlsx")
)

if (direct_export) {
  file.copy(
    from = paste0(delta_file, "_r", round_id, "_", output_tag, ".xlsx"),
    to = "T:/", overwrite = TRUE
  )
}

######################################################
#### prevalence for round 17 by vax status  ?? how to stratify
##### non-vax, vax 1 dose, vax 2 dose, vax 3 doses
######################################################

library(reaction)
library(progressr)

react_cache_environ <- new.env()
assign("current_round_number", 17, envir=react_cache_environ)
define_population(n_prior_rounds_selected = 0,
                  age_lower_bound = 0,
                  age_upper_bound = 150)

df_round <- with_progress(clean_linked_df("E:/dt20/linkedR17datANG.rds", progressr::progressor))

output_row <- NULL
temp_df <- NULL
mytable <- NULL
classic_df <- NULL

overall_prev_tab_exact <- NULL
classic_table_exact <- NULL

list_vax_status <- c("Unvaccinated", "Vaccinated - 1 dose", "Vaccinated - 2 doses", "Vaccinated - 3 doses")

for (vs in list_vax_status) {
  print(vs)

for (n in 1:length(symptoms_list)) {

  print(paste(n, symptoms_list[n], sympnames[n]))

  temp_df <- df_round

  temp_df <- as.data.frame(temp_df) %>%
    mutate(across(starts_with("sympt_any"),
                  ~recode(., `-92` = 0L,
                             `-91` = 0L,
                             `-77` = 0L))) %>%
    mutate(across(starts_with("sympt_any"),
                  ~replace_na(., 0)))

  temp_df <- temp_df[temp_df$vax_status == vs, ]

  temp_df <- as.data.frame(
    temp_df[
      with(temp_df,
           eval(parse(text=symptoms_list[n])) == 1 |
             eval(parse(text=symptoms_list[n])) == 0),])

  print(nrow(temp_df))

  # Unweighted
  overall_prev_tab_exact <- bind_rows(
    overall_prev(temp_df, method = "exact", outcome=symptoms_list[n]),
    .id = "Round") %>%
    mutate(Round = c(as.character(round_id)))

  overall_prev_tab_exact <- as.matrix(overall_prev_tab_exact)


  mytable <- cbind(
    sympnames[n],
    overall_prev_tab_exact[, 1, drop = FALSE],
    FormatCount(overall_prev_tab_exact[, 2:3, drop = FALSE]),
    FormatCI(FormatPrevalence(overall_prev_tab_exact[, 4:6, drop = FALSE]))
  )

  output_row <- rbind(output_row, mytable)

}

classic_df <- as.data.frame(df_round)
classic_df <- classic_df[, c("u_passcode", "sympt_cat")]
classic_df$sympt_cat <- classic_df$sympt_cat %>%
  recode("Classic COVID symptoms" = 1,
         "No symptoms" = 0,
         "Other symptoms" = 0,
         "NA" = 0 )

classic_table_exact <- as.matrix(
  overall_prev(classic_df, method = "exact", outcome="sympt_cat"))

classic_table <- cbind(
  "Classic symptoms",
  as.character(round_id),
  FormatCount(classic_table_exact[, 1:2, drop = FALSE]),
  FormatCI(FormatPrevalence(classic_table_exact[, 3:5, drop = FALSE]))
)

output_row <- rbind(classic_table, output_row)
colnames(output_row) <- c("Symptom", "Round", "Present", "Total COVID+ve", "Unweighted prevalence")


write.xlsx(output_row,
           colNames = TRUE,
           paste0(overall_file, "_", vs, "_r", round_id, "_", output_tag, ".xlsx")
)

if (direct_export) {
  file.copy(
    from = paste0(overall_file, "_r", round_id, "_", output_tag, ".xlsx"),
    to = "T:/", overwrite = TRUE
  )
}
}

  ## Prevalence stratified by covariates

  # Extracting covariate names
  # tmp <- read.xlsx(annot_file)
  # covs_names <- tmp[, 2]
  # names(covs_names) <- tmp[, 1]
  # covs_names <- covs_names[covs]
  #
  # # Removing unused variables
  # df_round <- df_round[, c(res_param, covs, weight_params)]
  #
  # # Recoding categorical variables
  # covs_to_recode <- getSheetNames(recoding_file)
  # covs_to_recode <- intersect(names(covs_names), covs_to_recode)
  # if (length(covs_to_recode) > 0) {
  #   for (i in 1:length(covs_to_recode)) {
  #     recoding <- read.xlsx(recoding_file, sheet = covs_to_recode[i])
  #     recoding[which(is.na(recoding[, 1])), 1] <- "NA"
  #     renaming <- recoding[, 2]
  #     names(renaming) <- recoding[, 1]
  #     x <- as.character(df_round[, covs_to_recode[i]])
  #     x[is.na(x)] <- "NA"
  #     x <- factor(x, levels = names(renaming), labels = renaming)
  #     df_round[, covs_to_recode[i]] <- x
  #   }
  # }
  #
  # # Recoding continuous to categorical
  # covs_to_recode <- getSheetNames(recoding_from_cont_file)
  # covs_to_recode <- intersect(names(covs_names), covs_to_recode)
  # if (length(covs_to_recode) > 0) {
  #   for (i in 1:length(covs_to_recode)) {
  #     recoding <- read.xlsx(recoding_from_cont_file, sheet = covs_to_recode[i])
  #     x <- as.numeric(df_round$age)
  #     x <- cut(x, breaks = c(min(x) - 10, recoding[, 1]), labels = recoding[, 2])
  #     df_round[, covs_to_recode[i]] <- x
  #   }
  # }
  #
  # # Specific recoding for r14/r15
  # if ("7+" %in% df_round$hh_size_cat) {
  #   df_round$hh_size_cat <- factor(df_round$hh_size_cat,
  #     levels = c(1:6, "7+"),
  #     labels = c(1:5, "6+", "6+")
  #   )
  # } else {
  #   df_round$hh_size_cat <- factor(df_round$hh_size_cat,
  #     levels = c(1:5, "6+"),
  #     labels = c(1:5, "6+")
  #   )
  # }
  #
  # # Make the prevalence tables for the above covariates using Vivi's code (unweighted)
  # system.time({
  #   mytable <- ExtractPrevalence(
  #     df_round = df_round,
  #     covs = covs, covs_names = covs_names,
  #     res_param = res_param, weighted = FALSE
  #   )
  # })
  #
  # if (weighted) {
  #   system.time({
  #     mytable_weighted <- ExtractPrevalence(
  #       df_round = df_round,
  #       covs = covs, covs_names = covs_names,
  #       res_param = res_param,
  #       weight_params = weight_params, weighted = TRUE
  #     )
  #   })
  #   tmp <- mytable
  #   tmp[rownames(mytable_weighted), 5] <- mytable_weighted[, 3]
  #   mytable <- cbind(mytable, tmp[, 5])
  # }
  #
  # # Loading template
  # wb <- loadWorkbook(template_file)
  #
  # # Checking consistency with template
  # template_required <- NULL
  # for (sheet_id in c(1, 2)) {
  #   print(sheet_id)
  #
  #   # Reading template
  #   mysheetname <- eval(parse(text = paste0("template_sheet", sheet_id)))
  #   mysheet <- readWorkbook(template_file,
  #     sheet = mysheetname
  #   )
  #
  #   # Removing legend
  #   mysheet <- mysheet[!grepl("^\\*", mysheet[, 1]), ]
  #   rownames(mysheet) <- paste0(
  #     ExtendList(mysheet[, 1]),
  #     "_",
  #     mysheet[, 2]
  #   )
  #   template_required <- c(template_required, rownames(mysheet))
  #
  #   # Checking that all variables are consistent
  #   if (any(!rownames(mysheet) %in% rownames(mytable))) {
  #     stop(paste0(
  #       "Inconsistencies in variable names: ",
  #       paste(rownames(mysheet)[!rownames(mysheet) %in% rownames(mytable)],
  #         collapse = " ; "
  #       )
  #     ))
  #   }
  #
  #   filled_sheet <- mysheet
  #   if (weighted) {
  #     filled_sheet[, 3:6] <- mytable[rownames(filled_sheet), 3:6]
  #   } else {
  #     filled_sheet[, 3:5] <- mytable[rownames(filled_sheet), 3:5]
  #   }
  #   print(filled_sheet)
  #
  #   # Checking that all observations are available for all variables
  #   mytotal <- sum(!is.na(df_round$estbinres))
  #   ids <- which(!duplicated(ExtendList(filled_sheet$Variable)))
  #   ids <- c(ids, nrow(filled_sheet) + 1)
  #   for (i in 1:(length(ids) - 1)) {
  #     mysum <- sum(as.numeric(gsub(",", "", filled_sheet$Total))[seq(ids[i], ids[i + 1] - 1)])
  #     if (mysum != mytotal) {
  #       stop(paste0(
  #         "Missing some observations for variable ",
  #         unique(na.exclude(filled_sheet$Variable))[i]
  #       ))
  #     }
  #   }
  #
  #   # Updating the sheet
  #   writeData(wb, sheet = mysheetname, x = filled_sheet, colNames = FALSE, startRow = 2)
  # }
  #
  # # Resizing the cells
  # for (sheet_id in c(1, 2)) {
  #   tmpsheet <- readWorkbook(template_file, sheet = sheet_id)
  #
  #   # Setting column widths
  #   removeColWidths(wb,
  #     sheet = sheet_id,
  #     cols = seq(3, ncol(filled_sheet))
  #   )
  #   setColWidths(wb,
  #     sheet = sheet_id,
  #     cols = seq(3, ncol(filled_sheet)),
  #     widths = column_widths,
  #     ignoreMergedCells = TRUE
  #   )
  #
  #   # Setting row heights
  #   removeRowHeights(wb,
  #     sheet = sheet_id,
  #     rows = seq(1, nrow(tmpsheet) + 1)
  #   )
  #   setRowHeights(wb,
  #     sheet = sheet_id,
  #     rows = seq(1, nrow(tmpsheet) + 1),
  #     heights = row_height
  #   )
  # }
  #
  # # Checking that all strata were in template
  # if (!all(rownames(mytable) %in% template_required)) {
  #   stop(paste0(
  #     "Stratifications not in template: ",
  #     paste(rownames(mytable)[!rownames(mytable) %in% template_required],
  #       collapse = " ; "
  #     )
  #   ))
  # }
  #
  # # Saving updated workbook
  # saveWorkbook(wb,
  #   file = paste0(output_file, "_r", round_id, "_", output_tag, ".xlsx"),
  #   overwrite = TRUE
  # )
  #
  # # Copying output to transfer folder
  # if (direct_export) {
  #   file.copy(
  #     from = paste0(output_file, "_r", round_id, "_", output_tag, ".xlsx"),
  #     to = "T:/", overwrite = TRUE
  #   )
  # }


