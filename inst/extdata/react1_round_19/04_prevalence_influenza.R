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
  "openxlsx", "ggvenn"
)
load_packages(pkgs)

# Source any functions from the local file
source("Scripts/functions/add_conf_ints.R")
source("Scripts/functions/make_tables.R")
source("Scripts/functions/overall_prev.R")
source("Scripts/functions/formatting_functions.R")


## Parametrisation

# Paths to files
data_file <- paste0("E:/Group/saved_objects/rep", round_id, ".rds")
output_file <- "Tables/Prevalence_influenza"
output_tag <- Sys.Date()

annot_file <- "Parameters/Variable_names.xlsx"
recoding_file <- "Parameters/Recoding.xlsx"
recoding_from_cont_file <- "Parameters/Recoding_from_continuous_flu_age.xlsx"

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

# Adding variable introduced in round 15
if (!"vax_status_noDate_v2" %in% colnames(df_round)) {
  df_round$vax_status_noDate_v2 <- df_round$vax_status_noDate
}

# Recoding influenza infection
for (mytest in c("influenzaa", "influenzab")){
  df_round[,mytest]=as.numeric(as.character(factor(df_round[,mytest], 
                                                   levels=c("negative", "positive"), 
                                                   labels=c(0,1))))
}
df_round$influenza=ifelse(df_round$influenzaa+df_round$influenzab>0, yes=1, no=0)

# Removing missing in estbinres/flu test results
df_round <- df_round %>%
  filter(!is.na(estbinres)) %>%
  mutate(group = "Overall")
df_round <- df_round %>%
  filter(!is.na(influenzaa)) %>%
  mutate(group = "Overall")
df_round <- df_round %>%
  filter(!is.na(influenzab)) %>%
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

# Extracting covariate names
tmp <- read.xlsx(annot_file)
covs_names <- tmp[, 2]
names(covs_names) <- tmp[, 1]
covs_names <- covs_names[covs]

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
    print(table(x))
    x[is.na(x)] <- "NA"
    x <- factor(x, levels = names(renaming), labels = renaming)
    print(table(x))
    df_round[, covs_to_recode[i]] <- x
  }
}

# Recoding continuous to categorical
covs_to_recode <- getSheetNames(recoding_from_cont_file)
covs_to_recode <- intersect(names(covs_names), covs_to_recode)
if (length(covs_to_recode) > 0) {
  for (i in 1:length(covs_to_recode)) {
    recoding <- read.xlsx(recoding_from_cont_file, sheet = covs_to_recode[i])
    x <- as.numeric(df_round$age)
    x <- cut(x, breaks = c(min(x) - 10, recoding[, 1]), labels = recoding[, 2])
    print(table(x))
    df_round[, covs_to_recode[i]] <- x
  }
}

# Creating combinations of outcomes
df_round$influenzaa=as.numeric(as.character(df_round$influenzaa))
df_round$influenzab=as.numeric(as.character(df_round$influenzab))
df_round$flua_only=ifelse((df_round$influenzaa==1)&(df_round$influenzab==0), yes=1, no=0)
df_round$flub_only=ifelse((df_round$influenzaa==0)&(df_round$influenzab==1), yes=1, no=0)
df_round$covid_and_flua=ifelse(df_round$estbinres+df_round$influenzaa==2, 
                               yes=1, no=0)
df_round$covid_and_flub=ifelse(df_round$estbinres+df_round$influenzab==2, 
                               yes=1, no=0)
df_round$flu=ifelse(df_round$influenzaa+df_round$influenzab>=1, 
                    yes=1, no=0)
df_round$covid_and_flu=ifelse(df_round$estbinres+df_round$flu==2, 
                              yes=1, no=0)
df_round$flua_and_flub=ifelse(df_round$influenzaa+df_round$influenzab==2, 
                              yes=1, no=0)
df_round$covid_and_flua_and_flub=ifelse(df_round$estbinres+df_round$influenzaa+df_round$influenzab==3, 
                                        yes=1, no=0)

write.table(cbind(df_round$u_passcode[which((df_round$estbinres==1)&(df_round$influenzaa==1)&(df_round$influenzab==1))]),
            quote=FALSE, "Tables/List_triple_infected.txt", row.names = FALSE, col.names = FALSE)
write.table(cbind(df_round$u_passcode[which((df_round$estbinres==1)&(df_round$influenzaa==0)&(df_round$influenzab==1))]),
            quote=FALSE, "Tables/List_double_infected_B.txt", row.names = FALSE, col.names = FALSE)

nrow(df_round)
table(df_round$age%in%c("05-18"))
table(df_round$fluvacc)
table(df_round$fluvacc, df_round$age%in%c("05-18"))

df_round_full=df_round
for (stratum in c("all", "Yes", "No", "Unknown")){
  if (stratum!="all"){
    df_round=df_round_full[which(df_round_full$fluvacc==stratum),]
  }
  print(stratum)
  print(nrow(df_round))
  print(sum(df_round$age%in%c("05-18")))
  print(sum(!df_round$age%in%c("05-18")))
  
  # pdf(paste0("Figures/Venn_diagram_vacc_",tolower(stratum),".pdf"))
  ggvenn(list(`COVID-19`=which(df_round$estbinres==1),
              `Influenza A`=which(df_round$influenzaa==1),
              `Influenza B`=which(df_round$influenzab==1)),
         fill_color = c("royalblue", "tomato", "forestgreen"),
         stroke_color = NA,
         text_size=7,
         set_name_size=8,
         show_percentage = FALSE)
  ggsave(paste0("Figures/Venn_diagram_vacc_",tolower(stratum),".pdf"))
  # dev.off()
  
  # pdf(paste0("Figures/Venn_diagram_vacc_",tolower(stratum),"_below_18.pdf"))
  ggvenn(list(`COVID-19`=which((df_round$estbinres==1)&(df_round$age%in%c("05-18"))),
              `Influenza A`=which((df_round$influenzaa==1)&(df_round$age%in%c("05-18"))),
              `Influenza B`=which((df_round$influenzab==1)&(df_round$age%in%c("05-18")))),
         fill_color = c("royalblue", "tomato", "forestgreen"),
         stroke_color = NA,
         text_size=7,
         set_name_size=8,
         show_percentage = FALSE)
  ggsave(paste0("Figures/Venn_diagram_vacc_",tolower(stratum),"_below_18.pdf"))
  # dev.off()
  
  # pdf(paste0("Figures/Venn_diagram_vacc_",tolower(stratum),"_above_18.pdf"))
  ggvenn(list(`COVID-19`=which((df_round$estbinres==1)&(!df_round$age%in%c("05-18"))),
              `Influenza A`=which((df_round$influenzaa==1)&(!df_round$age%in%c("05-18"))),
              `Influenza B`=which((df_round$influenzab==1)&(!df_round$age%in%c("05-18")))),
         fill_color = c("royalblue", "tomato", "forestgreen"),
         stroke_color = NA,
         text_size=7,
         set_name_size=8,
         show_percentage = FALSE)
  ggsave(paste0("Figures/Venn_diagram_vacc_",tolower(stratum),"_above_18.pdf"))
  # dev.off()
}
df_round=df_round_full


## Stratified by COVID symptoms

full_dataset=df_round
age_class=c("all_ages", "below_18", "above_18")
for (j in 1:3){
  if (j==1){
    ids=1:nrow(full_dataset)
  }
  if (j==2){
    ids=which(full_dataset$age%in%c("05-11","12-17"))
  } 
  if (j==3){
    ids=which(!full_dataset$age%in%c("05-11","12-17"))
  }
  df_round=full_dataset[ids,]
  
  for (k in 1:2){
    if (k==1){
      sympt_status_list=c("all", 
                          "Classic COVID symptoms", 
                          "Other symptoms",
                          "No symptoms",
                          "Unknown")
    } else {
      df_round$sympt_cat=factor(df_round$sympt_cat,
                                levels=c("Classic COVID symptoms",
                                         "Other symptoms",
                                         "No symptoms",
                                         "Unknown"),
                                labels=c("Any symptoms",
                                         "Any symptoms",
                                         "No symptoms",
                                         "Unknown"))
      
      sympt_status_list=c("all", 
                          "Any symptoms",
                          "No symptoms",
                          "Unknown")
    }
    
    outcome_vars=c("COVID-19",
                   "Influenza A",
                   "Influenza B",
                   "Influenza A or B",
                   "Influenza A and B",
                   "COVID-19 and influenza A",
                   "COVID-19 and influenza B",
                   "COVID-19 and influenza A or B",
                   "COVID-19 and influenza A and B")
    
    names(outcome_vars)=c("estbinres", 
                          "influenzaa",
                          "influenzab",
                          "flu",
                          "flua_and_flub",
                          "covid_and_flua",
                          "covid_and_flub",
                          "covid_and_flu",
                          "covid_and_flua_and_flub")
    
    myfulltable=NULL
    for (sympt_status in sympt_status_list){
      # Filtering based on symptoms status
      if (sympt_status=="all"){
        mydata=df_round
      } else {
        mydata=df_round[which(df_round$sympt_cat==sympt_status),]
      }
      
      mybigtable=NULL
      for (outcome_var in names(outcome_vars)){
        overall_prev_tab_exact <- bind_rows(
          overall_prev(mydata, method = "exact", outcome=outcome_var),
          .id = "Round"
        ) %>%
          mutate(Round = c(as.character(round_id)))
        overall_prev_tab_exact <- as.matrix(overall_prev_tab_exact)
        
        mytable <- cbind(
          overall_prev_tab_exact[, 1, drop = FALSE],
          FormatCount(overall_prev_tab_exact[, 2:3, drop = FALSE]),
          FormatCI(FormatPrevalence(overall_prev_tab_exact[, 4:6, drop = FALSE]))
        )
        colnames(mytable) <- c("Round", "Positive", "Total", "Unweighted prevalence")
        mytable=mytable[,-1,drop=FALSE]
        
        mybigtable=rbind(mybigtable, mytable)
      }
      rownames(mybigtable)=outcome_vars
      myfulltable=cbind(myfulltable, mybigtable)
    }
    
    # Saving updated workbook
    write.xlsx(as.data.frame(myfulltable), 
               col.names=TRUE, row.names=TRUE,
               file = paste0(output_file, "_r", round_id, "_", k, "_", age_class[j], "_", output_tag, ".xlsx"),
               overwrite = TRUE
    )
    
    # Copying output to transfer folder
    if (direct_export) {
      file.copy(
        from = paste0(output_file, "_r", round_id, "_", k, "_", age_class[j], "_", output_tag, ".xlsx"),
        to = "T:/", overwrite = TRUE
      )
    }
  }
}
