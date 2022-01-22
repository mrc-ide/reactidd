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
weighted <- FALSE # whether weighted prevalences should be included or not

# Copying output files directly to transfer folder
direct_export <- TRUE

# Choice of rounds
round_id <- 16
linked=FALSE

# Paths to files
overall_file <- "E:/Group/report/round16/Tables/Overall_prevalence"
output_file <- "E:/Group/report/round16/Tables/Prevalence"

output_tag <- Sys.Date()
annot_file <- "E:/Group/report/round16/Parameters/Variable_names.xlsx"
template_file <- "E:/Group/report/round16/Parameters/Table2a"
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
if (linked){
  covs <- c(
    "age",
    "gender_char",
    "ethnic_new_char",
    "covida",
    "covidcon_char", 
    "sympt_cat",
    "region",
    "imd_quintile",
    "hh_size_cat",
    "nchild2",
    "link_vax_status_booster14days"
  )
} else {
  covs <- c(
    "age",
    "gender_char",
    "ethnic_new_char",
    "covida",
    "covidcon_char", 
    "sympt_cat",
    "region",
    "urban",
    "imd_quintile",
    "hh_size_cat",
    "nchild2"
  ) 
}

# "urban",
# Defining the column widths / row heights
column_widths <- c(5.5, 5.5, 19)
if (weighted) {
  column_widths <- c(column_widths, 19)
}
row_height <- 15

# Path to the data
if (linked){
  data_file <- paste0("E:/dt20/linkedR", round_id, "datANG.rds")
} else {
  data_file <- paste0("E:/Group/saved_objects/rep", round_id, ".rds")
}

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
    x[is.na(x)] <- "NA"
    x <- factor(x, levels = names(renaming), labels = renaming)
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
    df_round[, covs_to_recode[i]] <- x
  }
}

# Specific recoding for r14/r15
if ("7+" %in% df_round$hh_size_cat) {
  df_round$hh_size_cat <- factor(df_round$hh_size_cat,
                                 levels = c(1:6, "7+"),
                                 labels = c(1:5, "6+", "6+")
  )
} else {
  df_round$hh_size_cat <- factor(df_round$hh_size_cat,
                                 levels = c(1:5, "6+"),
                                 labels = c(1:5, "6+")
  )
}

# Identifying delta and omicron cases
df_lineage=readRDS(paste0("E:/Group/saved_objects/rep",round_id,"_lineage.rds"))
ids_omicron=df_lineage$u_passcode[which(df_lineage$react_lineage=="BA.1")]
ids_delta=df_lineage$u_passcode[which(df_lineage$react_lineage!="BA.1")]
ids_positives=rownames(df_round)[which(df_round$estbinres==1)]
ids_unknown=ids_positives[!ids_positives%in%df_lineage$u_passcode]
df_round$lineage=rep(NA,nrow(df_round))
df_round[ids_omicron, "lineage"]="omicron"
df_round[ids_delta, "lineage"]="delta"
df_round[ids_unknown, "lineage"]="unknown"

# Extracting counts/proportions and apply chi-squared test
mytable=NULL
mytable_ci_delta=mytable_ci_omicron=mytable_ci_unknown=NULL
for (k in 1:length(covs)){
  conttable=table(df_round[,covs[k]], df_round[,"lineage"])
  conttable=conttable[apply(conttable,1,sum)!=0,]
  mynames=rownames(conttable)
  chi2=chisq.test(conttable)
  
  for (i in 1:nrow(conttable)){
    mytable_ci_delta=rbind(mytable_ci_delta, 
                           propCI(x=conttable[i,1], n=sum(conttable[,1]), method="wilson")[c("p","lower","upper")])
  }
  
  for (i in 1:nrow(conttable)){
    mytable_ci_omicron=rbind(mytable_ci_omicron, 
                             propCI(x=conttable[i,2], n=sum(conttable[,2]), method="wilson")[c("p","lower","upper")])
  }
  
  for (i in 1:nrow(conttable)){
    mytable_ci_unknown=rbind(mytable_ci_unknown, 
                             propCI(x=conttable[i,3], n=sum(conttable[,3]), method="wilson")[c("p","lower","upper")])
  }
  
  conttable=matrix(paste0(conttable, " (", 
                          formatC(prop.table(conttable,2)*100, format="f", digits=2),
                          "%)"), ncol=ncol(conttable))
  rownames(conttable)=mynames
  mytable=rbind(mytable, 
                cbind(conttable[,c(2,1,3)],
                      c(formatC(chi2$p.value, format="e", digits=2), 
                        rep(NA,nrow(conttable)-1))))
}

write.xlsx(mytable, row.names=TRUE, overwrite=TRUE,
           paste0("E:/Group/report/round16/Tables/Characteristics_omicron",
                  ifelse(linked, yes="_link", no=""),".xlsx"))

write.xlsx(mytable_ci_delta, row.names=TRUE, overwrite=TRUE,
           paste0("E:/Group/report/round16/Tables/Characteristics_omicron_ci_delta",
                  ifelse(linked, yes="_link", no=""),".xlsx"))

write.xlsx(mytable_ci_omicron, row.names=TRUE, overwrite=TRUE,
           paste0("E:/Group/report/round16/Tables/Characteristics_omicron_ci_omicron",
                  ifelse(linked, yes="_link", no=""),".xlsx"))

write.xlsx(mytable_ci_unknown, row.names=TRUE, overwrite=TRUE,
           paste0("E:/Group/report/round16/Tables/Characteristics_omicron_ci_unknown",
                  ifelse(linked, yes="_link", no=""),".xlsx"))

# x=df_round[,c("age","lineage")]
# x=na.exclude(x)
# ids=which((x$age>17)&(x$age<255))
# x=x[ids,]
# t.test(x$age[which(x$lineage=="delta")], y=x$age[which(x$lineage=="omicron")])
# mean(x$age[which(x$lineage=="delta")])
# mean(x$age[which(x$lineage=="omicron")])

