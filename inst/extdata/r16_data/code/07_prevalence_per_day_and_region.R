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
covs <- c(
  "gender_char", "age", "region",
  "work_new_alt", "ethnic_new_char",
  "hh_size_cat", "covidcon_char", "sympt_cat",
  "nchild2", "covida",
  "imd_quintile", "vax_status_noDate_v2",
  "d_comb"
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


## Overall prevalence

# Unweighted
overall_prev_tab_exact <- bind_rows(
  overall_prev(df_round, method = "exact"),
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

# Weighted
if (weighted) {
  dclus15g <- svydesign(id = ~id, strata = ~lacode, weights = ~wt_antigen, data = df_round, nest = TRUE)
  wt_prev_o_r15 <- svyby(~estbinres, by = ~group, design = dclus15g, FUN = svyciprop, vartype = "ci") %>% rename(level = group)
  
  wt_prev_tab <- bind_rows(wt_prev_o_r15,
                           .id = "round"
  ) %>%
    rename(
      wt_prev = estbinres,
      lower = ci_l,
      upper = ci_u
    ) %>%
    mutate(round = c("15"))
  
  mytable <- cbind(
    mytable,
    FormatCI(FormatPrevalence(wt_prev_tab[, 3:5, drop = FALSE]))
  )
  colnames(mytable)[ncol(mytable)] <- "Weighted prevalence"
}

write.xlsx(mytable,
           colNames = TRUE,
           paste0(overall_file, "_r", round_id, "_", output_tag, ".xlsx")
)

if (direct_export) {
  file.copy(
    from = paste0(overall_file, "_r", round_id, "_", output_tag, ".xlsx"),
    to = "T:/", overwrite = TRUE
  )
}

## Prevalence stratified by covariates

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

full_dataset=df_round
full_dataset$total=factor(rep(1,nrow(full_dataset)))

table(full_dataset$d_comb)
table(full_dataset$d_comb, full_dataset$estbinres)

plot(table(full_dataset$d_comb))

tab1=table(full_dataset$d_comb, full_dataset$region)
tab2=table(full_dataset$d_comb, full_dataset$region, full_dataset$estbinres)[,,1]
tab3=table(full_dataset$d_comb, full_dataset$region, full_dataset$estbinres)[,,2]
write.xlsx(tab1, "E:/Group/report/round16/Tables/Counts_full.xlsx")
write.xlsx(tab2, "E:/Group/report/round16/Tables/Counts_0.xlsx")
write.xlsx(tab3, "E:/Group/report/round16/Tables/Counts_1.xlsx")

mydates=unique(full_dataset$d_comb)
mydates=mydates[!is.na(mydates)]
mydates=sort(mydates)

wb=createWorkbook()
for (k in 1:length(mydates)){
  mydate=mydates[k]
  
  # Extracting data for corresponding date
  df_round=full_dataset[which(full_dataset$d_comb==mydate),]
  
  # Make the prevalence tables for the above covariates using Vivi's code (unweighted)
  system.time({
    mytable <- ExtractPrevalence(
      df_round = df_round,
      covs = "region", covs_names = "Region",
      res_param = res_param, weighted = weighted
    )
  })
  
  mytable=mytable[,2:ncol(mytable)]
  addWorksheet(wb, sheetName = k)
  writeData(wb, sheet=k, x=mytable)
}
saveWorkbook(wb, file="E:/Group/report/round16/Tables/Daily_prevalence.xlsx")

ids=readRDS("E:/Group/report/round16/Tables/omicron_cases.rds")
mytable=full_dataset[ids,c("age","gender_char","covida","covidcon_char","sympt_cat","region","nchild2")]

mytable=cbind(names(do.call(c,apply(mytable,2,table))),
paste0(do.call(c,apply(mytable,2,table)), " (",
formatC(do.call(c,apply(mytable,2,FUN=function(x){prop.table(table(x))}))*100,
        format="f", digits=2), "%)"))
write.xlsx(mytable, "E:/Group/report/round16/Tables/omicron.xlsx")

