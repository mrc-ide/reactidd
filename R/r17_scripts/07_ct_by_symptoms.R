rm(list = ls(all = TRUE))

round_id=17

# Setting working directory
setwd(paste0("E:/Group/report/round", round_id))

# Loading required packages
source("Scripts/functions/load_packages.R")
pkgs <- c(
  "prevalence", "mgcv", "mgcViz", "MASS", "dplyr",
  "tidyr", "forcats", "ggplot2", "qpcR", "survey", "reshape2",
  "openxlsx", "colorspace"
)
load_packages(pkgs)

# Source any functions from the local file
source("Scripts/functions/add_conf_ints.R")
source("Scripts/functions/make_tables.R")
source("Scripts/functions/overall_prev.R")
source("Scripts/functions/formatting_functions.R")

# Choice of round
for (round_id in 15:17) {
  ## Parametrisation

  # Paths to files
  data_file <- paste0("E:/Group/saved_objects/rep", round_id, ".rds")
  output_file <- "Tables/Logistic_models"

  output_tag <- Sys.Date()
  annot_file <- "Parameters/Variable_names.xlsx"
  template_file <- "Parameters/OR_table.xlsx"

  recoding_file <- "Parameters/Recoding.xlsx"
  recoding_from_cont_file <- "Parameters/Recoding_from_continuous.xlsx"

  # Copying output files directly to transfer folder
  direct_export <- TRUE

  # Variable for test results
  res_param <- "estbinres"

  # Variables for stratification
  covs <- c(
    "gender_char", "age", "region",
    "work_new_alt", "ethnic_new_char",
    "hh_size_cat", "nchild2", "imd_quintile",
    "sympt_cat", "ct1", "ct2"
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

  # Restricted to COVID cases
  df_round <- df_round[which(df_round$estbinres == 1), ]

  dir.create("Figures/ct_values_covid", showWarnings = FALSE)
  {
    pdf(paste0("Figures/ct_values_covid/ct_by_symptoms_r", round_id, ".pdf"),
        width = 20, height = 8)
    par(mfrow = c(1, 2), mar = c(5, 5, 7, 1))
    for (ct_outcome in c("ct1", "ct2")) {
      mylist <- list()

      ids <- which(df_round$sympt_cat == "No symptoms")
      mylist <- c(mylist, list(df_round[ids, ct_outcome]))

      ids <- which(df_round$sympt_cat %in% c("Classic COVID symptoms", "Other symptoms"))
      mylist <- c(mylist, list(df_round[ids, ct_outcome]))

      ids <- which(df_round$sympt_cat %in% c("Classic COVID symptoms"))
      mylist <- c(mylist, list(df_round[ids, ct_outcome]))

      ids <- which(df_round$sympt_cat %in% c("Other symptoms"))
      mylist <- c(mylist, list(df_round[ids, ct_outcome]))

      ids <- which(df_round$sympt_cat %in% c("Unknown"))
      mylist <- c(mylist, list(df_round[ids, ct_outcome]))

      N0 <- formatC(sapply(mylist, FUN = function(x) {
        sum(round(x, digits = 4) == 0)
      }),
      format = "f", digits = 0, big.mark = ","
      )
      mylist <- lapply(mylist, FUN = function(x) {
        x[which(round(x, digits = 4) != 0)]
      })

      names(mylist) <- c("No symptoms", "Any symptoms", "Classic COVID symptoms", "Other symptoms", "Unknown")
      names(mylist) <- paste0(
        names(mylist), "\n (N=",
        formatC(sapply(mylist, length), format = "f", digits = 0, big.mark = ","), ")"
      )

      mycolours <- lighten(c("darkred"), amount = 0.5)
      boxplot(mylist,
        range = 0,
        boxcol = "white", col = mycolours, ylim = c(0, 50),
        staplecol = mycolours, whiskcol = mycolours, lty = 1,
        las = 1, cex.axis = 1.5, cex.main = 2, xaxt = "n",
        main = "COVID-19 infection",
        ylab = paste0(
          "Ct values (",
          ifelse(ct_outcome == "ct1", yes = "N", no = "E"), " gene)"
        ),
        cex.lab = 2
      )

      set.seed(1)
      for (i in 1:length(mylist)) {
        if (length(mylist[[i]]) > 0) {
          points(i + ProportionalJitter(mylist[[i]]), mylist[[i]],
            pch = 19, cex = 0.5,
            col = "darkred"
          )
        }
      }

      axis(side = 1, at = 1:length(mylist), labels = NA)
      axis(
        side = 1, at = 1:length(mylist), labels = names(mylist),
        line = 2, cex.axis = 0.8, tick = FALSE
      )

      axis(
        side = 3, at = 1:length(mylist),
        labels = unlist(sapply(N0, FUN = function(x) {
          eval(parse(text = paste0("expression(N['Cp=0']*'=", x, "')")))
        })),
        line = -0.5, cex.axis = 1.5, tick = FALSE
      )
    }
    dev.off()
  }


  for (ct in c("ct1", "ct2")) {
    y <- df_round[, ct]
    mymodel <- lm(y ~ df_round$sympt_cat)
    write.xlsx(summary(mymodel)$coefficients,
      row.names = TRUE, overwrite = TRUE,
      paste0("Figures/ct_values_covid/Linear_model_", ct, "_r", round_id, ".xlsx")
    )
  }
}
