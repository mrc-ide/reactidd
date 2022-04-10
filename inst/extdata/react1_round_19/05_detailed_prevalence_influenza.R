rm(list = ls(all = TRUE))

# Choice of rounds
round_ids <- c(17)
round_id=round_ids[1]

# Setting working directory
setwd(paste0("E:/Group/report/round17"))

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

# Modelling options
weighted <- FALSE # whether weighted prevalences should be included or not

# Copying output files directly to transfer folder
direct_export <- TRUE

# Paths to files
overall_file <- "Tables/Overall_prevalence"
output_file <- "Tables/Prevalence"

output_tag <- Sys.Date()
annot_file <- "Parameters/Variable_names.xlsx"
template_file <- "Parameters/Table2a_flu_age"
template_sheet1 <- "Table2a"
template_sheet2 <- "Table2b"

recoding_file <- "Parameters/Recoding.xlsx"
recoding_from_cont_file <- "Parameters/Recoding_from_continuous_flu_age.xlsx"

# Fetching corresponding template
if (weighted) {
  template_file <- paste0(template_file, "_weighted.xlsx")
  output_file <- paste0(output_file, "_weighted")
} else {
  template_file <- paste0(template_file, "_unweighted.xlsx")
  output_file <- paste0(output_file, "_unweighted")
}

age_lowerbound=18

# Variable for test results
output_file_tmp=output_file
overall_file_tmp=overall_file
for (vaccinated in c("all","vaccinated", "unvaccinated")){
  for (res_param in c("influenzaa", "influenzab", "influenza", "influenzaab")){
    
    output_file <- paste0(output_file_tmp, "_", res_param, "_over_", age_lowerbound, "_", vaccinated)
    overall_file=paste0(overall_file_tmp, "_", res_param, "_over_", age_lowerbound, "_", vaccinated)
    
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
      "nchild2",
      "imd_quintile", "vax_status_noDate_v2"
    )
    
    # Defining the column widths / row heights
    column_widths <- c(5.5, 5.5, 19)
    if (weighted) {
      column_widths <- c(column_widths, 19)
    }
    row_height <- 15
    
    
    ## Computing the prevalences
    
    for (round_id in round_ids) {
      # Path to the data
      if (round_id==16){
        data_file <- paste0("E:/Group/saved_objects/rep", round_id, "_flu.rds")
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
      
      if (vaccinated=="vaccinated"){
        ids=which(df_round$fluvacc=="1")
        df_round=df_round[ids,]
      } 
      
      if (vaccinated=="unvaccinated"){
        ids=which(df_round$fluvacc%in%c("2","3"))
        df_round=df_round[ids,]
      }
      
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
      df_round$influenzaab=ifelse(df_round$influenzaa+df_round$influenzab==2, yes=1, no=0)
      
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
      
      # Filtering on age
      df_round=df_round[which(df_round$age>age_lowerbound),]
      
      ## Overall prevalence
      
      # Unweighted
      overall_prev_tab_exact <- bind_rows(
        overall_prev(df_round, method = "exact", outcome = res_param),
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
      
      # Make the prevalence tables for the above covariates using Vivi's code (unweighted)
      system.time({
        mytable <- ExtractPrevalence(
          df_round = df_round,
          covs = covs, covs_names = covs_names,
          res_param = res_param, weighted = FALSE
        )
      })
      if (weighted) {
        system.time({
          mytable_weighted <- ExtractPrevalence(
            df_round = df_round,
            covs = covs, covs_names = covs_names,
            res_param = res_param,
            weight_params = weight_params, weighted = TRUE
          )
        })
        tmp <- mytable
        tmp[rownames(mytable_weighted), 5] <- mytable_weighted[, 3]
        mytable <- cbind(mytable, tmp[, 5])
      }
      
      # Loading template
      wb <- loadWorkbook(template_file)
      
      # Checking consistency with template
      template_required <- NULL
      for (sheet_id in c(1, 2)) {
        print(sheet_id)
        
        # Reading template
        mysheetname <- eval(parse(text = paste0("template_sheet", sheet_id)))
        mysheet <- readWorkbook(template_file,
                                sheet = mysheetname
        )
        
        # Removing legend
        mysheet <- mysheet[!grepl("^\\*", mysheet[, 1]), ]
        rownames(mysheet) <- paste0(
          ExtendList(mysheet[, 1]),
          "_",
          mysheet[, 2]
        )
        template_required <- c(template_required, rownames(mysheet))
        
        # Checking that all variables are consistent
        if (any(!rownames(mysheet) %in% rownames(mytable))) {
          stop(paste0(
            "Inconsistencies in variable names: ",
            paste(rownames(mysheet)[!rownames(mysheet) %in% rownames(mytable)],
                  collapse = " ; "
            )
          ))
        }
        
        filled_sheet <- mysheet
        if (weighted) {
          filled_sheet[, 3:6] <- mytable[rownames(filled_sheet), 3:6]
        } else {
          filled_sheet[, 3:5] <- mytable[rownames(filled_sheet), 3:5]
        }
        print(filled_sheet)
        
        # Checking that all observations are available for all variables
        mytotal <- sum(!is.na(df_round[,res_param]))
        ids <- which(!duplicated(ExtendList(filled_sheet$Variable)))
        ids <- c(ids, nrow(filled_sheet) + 1)
        for (i in 1:(length(ids) - 1)) {
          mysum <- sum(as.numeric(gsub(",", "", filled_sheet$Total))[seq(ids[i], ids[i + 1] - 1)])
          if (mysum != mytotal) {
            stop(paste0(
              "Missing some observations for variable ",
              unique(na.exclude(filled_sheet$Variable))[i]
            ))
          }
        }
        
        # Updating the sheet
        writeData(wb, sheet = mysheetname, x = filled_sheet, colNames = FALSE, startRow = 2)
      }
      
      # Resizing the cells
      for (sheet_id in c(1, 2)) {
        tmpsheet <- readWorkbook(template_file, sheet = sheet_id)
        
        # Setting column widths
        removeColWidths(wb,
                        sheet = sheet_id,
                        cols = seq(3, ncol(filled_sheet))
        )
        setColWidths(wb,
                     sheet = sheet_id,
                     cols = seq(3, ncol(filled_sheet)),
                     widths = column_widths,
                     ignoreMergedCells = TRUE
        )
        
        # Setting row heights
        removeRowHeights(wb,
                         sheet = sheet_id,
                         rows = seq(1, nrow(tmpsheet) + 1)
        )
        setRowHeights(wb,
                      sheet = sheet_id,
                      rows = seq(1, nrow(tmpsheet) + 1),
                      heights = row_height
        )
      }
      
      # Checking that all strata were in template
      if (!all(rownames(mytable) %in% template_required)) {
        stop(paste0(
          "Stratifications not in template: ",
          paste(rownames(mytable)[!rownames(mytable) %in% template_required],
                collapse = " ; "
          )
        ))
      }
      
      # Saving updated workbook
      saveWorkbook(wb,
                   file = paste0(output_file, "_r", round_id, "_", output_tag, ".xlsx"),
                   overwrite = TRUE
      )
      
      # Copying output to transfer folder
      if (direct_export) {
        file.copy(
          from = paste0(output_file, "_r", round_id, "_", output_tag, ".xlsx"),
          to = "T:/", overwrite = TRUE
        )
      }
    }
  }
}
