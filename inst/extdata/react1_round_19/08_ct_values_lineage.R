rm(list = ls(all = TRUE))


# setup -------------------------------------------------------------------

# Choice of rounds
round_ids <- c(18)
round_id=round_ids[1]

# Setting working directory
setwd(paste0("E:/Group/report/round",round_id))

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


## Parametrisation

# Paths to files
data_file <- paste0("E:/Group/saved_objects/rep", round_id, "_lineage", ".rds")
vaccination_file <- paste0("E:/dt20/linkedR", round_id, "datANG.rds")
output_file <- "Tables/Prevalence_influenza"
output_tag <- Sys.Date()

annot_file <- "Parameters/Variable_names.xlsx"
recoding_file <- "Parameters/Recoding.xlsx"
recoding_from_cont_file <- "Parameters/Recoding_from_continuous.xlsx"

# Modelling options
weighted <- FALSE # whether weighted prevalences should be included or not

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
  "gender_char", "age", "region",
  "work_new_alt", "ethnic_new_char",
  "hh_size_cat", "covidcon_char", "sympt_cat",
  "nchild2",
  "imd_quintile", "vax_status_noDate_v2",
  "ct1", "ct2", "u_passcode",
  "res", "react_lineage"
)


# dataprep -----------------------------------------------------------------
## Loading and preparing the data



# Loading the data
df_round <- data.frame(readRDS(data_file))
vax_df <- data.frame(readRDS(vaccination_file))

rownames(df_round)=df_round$u_passcode

# Adding variable introduced in round 15
if (!"vax_status_noDate_v2" %in% colnames(df_round)) {
  df_round$vax_status_noDate_v2 <- df_round$vax_status_noDate
}

# Removing missing in estbinres/flu test results
df_round <- df_round %>%
  filter(!is.na(estbinres))

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

linkedCols <- c("u_passcode", "link_vaxdose",
"link_vax1date","link_vax2date","link_vax3date","link_vax1type"  ,
"link_vax2type","link_boostertype","unlinked","link_vax1toswab"  ,
"link_vaxrecenttoswab","link_vax3toswab","diff_t1_t2","diff_t2_t3"    ,
"link_vax_status_booster0days","link_vax_status_booster7days" , "link_vax_status_booster14days",
"link_v1time","link_v1tov2","link_v2tov3"      )

vax_df <- vax_df[, linkedCols]


df_round <- merge(df_round, vax_df, by = "u_passcode", all.x= TRUE)

df_round <- df_round %>%
  mutate(anySym = case_when(
    sympt_cat %in% c("Classic COVID symptoms", "Other symptoms") ~ "Any symptoms",
    sympt_cat == "No symptoms" ~ "No symptoms",
    sympt_cat == "Unknown" ~ "Unknown"
    ))


KtestVariants=function(mylist){
  pvalues=NULL

  pval1 = kruskal.test(x = c(mylist[[1]], mylist[[2]], mylist[[3]], mylist_1[[4]], recursive = TRUE),
                 g = rep(c(1,1,2,2),times = c(length(mylist_1[[1]]), length(mylist_1[[2]]), length(mylist_1[[3]]), length(mylist_1[[4]]))))$p.value

  pval2 = kruskal.test(x = c(mylist[[1]], mylist[[2]], mylist[[5]], mylist_1[[6]], recursive = TRUE),
                 g = rep(c(1,1,2,2),times = c(length(mylist_1[[1]]), length(mylist_1[[2]]), length(mylist_1[[5]]), length(mylist_1[[6]]))))$p.value

  pval3 = kruskal.test(x = c(mylist[[3]], mylist[[4]], mylist[[5]], mylist_1[[6]], recursive = TRUE),
                 g = rep(c(1,1,2,2),times = c(length(mylist_1[[3]]), length(mylist_1[[4]]), length(mylist_1[[5]]), length(mylist_1[[6]]))))$p.value

  pvalues = c(pval1, pval2, pval3)


  return(pvalues)
}


# Violinplots (Symptoms) -----------------------------------------------------------


mycolours <- rep(c("darkred", "navy"),3)
mycolours2 <- rep(c("darkorange", "navy"),3)

{
  pdf(paste0("E:/Group/report/round18/Figures/", "Ct_values_byLineage_classicSym",
             "_",  Sys.Date(), ".pdf"),
      height = 6,
      width = 9.708)

  par(mfrow=c(1,2))
  for (ct_outcome in c("ct1", "ct2")){

    mylist_1 <- list()

    for (sympStatus in c("Classic COVID symptoms","No symptoms")){
      ids <- which(df_round$sympt_cat == sympStatus &
                     df_round$react_lineage == "BA.1")
      mylist_1 <- c(mylist_1, list(df_round[ids, ct_outcome]))
    }

    for (sympStatus in c("Classic COVID symptoms","No symptoms")){
      ids <- which(df_round$sympt_cat == sympStatus &
                     df_round$react_lineage == "BA.1.1")
      mylist_1 <- c(mylist_1, list(df_round[ids, ct_outcome]))
    }

    for (sympStatus in c("Classic COVID symptoms","No symptoms")){
      ids <- which(df_round$sympt_cat == sympStatus &
                     df_round$react_lineage == "BA.2")
      mylist_1 <- c(mylist_1, list(df_round[ids, ct_outcome]))
    }

    N0 <- formatC(sapply(mylist_1, FUN = function(x) {
      sum(round(x, digits = 4) == 0)
    }),
    format = "f", digits = 0, big.mark = ","
    )
    mylist_1 <- lapply(mylist_1, FUN = function(x) {
      x[which(round(x, digits = 4) != 0)]
    })

    ViolinPlot(mylist = mylist_1, mycolours = mycolours, tick_length = 0.1)

    if (ct_outcome == "ct1"){
      mtext(expression(C[t]~Value~(N~Gene)), side = 2, cex = 1, line = 2)
    } else if (ct_outcome == "ct2") {
      mtext(expression(C[t]~Value~(E~Gene)), side = 2, cex = 1, line = 2)
    }

    print(KtestVariants(mylist_1))

    p_vals = NonParamTest(mylist_1)
    # mtext("Classic Symptoms vs Not Symptomatic", side = 3, cex = 1.5, line = 1)
    axis(side=1, at=c(1,2,3), labels=NA)
    axis(side=1, at=c(1,2,3), tick = FALSE, line = 1,
         labels = c(paste0("BA.1"  , "\n", round(p_vals[1],3)),
                    paste0("BA.1.1", "\n", if(round(p_vals[2],3)==0){"<0.001"}else{round(p_vals[2],3)}),
                    paste0("BA.2"  , "\n", round(p_vals[3],3))))
  }
  title(main=print("Classic Symptoms vs Not Symptomatic"),out=T, line=-2, cex=1.5)

  for (ct_outcome in c("ct1", "ct2")){

    mylist_1 <- list()



    for (sympStatus in c("Any symptoms","No symptoms")){
      ids <- which(df_round$anySym == sympStatus &
                     df_round$react_lineage == "BA.1")
      mylist_1 <- c(mylist_1, list(df_round[ids, ct_outcome]))
    }

    for (sympStatus in c("Any symptoms","No symptoms")){
      ids <- which(df_round$anySym == sympStatus &
                     df_round$react_lineage == "BA.1.1")
      mylist_1 <- c(mylist_1, list(df_round[ids, ct_outcome]))
    }

    for (sympStatus in c("Any symptoms","No symptoms")){
      ids <- which(df_round$anySym == sympStatus &
                     df_round$react_lineage == "BA.2")
      mylist_1 <- c(mylist_1, list(df_round[ids, ct_outcome]))
    }

    N0 <- formatC(sapply(mylist_1, FUN = function(x) {
      sum(round(x, digits = 4) == 0)
    }),
    format = "f", digits = 0, big.mark = ","
    )
    mylist_1 <- lapply(mylist_1, FUN = function(x) {
      x[which(round(x, digits = 4) != 0)]
    })

    ViolinPlot(mylist = mylist_1, mycolours = mycolours2, tick_length = 0.1)

    if (ct_outcome == "ct1"){
      mtext(expression(C[t]~Value~(N~Gene)), side = 2, cex = 1, line = 2)
    } else if (ct_outcome == "ct2") {
      mtext(expression(C[t]~Value~(E~Gene)), side = 2, cex = 1, line = 2)
    }

    print(KtestVariants(mylist_1))

    p_vals = NonParamTest(mylist_1)
    # mtext("Any Symptoms vs Not Symptomatic", side = 3, cex = 1.5, line = 1)
    axis(side=1, at=c(1,2,3), labels=NA)
    axis(side=1, at=c(1,2,3), tick = FALSE, line = 1,
         labels = c(paste0("BA.1"  , "\n", round(p_vals[1],3)),
                    paste0("BA.1.1", "\n",  if(round(p_vals[2],3)==0){"<0.001"}else{round(p_vals[2],3)}),
                    paste0("BA.2"  , "\n", round(p_vals[3],3))))
  }
  title(main=print("Any Symptoms vs Not Symptomatic"),out=T, line=-2, cex.main=1.5)

}

dev.off()
