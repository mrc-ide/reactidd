rm(list = ls(all = TRUE))
setwd("E:/Group/report/round16/Scripts/")

# Loading required packages
source("functions/load_packages.R")
pkgs <- c(
  "prevalence", "mgcv", "mgcViz", "MASS", "dplyr",
  "tidyr", "forcats", "ggplot2", "qpcR", "survey", "reshape2",
  "openxlsx", "colorspace"
)
load_packages(pkgs)

# Source any functions from the local file
source("functions/add_conf_ints.R")
source("functions/make_tables.R")
source("functions/overall_prev.R")
source("functions/formatting_functions.R")


## Parametrisation

# Paths to files
round_id <- 16
data_file <- paste0("E:/Group/saved_objects/rep", round_id, "_flu.rds")
output_file <- "E:/Group/report/round16/Tables/Prevalence_influenza"
output_tag <- Sys.Date()

annot_file <- "E:/Group/report/round16/Parameters/Variable_names.xlsx"
recoding_file <- "E:/Group/report/round16/Parameters/Recoding.xlsx"
recoding_from_cont_file <- "E:/Group/report/round16/Parameters/Recoding_from_continuous.xlsx"

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
  "gender_char", "age", "region",
  "work_new_alt", "ethnic_new_char",
  "hh_size_cat", "covidcon_char", "sympt_cat",
  "nchild2",
  "imd_quintile", "vax_status_noDate_v2",
  "influenzaa", "influenzab", "influenza",
  "ct1", "ct2",
  "res",
  "influenzaacpvalue", "influenzabcpvalue",
  "fluvacc"
)


## Loading and preparing the data

# Loading the data
df_round <- data.frame(readRDS(data_file))
rownames(df_round)=df_round$u_passcode

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
df_round$influenzaa <- as.numeric(as.character(df_round$influenzaa))
df_round$influenzab <- as.numeric(as.character(df_round$influenzab))
df_round$covid_and_flua <- ifelse(df_round$estbinres + df_round$influenzaa == 2,
                                  yes = 1, no = 0
)
df_round$covid_and_flub <- ifelse(df_round$estbinres + df_round$influenzab == 2,
                                  yes = 1, no = 0
)
df_round$flu <- ifelse(df_round$influenzaa + df_round$influenzab >= 1,
                       yes = 1, no = 0
)
df_round$covid_and_flu <- ifelse(df_round$estbinres + df_round$flu == 2,
                                 yes = 1, no = 0
)
df_round$flua_and_flub <- ifelse(df_round$influenzaa + df_round$influenzab == 2,
                                 yes = 1, no = 0
)
df_round$covid_and_flua_and_flub <- ifelse(df_round$estbinres + df_round$influenzaa + df_round$influenzab == 3,
                                           yes = 1, no = 0
)

df_round$any_sympt=factor(df_round$sympt_cat,
                          levels=c("Classic COVID symptoms",
                                   "Other symptoms",
                                   "No symptoms",
                                   "Unknown"),
                          labels=c("S+","S+","S-","S?"))

df_round$fluvacc=factor(df_round$fluvacc, 
                        levels=c("Yes", "No", "Unknown"),
                        labels=c("V+", "V-", "V?"))

df_round$vacc_and_sympt=paste0(as.character(df_round$fluvacc),
                               "/",
                               as.character(df_round$any_sympt))

df_round$vacc_and_sympt=factor(df_round$vacc_and_sympt, 
                               levels=c("V+/S-", "V+/S+", "V+/S?",
                                        "V-/S-", "V-/S+", "V-/S?",
                                        "V?/S-", "V?/S+", "V?/S?"))

df_round$age_binary=ifelse(df_round$age%in%c("05-12", "13-17"),
                           yes="<18", no="18+")

# for (outcome in c("influenzaa", "influenzab", "influenza")){
#   table(df_round[,outcome], df_round$vacc_and_sympt, df_round$age_binary)
# }
mytable=rbind(table(df_round[,outcome], df_round$vacc_and_sympt),
              table(df_round[,outcome], df_round$vacc_and_sympt, df_round$age_binary)[,,1],
              table(df_round[,outcome], df_round$vacc_and_sympt, df_round$age_binary)[,,2])
write.xlsx(as.data.frame(mytable), 
           row.names=TRUE,
           paste0("../Tables/Contingency_table_influenza_vacc_and_sympt_", Sys.Date(), ".xlsx"))

dir.create("../Figures", showWarnings = FALSE)

{
  pdf(paste0("../Figures/Ct_values_infection_covid_r", round_id, "_", Sys.Date(), ".pdf"),
      width = 12, height = 7
  )
  par(mfrow = c(1, 2), mar = c(5, 5, 7, 1))
  for (ct_outcome in c("ct1", "ct2")) {
    mylist <- list()
    
    # Negative to flu
    ids <- which(df_round$estbinres == 0)
    mylist <- c(mylist, list(df_round[ids, ct_outcome]))
    
    # Negative to flu
    ids <- which(df_round$estbinres == 1)
    mylist <- c(mylist, list(df_round[ids, ct_outcome]))
    
    N0 <- formatC(sapply(mylist, FUN = function(x) {
      sum(round(x, digits = 4) == 0)
    }),
    format = "f", digits = 0, big.mark = ","
    )
    mylist <- lapply(mylist, FUN = function(x) {
      x[which(round(x, digits = 4) != 0)]
    })
    
    names(mylist) <- c("Negative", "Positive")
    names(mylist) <- paste0(
      names(mylist), "\n (N=",
      formatC(sapply(mylist, length), format = "f", digits = 0, big.mark = ","), ")"
    )
    
    mycolours <- lighten(c("grey30", "darkred"), amount = 0.5)
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
    
    # set.seed(1)
    # stripchart(mylist,
    #            pch=19, cex=0.5,
    #            col="darkred",
    #            method="jitter",
    #            vertical=TRUE,
    #            add=TRUE)
    
    axis(side = 1, at = 1:length(mylist), labels = NA)
    axis(
      side = 1, at = 1:length(mylist), labels = names(mylist),
      line = 2, cex.axis = 1.5, tick = FALSE
    )
    
    axis(
      side = 3, at = 1:length(mylist),
      # labels = paste0("N0=", N0),
      labels=unlist(sapply(N0, FUN=function(x){eval(parse(text=paste0("expression(N['Cp=0']*'=",x, "')")))})),
      line = -0.5, cex.axis = 1.5, tick = FALSE
    )
  }
  dev.off()
}


{
  pdf(paste0("../Figures/Ct_values_infection_detection_", round_id, "_", Sys.Date(), ".pdf"), width = 7, height = 7)
  par(mar = c(5, 5, 5, 1))
  plot(df_round$ct1, df_round$ct2,
       pch = 19, cex = 0.8, las = 1, cex.lab = 1.5,
       xlab = "Ct value (N gene)",
       ylab = "Ct value (E gene)",
       main = "COVID-19 infection",
       cex.main = 1.5,
       col = adjustcolor(ifelse(df_round$estbinres, yes = "red", no = "darkgrey"),
                         alpha.f = 0.5
       )
  )
  ids <- which(df_round$res == "Not Detected")
  print(table(df_round$estbinres[ids]))
  points(df_round$ct1[ids], df_round$ct2[ids],
         pch = 19, cex = 0.8,
         col = ifelse(df_round$estbinres[ids], yes = "red", no = "darkgrey")
  )
  points(df_round$ct1[ids], df_round$ct2[ids])
  legend("top",
         pch = c(19, 19, 21), pt.cex = c(0.8, 0.8, 1),
         col = c("grey", "red", "black"), pt.bg = c(NA, NA, "grey"),
         legend = c("Swab negative", "Swab positive", "Swab negative / undetected")
  )
  dev.off()
}


if (as.numeric(round_id) > 15) {{
  pdf(paste0("../Figures/Ct_values_infection_", Sys.Date(), ".pdf"), width = 15, height = 6)
  par(mfrow = c(1, 3), mar = c(5, 5, 5, 1))
  plot(df_round$ct1, df_round$ct2,
       pch = 19, cex = 0.8, las = 1, cex.lab = 1.5,
       xlab = "Ct value (N gene)",
       ylab = "Ct value (E gene)",
       main = "COVID-19 infection",
       cex.main = 1.5,
       col = adjustcolor(ifelse(df_round$estbinres, yes = "red", no = "darkgrey"),
                         alpha.f = 0.5
       )
  )
  ids <- which(df_round$res == "Not Detected")
  points(df_round$ct1[ids], df_round$ct2[ids])
  
  plot(df_round$influenzaacpvalue, df_round$influenzabcpvalue,
       pch = 19, cex = 0.8, las = 1, cex.lab = 1.5,
       xlab = "Cp value (influenza A)",
       ylab = "Cp value (influenza B)",
       main = "Influenza A infection",
       cex.main = 1.5,
       col = adjustcolor(ifelse(df_round$influenzaa, yes = "red", no = "darkgrey"),
                         alpha.f = 0.5
       )
  )
  
  plot(df_round$influenzaacpvalue, df_round$influenzabcpvalue,
       pch = 19, cex = 0.8, las = 1, cex.lab = 1.5,
       xlab = "Cp value (influenza A)",
       ylab = "Cp value (influenza B)",
       main = "Influenza B infection",
       cex.main = 1.5,
       col = adjustcolor(ifelse(df_round$influenzab, yes = "red", no = "darkgrey"),
                         alpha.f = 0.5
       )
  )
  dev.off()
}
  
  {
    pdf(paste0("../Figures/Ct_values_infection_influenza_A_B_sep_", Sys.Date(), ".pdf"), width = 12, height = 7)
    par(mfrow = c(1, 2), mar = c(5, 5, 7, 1))
    for (ct_outcome in c("influenzaacpvalue", "influenzabcpvalue")) {
      mylist <- list()
      
      # Negative to flu
      ids <- which(df_round[, gsub("cpvalue", "", ct_outcome)] == 0)
      mylist <- c(mylist, list(df_round[ids, ct_outcome]))
      
      # Positive to flu
      ids <- which(df_round[, gsub("cpvalue", "", ct_outcome)] == 1)
      mylist <- c(mylist, list(df_round[ids, ct_outcome]))
      
      N0 <- formatC(sapply(mylist, FUN = function(x) {
        sum(round(x, digits = 4) == 0)
      }),
      format = "f", digits = 0, big.mark = ","
      )
      mylist <- lapply(mylist, FUN = function(x) {
        x[which(round(x, digits = 4) != 0)]
      })
      
      names(mylist) <- c("Negative", "Positive")
      names(mylist) <- paste0(
        names(mylist), "\n (N=",
        formatC(sapply(mylist, length), format = "f", digits = 0, big.mark = ","), ")"
      )
      
      mycolours <- lighten(c("grey30", "darkred"), amount = 0.5)
      
      boxplot(mylist,
              range = 0,
              boxcol = "white", col = mycolours,
              staplecol = mycolours, whiskcol = mycolours, lty = 1,
              las = 1, cex.axis = 1.5, cex.main = 2, xaxt = "n", ylim = c(0, 50),
              main = paste0(
                "Influenza ",
                ifelse(ct_outcome == "influenzaacpvalue", yes = "A", no = "B"), " infection"
              ),
              ylab = paste0(
                "Cp values (influenza ",
                ifelse(ct_outcome == "influenzaacpvalue", yes = "A", no = "B"), ")"
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
        line = 2, cex.axis = 1.5, tick = FALSE
      )
      
      axis(
        side = 3, at = 1:length(mylist),
        # labels = paste0("N0=", N0),
        labels=unlist(sapply(N0, FUN=function(x){eval(parse(text=paste0("expression(N['Cp=0']*'=",x, "')")))})),
        line = -0.5, cex.axis = 1.5, tick = FALSE
      )
      
    }
    dev.off()
  }
  
  {
    pdf(paste0("../Figures/Ct_values_infection_influenza_A_B_", Sys.Date(), ".pdf"), width = 15, height = 7)
    mycolours <- "tomato"
    par(mfrow = c(1, 2), mar = c(5, 5, 5, 1))
    for (ct_outcome in c("influenzaacpvalue", "influenzabcpvalue")) {
      mylist <- list()
      
      # Flu A-/B-
      ids <- which((df_round[, "influenzaa"] == 0) & (df_round[, "influenzab"] == 0))
      mylist <- c(mylist, list(df_round[ids, ct_outcome]))
      
      # Flu A+/B-
      ids <- which((df_round[, "influenzaa"] == 1) & (df_round[, "influenzab"] == 0))
      mylist <- c(mylist, list(df_round[ids, ct_outcome]))
      
      # Flu A-/B+
      ids <- which((df_round[, "influenzaa"] == 0) & (df_round[, "influenzab"] == 1))
      mylist <- c(mylist, list(df_round[ids, ct_outcome]))
      
      # Flu A and B
      ids <- which((df_round[, "influenzaa"] == 1) & (df_round[, "influenzab"] == 1))
      mylist <- c(mylist, list(df_round[ids, ct_outcome]))
      
      N0 <- formatC(sapply(mylist, FUN = function(x) {
        sum(round(x, digits = 4) == 0)
      }),
      format = "f", digits = 0, big.mark = ","
      )
      mylist <- lapply(mylist, FUN = function(x) {
        x[which(round(x, digits = 4) != 0)]
      })
      
      names(mylist) <- c("A-/B-", "A+/B-", "A-/B+", "A+/B+")
      names(mylist) <- paste0(
        names(mylist), "\n (N=",
        formatC(sapply(mylist, length), format = "f", digits = 0, big.mark = ","), ")"
      )
      
      boxplot(mylist,
              range = 0,
              boxcol = "white", col = mycolours,
              staplecol = mycolours, whiskcol = mycolours, lty = 1,
              las = 1, cex.axis = 1.5, cex.main = 2, xaxt = "n", ylim = c(0, 50),
              ylab = paste0(
                "Cp values (influenza ",
                ifelse(ct_outcome == "influenzaacpvalue", yes = "A", no = "B"), ")"
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
        line = 2, cex.axis = 1.5, tick = FALSE
      )
      
      axis(
        side = 3, at = 1:length(mylist),
        # labels = paste0("N0=", N0),
        labels=unlist(sapply(N0, FUN=function(x){eval(parse(text=paste0("expression(N['Cp=0']*'=",x, "')")))})),
        line = -0.5, cex.axis = 1.5, tick = FALSE
      )
    }
    dev.off()
  }
  
  
  {
    pdf(paste0("../Figures/Ct_values_infections_covid_flu_1_", Sys.Date(), ".pdf"), width = 18, height = 7)
    par(mfrow = c(1, 2), mar = c(5, 5, 5, 1))
    mycolours <- "tomato"
    for (ct_outcome in c("ct1", "ct2")) {
      mylist <- list()
      
      # Negative to COVID-19 and flu
      ids <- which((df_round$estbinres == 0) & (df_round$influenzaa == 0) & (df_round$influenzab == 0))
      mylist <- c(mylist, list(df_round[ids, ct_outcome]))
      
      # Positive to COVID-19 and negative to flu
      ids <- which((df_round$estbinres == 1) & (df_round$influenzaa == 0) & (df_round$influenzab == 0))
      mylist <- c(mylist, list(df_round[ids, ct_outcome]))
      
      # Positive to COVID-19 and positive to flu
      ids <- which((df_round$estbinres == 1) & (df_round$influenzaa == 1) & (df_round$influenzab == 1))
      mylist <- c(mylist, list(df_round[ids, ct_outcome]))
      
      # Negative to COVID-19 and positive to flu
      ids <- which((df_round$estbinres == 0) & (df_round$influenzaa == 1) & (df_round$influenzab == 1))
      mylist <- c(mylist, list(df_round[ids, ct_outcome]))
      
      N0 <- formatC(sapply(mylist, FUN = function(x) {
        sum(round(x, digits = 4) == 0)
      }),
      format = "f", digits = 0, big.mark = ","
      )
      mylist <- lapply(mylist, FUN = function(x) {
        x[which(round(x, digits = 4) != 0)]
      })
      
      sapply(mylist, length)
      names(mylist) <- c("C-/F-", "C+/F-", "C+/F+", "C-/F+")
      names(mylist) <- paste0(
        names(mylist), "\n (N=",
        formatC(sapply(mylist, length), format = "f", digits = 0, big.mark = ","), ")"
      )
      
      boxplot(mylist,
              range = 0,
              boxcol = "white", col = mycolours,
              staplecol = mycolours, whiskcol = mycolours, lty = 1,
              las = 1, cex.axis = 1.5, cex.main = 2,
              xaxt = "n", ylim = c(0, 50),
              ylab = ifelse(ct_outcome %in% c("ct1", "ct2"),
                            yes = ifelse(ct_outcome == "ct1", yes = "Ct values (N gene)", no = "Ct values (E gene)"),
                            no = ifelse(ct_outcome == "influenzaacpvalue", yes = "Cp values (influenza A)", no = "Cp values (influenza B)")
              ),
              cex.lab = 2, axis.line = 2
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
        line = 2, cex.axis = 1.5, tick = FALSE
      )
      
      axis(
        side = 3, at = 1:length(mylist),
        # labels = paste0("N0=", N0),
        labels=unlist(sapply(N0, FUN=function(x){eval(parse(text=paste0("expression(N['Cp=0']*'=",x, "')")))})),
        line = -0.5, cex.axis = 1.5, tick = FALSE
      )
    }
    dev.off()
  }
  
  
  {
    pdf(paste0("../Figures/Ct_values_infections_covid_flu_2_", Sys.Date(), ".pdf"), width = 18, height = 7)
    par(mfrow = c(1, 2), mar = c(5, 5, 5, 1))
    mycolours <- "tomato"
    for (ct_outcome in c("influenzaacpvalue", "influenzabcpvalue")) {
      mylist <- list()
      
      # Negative to COVID-19 and flu
      ids <- which((df_round$estbinres == 0) & (df_round$influenzaa == 0) & (df_round$influenzab == 0))
      mylist <- c(mylist, list(df_round[ids, ct_outcome]))
      
      # Positive to COVID-19 and negative to flu
      ids <- which((df_round$estbinres == 1) & (df_round$influenzaa == 0) & (df_round$influenzab == 0))
      mylist <- c(mylist, list(df_round[ids, ct_outcome]))
      
      # Positive to COVID-19 and positive to flu
      ids <- which((df_round$estbinres == 1) & (df_round$influenzaa == 1) & (df_round$influenzab == 1))
      mylist <- c(mylist, list(df_round[ids, ct_outcome]))
      
      # Negative to COVID-19 and positive to flu
      ids <- which((df_round$estbinres == 0) & (df_round$influenzaa == 1) & (df_round$influenzab == 1))
      mylist <- c(mylist, list(df_round[ids, ct_outcome]))
      
      N0 <- formatC(sapply(mylist, FUN = function(x) {
        sum(round(x, digits = 4) == 0)
      }),
      format = "f", digits = 0, big.mark = ","
      )
      mylist <- lapply(mylist, FUN = function(x) {
        x[which(round(x, digits = 4) != 0)]
      })
      
      sapply(mylist, length)
      names(mylist) <- c("C-/F-", "C+/F-", "C+/F+", "C-/F+")
      names(mylist) <- paste0(
        names(mylist), "\n (N=",
        formatC(sapply(mylist, length), format = "f", digits = 0, big.mark = ","), ")"
      )
      
      boxplot(mylist,
              range = 0,
              boxcol = "white", col = mycolours,
              staplecol = mycolours, whiskcol = mycolours, lty = 1,
              las = 1, cex.axis = 1.5, cex.main = 2,
              xaxt = "n", ylim = c(0, 50),
              ylab = ifelse(ct_outcome %in% c("ct1", "ct2"),
                            yes = ifelse(ct_outcome == "ct1", yes = "Ct values (N gene)", no = "Ct values (E gene)"),
                            no = ifelse(ct_outcome == "influenzaacpvalue", yes = "Cp values (influenza A)", no = "Cp values (influenza B)")
              ),
              cex.lab = 2, axis.line = 2
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
        line = 2, cex.axis = 1.5, tick = FALSE
      )
      
      axis(
        side = 3, at = 1:length(mylist),
        # labels = paste0("N0=", N0),
        labels=unlist(sapply(N0, FUN=function(x){eval(parse(text=paste0("expression(N['Cp=0']*'=",x, "')")))})),
        line = -0.5, cex.axis = 1.5, tick = FALSE
      )
    }
    dev.off()
  }
  
  
  {
    pdf(paste0("../Figures/Ct_values_infections_covid_flu_A_B_1_", Sys.Date(), ".pdf"), width = 20, height = 7)
    par(mfrow = c(1, 2), mar = c(5, 5, 5, 1))
    mycolours <- "tomato"
    for (ct_outcome in c("ct1", "ct2")) {
      mylist <- list()
      
      # C-/A-/B-
      ids <- which((df_round$estbinres == 0) & (df_round$influenzaa == 0) & (df_round$influenzab == 0))
      mylist <- c(mylist, list(df_round[ids, ct_outcome]))
      
      # C+/A-/B-
      ids <- which((df_round$estbinres == 1) & (df_round$influenzaa == 0) & (df_round$influenzab == 0))
      mylist <- c(mylist, list(df_round[ids, ct_outcome]))
      
      # C+/A+/B-
      ids <- which((df_round$estbinres == 1) & (df_round$influenzaa == 1) & (df_round$influenzab == 0))
      mylist <- c(mylist, list(df_round[ids, ct_outcome]))
      
      # C+/A-/B+
      ids <- which((df_round$estbinres == 1) & (df_round$influenzaa == 0) & (df_round$influenzab == 1))
      mylist <- c(mylist, list(df_round[ids, ct_outcome]))
      
      # C+/A+/B+
      ids <- which((df_round$estbinres == 1) & (df_round$influenzaa == 1) & (df_round$influenzab == 1))
      mylist <- c(mylist, list(df_round[ids, ct_outcome]))
      
      # C-/A+/B-
      ids <- which((df_round$estbinres == 0) & (df_round$influenzaa == 1) & (df_round$influenzab == 0))
      mylist <- c(mylist, list(df_round[ids, ct_outcome]))
      
      # C-/A+/B+
      ids <- which((df_round$estbinres == 0) & (df_round$influenzaa == 1) & (df_round$influenzab == 1))
      mylist <- c(mylist, list(df_round[ids, ct_outcome]))
      
      # C-/A-/B+
      ids <- which((df_round$estbinres == 0) & (df_round$influenzaa == 0) & (df_round$influenzab == 1))
      mylist <- c(mylist, list(df_round[ids, ct_outcome]))
      
      N0 <- formatC(sapply(mylist, FUN = function(x) {
        sum(round(x, digits = 4) == 0)
      }),
      format = "f", digits = 0, big.mark = ","
      )
      mylist <- lapply(mylist, FUN = function(x) {
        x[which(round(x, digits = 4) != 0)]
      })
      
      sapply(mylist, length)
      names(mylist) <- c(
        "C-/A-/B-", "C+/A-/B-", "C+/A+/B-", "C+/A-/B+",
        "C+/A+/B+",
        "C-/A+/B-", "C-/A+/B+", "C-/A-/B+"
      )
      names(mylist) <- paste0(
        names(mylist), "\n (N=",
        formatC(sapply(mylist, length), format = "f", digits = 0, big.mark = ","), ")"
      )
      
      boxplot(mylist,
              range = 0,
              boxcol = "white", col = mycolours,
              staplecol = mycolours, whiskcol = mycolours, lty = 1,
              las = 1, cex.axis = 1.5, cex.main = 2, ylim = c(0, 50),
              xaxt = "n",
              ylab = ifelse(ct_outcome %in% c("ct1", "ct2"),
                            yes = ifelse(ct_outcome == "ct1", yes = "Ct values (N gene)", no = "Ct values (E gene)"),
                            no = ifelse(ct_outcome == "influenzaacpvalue", yes = "Cp values (influenza A)", no = "Cp values (influenza B)")
              ),
              cex.lab = 2, axis.line = 2
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
        line = 2, cex.axis = 1, tick = FALSE
      )
      
      axis(
        side = 3, at = 1:length(mylist),
        # labels = paste0("N0=", N0),
        labels=unlist(sapply(N0, FUN=function(x){eval(parse(text=paste0("expression(N['Cp=0']*'=",x, "')")))})),
        line = -0.5, cex.axis = 1, tick = FALSE
      )
    }
    dev.off()
  }
  
  
  {
    pdf(paste0("../Figures/Ct_values_infections_covid_flu_A_B_2_", Sys.Date(), ".pdf"), width = 20, height = 7)
    par(mfrow = c(1, 2), mar = c(5, 5, 5, 1))
    mycolours <- "tomato"
    for (ct_outcome in c("influenzaacpvalue", "influenzabcpvalue")) {
      mylist <- list()
      
      # C-/A-/B-
      ids <- which((df_round$estbinres == 0) & (df_round$influenzaa == 0) & (df_round$influenzab == 0))
      mylist <- c(mylist, list(df_round[ids, ct_outcome]))
      
      # C+/A-/B-
      ids <- which((df_round$estbinres == 1) & (df_round$influenzaa == 0) & (df_round$influenzab == 0))
      mylist <- c(mylist, list(df_round[ids, ct_outcome]))
      
      # C+/A+/B-
      ids <- which((df_round$estbinres == 1) & (df_round$influenzaa == 1) & (df_round$influenzab == 0))
      mylist <- c(mylist, list(df_round[ids, ct_outcome]))
      
      # C+/A-/B+
      ids <- which((df_round$estbinres == 1) & (df_round$influenzaa == 0) & (df_round$influenzab == 1))
      mylist <- c(mylist, list(df_round[ids, ct_outcome]))
      
      # C+/A+/B+
      ids <- which((df_round$estbinres == 1) & (df_round$influenzaa == 1) & (df_round$influenzab == 1))
      mylist <- c(mylist, list(df_round[ids, ct_outcome]))
      
      # C-/A+/B-
      ids <- which((df_round$estbinres == 0) & (df_round$influenzaa == 1) & (df_round$influenzab == 0))
      mylist <- c(mylist, list(df_round[ids, ct_outcome]))
      
      # C-/A+/B+
      ids <- which((df_round$estbinres == 0) & (df_round$influenzaa == 1) & (df_round$influenzab == 1))
      mylist <- c(mylist, list(df_round[ids, ct_outcome]))
      
      # C-/A-/B+
      ids <- which((df_round$estbinres == 0) & (df_round$influenzaa == 0) & (df_round$influenzab == 1))
      mylist <- c(mylist, list(df_round[ids, ct_outcome]))
      
      N0 <- formatC(sapply(mylist, FUN = function(x) {
        sum(round(x, digits = 4) == 0)
      }),
      format = "f", digits = 0, big.mark = ","
      )
      mylist <- lapply(mylist, FUN = function(x) {
        x[which(round(x, digits = 4) != 0)]
      })
      
      sapply(mylist, length)
      names(mylist) <- c(
        "C-/A-/B-", "C+/A-/B-", "C+/A+/B-", "C+/A-/B+",
        "C+/A+/B+",
        "C-/A+/B-", "C-/A+/B+", "C-/A-/B+"
      )
      names(mylist) <- paste0(
        names(mylist), "\n (N=",
        formatC(sapply(mylist, length), format = "f", digits = 0, big.mark = ","), ")"
      )
      
      boxplot(mylist,
              range = 0,
              boxcol = "white", col = mycolours,
              staplecol = mycolours, whiskcol = mycolours, lty = 1,
              las = 1, cex.axis = 1.5, cex.main = 2, ylim = c(0, 50),
              xaxt = "n",
              ylab = ifelse(ct_outcome %in% c("ct1", "ct2"),
                            yes = ifelse(ct_outcome == "ct1", yes = "Ct values (N gene)", no = "Ct values (E gene)"),
                            no = ifelse(ct_outcome == "influenzaacpvalue", yes = "Cp values (influenza A)", no = "Cp values (influenza B)")
              ),
              cex.lab = 2, axis.line = 2
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
        line = 2, cex.axis = 1, tick = FALSE
      )
      
      axis(
        side = 3, at = 1:length(mylist),
        # labels = paste0("N0=", N0),
        labels=unlist(sapply(N0, FUN=function(x){eval(parse(text=paste0("expression(N['Cp=0']*'=",x, "')")))})),
        line = -0.5, cex.axis = 1, tick = FALSE
      )
    }
    dev.off()
  }
  
  
  outcome="influenza"
  {
    pdf(paste0("../Figures/Ct_values_infection_", outcome, "_", Sys.Date(), ".pdf"), width = 12, height = 7)
    par(mfrow = c(1, 2), mar = c(5, 5, 3, 1))
    for (ct_outcome in c("influenzaacpvalue", "influenzabcpvalue")) {
      mylist <- list()
      
      # Negative to flu
      ids <- which(df_round[, outcome] == 0)
      mylist <- c(mylist, list(df_round[ids, ct_outcome]))
      
      # Negative to flu
      ids <- which(df_round[, outcome] == 1)
      mylist <- c(mylist, list(df_round[ids, ct_outcome]))
      
      N0 <- formatC(sapply(mylist, FUN = function(x) {
        sum(round(x, digits = 4) == 0)
      }),
      format = "f", digits = 0, big.mark = ","
      )
      mylist <- lapply(mylist, FUN = function(x) {
        x[which(round(x, digits = 4) != 0)]
      })
      
      names(mylist) <- c("Negative", "Positive")
      names(mylist) <- paste0(
        names(mylist), "\n (N=",
        formatC(sapply(mylist, length), format = "f", digits = 0, big.mark = ","), ")"
      )
      
      mycolours <- lighten(c("grey30", "darkred"), amount = 0.5)
      
      boxplot(mylist,
              range = 0,
              boxcol = "white", col = mycolours,
              staplecol = mycolours, whiskcol = mycolours, lty = 1,
              las = 1, cex.axis = 1.5, cex.main = 2, xaxt = "n", ylim = c(0, 50),
              ylab = paste0(
                "Cp values (influenza ",
                ifelse(ct_outcome == "influenzaacpvalue", yes = "A", no = "B"), ")"
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
        line = 2, cex.axis = 1.5, tick = FALSE
      )
      
      axis(
        side = 3, at = 1:length(mylist),
        # labels = paste0("N0=", N0),
        labels=unlist(sapply(N0, FUN=function(x){eval(parse(text=paste0("expression(N['Cp=0']*'=",x, "')")))})),
        line = -0.5, cex.axis = 1.5, tick = FALSE
      )
    }
    dev.off()
  }
  
  {
    pdf(paste0("../Figures/Ct_values_infection_", outcome, "_by_age_", Sys.Date(), ".pdf"), width = 15, height = 7)
    par(mfrow = c(1, 2), mar = c(5, 5, 3, 1))
    for (ct_outcome in c("influenzaacpvalue", "influenzabcpvalue")) {
      mylist <- list()
      
      # Negative to flu and <18
      ids <- which((df_round[, outcome] == 0)&(df_round$age%in%c("05-12", "13-17")))
      mylist <- c(mylist, list(df_round[ids, ct_outcome]))
      
      # Positive to flu and <18
      ids <- which((df_round[, outcome] == 1)&(df_round$age%in%c("05-12", "13-17")))
      mylist <- c(mylist, list(df_round[ids, ct_outcome]))
      
      # Negative to flu and 18+
      ids <- which((df_round[, outcome] == 0)&(!df_round$age%in%c("05-12", "13-17")))
      mylist <- c(mylist, list(df_round[ids, ct_outcome]))
      
      # Positive to flu and 18+
      ids <- which((df_round[, outcome] == 1)&(!df_round$age%in%c("05-12", "13-17")))
      mylist <- c(mylist, list(df_round[ids, ct_outcome]))
      
      N0 <- formatC(sapply(mylist, FUN = function(x) {
        sum(round(x, digits = 4) == 0)
      }),
      format = "f", digits = 0, big.mark = ","
      )
      mylist <- lapply(mylist, FUN = function(x) {
        x[which(round(x, digits = 4) != 0)]
      })
      
      names(mylist) <- c("F- (<18)", "F+ (<18)", "F- (18+)", "F+ (18+)")
      names(mylist) <- paste0(
        names(mylist), "\n (N=",
        formatC(sapply(mylist, length), format = "f", digits = 0, big.mark = ","), ")"
      )
      
      mycolours <- lighten(c("grey30", "darkred"), amount = 0.5)
      
      boxplot(mylist,
              range = 0,
              boxcol = "white", col = mycolours,
              staplecol = mycolours, whiskcol = mycolours, lty = 1,
              las = 1, cex.axis = 1.5, cex.main = 2, xaxt = "n", ylim = c(0, 50),
              ylab = paste0(
                "Cp values (influenza ",
                ifelse(ct_outcome == "influenzaacpvalue", yes = "A", no = "B"), ")"
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
        line = 2, cex.axis = 1.5, tick = FALSE
      )
      
      axis(
        side = 3, at = 1:length(mylist),
        labels=unlist(sapply(N0, FUN=function(x){eval(parse(text=paste0("expression(N['Cp=0']*'=",x, "')")))})),
        line = -0.5, cex.axis = 1.5, tick = FALSE
      )
    }
    dev.off()
  }
  
  {
    pdf(paste0("../Figures/Ct_values_infection_", outcome, "_by_symptom_", Sys.Date(), ".pdf"), width = 20, height = 10)
    par(mfrow = c(2, 1), mar = c(5, 5, 3, 1))
    for (ct_outcome in c("influenzaacpvalue", "influenzabcpvalue")) {
      mylist <- list()
      
      # Negative to flu 
      for (id in 1:length(levels(df_round$sympt_cat))){
        ids <- which((df_round[, outcome] == 0)&(df_round$sympt_cat==levels(df_round$sympt_cat)[id]))
        mylist <- c(mylist, list(df_round[ids, ct_outcome]))
      }
      
      # Positive to flu 
      for (id in 1:length(levels(df_round$sympt_cat))){
        ids <- which((df_round[, outcome] == 1)&(df_round$sympt_cat==levels(df_round$sympt_cat)[id]))
        mylist <- c(mylist, list(df_round[ids, ct_outcome]))
      }
      
      N0 <- formatC(sapply(mylist, FUN = function(x) {
        sum(round(x, digits = 4) == 0)
      }),
      format = "f", digits = 0, big.mark = ","
      )
      mylist <- lapply(mylist, FUN = function(x) {
        x[which(round(x, digits = 4) != 0)]
      })
      
      names(mylist)=c(paste0("F- (", levels(df_round$sympt_cat), ")"), 
                      paste0("F+ (", levels(df_round$sympt_cat), ")"))
      names(mylist) <- paste0(
        names(mylist), "\n (N=",
        formatC(sapply(mylist, length), format = "f", digits = 0, big.mark = ","), ")"
      )
      
      mycolours <- lighten(c("grey30", "darkred"), amount = 0.5)
      
      boxplot(mylist,
              range = 0,
              boxcol = "white", col = mycolours,
              staplecol = mycolours, whiskcol = mycolours, lty = 1,
              las = 1, cex.axis = 1.5, cex.main = 2, xaxt = "n", ylim = c(0, 50),
              ylab = paste0(
                "Cp values (influenza ",
                ifelse(ct_outcome == "influenzaacpvalue", yes = "A", no = "B"), ")"
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
        line = 2, cex.axis = 1, tick = FALSE
      )
      
      axis(
        side = 3, at = 1:length(mylist),
        labels=unlist(sapply(N0, FUN=function(x){eval(parse(text=paste0("expression(N['Cp=0']*'=",x, "')")))})),
        line = -0.5, cex.axis = 1.5, tick = FALSE
      )
    }
    dev.off()
  }
  
  # df_round_full=df_round
  # 
  # for (age_class in c("all", "below_18", "above_18")){
  #   if (age_class=="all"){
  #     df_round=df_round_full
  #   } else {
  #     if (age_class=="below_18"){
  #       df_round=df_round_full[which(df_round_full$age%in%c("05-12", "13-17")),]
  #     } else {
  #       df_round=df_round_full[which(!df_round_full$age%in%c("05-12", "13-17")),]
  #     }
  #   }
  #   
  #   {
  #     pdf(paste0("../Figures/Ct_values_infection_", outcome, "_by_symptom_", age_class, "_", Sys.Date(), ".pdf"), width = 20, height = 10)
  #     par(mfrow = c(2, 1), mar = c(5, 5, 3, 1))
  #     for (ct_outcome in c("influenzaacpvalue", "influenzabcpvalue")) {
  #       mylist <- list()
  #       
  #       
  #       # Negative to flu 
  #       for (id in 1:length(levels(df_round$sympt_cat))){
  #         ids <- which((df_round[, outcome] == 0)&(df_round$sympt_cat==levels(df_round$sympt_cat)[id]))
  #         mylist <- c(mylist, list(df_round[ids, ct_outcome]))
  #       }
  #       
  #       # Positive to flu 
  #       for (id in 1:length(levels(df_round$sympt_cat))){
  #         ids <- which((df_round[, outcome] == 1)&(df_round$sympt_cat==levels(df_round$sympt_cat)[id]))
  #         mylist <- c(mylist, list(df_round[ids, ct_outcome]))
  #       }
  #       
  #       N0 <- formatC(sapply(mylist, FUN = function(x) {
  #         sum(round(x, digits = 4) == 0)
  #       }),
  #       format = "f", digits = 0, big.mark = ","
  #       )
  #       mylist <- lapply(mylist, FUN = function(x) {
  #         x[which(round(x, digits = 4) != 0)]
  #       })
  #       
  #       names(mylist)=c(paste0("F- (", levels(df_round$sympt_cat), ")"), 
  #                       paste0("F+ (", levels(df_round$sympt_cat), ")"))
  #       names(mylist) <- paste0(
  #         names(mylist), "\n (N=",
  #         formatC(sapply(mylist, length), format = "f", digits = 0, big.mark = ","), ")"
  #       )
  #       
  #       mycolours <- lighten(c("grey30", "darkred"), amount = 0.5)
  #       
  #       boxplot(mylist,
  #               range = 0,
  #               boxcol = "white", col = mycolours,
  #               staplecol = mycolours, whiskcol = mycolours, lty = 1,
  #               las = 1, cex.axis = 1.5, cex.main = 2, xaxt = "n", ylim = c(0, 50),
  #               ylab = paste0(
  #                 "Cp values (influenza ",
  #                 ifelse(ct_outcome == "influenzaacpvalue", yes = "A", no = "B"), ")"
  #               ),
  #               cex.lab = 2
  #       )
  #       
  #       set.seed(1)
  #       for (i in 1:length(mylist)) {
  #         points(i + ProportionalJitter(mylist[[i]]), mylist[[i]],
  #                pch = 19, cex = 0.5,
  #                col = "darkred"
  #         )
  #       }
  #       
  #       axis(side = 1, at = 1:length(mylist), labels = NA)
  #       axis(
  #         side = 1, at = 1:length(mylist), labels = names(mylist),
  #         line = 2, cex.axis = 1, tick = FALSE
  #       )
  #       
  #       axis(
  #         side = 3, at = 1:length(mylist),
  #         labels=unlist(sapply(N0, FUN=function(x){eval(parse(text=paste0("expression(N['Cp=0']*'=",x, "')")))})),
  #         line = -0.5, cex.axis = 1.5, tick = FALSE
  #       )
  #     }
  #     dev.off()
  #   }
  # }
  
  
  df_round_full=df_round
  
  for (strat_var in c("fluvacc", "any_sympt", "vacc_and_sympt")){
    
    for (age_class in c("all", "below_18", "above_18")){
      if (age_class=="all"){
        df_round=df_round_full
      } else {
        if (age_class=="below_18"){
          df_round=df_round_full[which(df_round_full$age%in%c("05-11", "12-17")),]
        } else {
          df_round=df_round_full[which(!df_round_full$age%in%c("05-11", "12-17")),]
        }
      }
      
      {
        pdf(paste0("../Figures/Ct_values_infection_", outcome, "_by_",strat_var,"_", age_class, "_", Sys.Date(), ".pdf"), width = 20, height = 10)
        par(mfrow = c(2, 1), mar = c(5, 5, 3, 1))
        for (ct_outcome in c("influenzaacpvalue", "influenzabcpvalue")) {
          mylist <- list()
          
          for (id in 1:length(levels(df_round[,strat_var]))){
            # Negative to flu 
            ids <- which((df_round[, outcome] == 0)&(df_round[,strat_var]==levels(df_round[,strat_var])[id]))
            mylist <- c(mylist, list(df_round[ids, ct_outcome]))
            
            # Positive to flu 
            ids <- which((df_round[, outcome] == 1)&(df_round[,strat_var]==levels(df_round[,strat_var])[id]))
            mylist <- c(mylist, list(df_round[ids, ct_outcome]))
          }
          
          N0 <- formatC(sapply(mylist, FUN = function(x) {
            sum(round(x, digits = 4) == 0)
          }),
          format = "f", digits = 0, big.mark = ","
          )
          mylist <- lapply(mylist, FUN = function(x) {
            x[which(round(x, digits = 4) != 0)]
          })
          
          names(mylist)=paste0(c("F-/","F+/"), rep(levels(df_round[,strat_var]),each=2))
          names(mylist) <- paste0(
            names(mylist), "\n (N=",
            formatC(sapply(mylist, length), format = "f", digits = 0, big.mark = ","), ")"
          )
          
          mycolours <- lighten(c("grey30", "darkred"), amount = 0.5)
          
          boxplot(mylist,
                  range = 0,
                  boxcol = "white", col = mycolours,
                  staplecol = mycolours, whiskcol = mycolours, lty = 1,
                  las = 1, cex.axis = 1.5, cex.main = 2, xaxt = "n", ylim = c(0, 50),
                  ylab = paste0(
                    "Cp values (influenza ",
                    ifelse(ct_outcome == "influenzaacpvalue", yes = "A", no = "B"), ")"
                  ),
                  cex.lab = 2
          )
          
          set.seed(1)
          for (i in 1:length(mylist)) {
            if (length(mylist[[i]])>0){
              points(i + ProportionalJitter(mylist[[i]]), mylist[[i]],
                     pch = 19, cex = 0.5,
                     col = "darkred"
              )
            }
          }
          
          axis(side = 1, at = 1:length(mylist), labels = NA)
          axis(
            side = 1, at = 1:length(mylist), labels = names(mylist),
            line = 2, cex.axis = 1, tick = FALSE
          )
          
          axis(
            side = 3, at = 1:length(mylist),
            labels=unlist(sapply(N0, FUN=function(x){eval(parse(text=paste0("expression(N['Cp=0']*'=",x, "')")))})),
            line = -0.5, cex.axis = 1, tick = FALSE
          )
        }
        dev.off()
      }
    }
  }
}

# for (sheet in 1:2){
#   {
#     pdf(paste0("../Figures/Ct_values_infection_influenza_A_B_sep_sheet_",sheet,"_", Sys.Date(), ".pdf"), width = 12, height = 7)
#     par(mfrow = c(1, 2), mar = c(5, 5, 7, 1))
#     for (ct_outcome in c("influenzaacpvalue", "influenzabcpvalue")) {
#       mylist <- list()
#       
#       df_round=df_round_full
# 
# 
#       mydata=read.xlsx("../Data/FLU Confirmation Results 2021 01 03.xlsx", sheet = sheet)
#       myu_passcode=mydata$Sample.Name
#       myu_passcode=gsub("UK", "", gsub(" \\(.*", "", myu_passcode))
#       df_round=df_round[myu_passcode,]
#       
#       # Negative to flu
#       ids <- which(df_round[, gsub("cpvalue", "", ct_outcome)] == 0)
#       mylist <- c(mylist, list(df_round[ids, ct_outcome]))
#       
#       # Positive to flu
#       ids <- which(df_round[, gsub("cpvalue", "", ct_outcome)] == 1)
#       mylist <- c(mylist, list(df_round[ids, ct_outcome]))
#       
#       N0 <- formatC(sapply(mylist, FUN = function(x) {
#         sum(round(x, digits = 4) == 0)
#       }),
#       format = "f", digits = 0, big.mark = ","
#       )
#       mylist <- lapply(mylist, FUN = function(x) {
#         x[which(round(x, digits = 4) != 0)]
#       })
#       
#       names(mylist) <- c("Negative", "Positive")
#       names(mylist) <- paste0(
#         names(mylist), "\n (N=",
#         formatC(sapply(mylist, length), format = "f", digits = 0, big.mark = ","), ")"
#       )
#       
#       mycolours <- lighten(c("grey30", "darkred"), amount = 0.5)
#       
#       plot(NULL, xlim=c(0.5,2.5), ylim=c(0, 50), xlab="", ylab="", xaxt="n", yaxt="n",
#            panel.first=c(abline(h=seq(0,50), lty=3, col="grey"),
#                          abline(h=seq(0,50,by=5), lty=1, col="grey")))
#       boxplot(mylist,
#               range = 0,
#               boxcol = "white", col = mycolours,
#               staplecol = mycolours, whiskcol = mycolours, lty = 1,
#               las = 1, cex.axis = 1.5, cex.main = 2, xaxt = "n", ylim = c(0, 50),
#               main = paste0(
#                 "Influenza ",
#                 ifelse(ct_outcome == "influenzaacpvalue", yes = "A", no = "B"), " infection"
#               ),
#               ylab = paste0(
#                 "Cp values (influenza ",
#                 ifelse(ct_outcome == "influenzaacpvalue", yes = "A", no = "B"), ")"
#               ),
#               cex.lab = 2, add = TRUE
#       )
#       
#       set.seed(1)
#       for (i in 1:length(mylist)) {
#         if (length(mylist[[i]]) > 0) {
#           points(i + ProportionalJitter(mylist[[i]]), mylist[[i]],
#                  pch = 19, cex = 0.5,
#                  col = "darkred"
#           )
#         }
#       }
#       
#       axis(side = 1, at = 1:length(mylist), labels = NA)
#       axis(
#         side = 1, at = 1:length(mylist), labels = names(mylist),
#         line = 2, cex.axis = 1.5, tick = FALSE
#       )
#       
#       axis(
#         side = 3, at = 1:length(mylist),
#         # labels = paste0("N0=", N0),
#         labels=unlist(sapply(N0, FUN=function(x){eval(parse(text=paste0("expression(N['Cp=0']*'=",x, "')")))})),
#         line = -0.5, cex.axis = 1.5, tick = FALSE
#       )
#       
#     }
#     dev.off()
#   }
# }


# Copying output to transfer folder
myfiles=list.files("../Figures", pattern = as.character(Sys.Date()))
if (direct_export) {
  for (i in 1:length(myfiles)){
    file.copy(
      from = paste0("../Figures/", myfiles[i]),
      to = "T:/", overwrite = TRUE
    )
  }
}

