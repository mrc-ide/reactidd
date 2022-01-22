# First clear the environment of variables
rm(list = ls(all = TRUE))
setwd("E:/Group/report/round15/Tmp_scripts/")

source("functions/load_packages.R")

# Loading required packages
pkgs <- c(
  "magrittr",
  "dplyr",
  "openxlsx"
)
load_packages(pkgs)

# Source any functions from the local file
source("functions/add_conf_ints.R")
source("functions/make_tables.R")
source("functions/overall_prev.R")
source("functions/formatting_functions.R")

# Choice of rounds
round_id=c(16)
print(round_id)

# Paths to files
output_file <- "E:/Group/report/round16/Tables/Table_one_triple_infected_vs_healthy_"
recoding_file <- "E:/Group/report/round16/Parameters/Recoding.xlsx"
recoding_from_cont_file <- "E:/Group/report/round16/Parameters/Recoding_from_continuous.xlsx"
output_tag <- Sys.Date()
annot_file <- "E:/Group/report/round16/Parameters/Variable_names.xlsx"

# Variable for test results
res_param <- c("estbinres", "influenzaa", "influenzab")

# Variable with participant ID
id_param="u_passcode"

# Variables for stratification
covs <- c(
  "gender_char", 
  "age",
  "region",
  "ethnic_new_char",
  "work_new_alt",
  "hh_size_cat", 
  "nchild",
  "imd_quintile", 
  # "covidcon_char", 
  "sympt_cat",
  "covida",
  "indoor",
  "outdoor",
  "shield2",
  "immunocompromised", "organ_specific", "psychological",
  "vax_status"
)

# Age bounds
age_lower_bound <- 0
age_upper_bound <- 200
# age_lower_bound <- 11
# age_upper_bound <- 18


## Loading and preparing the data

data_file <- paste0("E:/Group/saved_objects/rep", round_id, ".rds")

df_round <- data.frame(readRDS(data_file))

# Adding variable introduced in round 15
if (!"vax_status_noDate_v2" %in% colnames(df_round)) {
  df_round$vax_status_noDate_v2 <- df_round$vax_status_noDate
}

# Removing missing in estbinres
df_round <- df_round %>% filter(!is.na(estbinres)) %>% mutate(group = "Overall")

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

# Recoding indoor activity
ids=colnames(df_round)[grep("indoor_", colnames(df_round))]
ids=ids[!ids%in%c("indoor_7","indoor_8")]
indoor_yes=indoor_unknown=rep(0, nrow(df_round))
for (i in 1:length(ids)){
  table(df_round[,ids[i]])
  indoor_yes=indoor_yes+ifelse(df_round[,ids[i]]==1, yes=1, no=0)
  indoor_unknown=indoor_unknown+ifelse(!df_round[,ids[i]]%in%c(0,1), yes=1, no=0)
  # indoor_unknown=indoor_unknown+ifelse(df_round[,ids[i]]==-77, yes=1, no=0)
}
indoor_yes=ifelse(indoor_yes>0, yes=1, no=0)
indoor_unknown=ifelse(indoor_unknown>0, yes=1, no=0)

indoor=rep("No", nrow(df_round))
indoor[which(indoor_unknown==1)]="Unknown"
indoor[which(indoor_yes==1)]="Yes"

df_round$indoor=factor(indoor, levels=c("No","Yes","Unknown"))

# Recoding outdoor activity
ids=colnames(df_round)[grep("outdoor_", colnames(df_round))]
ids=ids[!ids%in%c("outdoor_7","outdoor_8")]
indoor_yes=indoor_unknown=rep(0, nrow(df_round))
for (i in 1:length(ids)){
  table(df_round[,ids[i]])
  indoor_yes=indoor_yes+ifelse(df_round[,ids[i]]==1, yes=1, no=0)
  indoor_unknown=indoor_unknown+ifelse(!df_round[,ids[i]]%in%c(0,1), yes=1, no=0)
  # indoor_unknown=indoor_unknown+ifelse(df_round[,ids[i]]==-77, yes=1, no=0)
}
indoor_yes=ifelse(indoor_yes>0, yes=1, no=0)
indoor_unknown=ifelse(indoor_unknown>0, yes=1, no=0)

outdoor=rep("No", nrow(df_round))
outdoor[which(indoor_unknown==1)]="Unknown"
outdoor[which(indoor_yes==1)]="Yes"

df_round$outdoor=factor(outdoor, levels=c("No","Yes","Unknown"))

# Recoding comorbidities
ids=paste0("healtha_", c(1, 14, 12, 
                         3, 6, 7, 8, 9,
                         15, 16, 17))
comorb_unknown=rep(0, nrow(df_round))
for (i in 1:length(ids)){
  table(df_round[,ids[i]])
  comorb_unknown=comorb_unknown+ifelse(!df_round[,ids[i]]%in%c(0,1), yes=1, no=0)
  # indoor_unknown=indoor_unknown+ifelse(df_round[,ids[i]]==-77, yes=1, no=0)
}
comorb_unknown=ifelse(comorb_unknown>0, yes=1, no=0)

comorbidities=rep("No", nrow(df_round))
comorbidities[which(comorb_unknown==1)]="Unknown"

mynames=c("immunocompromised", "organ_specific", "psychological")
mylist=list(c(1, 14, 12),
            c(3, 6, 7, 8, 9),
            c(15, 16, 17))
for (l in 1:length(mylist)){
  tmp=rep(0, nrow(df_round))
  ids=paste0("healtha_", mylist[[l]])
  for (i in 1:length(ids)){
    table(df_round[,ids[i]])
    tmp=tmp+ifelse(df_round[,ids[i]]==1, yes=1, no=0)
  }
  tmp=ifelse(tmp>0, yes=1, no=0)
  
  tmpcomorb=rep("No", nrow(df_round))
  tmpcomorb[which(comorb_unknown==1)]="Unknown"
  tmpcomorb[which(tmp==1)]="Yes"
  
  assign(mynames[l], tmpcomorb)
}
df_round$immunocompromised=factor(immunocompromised, levels=c("No","Yes","Unknown"))
df_round$organ_specific=factor(organ_specific, levels=c("No","Yes","Unknown"))
df_round$psychological=factor(psychological, levels=c("No","Yes","Unknown"))

# Restricting to 18-64
df_round=filter(df_round, 
                age > as.numeric(age_lower_bound) & 
                  age < as.numeric(age_upper_bound))

# Extracting covariate names
tmp <- read.xlsx(annot_file)
covs_names <- tmp[, 2]
names(covs_names) <- tmp[, 1]
covs_names <- covs_names[covs]

# Removing unused variables
df_round <- df_round[, c(id_param, res_param, covs)]

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
  
  # # Removing categories with less than 10 observations
  # if (any(table(x) < 10)) {
  #   toremove <- names(table(x))[which(table(x) < 10)]
  #   print(paste0("Excluding category ", renaming[toremove]))
  #   x[which(x == toremove)] <- NA
  #   renaming <- renaming[!names(renaming) %in% toremove]
  # }
  
  x <- factor(x, levels = names(renaming), labels = renaming)
  print(table(x, useNA = "always"))
  
  # Defining the reference level
  x <- relevel(x, ref = recoding[which(recoding[, 3] == "*"), 2])
  
  # # Removing unused levels
  # x=DropLevels(x)
  
  df_round[, covs_to_recode[i]] <- x
}

# Recoding continuous to categorical
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

df_round$influenzaa=as.numeric(as.character(factor(df_round$influenzaa, 
                                                   levels=c("negative","positive"),
                                                   labels=c(0,1))))

df_round$influenzab=as.numeric(as.character(factor(df_round$influenzab, 
                                                   levels=c("negative","positive"),
                                                   labels=c(0,1))))

df_round$triple_vs_healthy=factor(df_round$estbinres+df_round$influenzaa+df_round$influenzab,
                                  levels=c(0,3), labels=c("Healthy", "Triple"))

df_round$triple_vs_covid_only=rep(NA, nrow(df_round))
df_round$triple_vs_covid_only[which((df_round$estbinres==1)&(df_round$influenzaa==1)&(df_round$influenzab==1))]=1
df_round$triple_vs_covid_only[which((df_round$estbinres==1)&(df_round$influenzaa==0)&(df_round$influenzab==0))]=0

# Creating the round ID
df_round$round=rep(paste0("Round ", round_id), nrow(df_round))
mydata=df_round

# Using vaccinated (2 doses) vs not vaccinated 
mydata$linked=mydata$triple_vs_healthy


mydata=mydata[which(!is.na(as.character(mydata$linked))),]
table(mydata$linked, useNA="always")

mytable=NULL
for (k in 1:(length(covs)-1)){
  covariate=covs[k]
  
  # Full population
  tmp_full=table(mydata[,covs[k]])
  tmp_full_prop=FormatPrevalence(prop.table(tmp_full))
  
  # Stratified by linkage status
  tmp=table(mydata[,covs[k]], mydata[,"linked"])
  tmp_prop=FormatPrevalence(prop.table(tmp, margin=2))
  chi2=chisq.test(x=tmp[which(apply(tmp,1,sum)>10),])
  pval=formatC(chi2$p.value, format="e", digits=2)
  
  mytable_init=cbind(
    c(covs_names[k], rep(NA, nrow(tmp)-1)),
    rownames(tmp), 
    FormatCount(tmp_full), tmp_full_prop,
    FormatCount(tmp[,1]), tmp_prop[,1], 
    FormatCount(tmp[,2]), tmp_prop[,2],
    c(pval, rep(NA, nrow(tmp)-1)))
  
  # Logistic regression
  for (subset in c("full")){
    if (subset=="full"){
      tmpdata=mydata
    } else {
      tmpdata=mydata[which(mydata$round==paste0("Round ", subset)),]
    }
    
    tmpcont=apply(table(tmpdata[,covs[k]], tmpdata[,"linked"]),1,min)
    thr=2
    if (any(tmpcont<thr)){
      tmpdata[,covs[k]]=factor(tmpdata[,covs[k]], 
                               levels=levels(tmpdata[,covs[k]])[which(tmpcont>=thr)])
      tmpdata=tmpdata[which(!is.na(tmpdata[,covs[k]])),]
    }
    
    myformula=paste0("linked~",covs[k])
    mymodel=glm(as.formula(myformula), data=tmpdata, 
                family=binomial(link = "logit"))
    mymodel0=glm(linked~1, data=tmpdata, 
                 family=binomial(link = "logit"))
    coefs <- exp(mymodel$coefficients)
    coefs_ci <- exp(confint.default(mymodel))
    tmp_coef <- cbind(
      coefs[grep(paste0("^", covariate), names(coefs))],
      coefs_ci[grep(paste0("^", covariate), rownames(coefs_ci)), , drop = FALSE]
    )
    base_model <- FormatCI(FormatOR(tmp_coef))
    tmpnames <- names(coefs)[grep(paste0("^", covariate), names(coefs))]
    tmpnames <- gsub(covariate, "", tmpnames)
    
    myref <- levels(tmpdata[, covariate])[!levels(tmpdata[, covariate]) %in% tmpnames]
    base_model <- rbind("Ref", base_model)[,1]
    tmpnames <- c(myref, tmpnames)
    names(base_model)=tmpnames
    
    myanova=formatC(anova(mymodel, mymodel0, test="LRT")$`Pr(>Chi)`[2],
                    format="e", digits=2)
    
    assign(paste0("base_model_", subset), base_model)
    assign(paste0("myanova_", subset), myanova)
    
    mytable_init=cbind(mytable_init, 
                       base_model[rownames(tmp)],
                       c(myanova, rep(NA, nrow(tmp)-1)))
  }
  
  mytable=rbind(mytable, mytable_init)
}

# Including in file name
output_file=paste0(output_file, 
                   "_r", round_id, 
                   "_", age_lower_bound, 
                   "_", age_upper_bound)

# Saving the table
write.xlsx(mytable, rowNames=FALSE, colNames=TRUE, 
           paste0(output_file,".xlsx"))
file.copy(from=paste0(output_file,".xlsx"), 
          to="T:/", overwrite = TRUE)
