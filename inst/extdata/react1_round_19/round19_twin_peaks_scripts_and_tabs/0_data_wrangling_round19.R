# Data Wrangling of REACT-1-19 data


## Preamble and script setup

# First clear the environment of variables
rm(list=ls(all=TRUE))
setwd("E:/Group/")

# Source functions
source("E:/Group/functions/load_packages.R")
source("E:/Group/functions/data_cleaning_functions.R")
source("E:/Group/functions/create_react_week_number.R")

# Pull in packages needed
pkgs <- c("knitr", "dplyr","tidyr", "readr", "tidyverse", "janitor", 'data.table')
load_packages(pkgs)

# Paths to files
dirData <- "E:/Incoming _data/AG/Antigen19/"
fnData <- paste(
  dirData,"AntigenR19_FINAL File_INTUSE_03-APR-22.csv", sep="")
dirOutput="./saved_objects"
fnOutput="rep19.rds"

## Load data and add variables to working dataframe
# Load up the data.
dfDatRaw=fread(fnData, data.table = FALSE, na='')

# Look at the size of the dataframe and all the column names
dim(dfDatRaw)
names(dfDatRaw)
table(dfDatRaw$FinalResult,useNA = "always")



#' ## Create a dataframe with all covariates of interest
#' Make dataframe by grabbing the desired covariates of interest from the
#' raw data. Further manipulation into the final, clean data set will
#' happen below. Use standard data frame creation first.

#########################
# clean column names
########################
dfY <- dfDatRaw %>%
  clean_names()

dfRes <- dfY %>%
  dplyr::select(
    # unique codes for participants
    u_passcode, date_last_modified,
    # Substudy
    #postcourier,
    # weights
    wt_antigen,
    #DOB
    #r4_day, r4_month, r4_year,#r4_dk,
    # location
    region, u_lacode, u_lsoa, u_output_area, u_region, u_postcode,
    # u_utla_code,
    # individual characteritics (covariates of interest)
    u_gender, u_age, ethnic, smokenow, smokecig, vapnow, smokevap, empl, educ,
    #hosp01:hosp05,
    nadults, nadults1,
    nchild, nchild1, adultage_1:adultage_10, adultage1_1:adultage1_10,
    childage_1:childage_14, childage1_1:childage1_14,
    # work type
    worktyp1_1:worktyp1_7,
    worktyp2_1:worktyp2_2,worktyp2_10:worktyp2_12,worktyp2_3:worktyp2_7,worktyp2_8:worktyp2_9,
    worktyp3_1:worktyp3_7,
    worktyp4_1,worktyp4_10:worktyp4_12, worktyp4_2:worktyp4_7,worktyp4_8:worktyp4_9,
    # symptoms
    sympt_any1_1:sympt_any3_10, symptoth, #symptnowa,
    symptnowaw_1:symptnowaw_28, symptnowaw_27, symptfirst_1:symptfirst_27,
    symptst, symptlast_1:symptlast_28, symptfn, #cursympt,
    brediff, seekmed, feelun,
    # contacts
    covidcon_1:covidcon_3,
    contact1, contact4_1, contact4_2, contact4_3, contact4_4, #contact4_5,
    contact5, #carehome, #caretype, perscare,
    covidconnum,
    # deprivation
    imd_decile,
    # university
    campus2, edtype,
    contact1,# school,
    #indoor_1:outdoor_8,
    # contact1:contact5, carehome, caretype, perscare,
    # vaccination
    vaccine3, vaccdose, vaccinefirst, vaccinesecond, vaccinetype, vaccinetype_1:vaccinetypesym_4, 
    vaccineapp2, vaccineaccept,
    vaccine3sym, vaccdosesym, vaccinefirstsym, vaccinesecondsym,
    vaccinetypesym,
    #Demographics
    #dwelltyp, workstudypers1_1, workstudypers1_2,workstudypers1_3,workstudypers1_4,workstudypers1_5,
    #workstudypers2, studstay, edtype, campus2,
    # Abroad questions
    abroad, countryvisit1, countryvisit2,
    # pre-existing health conditions
    healtha_1:healtha_18, shield1, shield2,
    #leave1,
    leave2_1:leave2_3, leave2_5:leave2_10,# leave2_11_oth,
    #transp_1:transp_13,
    facecov,indmask,# outmask,
    #indoor_1:indoor_9, indoor_7:indoor_8, outdoor_1:outdoor_8,
    #bubble2, bubblenum, contactbub4_1:contactbub4_4, contactbub4_dk,
    # behava1:behava17, didn't find varibales in the round 8
    # results
    e_gene_c_tvalue, n_gene_c_tvalue, result, final_result, swabdate_new,
    #collected_date_time,
    swabdatepost, 
    # date_eurofinreceived,
    collecteddate,
    # flu variable
    fl_ureportdate, fl_utestdate, fluvacc, fluvaccdate,
    influenzaa, influenzaacpvalue, influenzab, influenzabcpvalue,
    date_flu_last_modified,
    #rm_postal_collection,r_mdaysinpost,
    # shedding data
    #assigned_shed,
    #shed_id_1, n_gene_c_tvalue_shed1, e_gene_shed1, e_gene_c_tvalue_shed1,
    #collected_date_time_shed1, result_shed1, final_result_shed1,
    #shed_id_2, n_gene_c_tvalue_shed2, e_gene_shed2, e_gene_c_tvalue_shed2,
    #collected_date_time_shed2, result_shed2, final_result_shed1
    # prior antigen test
    #agprev2:agprev4, agprevtype_1:agprevtype_3, agprevnum,
    covida, covidb_new
  ) %>%
  rename( res = result,
          gender = u_gender,
          region_id = region,
          lacode = u_lacode,
          lsoa = u_lsoa,
          postcode = u_postcode,
          output_area = u_output_area,
          age = u_age,
          ct1 = n_gene_c_tvalue,
          ct2= e_gene_c_tvalue,
          d_swab = swabdate_new,
          d_post = swabdatepost,
          d_rmcol = collecteddate,
          # d_euro = date_eurofinreceived,
          #rm_dip = r_mdaysinpost,
          vaccinated = vaccine3,
          #vaccine_doses = vaccdose,
          vaccine_first = vaccinefirst,
          vaccine_second = vaccinesecond,
          vaccine_type = vaccinetype,
          face_cov = facecov,
          # shedding
          #n_gene_ct_value_shed1 = n_gene_c_tvalue_shed1,
          #n_gene_ct_value_shed2 = n_gene_c_tvalue_shed2,
          #e_gene_ct_value_shed1 = e_gene_c_tvalue_shed1,
          #e_gene_ct_value_shed2 = e_gene_c_tvalue_shed2
  )

######################################################################
# Dates
######################################################################

dfRes <- dfRes %>%
  mutate(d_swab = as.Date(d_swab, format = "%m/%d/%Y"),
         d_post = as.Date(d_post, format = "%m/%d/%Y"),
         #d_euro = as.Date(d_euro, format = "%m/%d/%Y"),
         d_rmcol = as.Date(d_rmcol, format = "%m/%d/%Y"),
         date_last_modified = as.Date(date_last_modified, format = "%m/%d/%Y"))


######################################################################
# Date of Birth
######################################################################

#dfRes <- dfRes %>%
#  mutate(d_birth = ifelse(r4_dk==1, NA,
#                          paste(r4_year,"-",r4_month,"-",r4_day, sep="")),
#         d_birth = as.Date(d_birth))

######################################################################
# Results
######################################################################
dfRes <- dfRes %>%
  mutate(ct1 = as.numeric(ct1),
         ct2 = as.numeric(ct2))

dfRes <- dfRes %>%
  mutate(# convert VOID result to NA
    res = ifelse( res == "Void", NA, res),
    estbinres = ifelse(is.na(res), NA,
                       ifelse(res!="Detected", 0,
                              ifelse(ct1 > 0 & ct2 > 0, 1,
                                     ifelse(ct1 < 37 & ct1 > 0, 1, 0)))),
    estbinres35 = ifelse(is.na(res), NA,
                         ifelse(res!="Detected", 0,
                                ifelse(ct1 > 0 & ct2 > 0, 1,
                                       ifelse(ct1 < 35 & ct1 > 0, 1, 0)))),
    estbinres33 = ifelse(is.na(res), NA,
                         ifelse(res!="Detected", 0,
                                ifelse(ct1 > 0 & ct2 > 0, 1,
                                       ifelse(ct1 < 33 & ct1 > 0, 1, 0)))),
    estbinres30 = ifelse(is.na(res), NA,
                         ifelse(res!="Detected", 0,
                                ifelse(ct1 > 0 & ct2 > 0, 1,
                                       ifelse(ct1 < 30 & ct1 > 0, 1, 0)))),
    estbinres27 = ifelse(is.na(res), NA,
                         ifelse(res!="Detected", 0,
                                ifelse(ct1 > 0 & ct2 > 0, 1,
                                       ifelse(ct1 < 27 & ct1 > 0, 1, 0)))),
    res_ct12_gt0 = ifelse(is.na(res), NA,
                          ifelse(res != "Detected", 0,
                                 ifelse(ct1 > 0 & ct2 > 0, 1, 0))))


######################################################################
# number of adults and children in a household
######################################################################
dfRes <- dfRes %>%
  mutate(nadults = ifelse(nadults < 0, NA, nadults),
         nadults1 = ifelse(nadults1 < 0, NA, nadults1),
         nchild = ifelse(nchild < 0, NA, nchild),
         nchild1 = ifelse(nchild1 < 0, NA, nchild1),
         nchild2 = ifelse(nchild==0, 0, "1+")) %>%
  mutate(nchild2 = ifelse(is.na(nchild2), "NA", nchild2))


######################################################################
# Covid Contact
######################################################################
dfRes <- dfRes %>%
  mutate(covidcon = case_when(
    covidcon_1 == 1 ~ 1,
    covidcon_2 == 1 ~ 2,
    covidcon_3 == 1 ~ 3
  ))


######################################################################
# Deprivation
######################################################################
dfRes <- dfRes %>%
  mutate(imd_quintile = case_when(
    imd_decile %in% c(1,2) ~ 1,
    imd_decile %in% c(3,4) ~ 2,
    imd_decile %in% c(5,6) ~ 3,
    imd_decile %in% c(7,8) ~ 4,
    imd_decile %in% c(9,10) ~ 5
  ))


######################################################################
# Smoke
######################################################################
dfRes <- dfRes %>%
  mutate(smokenow = ifelse(smokenow < 0 | smokenow == 3, NA, smokenow))

######################################################################
# Gender
######################################################################
dfRes <- dfRes %>%
  mutate(gender = ifelse(gender %in% c(9, -92), NA, gender))

######################################################################
# Shield
######################################################################
dfRes <- dfRes %>%
  mutate(shield1 = ifelse(shield1 < 0, NA, shield1),
         shield2 = ifelse(shield2 < 0, NA, shield2))

######################################################################
# the number of contacts you made yesterday (not inclduing your household)
######################################################################
dfRes <- dfRes %>%
  mutate(contact1 = ifelse(contact1 < 0, NA, contact1))


######################################################################
# Have you/your children gone back to school yet in person?
######################################################################
#dfRes <- dfRes %>%
#  mutate(school = case_when(
#    school < 0 ~ NA_character_,
#    school == 1 ~ "Yes",
#    school %in% c(2,3) ~ "No",
#    school == 4 ~ NA_character_))


######################################################################
# Add id variable
######################################################################
dfRes$id <- seq(1, nrow(dfRes))


######################################################################
#################################################
# Date variables
######################################################################
### OLIVER TO UPDATE XXXXXXXX!!!!!!!!!!!!!!!!!!!
## SHOULD BE GETTING RM date to replace d_euro -1 later

dfRes <- mutate(dfRes, d_comb = as.Date(ifelse(is.na(d_swab)==TRUE, d_rmcol,
                                               ifelse(d_swab<as.Date("2022-03-08"), d_rmcol,
                                                      ifelse(d_swab>as.Date("2022-05-11"), d_rmcol, d_swab))), origin="1970-01-01"))

dfRes <- mutate(dfRes, d_comb = as.Date(ifelse(d_comb<as.Date("2022-03-08"), NA,d_comb), origin="1970-01-01"))


########################################################################
# Region
########################################################################
#' create a new region (character) variable
dfRes <- get_region_name(df = dfRes,
                         numeric_region_column = "region_id")

####################################################################
# Age
####################################################################
# create age groups: 5-12, 13-17, 18-24, 25-34, ..., 65+
age_breaks <- c(4,12,17,24,34,44,54,64,max(dfRes$age))
dfRes$age_group <- cut(dfRes$age,
                       breaks = age_breaks,
                       right = TRUE)

## add two more age groups
age_breaks2 <- c(4,12,17,24,34,44,54,64,74, max(dfRes$age))
dfRes$age_group2 <- cut(dfRes$age,
                        breaks = age_breaks2,
                        right = TRUE)

#' split age into thirds 5 - 35, 36 - 56, 57 +
dfRes <- dfRes %>%
  mutate(age_group_tert = factor(ifelse(age <= 35, "Young",
                                        ifelse( age < 57, "Middle", "Old")),
                                 levels = c("Young", "Middle", "Old"))
  )

##########################################################################
## Symptoms
##########################################################################
# mutate symptom variables
neg_to_na <- function(x) ifelse(x < 0, NA, x)

dfRes <- dfRes %>%
  mutate_at(vars(matches("sympt_any")), neg_to_na) %>%
  mutate_at(vars(matches("symptnowa")), neg_to_na) %>%
  mutate(symptnow_bin = ifelse( feelun == 2 |
                                  (sympt_any1_5 == 1 & sympt_any2_14 == 1 & sympt_any3_10 == 1), 0, 1),
         # Create a categorical variable for classic symptoms/other symptoms/no symptoms in last week
         sympt_cat = factor(ifelse(symptnow_bin == 0, "No symptoms",
                                   ifelse(symptnow_bin == 1 &
                                            (sympt_any1_1 == 1 | sympt_any1_2 == 1 | sympt_any1_3 == 1 | sympt_any1_4 == 1), "Classic COVID symptoms",
                                          "Other symptoms")),levels = c("Classic COVID symptoms","Other symptoms","No symptoms")),
         res_asympt = ifelse(is.na(estbinres), NA,
                             ifelse(symptnow_bin == 0 & estbinres == 1, 1, 0))
  ) %>%
  mutate(sympt_cat = ifelse(is.na(sympt_cat), "NA", as.character(sympt_cat))) %>%
  mutate(sympt_cat = factor(sympt_cat, levels = c("Classic COVID symptoms","Other symptoms","No symptoms", "NA")))

# %>%
# activities variables
#mutate_at(vars(matches("indoor")), neg_to_na) %>%
#mutate_at(vars(matches("outdoor")), neg_to_na) %>%
#mutate_at(vars(matches("transp")), neg_to_na) %>%
#mutate_at(vars(matches("leave2")), neg_to_na) %>%
#mutate(
#  transp_single = ifelse(transp_1 == 1 | transp_2 == 1 | transp_3 == 1, 1,
#                         ifelse(is.na(transp_1) & is.na(transp_2) & is.na(transp_3), NA, 0)),
#  transp_car = ifelse(transp_4 == 1 | transp_5 == 1 | transp_6 == 1 | transp_7 == 1 | transp_8 == 1, 1,
#                      ifelse(is.na(transp_4) & is.na(transp_5) & is.na(transp_6) & is.na(transp_7) & is.na(transp_8), NA, 0)),
#  transp_mass = ifelse(transp_9 == 1 | transp_10 == 1 | transp_11 == 1 | transp_12 == 1, 1,
#                       ifelse(is.na(transp_9) & is.na(transp_10) & is.na(transp_11) & is.na(transp_12), NA, 0))
#)

####################################################################
#' #' Create new symptom now (symptoms in last week) variables
#' #' Each new variable will have three levels:
#' #' * 0 = no symptoms
#' #' * 1 = yes, has symptom i in last week
#' #' * 2 = yes, has symptoms other than symptom i in last week
index <- which(startsWith(names(dfRes),"symptnowaw_")) # index of symptnowaw columns in dfRes
n_col <- ncol(dfRes)

for (i in 1:length(index)){
  dfRes[ , n_col + i ] <- ifelse(dfRes$symptnow_bin == 0, 0,
                                 ifelse(dfRes[,index[i]] == 1, 1, 2))
  names(dfRes)[ n_col + i ] <- paste0("symptnow_new_", i)
}

# Vivi's code - create last week symptoms variable
# sympt_cat2 is for people still have any symptoms in last week, then classify their symptoms occurred in last week
# sympt_cat3 is for people still have any symptoms in last week, then classify their symptoms occurred in last month
# dfRes <- dfRes %>%
#   mutate(sympt_lastweek_bin = ifelse( feelun == 2 |
#                                         (sympt_any1_5 == 1 & sympt_any2_14 == 1 & sympt_any3_10 == 1) |
#                                         symptnowaw_28 == 1, 0, 1),
#          sympt_cat2 = factor(ifelse(sympt_lastweek_bin == 0, "No symptoms",
#                                     ifelse(sympt_lastweek_bin == 1 &
#                                              (symptnowaw_1 == 1 | symptnowaw_2 == 1 | symptnowaw_3 == 1 | symptnowaw_4 == 1), "Classic COVID symptoms",
#                                            "Other symptoms")),levels = c("Classic COVID symptoms","Other symptoms","No symptoms")),
#          sympt_cat3 = factor(ifelse(sympt_lastweek_bin == 0, "No symptoms",
#                                     ifelse(sympt_lastweek_bin == 1 &
#                                              (sympt_any1_1 == 1 | sympt_any1_2 == 1 | sympt_any1_3 == 1 | sympt_any1_4 == 1), "Classic COVID symptoms",
#                                            "Other symptoms")),levels = c("Classic COVID symptoms","Other symptoms","No symptoms"))) %>%
#   mutate(sympt_cat = ifelse(is.na(sympt_cat), "NA", as.character(sympt_cat)),
#          sympt_cat2 = ifelse(is.na(sympt_cat2), "NA", as.character(sympt_cat2)),
#          sympt_cat3 = ifelse(is.na(sympt_cat3), "NA", as.character(sympt_cat3))) %>%
#   mutate(sympt_cat = factor(sympt_cat, levels = c("Classic COVID symptoms","Other symptoms","No symptoms", "NA")),
#          sympt_cat2 = factor(sympt_cat2, levels = c("Classic COVID symptoms","Other symptoms","No symptoms", "NA")),
#          sympt_cat3 = factor(sympt_cat3, levels = c("Classic COVID symptoms","Other symptoms","No symptoms", "NA")))

names(dfRes)

########################################################################
### Household size/ composition
######################################################################
# To get household size, we need to sum NADULTS and NCHILD for respondents
# aged 5-12 or 18+ and NADULTS1 and NCHILD1 for respondents aged 13-17.
# Create new household size variable
dfRes$hh_size <- apply(dfRes[,c("nadults", "nadults1", "nchild", "nchild1")],1,
                       sum, na.rm = TRUE)
dfRes <- dfRes %>%
  mutate(hh_size = ifelse(hh_size <= 0 , NA, hh_size),
         hh_size_cat = as.factor(ifelse(hh_size >= 6, "6+", hh_size)),
         hh_size_tert = factor(ifelse(hh_size %in% c(1,2), "1-2",
                                      ifelse(hh_size %in% c(3:5), "3-5", "6+")),
                               levels = c("1-2", "3-5", "6+")),
         # #hh_any_school_age = ifelse((nchild > 0 | nchild1 > 0) & , 1, 0),
         # hh_primary_only = ifelse(nchild > 0 | nchild1 > 0,
         #                          ifelse(childage1_1:childage1_10 > 5, 1, 0)),
  )

#' ### Create new household variables based on school aged children
# Sorry the code is messy and long!
# households with primary school aged children only
test_dat <- dfRes %>%
  dplyr::select(id, childage_1:childage1_14) %>%
  filter_at(vars(starts_with("childage")), any_vars(. > 4)) %>%
  filter_at(vars(starts_with("childage")), all_vars(. < 12))

dfRes$household_primary_only <- 0
dfRes$household_primary_only[test_dat$id] <- 1

# households with secondary school aged children only
neg_to_pos <- function(x) ifelse(x == -91, 91, x)

test_dat <- dfRes %>%
  dplyr::select(id, childage_1:childage1_14) %>%
  # change -91 to 91 for filtering
  mutate_at(vars(starts_with("childage")), neg_to_pos) %>%
  filter_at(vars(starts_with("childage")), all_vars(. > 11)) %>%
  filter_at(vars(starts_with("childage")), any_vars(. < 91))

dfRes$household_secondary_only <- 0
dfRes$household_secondary_only[test_dat$id] <- 1

# households with any school aged children
test_dat <- dfRes %>%
  dplyr::select(id, childage_1:childage1_14) %>%
  # change -91 to 91 for filtering
  mutate_at(vars(starts_with("childage")), neg_to_pos) %>%
  filter_at(vars(starts_with("childage")), all_vars(. > 4)) %>%
  filter_at(vars(starts_with("childage")), any_vars(. < 91))

dfRes$household_any_school_aged <- 0
dfRes$household_any_school_aged[test_dat$id] <- 1

####################################################################
### Visited a hospital in the last 2 weeks
####################################################################
# Convert hosp01:hosp05 into a single variable
hosp_dat <- dfRes %>%
  dplyr::select(id, hosp01:hosp05) %>%
  mutate_at(vars(matches("hosp")), neg_to_na) %>%
  group_by(id) %>%
  gather(hosp, val, -id) %>%
  filter(val == 1) %>%
  dplyr::select(-val) %>%
  arrange(id) %>%
  mutate(hosp = as.numeric(gsub("hosp0","", hosp))) %>%
  summarise(min = min(hosp)) %>%
  rename("hosp" = "min")

#' Bind to dfRes  and drop original variables
dfRes <- left_join(dfRes, hosp_dat, by = "id")

#####################################################################
### Ethnicity
#####################################################################
# We'll create broad categories for the ethnic variable:
# * 1 = white
# * 2 = mixed
# * 3 = asian / asian british
# * 4 = black / african / caribbean / black british
# * 5 = other

dfRes$ethnic_new <- ifelse(dfRes$ethnic %in% c(1:4), "white",
                           ifelse(dfRes$ethnic %in% c(5:8), "mixed",
                                  ifelse(dfRes$ethnic %in% c(9:13), "asian",
                                         ifelse(dfRes$ethnic %in% c(14:16), "black",
                                                ifelse(dfRes$ethnic %in% c(17,18),"other", NA)
                                         )
                                  )
                           )
)

#####################################################################
### Detail ethnicity
#####################################################################
dfRes <- dfRes %>%
  mutate(ethnic19_char = ifelse(ethnic==1,"English/Welsh/Scottish/Northern Irish/British",
                                ifelse(ethnic==2,"Irish",
                                       ifelse(ethnic==3,"Gypsy or Irish Traveller",
                                              ifelse(ethnic==4,"Other white background",
                                                     ifelse(ethnic==5,"White and Black Caribbean",
                                                            ifelse(ethnic==6,"White and Black African",
                                                                   ifelse(ethnic==7,"White and Asian",
                                                                          ifelse(ethnic==8,"Other Mixed/Multiple ethnic background",
                                                                                 ifelse(ethnic==9,"Indian",
                                                                                        ifelse(ethnic==10,"Pakistani",
                                                                                               ifelse(ethnic==11,"Bangladeshi",
                                                                                                      ifelse(ethnic==12,"Chinese",
                                                                                                             ifelse(ethnic==13,"Other Asian background",
                                                                                                                    ifelse(ethnic==14,"African",
                                                                                                                           ifelse(ethnic==15,"Caribbean",
                                                                                                                                  ifelse(ethnic==16,"Other Black/African/Caribbean background",
                                                                                                                                         ifelse(ethnic==17,"Arab",
                                                                                                                                                ifelse(ethnic==18|ethnic==20,"Other ethnic group",
                                                                                                                                                       ifelse(ethnic==19,"Prefer not to say", "NA")))))))))))))))))))) %>%
  mutate(ethnic19_char = as.factor(ethnic19_char))




######################################################################
# Work Type (Key worker status)
######################################################################
work_dat1 <- dfRes %>%
  dplyr::select(id, worktyp1_1:worktyp1_7) %>%
  mutate_at(vars(matches("worktyp1_")), neg_to_na) %>%
  group_by(id) %>%
  gather(work_type, val, -id) %>%
  filter(val == 1) %>%
  dplyr::select(-val) %>%
  arrange(id) %>%
  mutate(work_type = as.numeric(gsub("worktyp1_","", work_type))) %>%
  summarise(min = min(work_type)) %>%
  rename("work_type1" = "min")

#' Bind to dfRes  and drop original variables
dfRes <- left_join(dfRes, work_dat1, by = "id") %>%
  dplyr::select(-worktyp1_1, -worktyp1_2, -worktyp1_3, -worktyp1_4,
                -worktyp1_5, -worktyp1_6, -worktyp1_7)

###################################
# detail work type 1
##################################

dfRes <- dfRes %>%
  mutate(work_type1_char = ifelse(work_type1==1, "Health care workers with direct patient contact",
                                  ifelse(work_type1==2, "Health care workers with no patient contact",
                                         ifelse(work_type1==3,"Care home workers with direct contact with clients",
                                                ifelse(work_type1==4,"Care home workers without contact with clients",
                                                       ifelse(work_type1==5, "Other essential/key workers",
                                                              ifelse(work_type1==6,"None of these",
                                                                     ifelse(work_type1==7,"Don't know","WT1 - NA")))))))) %>%
  mutate(work_type1_char = ifelse(is.na(work_type1_char),"WT1 - NA",work_type1_char)) %>%
  mutate(work_type1_char = as.factor(work_type1_char))




#' Work type 2 - vivi
#' Only answered if work_type1 = 5, 6, 7
work_dat2 <- dfRes %>%
  dplyr::select(id, worktyp2_1:worktyp2_7,worktyp2_8:worktyp2_11) %>%
  mutate_at(vars(matches("worktyp2_")), neg_to_na) %>%
  group_by(id) %>%
  gather(work_type2, val, -id) %>%
  filter(val == 1) %>%
  dplyr::select(-val) %>%
  arrange(id) %>%
  mutate(work_type2 = as.numeric(gsub("worktyp2_","", work_type2))) %>%
  summarise(min = min(work_type2)) %>%
  rename("work_type2" = "min")

#' Bind to dfRes  and drop original variables
dfRes <- left_join(dfRes, work_dat2, by = "id") %>%
  dplyr::select(-worktyp2_1, -worktyp2_2, -worktyp2_3, -worktyp2_4,
                -worktyp2_5, -worktyp2_6, -worktyp2_7,-worktyp2_8,
                -worktyp2_10,-worktyp2_11)

###################################
# detail work type 3
##################################
work_dat3 <- dfRes %>%
  dplyr::select(id, worktyp3_1:worktyp3_7) %>%
  mutate_at(vars(matches("worktyp3_")), neg_to_na) %>%
  group_by(id) %>%
  gather(work_type, val, -id) %>%
  filter(val == 1) %>%
  dplyr::select(-val) %>%
  arrange(id) %>%
  mutate(work_type = as.numeric(gsub("worktyp3_","", work_type))) %>%
  summarise(min = min(work_type)) %>%
  rename("work_type3" = "min")

#' Bind to dfRes  and drop original variables
dfRes <- left_join(dfRes, work_dat3, by = "id") %>%
  dplyr::select(-worktyp3_1, -worktyp3_2, -worktyp3_3, -worktyp3_4,
                -worktyp3_5, -worktyp3_6, -worktyp3_7)


dfRes <- dfRes %>%
  mutate(work_type3_char = ifelse(work_type1==1|work_type1==2|work_type1==3|work_type1==4|work_type1==5, ifelse(work_type3==1,"Health care workers with direct patient contact",
                                                                                                                ifelse(work_type3==2,"Healthcare workers with no patient contact",
                                                                                                                       ifelse(work_type3==3,"Care home workers with direct contact with clients",
                                                                                                                              ifelse(work_type3==4,"Care home workers without contact with clients",
                                                                                                                                     ifelse(work_type3==5,"Other essential/key workers",
                                                                                                                                            ifelse(work_type3==6,"None of these",
                                                                                                                                                   ifelse(work_type3==7,"Don't know","WT3 - NA"))))))),"WT3_No")) %>%
  mutate(work_type3_char = ifelse(is.na(work_type3_char),"WT3 - NA",work_type3_char)) %>%
  mutate(work_type3_char = as.factor(work_type3_char))


#' Quick check for duplicate ids.
length(dfRes$id) - length(unique(dfRes$id))
dim(dfDatRaw)[1] - length(dfRes$id)
dim(dfRes)[1]
#' Should be 0.
#'
#' ### Combine work type 1 with empl into single variable
#' New variable is called work_new with the following levels:
#' 1. health care worker with direct contact with patients
#' 2. health care worker without direct contact with patients
#' 3. care home worker with direct contact with clients
#' 4. care home worker without direct contact with clients
#' 5. essential/key worker
#' 6. other worker (empl %in% (full-time, part-time, self-employed))
#' 7. not full-time, part-time, self-employed
#' 8. unknown
dfWorkNew <- dfRes %>%
  dplyr::select(id, work_type1, empl) %>%
  mutate(work_new = work_type1,
         # change unknown work status from 7 to 8
         work_new = ifelse(work_new == 7, 8, work_new),
         # people who are NOT full-time, part-time, or self-employed
         work_new = ifelse(empl %in% 4:11, 7, work_new),
         work_new = ifelse(work_new == 8, NA, work_new)
  )

#' Join new column to dfRes
dfWorkNew <- dfWorkNew %>%
  dplyr::select(-work_type1, -empl)

dfRes <- left_join(dfRes, dfWorkNew, by = "id")


#' Collapse categories of work_new, so there are no zeros
dfRes <- dfRes %>%
  mutate(work_new_alt = case_when(
    work_new %in% 1:4 ~ "HCW/CHW",
    work_new == 5 ~ "Key worker (other)",
    work_new == 6 ~ "Other worker",
    work_new == 7 ~ "Not FT, PT, SE"
  ),
  work_new_alt_ext = case_when(
    work_new_alt == "HCW/CHW" ~ "HCW/CHW",
    work_new_alt == "Key worker (other)" & work_type2 %in% 1:7 ~ "Key worker (other) - public facing",
    work_new_alt == "Key worker (other)" & work_type2 %in% c(8,9) ~ "Key worker (other) - not public facing",
    #work_new_alt == "Key worker (other)" & work_type2 == 9 ~ "Key worker (other) - not outside home",
    work_new_alt == "Other worker" & work_type2 %in% 1:7 ~ "Other worker - public facing",
    work_new_alt == "Other worker" & work_type2 %in% c(8,9) ~ "Other worker - not public facing",
    #work_new_alt == "Other worker" & work_type2 == 9 ~ "Other worker - not outside home",
    work_new_alt == "Not FT, PT, SE" ~ "Not FT, PT, SE"
  ),
  covidcon_char = case_when(
    covidcon == 1 ~ "Yes, contact with confirmed/tested COVID-19 case",
    covidcon == 2 ~ "Yes, contact with suspected COVID-19 case",
    covidcon == 3 ~ "No"
  )
  )

dfRes <- dfRes %>%
  mutate(work_new_alt_ext = ifelse(is.na(work_new_alt_ext),"NA",work_new_alt_ext),
         work_new_alt = ifelse(is.na(work_new_alt),"NA",work_new_alt),
         covidcon_char = ifelse(is.na(covidcon_char),"NA",covidcon_char)) %>%
  mutate(work_new_alt_ext = factor(work_new_alt_ext,
                                   levels = c("HCW/CHW",
                                              "Key worker (other) - public facing",
                                              "Key worker (other) - not public facing",
                                              "Other worker - public facing",
                                              "Other worker - not public facing", "Not FT, PT, SE","NA")),
         work_new_alt = factor(work_new_alt,
                               levels = c("HCW/CHW", "Key worker (other)",
                                          "Other worker", "Not FT, PT, SE","NA")),
         covidcon_char = factor(covidcon_char,
                                levels = c("Yes, contact with confirmed/tested COVID-19 case",
                                           "Yes, contact with suspected COVID-19 case",
                                           "No", "NA"))
  )



#' Create new categorical Re-label and re-order factor levels
dfRes <- dfRes %>%
  mutate(gender_char = factor(recode(gender, `1`= "Male", `2` = "Female"),
                              levels = c("Male", "Female")),
         age_group_char = case_when(age %in% c(5:12) ~ "5-12",
                                    age %in% c(13:17) ~ "13-17",
                                    age %in% c(18:24) ~ "18-24",
                                    age %in% c(25:34) ~ "25-34",
                                    age %in% c(35:44) ~ "35-44",
                                    age %in% c(45:54) ~ "45-54",
                                    age %in% c(55:64) ~ "55-64",
                                    age >= 65 ~ "65+"),
         age_group_char2 = case_when(age %in% c(5:12) ~ "5-12",
                                     age %in% c(13:17) ~ "13-17",
                                     age %in% c(18:24) ~ "18-24",
                                     age %in% c(25:34) ~ "25-34",
                                     age %in% c(35:44) ~ "35-44",
                                     age %in% c(45:54) ~ "45-54",
                                     age %in% c(55:64) ~ "55-64",
                                     age %in% c(65:74) ~ "65-74",
                                     age >= 75 ~ "75+"),
         age_group_char3 = case_when(age %in% c(5:11) ~ "5-11",
                                     age %in% c(12:17) ~ "12-17",
                                     age %in% c(18:24) ~ "18-24",
                                     age %in% c(25:34) ~ "25-34",
                                     age %in% c(35:44) ~ "35-44",
                                     age %in% c(45:54) ~ "45-54",
                                     age %in% c(55:64) ~ "55-64",
                                     age %in% c(65:74) ~ "65-74",
                                     age >= 75 ~ "75+"),
         region = factor(region, levels = c("North East", "North West", "Yorkshire and The Humber",
                                            "East Midlands","West Midlands","East of England",
                                            "London", "South East", "South West")),
         work_new_char = factor(recode(work_new,
                                       `1` = "Health care worker with direct contact with patients",
                                       `2` = "Health care worker without direct contact with patients",
                                       `3` = "Care home worker with direct contact with clients",
                                       `4` = "Care home worker without direct contact with clients",
                                       `5` = "Other essential/key worker",
                                       `6` = "Other worker",
                                       `7` = "Not full-time, part-time, or self-employed"
         ),
         levels = c("Health care worker with direct contact with patients",
                    "Health care worker without direct contact with patients",
                    "Care home worker with direct contact with clients",
                    "Care home worker without direct contact with clients",
                    "Other essential/key worker", "Other worker",
                    "Not full-time, part-time, or self-employed")),
         ethnic_new_char = factor(recode(ethnic_new, white = "White", asian = "Asian / Asian British",
                                         black = "Black / African / Caribbean / Black British",
                                         mixed = "Mixed", other = "Other"),
                                  levels = c("Asian / Asian British",
                                             "Black / African / Caribbean / Black British",
                                             "Mixed", "Other","White")),
         # smokenow_char = factor(recode(smokenow, `1` = "Yes", `2` = "No"),
         #                        levels = c("Yes", "No")),
         hosp_char = factor(recode(hosp,
                                   `1` = "Yes, I have",
                                   `2` = "Yes, my child has",
                                   `3` = "Yes, someone in my household has",
                                   `4` = "No",
                                   `5` = "Don't Know"),
                            levels = c("Yes, I have", "Yes, my child has", "Yes, someone in my household has", "No", "Don't Know"))
         #school = factor(school, levels = c("Yes", "No"))
         # leave1_char = factor(recode(leave1, `1` = "Yes", `2` = "No"),
         #                      levels = c("Yes", "No"))

  ) %>%
  mutate(gender_char = ifelse(is.na(gender_char),"NA",as.character(gender_char)),
         age_group_char3 = factor(age_group_char3,levels = c("5-11", "12-17", "18-24", "25-34",
                                                             "35-44", "45-54", "55-64", "65-74","75+")),
         # age_group_char = ifelse(is.na(age_group_char),"NA",age_group_char),
         region = ifelse(is.na(region),"NA",as.character(region)),
         ethnic_new_char = ifelse(is.na(ethnic_new_char),"NA",as.character(ethnic_new_char)),
         hh_size_cat = ifelse(is.na(hh_size_cat),"NA",as.character(hh_size_cat)),
         hosp_char = ifelse(is.na(hosp_char),"NA",as.character(hosp_char)),
         imd_quintile = ifelse(is.na(imd_quintile),"NA",as.character(imd_quintile))) %>%
  mutate(gender_char = factor(gender_char, levels = c("Male", "Female", "NA")),
         age_group_char = as.factor(age_group_char),
         region = factor(region, levels = c("North East", "North West", "Yorkshire and The Humber",
                                            "East Midlands","West Midlands","East of England",
                                            "London", "South East", "South West")),
         work_new_char = factor(work_new_char,levels = c("Health care worker with direct contact with patients",
                                                         "Health care worker without direct contact with patients",
                                                         "Care home worker with direct contact with clients",
                                                         "Care home worker without direct contact with clients",
                                                         "Other essential/key worker", "Other worker",
                                                         "Not full-time, part-time, or self-employed","NA")),
         # work_new_alt_ext = as.factor(work_new_alt_ext),
         ethnic_new_char = factor(ethnic_new_char, levels = c("Asian / Asian British",
                                                              "Black / African / Caribbean / Black British",
                                                              "Mixed", "Other","White","NA")),
         hh_size_cat = as.factor(hh_size_cat),
         # smokenow_char = factor(smokenow_char, levels = c("Yes", "No","NA")),
         hosp_char = factor(hosp_char,
                            levels = c("Yes, I have", "Yes, my child has", "Yes, someone in my household has", "No", "Don't Know","NA")),
         imd_quintile = as.factor(imd_quintile)
  )



#######################################################################
# Assign week number to the data
#######################################################################

end_date <- max(dfRes$d_comb, na.rm = TRUE)

dfRes <- dfRes %>%
  left_join(create_react_week_number(end_date = end_date),
            by = c("d_comb" = "date"))


#################################################################
# vaccination questions
#################################################################
# make the columns all nice and tidy

dfRes <- dfRes %>%
  mutate(vaccinated = case_when(
    vaccinated %in% c(-77, -91, -92) ~ NA_character_,
    vaccinated == 0 ~ NA_character_,
    vaccinated == 1 ~ "yes",
    vaccinated == 2 ~ "no",
    vaccinated == 3 ~ "unsure if received as part of trial"
  ),
  vaccine_doses = case_when(
    vaccdose %in% c(-77, -91, -92) ~ NA_character_,
    vaccdose == 0 ~ NA_character_,
    vaccdose == 1 ~ "yes",
    vaccdose == 2 ~ "no"),
  vaccine_first = as.Date(vaccine_first, format = "%m/%d/%Y"),
  vaccine_second = as.Date(vaccine_second, format = "%m/%d/%Y"),
  vaccinefirstsym = as.Date(vaccinefirstsym, format = "%m/%d/%Y"),
  vaccinesecondsym = as.Date(vaccinesecondsym, format = "%m/%d/%Y"),
  had_vaccine_first = ifelse(is.na(vaccine_first), FALSE, TRUE),
  had_vaccine_second = ifelse(is.na(vaccine_second), FALSE, TRUE),
  # vaccine_type = case_when(
  #   vaccine_type %in% c(-77, -91, -92) ~ NA_character_,
  #   vaccine_type == 1 ~ "Pfizer/BioNTtech",
  #   vaccine_type == 2 ~ "AstraZeneca/Oxford",
  #   vaccine_type == 3 ~ "Moderna",
  #   vaccine_type == 4 ~ "don't know"
  # )
  )


###############################
# Vaccination - David's code, edited by Vivi
#############################

dfRes <- dfRes %>%
  mutate(vaccine_first_sym = as.Date(vaccinefirstsym, "%m/%d/%Y"),
         vaccine_to_swab_time = d_comb-vaccinefirstsym,
         vaccine_first_num = as.Date(vaccine_first, "%m/%d/%Y"),
         vaccine_to_swab_time = ifelse(vaccinefirstsym>0, d_comb - vaccinefirstsym, d_comb-vaccine_first_num),
         over70 = ifelse(age>=70,1,0),
         # vax_type (type = pfizer, oxford,moderna, other) here
         # refers to people who only received one does prior to the swab
         vax_pfizer = ifelse(vaccine_type==1 | vaccinetypesym == 1,1,0),# & vaccdosesym==1 & vaccine_to_swab_time>0,1,0),
         vax_oxford = ifelse(vaccine_type==2 | vaccinetypesym == 2,1,0),# & vaccdosesym==1 & vaccine_to_swab_time>0,1,0),
         vax_moderna = ifelse(vaccine_type==3 | vaccinetypesym == 3,1,0),# & vaccdosesym==1 & vaccine_to_swab_time>0,1,0),
         vax_other = ifelse(vaccine_type==4 | vaccinetypesym == 4,1,0),#  & vaccdosesym==1 & vaccine_to_swab_time>0,1,0),
         vax_any = ifelse(vax_pfizer+vax_oxford+vax_moderna+vax_other>0,1,0),
         vax_double = ifelse(max(vaccdose,vaccdosesym)>1,1,0))
#unvax = ifelse(vaccdosesym<0 | is.na(vaccdosesym)==TRUE | vaccine_to_swab_time<=0,1,0)) #vaccdosessym never NA)

# People who received 2 does prior to the swab or received the first does >=14 days prior to the swab,
# will be treated as vaccinated
# Note that there there people receievd more than two does vaccdosesym ==3,
# we ignore them here and David will check against linkage data later
# Also, combine separate vax types into one single variable

#Variable names (notrename commented out)
#vaccinated = vaccine3,
#vaccine_doses = vaccdose,
#vaccine_first = vaccinefirst,
#vaccine_second = vaccinesecond,
#vaccine_type = vaccinetype,

dfRes <- dfRes %>%
  mutate(vax_status = ifelse(vaccine_to_swab_time > 0 | vaccdosesym > 0, "Yes", "No"),
         vax_type_char = case_when(vax_pfizer == 1 ~ "pfizer",
                                   vax_oxford == 1 ~ "oxford",
                                   vax_moderna == 1 ~ "moderna",
                                   vax_other == 1 ~ "other")) %>%
  mutate(vax_status = ifelse(is.na(vax_status), "NA", vax_status),
         vax_type_char = ifelse(is.na(vax_type_char), "NA", vax_type_char)) %>%
  mutate(vax_status = factor(vax_status, levels = c("No", "Yes", "NA")),
         vax_type_char = factor(vax_type_char, levels = c("pfizer", "oxford", "moderna", "other", "NA"))) %>%
  #mutate(vaccine3sym = as.numeric(vaccine3sym),
  #vacine_doses = as.numeric(vaccine_doses),
  #vaccindosessym = as.numeric(vaccinedosessym))
  #mutate(vaccdose = as.numeric(vaccine_doses)) %>%
  mutate(vaccine_second_sym = as.Date(vaccinesecondsym, "%m/%d/%Y"),
         time_after_2dose = ifelse(vaccinesecondsym>0, d_comb - vaccinesecondsym, d_comb-vaccine_second_sym)) %>%
  mutate(vax_status_number = case_when((vaccine3sym==1 & vaccdose==1 & vaccdosesym<=1) ~ 1,
                                       (vaccine3sym==1 & vaccdosesym==1) ~ 1,
                                       (vaccine3sym==1 & vaccdose==2 & vaccdosesym<=2) ~ 2,
                                       (vaccine3sym==1 & vaccdosesym==2) ~ 2,
                                       (vaccine3sym==1 & vaccdose>=3 & vaccdosesym<=3) ~ 3,
                                       (vaccine3sym==1 & vaccdosesym>=3) ~ 3,
                                       (vaccine3sym==1 & vaccdose<=0 & vaccdosesym<=0) ~ -1,#<=100
                                       vaccine3sym==2 ~ 0),
         vax_status_cat = case_when(vax_status_number == 0 ~ "Not vaccinated",
                                    vax_status_number == 1 & vaccine_to_swab_time >= 1 & vaccine_to_swab_time < 14 ~ "Not vaccinated",
                                    vax_status_number == 1 & vaccine_to_swab_time >= 14 ~ "One does",
                                    vax_status_number %in% c(2,3) & time_after_2dose < 14 ~ "One does",
                                    vax_status_number %in% c(2,3) & time_after_2dose >= 14 ~ "Two does",
                                    vax_status_number == -1 ~ "Unknown does",
                                    is.na(vax_status_number) == TRUE ~ "NA"),
         vax_status_noDate = case_when(vax_status_number == 0 ~ "Not vaccinated",
                                       vax_status_number == 1 ~ "One does",
                                       vax_status_number %in% c(2,3) ~ "Two does",
                                       vax_status_number == -1 ~ "Unknown does",
                                       is.na(vax_status_number) == TRUE ~ "NA"),
         vax_status_noDate_v2 = case_when(vax_status_number == 0 ~ "Not vaccinated",
                                          vax_status_number == 1 ~ "One does",
                                          vax_status_number == 2 ~ "Two does",
                                          vax_status_number == 3 ~ "Three does",
                                          vax_status_number == -1 ~ "Unknown does",
                                          is.na(vax_status_number) == TRUE ~ "NA")) %>%
  #age<16 ~ 0)) %>%
  mutate(vax_wane = case_when(vax_status_cat == "Not vaccinated" ~ "Unvaccinated",
                              vax_status_cat == "One does" ~ "1 dose",
                              vax_status_cat == "Two does" & time_after_2dose <= 90 &  time_after_2dose >= 14 ~ "2 dose < 3 months",
                              vax_status_cat == "Two does" & time_after_2dose > 90 &  time_after_2dose <= 180 ~ "2 dose 3-6 months",
                              vax_status_cat == "Two does" & time_after_2dose > 180 &  time_after_2dose <= 365.25 ~ "2 dose > 6 months"))


# fix all of the broken dates - dates will need manually updating each time the
# data changes

# choose a sensible most recent date for vaccination - I am choosing the date in
# the raw data file

max_vaccination_date <- as.Date(max(dfRes$d_comb, na.rm = TRUE))

# I cannot get ifelse to keep things in date format, hence repeatedly having to
#  convert columns back to date type

# assume that vaccination during Jan 2020 should really by Jan 2021.
# change any dates after 2020-12-01 to have year 2021, then set any that occur
#  after max_vaccination_date to NA.

# do this for vaccine_first and vaccine_second, with additional restriction that
#  vaccination_second must come after vaccination_first
#
# dfRes <- dfRes %>%
#   mutate(vaccine_first = ifelse( vaccinated == "yes" & vaccine_first < as.Date("2020-12-01"),
#                                  format(vaccine_first, "2021-%m-%d"),
#                                  format(vaccine_first, "%Y-%m-%d")),
#          vaccine_first = as.Date(vaccine_first),
#          vaccine_first = if_else( vaccinated == "yes" & vaccine_first < max_vaccination_date,
#                                   format(vaccine_first, "%Y-%m-%d"),
#                                   NA_character_),
#          vaccine_first = as.Date(vaccine_first),
#
#   ) %>%
#   mutate(vaccine_second = ifelse( vaccinated == "yes" & vaccine_second < as.Date("2020-12-01"),
#                                   format(vaccine_second, "2021-%m-%d"),
#                                   format(vaccine_second, "%Y-%m-%d")),
#          vaccine_second = as.Date(vaccine_second),
#          vaccine_second = ifelse( vaccinated == "yes" & vaccine_second < max_vaccination_date,
#                                   format(vaccine_second, "%Y-%m-%d"),
#                                   NA_character_),
#          vaccine_second = as.Date(vaccine_second),
#          vaccine_second = ifelse( format(vaccine_second, "%Y-%m-%d") > format(vaccine_first, "%Y-%m-%d"),
#                                   format(vaccine_second, "%Y-%m-%d"),
#                                   NA_character_),
#          vaccine_second = as.Date(vaccine_second)
#   )
#
# dfRes <- dfRes %>%
#   mutate(vaccinefirstsym = ifelse( vaccinated == "yes" & vaccinefirstsym < as.Date("2020-12-01"),
#                                  format(vaccinefirstsym, "2021-%m-%d"),
#                                  format(vaccinefirstsym, "%Y-%m-%d")),
#          vaccinefirstsym = as.Date(vaccinefirstsym),
#          vaccinefirstsym = ifelse( vaccinated == "yes" & vaccinefirstsym < max_vaccination_date,
#                                   format(vaccinefirstsym, "%Y-%m-%d"),
#                                   NA_character_),
#          vaccinefirstsym = as.Date(vaccinefirstsym),
#
#   ) %>%
#   mutate(vaccinesecondsym = ifelse( vaccinated == "yes" & vaccinesecondsym < as.Date("2020-12-01"),
#                                    format(vaccinesecondsym, "2021-%m-%d"),
#                                    format(vaccinesecondsym, "%Y-%m-%d")),
#          vaccinesecondsym = as.Date(vaccinesecondsym),
#          vaccinesecondsym = if_else( vaccinated == "yes" & vaccinesecondsym < max_vaccination_date,
#                                     format(vaccinesecondsym, "%Y-%m-%d"),
#                                     NA_character_),
#          vaccinesecondsym = as.Date(vaccinesecondsym),
#
#   )
#
#
# # variable to indicate if the vaccine was before or after swab
# dfRes <- dfRes %>%
#   mutate( vaccine_first_before_swab = case_when(
#     vaccine_first < d_comb ~ "yes",
#     vaccine_first > d_comb ~ "no",
#     vaccine_first == d_comb ~ "same day")
#   )
#
# # create a new varibale called vacc_status:
# # if vaccinated == "yes" and vaccine_first_before_swab == "yes", --> vacc_status == "yes"
# # if vaccinated == "no" or vaccinated == "yes" and vaccine_first_before_swab == "no", --> vacc_status == "no"
# dfRes <- dfRes %>%
#   mutate(vacc_status = ifelse(vaccinated=="yes" & vaccine_first_before_swab=="yes", "Yes",
#                               ifelse(vaccinated=="no"|(vaccinated=="yes"&vaccine_first_before_swab=="no"), "No", NA)))

# ##################################################################
# # Type of house
# ##################################################################
# dfRes <- mutate(dfRes, dwelling = ifelse(dwelltyp==1,"House/Bungalow",
#                                          ifelse(dwelltyp==2,"Flat/Apartment",
#                                                 ifelse(dwelltyp==3,"Hostel",
#                                                        ifelse(dwelltyp==4,"caravan/mobile home",
#                                                               ifelse(dwelltyp==5,"Sheltered House",
#                                                                      ifelse(dwelltyp==6, "Homeless",
#                                                                             ifelse(dwelltyp==7, "Other",
#                                                                                    ifelse(dwelltyp==8, "Prefer not to say", "NA")))))))))
# ##################################################################
# # Household member work at school
# ##################################################################
# dfRes <- mutate(dfRes, hhmem_work_preschool = ifelse(workstudypers1_1==1, "yes",
#                                                      ifelse(workstudypers1_1==0,"no",
#                                                             "NA")))
# dfRes <- mutate(dfRes, hhmem_work_prischool = ifelse(workstudypers1_2==1, "yes",
#                                                      ifelse(workstudypers1_2==0,"no",
#                                                             NA)))
# dfRes <- mutate(dfRes, hhmem_work_secschool = ifelse(workstudypers1_3==1, "yes",
#                                                      ifelse(workstudypers1_3==0,"no",
#                                                             NA)))
# dfRes <- mutate(dfRes, hhmem_work_coluni = ifelse(workstudypers1_4==1, "yes",
#                                                   ifelse(workstudypers1_4==0,"no",
#                                                          NA)))
# dfRes <- mutate(dfRes, hhmem_work_preprisecschoolcoluni = ifelse(workstudypers1_1==0, "yes",
#                                                                  ifelse(workstudypers1_1==1,"no",
#                                                                         NA)))
# dfRes <- mutate(dfRes, workstudypers1_char = ifelse(workstudypers1_1==1, "Pre_school",
#                                                     ifelse(workstudypers1_2 ==1, "Primary school",
#                                                            ifelse(workstudypers1_3==1, "Secondary school",
#                                                                   ifelse(workstudypers1_4==1, "College/university",
#                                                                          ifelse(workstudypers1_5==1, "None of these", "NA"))))))
#
# #####################################################################
# # Student questions
# #####################################################################
#
# dfRes <- mutate(dfRes, hhmem_attend_schooluni = ifelse(workstudypers2==1, "Yes",
#                                                        ifelse(workstudypers2==2,"No",
#                                                               ifelse(workstudypers2==3, "Don't know", "NA"))))
#
# dfRes <- mutate(dfRes, christmas_study_return = ifelse(studstay==1, "Yes",
#                                                        ifelse(studstay==2,"No",
#                                                               ifelse(studstay==3, "Prefer not to say", "NA"))))
#
#
# dfRes <- mutate(dfRes, study_where = ifelse(edtype==1, "Further Education/Vocational Training College",
#                                             ifelse(edtype==2,"Uni undergrad",
#                                                    ifelse(edtype==3,"Uni postgrad",
#                                                           ifelse(edtype==4,"Another type of institution",
#                                                                  ifelse(edtype==5, "Don't know", "NA"))))))
#
# dfRes <- mutate(dfRes, student_curr_living = ifelse(campus2==1, "Uni halls",
#                                                     ifelse(campus2==2,"Private student halls",
#                                                            ifelse(campus2 ==3, "Private rented house/flat with other students",
#                                                                   ifelse(campus2==4, "private rented house/flat with NO other students",
#                                                                          ifelse(campus2==5, "Own home which you own",
#                                                                                 ifelse(campus2==6, "Parent/Guardians home",
#                                                                                        ifelse(campus2==7, "Other", "NA"))))))))
#
# ################################################################################
# # Travel abroad questions
# ################################################################################
dfRes <- mutate(dfRes, abroad_2weeksprior = ifelse(abroad==1, "Yes",
                                                   ifelse(abroad==2,"No",
                                                          "NA")))
dfRes <- mutate(dfRes, abroad_country1 = ifelse(abroad==1 & countryvisit1!=-66, countryvisit1, NA))
dfRes <- mutate(dfRes, abroad_country2 = ifelse(abroad==1 & countryvisit2!=-66, countryvisit2, NA))

# ####################################################################
# # Pre-existing health conditions - Health
# ####################################################################
# dfRes <- dfRes %>%
#   mutate(healtha_char = ifelse(healtha_1==1,"Organ transplant recipient",
#                                ifelse(healtha_2==1,"Diabetes",
#                                       ifelse(healtha_3==1,"Heart disease/heart problems",
#                                              ifelse(healtha_4==1,"Hypertension",
#                                                     ifelse(healtha_6==1,"Stroke",
#                                                            ifelse(healtha_7==1,"Kindey disease",
#                                                                   ifelse(healtha_8==1,"Liver disease",
#                                                                          ifelse(healtha_9==1,"Anaemia",
#                                                                                 ifelse(healtha_10==1,"Asthma",
#                                                                                        ifelse(healtha_11,"Other lung condition",
#                                                                                               ifelse(healtha_12==1,"Cancer",
#                                                                                                      ifelse(healtha_13==1,"Condition affecting the brain and nerves",
#                                                                                                             ifelse(healtha_14==1,"A weakened immune system/reduced ability to deal with infections",
#                                                                                                                    ifelse(healtha_15==1,"Depression",
#                                                                                                                           ifelse(healtha_16==1,"Anxiety",
#                                                                                                                                  ifelse(healtha_17==1,"Psychiatric disorder",
#                                                                                                                                         ifelse(healtha_18==1,"None of these", "NA")))))))))))))))))
#   )
#
# ####################################################################
# # Pre-existing health conditions - Shield 1
# ####################################################################
# dfRes <- dfRes %>%
#   mutate(shield1_char = case_when(
#     shield1 == 1 ~ "Yes",
#     shield1 == 2 ~ "No",
#     shield1 == 3 ~ "Don't know",
#     is.na(shield1) ~ "NA")
#     )
#
# ####################################################################
# # Pre-existing health conditions - Shield 2
# ####################################################################
# dfRes <- dfRes %>%
#   mutate(shield2_char = case_when(
#     shield2 == 1 ~ "Yes",
#     shield2 == 2 ~ "No",
#     is.na(shield2) ~ "NA")
#     )
#
# ####################################################################
# # Pre-existing health conditions - Transportation
# ####################################################################
# dfRes <- dfRes %>%
#   mutate(
#     transp_char = ifelse(transp_single==1, "transp_single",
#                          ifelse(transp_car==1, "transp_car",
#                                 ifelse(transp_mass==1,"transp_mass", "NA")))
#     )
#
# ####################################################################
# # Pre-existing health conditions - Face cover
# ####################################################################
# dfRes <- dfRes %>%
#   mutate(face_cov_char = case_when(
#     face_cov %in% c(-77, -02) ~ "NA",
#     face_cov == 1 ~ "No",
#     face_cov == 2 ~ "Yes, at work/school only",
#     face_cov == 3 ~ "Yes, in other situations only",
#     face_cov == 4 ~ "Yes, usually both at work/school and in other situations",
#     face_cov == 5 ~ "My/Their face is already covered for other reasons")
#     )
#
# ####################################################################
# # Pre-existing health conditions - Indoor mask
# ####################################################################
# dfRes <- dfRes %>%
#   mutate(indmask_char = case_when(
#     indmask %in% c(-77, -91, -92) ~ "NA",
#     indmask == 1 ~ "All of the time",
#     indmask == 2 ~ "Some of the time",
#     indmask == 3 ~ "Hardly ever",
#     indmask == 4 ~ "Never",
#     indmask == 5 ~ "Don't know")
#     )
#
# ####################################################################
# # Pre-existing health conditions - Outdoor mask
# ####################################################################
# dfRes <- dfRes %>%
#   mutate(outmask_char = case_when(
#     outmask %in% c(-77, -91, -92) ~ "NA",
#     outmask == 1 ~ "All of the time",
#     outmask == 2 ~ "Some of the time",
#     outmask == 3 ~ "Hardly ever",
#     outmask == 4 ~ "Never",
#     outmask == 5 ~ "Don't know")
#     )
#
# ####################################################################
# # Pre-existing health conditions - Bubble2 (Is your household currently in a bubble of any kind)
# ####################################################################
# dfRes <- dfRes %>%
#   mutate(bubble2_char = case_when(
#     bubble2 %in% c(-77,-92,0) ~ "NA",
#     bubble2 == 1 ~ "Yes",
#     bubble2 == 2 ~ "No",
#     bubble2 == 3 ~ "Don't know")
#     )
# ####################################################################
# # Pre-existing health conditions - Bubble num
# ####################################################################
# dfRes <- dfRes %>%
#   mutate(bubblenum_char = case_when(
#     bubblenum %in% c(-77,-66,0) ~ "NA",
#     bubblenum %in% c(1,2,3) ~ "1-3",
#     bubblenum %in% c(4,5,6) ~ "4-5",
#     bubblenum %in% c(7,8.9) ~ "7-9",
#     bubblenum %in% c(10:50) ~ "9+")
#     )
#
# ####################################################################
# # Pre-existing health conditions - convert variables to factor
# ####################################################################
# dfRes <- dfRes %>%
#   mutate(healtha_char = as.factor(healtha_char),
#          face_cov_char = as.factor(face_cov_char),
#          indmask_char = as.factor(indmask_char),
#          outmask_char = as.factor(outmask_char),
#          bubble2_char = as.factor(bubble2_char),
#          bubblenum_char = as.factor(bubblenum_char),
#          shield1_char = as.factor(shield1_char),
#          shield2_char = as.factor(shield2_char),
#          transp_char = as.factor(transp_char))
#
# #####################################################################
# # shedding data
# #####################################################################
#
# # shedding 1
# dfRes <- dfRes %>%
# mutate(collected_date_shed1 = as.Date(collected_date_time_shed1, format = "%m/%d/%Y")) %>%
#   dplyr::select(!(collected_date_time_shed1)) %>%
#   mutate(result_shed1 = ifelse(result_shed1 == "Void", NA, result_shed1),
#          estbinres_shed1 = ifelse(is.na(result_shed1), NA,
#                                   ifelse(result_shed1 != "Detected", 0,
#                                          ifelse(n_gene_ct_value_shed1 > 0 & e_gene_ct_value_shed1 > 0, 1,
#                                                 ifelse(n_gene_ct_value_shed1 < 37 & n_gene_ct_value_shed1 > 0, 1, 0)))))
#
# # shedding 2
# dfRes <- dfRes %>%
#   mutate(collected_date_shed2 = as.Date(collected_date_time_shed2, format = "%m/%d/%Y")) %>%
#   dplyr::select(!(collected_date_time_shed2)) %>%
#   mutate(result_shed2 = ifelse(result_shed2 == "Void", NA, result_shed2),
#          estbinres_shed2 = ifelse(is.na(result_shed2), NA,
#                                   ifelse(result_shed2 != "Detected", 0,
#                                          ifelse(n_gene_ct_value_shed2 > 0 & e_gene_ct_value_shed2 > 0, 1,
#                                                 ifelse(n_gene_ct_value_shed2 < 37 & n_gene_ct_value_shed2 > 0, 1, 0)))))

#### prior infection ####
# dfRes <- dfRes %>%
#   mutate(prior_inf = ifelse(agprev1 == 1,
#                             ifelse((agprevtype_1==1 & agprevtype_2 == 1) | agprevnum == 2, case_when((agprev3 == 1)| ((agprev3 %in% c(2,3,4))&agprev4==1)~ "Prior infection",
#                                                                                                       (agprev3 %in% c(2,3,4))&agprev4==2 ~ "No prior infection"), "Unkown"), "Unkown")
#   ) %>%
#   mutate(prior_inf = ifelse(is.na(prior_inf), "Unkown", prior_inf))


dfRes <- dfRes %>%
  mutate(covidb_new = as.Date(covidb_new, "%m/%d/%Y"),
         infection_to_swab =  d_comb - covidb_new) %>%
  mutate(prior_inf = case_when(covida == 1 & (infection_to_swab >28) ~ "Previous infection greater than 28 days",
                               covida == 1 & (infection_to_swab <= 28) ~ "Previous infection within 28 days",
                               covida == 1 & is.na(infection_to_swab) ~ "Previous infection with unknow time",
                               covida %in% c(2,3) ~ "Suspected previous infection",
                               covida == 4 ~ "No previous infection")
  ) %>%
  mutate(prior_inf = ifelse(is.na(prior_inf), "Unkown", prior_inf))

#####################################################################
# Living in urban area: yes or no
#####################################################################
bridge <- fread("E:/Group/report/round16/Parameters/urban_vs_rural.csv", data.table=FALSE)

bridge$postcode <- gsub(" ", "", bridge$postcode)
dfRes$postcode <- gsub(" ", "", dfRes$postcode)

# table(!df_round$postcode%in%bridge$postcode)

mapping <- bridge$URBAN
names(mapping) <- bridge$postcode

dfRes$urban <- mapping[dfRes$postcode]
dfRes <- dfRes %>%
  mutate(urban_char = case_when(urban == 1 ~ "Yes",
                                urban == 0 ~ "No",
                                is.na(urban) == TRUE ~ "urban - unknown"))

#####################################################################
# Time delta (of last positive PCR to current REACT round swab)
#####################################################################
# edit : addendum 8 Feb 22, David Tang

dfRes <- dfRes %>%
  mutate(d_swab = as.Date(d_swab, "%Y-%m-%d")) %>%
  mutate(agprev5_new = as.Date(agprev5_new, "%m/%d/%Y")) %>%
  mutate(timedelta_prevag5_swab = as.integer(difftime(d_swab, agprev5_new, units = "days")))

#####################################################################
# Re-level factors
#####################################################################
dfRes <- within(dfRes, gender_char <- relevel(gender_char,ref="Male"))
dfRes <- within(dfRes, age_group_char <- relevel(age_group_char,ref="35-44"))
dfRes <- within(dfRes, age_group_char2 <- relevel(age_group_char2,ref="35-44"))
dfRes <- within(dfRes, region <- relevel(region,ref= "South East"))
dfRes <- within(dfRes, work_new_alt <- relevel(work_new_alt,ref="Other worker"))
dfRes <- within(dfRes, work_type1_char <- relevel(work_type1_char,ref="Other essential/key workers"))
dfRes <- within(dfRes, work_type3_char <- relevel(work_type3_char,ref="Other essential/key workers"))
dfRes <- within(dfRes, work_new_alt_ext <- relevel(work_new_alt_ext,ref="Other worker - not public facing"))
dfRes <- within(dfRes, ethnic_new_char <- relevel(ethnic_new_char, ref = "White"))
dfRes <- within(dfRes, ethnic19_char <- relevel(ethnic19_char, ref = "English/Welsh/Scottish/Northern Irish/British"))
# dfRes <- within(dfRes, smokenow_char <- relevel(smokenow_char,ref="No"))
dfRes <- within(dfRes, hh_size_cat <- relevel(hh_size_cat, ref="1"))
dfRes <- within(dfRes, hosp_char <- relevel(hosp_char, ref="No"))
dfRes <- within(dfRes, covidcon_char <- relevel(covidcon_char, ref="No"))
dfRes <- within(dfRes, age_group_tert <- relevel(age_group_tert, ref="Middle"))
dfRes <- within(dfRes, hh_size_tert <- relevel(hh_size_tert, ref="1-2"))
#dfRes <- within(dfRes, school <- relevel(school,ref="No"))
# dfRes <- within(dfRes, leave1_char <- relevel(leave1_char,ref="No"))
dfRes <- within(dfRes, sympt_cat <- relevel(sympt_cat, ref="No symptoms"))

# dfRes <- within(dfRes, dwelling  <- relevel(as.factor(dwelling),ref="Prefer not to say"))
# dfRes <- within(dfRes, workstudypers1_char <- relevel(as.factor(workstudypers1_char),ref="None of these"))
# dfRes <- within(dfRes, hhmem_attend_schooluni <- relevel(as.factor(hhmem_attend_schooluni),ref="No"))
# dfRes <- within(dfRes, christmas_study_return <- relevel(as.factor(christmas_study_return),ref="No"))
# dfRes <- within(dfRes, study_where <- relevel(as.factor(study_where),ref="Another type of institution"))
# dfRes <- within(dfRes, student_curr_living <- relevel(as.factor(student_curr_living),ref="Other"))
# dfRes <- within(dfRes, abroad_2weeksprior <- relevel(as.factor(abroad_2weeksprior),ref="No"))
# dfRes <- within(dfRes, healtha_char <- relevel(healtha_char,ref="None of these"))
# dfRes <- within(dfRes, face_cov_char <- relevel(face_cov_char,ref="No"))
# dfRes <- within(dfRes, indmask_char <- relevel(indmask_char,ref="Never"))
# dfRes <- within(dfRes, outmask_char <- relevel(outmask_char,ref="Never"))
# dfRes <- within(dfRes, bubble2_char <- relevel(bubble2_char,ref="No"))
# dfRes <- within(dfRes, bubblenum_char <- relevel(bubblenum_char,ref="1-3"))
# dfRes <- within(dfRes, shield1_char <- relevel(shield1_char,ref="Don't know"))
# dfRes <- within(dfRes, shield2_char <- relevel(shield2_char,ref="No"))
# dfRes <- within(dfRes, transp_char <- relevel(transp_char,ref="transp_single"))


##################################################################
# save out the files: all data and interim sets
##################################################################

# final_date_in_a <- as.Date("2021-10-29") # Edit when we create rep10a and rep10b
#
# dfRes <- dfRes %>%
#   mutate(interim_label = if_else(is.na(d_comb), date_last_modified, d_comb),
#          interim_label = if_else(interim_label <= final_date_in_a, "a", "b"))
#
# dfRes_a <- dfRes %>%
#   filter(interim_label == "a")
#
# dfRes_b <- dfRes %>%
#   filter(interim_label == "b")

# FALSE positive corrected below
dfRes[dfRes$u_passcode%in%"KKDQ79YX",]$estbinres <- 0
# dfRes[dfRes$u_passcode%in%"KKDQ79YX",]$estbinres
# Save data wrangled object
#saveRDS(dfRes, file="./saved_objects/rep15.rds")
saveRDS(dfRes, file=paste0(dirOutput, "/", fnOutput))
# saveRDS(dfRes_a, file="./saved_objects/rep15a.rds")
# saveRDS(dfRes_b, file="./saved_objects/rep15b.rds")


