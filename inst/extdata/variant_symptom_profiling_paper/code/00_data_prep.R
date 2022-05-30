

# Import data -------------------------------------------------------------

dfRes <- readRDS("E:/Group/react2_study5/saved_objects/react_1_saved_objects/react1_r1_to_r19.rds")

table(dfRes$round, dfRes$vax_status_number, exclude = "none")

table(dfRes$sympt_cat, exclude = "none")

### define lists of variable names
covid_yesnos <- paste0("symptnowaw_",1:26)
# covid_yesnos_1month <- dfRes %>% select(sympt_any1_1:sympt_any3_9) %>% colnames()
# covid_first <- dfRes %>% select(symptfirst_1:symptfirst_26) %>% colnames()


### clean up the data
dfRes <- dfRes %>% 
  mutate_at(all_of(covid_yesnos),binaryCleaner_1_0)  %>% 
  # mutate_at(all_of(covid_yesnos_1month),binaryCleaner_1_0)  %>% 
  # mutate_at(all_of(covid_first),binaryCleaner_1_0)  %>% 
  mutate_at(c("gender","age", "ethnic", "smokecig", 
              "empl", "educ","covidcon"),
            continuousCleaner)  %>% 
  mutate(age_group_pred = case_when(age <18 ~ 1,
                                    age <=54 ~ 2,
                                    age >54  ~ 3,
                                    TRUE ~NA_real_))


dfRes <- dfRes %>% 
  mutate_at(all_of(covid_yesnos),binaryCleaner_1_0)  %>% 
  # mutate_at(all_of(covid_yesnos_1month),binaryCleaner_1_0)  %>% 
  # mutate_at(all_of(covid_first),binaryCleaner_1_0)  %>% 
  mutate_at(c("gender","age", "ethnic", "smokecig", 
              "empl", "educ","covidcon"),
            continuousCleaner)  %>% 
  mutate(age_group_pred = case_when(age <18 ~ 1,
                                    age <=54 ~ 2,
                                    age >54  ~ 3,
                                    TRUE ~NA_real_),
         vax_status_number=case_when(vax_status_number <0 ~ NA_real_,
                                     T ~ vax_status_number)) %>% 
  mutate(vax_one_or_more_dose = case_when(vax_status_number >=1 ~ 1,
                                          TRUE ~0),
         vax_two_dose = case_when(vax_status_number >=2 ~ 1,
                                  TRUE ~0),
         vaccine_status_named = case_when(vax_status_number ==1 ~ "One dose",
                                          vax_status_number ==2 ~ "Two doses",
                                          vax_status_number ==3 ~ "Three doses",
                                          vax_status_number ==0 ~ "Unvaccinated",
                                          TRUE ~ "Vaccine status unknown/unclear")) %>% 
  mutate(smokever = case_when(smokenow == 1 ~ "Current smoker",
                              smokecig == 1 & smokenow == 2 ~ "Former smoker",
                              smokecig == 2 ~ "Never smoker",
                              smokecig == 3 | smokenow == 3 ~"Prefer not to say",
                              TRUE ~ NA_character_),
         variant_named=case_when(react_lineage=="B.1.1.7" ~ "Alpha",
                                 react_lineage=="B.1.617.2" ~ "Delta",
                                 react_lineage=="Wildtype" ~ "Wild type",
                                 react_lineage %in% c("BA.1","BA.1.1","BA.2",
                                                      "B.A.1","B.A.1.1","B.A.2") ~ "Omicron",
                                 !is.na(react_lineage) ~ "Other",
                                 estbinres==1 ~ "Not typed",
                                 TRUE ~ NA_character_),
         round=as.factor(round),
         all_participants = "All participants")


# dfRes$vax_status_number %>% table(exclude="none")
# 
# 
# table(dfRes$vax_status,dfRes$vax_status_number,exclude="none")
# 

### Add a variable for any one of the 26 reported symptoms
dfRes$symptomatic <- as.numeric(rowSums(dfRes[,covid_yesnos], na.rm=T) > 0)

### Add variable for one of four classical symptoms
dfRes$one_of_four <- as.numeric(rowSums(dfRes[, covid_yesnos[c(1:4)]], na.rm=T) > 0)

### get symptomatic pos
dfRes$symp_pos <- case_when(dfRes$symptomatic ==1 & dfRes$estbinres == 1 ~ 1,
                            T ~0)

## covid hospitality
dfRes <- dfRes %>% mutate(seekmed_cat =case_when(seekmed==1 ~ "Yes",
                                        seekmed==2 ~ "No",
                                        TRUE ~NA_character_),
                 kindmed_cat = case_when(kindmed_9==1 ~ "ICU",
                                         kindmed_8==1 ~ "Hospital admission",
                                         kindmed_10==1|
                                           kindmed_6==1| 
                                           kindmed_5==1 |
                                           kindmed_4==1 |
                                           kindmed_3==1~ "GP/Walk-in/A&E/consultation",
                                         kindmed_2==1 | kindmed_1==1 ~ "111 / pharmacy",
                                         seekmed==1 ~ "Other",
                                         TRUE ~ NA_character_),
                 hospitalised=case_when(kindmed_9==1 | kindmed_8==1~ "Yes",
                                        TRUE ~ "No")) %>% 
  mutate(variant_inferred = factor(case_when(round %in% c(2:7) ~ "Rounds 2-7 (Wild type)",
                                             round %in% c(8:10) ~ "Rounds 8-10 (Alpha)",
                                             round %in% c(12:15) ~ "Rounds 12-15 (Delta)",
                                             round %in% c(17,19) ~ "Rounds 17-19 (Omicron)",
                                             TRUE ~NA_character_), levels = c("Rounds 2-7 (Wild type)",
                                                                              "Rounds 8-10 (Alpha)",
                                                                              "Rounds 12-15 (Delta)",
                                                                              "Rounds 17-19 (Omicron)")),
         variant_inferred_detail = factor(case_when(round %in% c(2:7) ~ "Rounds 2-7 (Wild type)",
                                             round %in% c(8:10) ~ "Rounds 8-10 (Alpha)",
                                             round %in% c(12:15) ~ "Rounds 12-15 (Delta)",
                                             grepl("BA.1",react_lineage) & round %in% c(16,17,18,19) ~ "BA.1 (Omicron)",
                                             grepl("BA.2",react_lineage) & round %in% c(16,17,18,19) ~ "BA.2 (Omicron)",
                                             TRUE ~NA_character_), levels = c("Rounds 2-7 (Wild type)",
                                                                              "Rounds 8-10 (Alpha)",
                                                                              "Rounds 12-15 (Delta)",
                                                                              "BA.1 (Omicron)",
                                                                              "BA.2 (Omicron)")),
         ) %>% 
  ### Vaccination status
  mutate(doublevaxxed=case_when(vax_status_number==2 ~1,
                                TRUE ~0),
         boosted=case_when(vax_status_number==3 ~1,
                           TRUE ~0),
         doublevaxxed_or_boosted=case_when(vax_status_number>=2 ~1,
                                           TRUE ~0))




### Replace NA with 0 in symptoms
dfRes[covid_yesnos][is.na(dfRes[covid_yesnos])] <- 0


### Create new symptomatic 26 var
dfRes$symptomatic_26 = as.numeric(rowSums(dfRes[,covid_yesnos]==1,na.rm = T)>0)

#### FIrst symptoms ###

# add first symptoms to covnamelist
cov_name_list[sympnames_type_df$first_symptom_code[1:26]] <- sympnames_type_df$symptom[1:26]

# Sort out symptoms
dfRes <- dfRes %>% 
  mutate_at(all_of(sympnames_type_df$first_symptom_code[1:26]),binaryCleaner_1_0) 


### Replace NA with 0 in symptoms
dfRes[sympnames_type_df$first_symptom_code[1:26]][is.na(dfRes[sympnames_type_df$first_symptom_code[1:26]])] <- 0

# Tidy up time since sympton onset
dfRes <- dfRes %>% mutate(days_since_symptom_onset = case_when(symptst <0 ~ NA_real_,
                                                               symptst >11 ~ NA_real_,
                                                               T ~ symptst))
cov_name_list$days_since_symptom_onset="Days since onset of first symptom"




# Data filtering ----------------------------------------------------------

# make exclusions and calculate numbers after each exclusion
n_all = nrow(dfRes %>% filter(round %in%c(2:10,13:19)))

dfRes <- dfRes %>% filter(age >=18, age <=200)
n_exclude_age_18_64 = nrow(dfRes%>% filter(round %in%c(2:10,13:19)))

dfRes=dfRes %>% filter(!is.na(estbinres))
n_exclude_missing_res = nrow(dfRes%>% filter(round %in%c(2:10,13:19)))

dfRes <- dfRes %>% filter(!is.na(gender), !is.na(age)) 
n_exclude_missing_age_sex = nrow(dfRes%>% filter(round %in%c(2:10,13:19)))

dfRes <- dfRes %>% filter(sympt_cat != "NA" & !is.na(sympt_cat)) 
n_exclude_unknown_symptoms = nrow(dfRes%>% filter(round %in%c(2:10,13:19)))

# write file summarising exclusions
fileConn <- file(paste0(outpath,"exclusions_summary.txt"))
writeLines(c(paste0("Starting population was ", n_all),
             paste0(n_all-n_exclude_age_18_64, 
                    " were excluded because <18"),
             paste0("A further ",n_exclude_age_18_64-n_exclude_missing_res, " were excluded because of missing test results"),
             paste0("A further ",n_exclude_missing_res - n_exclude_missing_age_sex, 
                    " were excluded because of missing age or sex data"),
             paste0("A further ",n_exclude_missing_age_sex - n_exclude_unknown_symptoms, 
                    " were excluded because of unknown symptom data"),
             paste0("The final study population, after exclusions, was ", nrow(dfRes%>% filter(round %in%c(2:10,13:19))))
             ), fileConn)
close(fileConn)



# Test train split --------------------------------------------------------

### Create test-train split
set.seed(123)
ttsplit <- caret::createDataPartition( y = dfRes[dfRes$symptomatic == 1,]$estbinres,times = 1,p = 0.7,list = F)

### add variable for ttsplot
dfRes$test_train_split=NA_character_
dfRes[dfRes$symptomatic == 1,]$test_train_split <- case_when(as.numeric(rownames(dfRes[dfRes$symptomatic == 1,])) %in% ttsplit ~ "train",
                                                             TRUE~ "test")
# dfRes[dfRes$symptomatic == 0,]$test_train_split = NA_char
table(dfRes$test_train_split, dfRes$symptomatic, exclude = "none")


# saveRDS(object = dfRes, file = "E:/Group/react2_study5/saved_objects/react_1_saved_objects/delta_symptom_prediction_data.rds")

