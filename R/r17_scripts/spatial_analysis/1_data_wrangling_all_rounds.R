library(dplyr)
library(sf)

setwd("E:/Group/react1_spatial_analysis")

# load postcode data for joining
postcode_locs <- read.csv("../aux_data/ukpostcodes.csv")
postcode_locs$postcode=gsub(" ", "", postcode_locs$postcode)

# load wrangled data
round1_res <- readRDS("../saved_objects/rep1.rds")
round2_res <- readRDS("../saved_objects/rep2.rds")
round3_res <- readRDS("../saved_objects/rep3.rds")
round4_res <- readRDS("../saved_objects/rep4.rds")
round5_res <- readRDS("../saved_objects/rep5.rds")
round6_res <- readRDS("../saved_objects/rep6.rds")
round6a_res <- readRDS("../saved_objects/rep6_interim.rds")
round6b_res <- readRDS("../saved_objects/rep6_halfTerm.rds")
round7_res <- readRDS("../saved_objects/rep7.rds")
round7a_res <- readRDS("../saved_objects/rep7a.rds")
round7b_res <- readRDS("../saved_objects/rep7b.rds")
round8_res <- readRDS("../saved_objects/rep8.rds")
round8a_res <- readRDS("../saved_objects/rep8a.rds")
round8b_res <- readRDS("../saved_objects/rep8b.rds")
round9_res <- readRDS("../saved_objects/rep9.rds")
round9a_res <- readRDS("../saved_objects/rep9a.rds")
round9b_res <- readRDS("../saved_objects/rep9b.rds")
round10_res <- readRDS("../saved_objects/rep10.rds")
round11_res <- readRDS("../saved_objects/rep11.rds")
round12_res <- readRDS("../saved_objects/rep12.rds")
round13_res <- readRDS("../saved_objects/rep13.rds")
round13a_res <- readRDS("../saved_objects/rep13a.rds")
round13b_res <- readRDS("../saved_objects/rep13b.rds")
round14_res <- readRDS("../saved_objects/rep14.rds")
round15_res <- readRDS("../saved_objects/rep15.rds")
round16_res <- readRDS("../saved_objects/rep16.rds")
round17_res <- readRDS("../saved_objects/rep17.rds")

# choose relevant columns
names_to_keep <- c(
  "id",
  "estbinres",
  "region",
  "postcode",
  "lacode",
  "d_comb",
  "ct1",
  "ct2"
)

for (round_id in c(1:6, "6a", "6b", 
                   7, "7a", "7b", 
                   8, "8a", "8b", 
                   9, "9a", "9b",
                   10:13, "13a", "13b", 14:17)){
  print(round_id)
  
  # Re-formatting postcode for consistency across rounds
  mydata=eval(parse(text=paste0("round", round_id, "_res")))
  mydata$postcode=gsub(" ", "", mydata$postcode)
  
  # Creating a column 'round'
  mydata <- mydata %>%
    dplyr::select(all_of(names_to_keep)) %>%
    mutate(round = rep(paste0("round_", round_id), length = nrow(mydata)))
  
  assign(paste0("r", round_id), mydata)
}

# # tidy up data for each round, creating a column 'round'
# r1 <- round1_res %>%
#   dplyr::select(all_of(names_to_keep)) %>%
#   mutate(round = rep("round_1", length = nrow(round1_res)))
# 
# 
# r2 <- round2_res %>%
#   dplyr::select(all_of(names_to_keep)) %>%
#   mutate(round = rep("round_2", length = nrow(round2_res)))
# 
# r3 <- round3_res %>%
#   dplyr::select(all_of(names_to_keep)) %>%
#   mutate(round = rep("round_3", length = nrow(round3_res)))
# 
# r4 <- round4_res %>%
#   dplyr::select(all_of(names_to_keep)) %>%
#   mutate(round = rep("round_4", length = nrow(round4_res)))
# 
# r5 <- round5_res %>%
#   dplyr::select(all_of(names_to_keep)) %>%
#   mutate(round = rep("round_5", length = nrow(round5_res)))
# 
# r6 <- round6_res %>%
#   dplyr::select(all_of(names_to_keep)) %>%
#   mutate(round = rep("round_6", length = nrow(round6_res)))
# 
# r6a <- round6a_res %>%
#   dplyr::select(all_of(names_to_keep)) %>%
#   mutate(round = rep("round_6a", length = nrow(round6a_res)))
# 
# r6b <- round6b_res %>%
#   dplyr::select(all_of(names_to_keep)) %>%
#   mutate(round = rep("round_6b", length = nrow(round6b_res)))
# 
# r7 <- round7_res %>%
#   dplyr::select(all_of(names_to_keep)) %>%
#   mutate(round = rep("round_7", length = nrow(round7_res)))
# 
# r7a <- round7a_res %>%
#   dplyr::select(all_of(names_to_keep)) %>%
#   mutate(round = rep("round_7a", length = nrow(round7a_res)))
# 
# r7b <- round7b_res %>%
#   dplyr::select(all_of(names_to_keep)) %>%
#   mutate(round = rep("round_7b", length = nrow(round7b_res)))
# 
# r8 <- round8_res %>%
#   dplyr::select(all_of(names_to_keep)) %>%
#   mutate(round = rep("round_8", length = nrow(round8_res)))
# 
# r8a <- round8a_res %>%
#   dplyr::select(all_of(names_to_keep)) %>%
#   mutate(round = rep("round_8a", length = nrow(round8a_res)))
# 
# r8b <- round8b_res %>%
#   dplyr::select(all_of(names_to_keep)) %>%
#   mutate(round = rep("round_8b", length = nrow(round8b_res)))
# 
# r9 <- round9_res %>%
#   dplyr::select(all_of(names_to_keep)) %>%
#   mutate(round = rep("round_9", length = nrow(round9_res)))
# 
# r9a <- round9a_res %>%
#   dplyr::select(all_of(names_to_keep)) %>%
#   mutate(round = rep("round_9a", length = nrow(round9a_res)))
# 
# r9b <- round9b_res %>%
#   dplyr::select(all_of(names_to_keep)) %>%
#   mutate(round = rep("round_9b", length = nrow(round9b_res)))
# 
# r10 <- round10_res %>%
#   dplyr::select(all_of(names_to_keep)) %>%
#   mutate(round = rep("round_10", length = nrow(round10_res)))
# 
# r11 <- round11_res %>%
#   dplyr::select(all_of(names_to_keep)) %>%
#   mutate(round = rep("round_11", length = nrow(round11_res)))
# 
# r12 <- round12_res %>%
#   dplyr::select(all_of(names_to_keep)) %>%
#   mutate(round = rep("round_12", length = nrow(round12_res)))
# 
# r13 <- round13_res %>%
#   dplyr::select(all_of(names_to_keep)) %>%
#   mutate(round = rep("round_13", length = nrow(round13_res)))
# 
# r13a <- round13a_res %>%
#   dplyr::select(all_of(names_to_keep)) %>%
#   mutate(round = rep("round_13a", length = nrow(round13a_res)))
# 
# r13b <- round13b_res %>%
#   dplyr::select(all_of(names_to_keep)) %>%
#   mutate(round = rep("round_13b", length = nrow(round13b_res)))
# 
# r14 <- round14_res %>%
#   dplyr::select(all_of(names_to_keep)) %>%
#   mutate(round = rep("round_14", length = nrow(round14_res)))
# 
# r15 <- round15_res %>%
#   dplyr::select(all_of(names_to_keep)) %>%
#   mutate(round = rep("round_15", length = nrow(round15_res)))
# 
# r16 <- round16_res %>%
#   dplyr::select(all_of(names_to_keep)) %>%
#   mutate(round = rep("round_16", length = nrow(round16_res)))
# 
# r17 <- round17_res %>%
#   dplyr::select(all_of(names_to_keep)) %>%
#   mutate(round = rep("round_17", length = nrow(round17_res)))

# bind together the datasets and rename id column to 'old id'
all_rounds <- bind_rows(
  r1, r2, r3, r4, r5, r6, r6a, r6b, r7, r7a, r7b, r8, r8a,
  r8b, r9, r9a, r9b, r10, r11, r12, r13, r13a, r13b, r14, r15, r16, r17
) %>%
  rename(old_id = id)

# join the postcode data
all_rounds <- left_join(all_rounds,
  postcode_locs %>% dplyr::select(-c(id)),
  by = "postcode"
)

#################################################################
# get local authority district names
#################################################################
lad_shp <- st_read("E:/Group/aux_data/Local_Authority_Districts__December_2019__Boundaries_UK_BGC-shp/Local_Authority_Districts__December_2019__Boundaries_UK_BGC.shp")

lads <- lad_shp %>%
  st_set_geometry(NULL) %>%
  dplyr::select(lad19cd, lad19nm) %>%
  rename(
    lacode = lad19cd,
    laname = lad19nm
  )

# join the ltla names data
all_rounds <- left_join(all_rounds,
  lads,
  by = "lacode"
)

#################################################################
# remove and track missing data rows
#################################################################

# find the NAs in estbinres, latitude, longitude and remove from the data
na_rows <- which(is.na(all_rounds$latitude) | is.na(all_rounds$longitude) | is.na(all_rounds$estbinres))

# remove a dodgy result that has latitude 100
incorrect_rows <- which(all_rounds$latitude >= 99)

# remove the NAs and incorrect results, but save them as a different data frame
rm_rows <- c(na_rows, incorrect_rows)

dat_removed <- all_rounds[rm_rows, ]
all_rounds <- all_rounds[-rm_rows, ]

# # create a new id columns to ensure unique id
# all_rounds$id <- seq(1, nrow(all_rounds), by = 1)

# ensure region column is character not factor
all_rounds$region <- as.character(all_rounds$region)


#################################################################
# include jittered lat and long
#################################################################

all_rounds <- all_rounds %>%
  mutate(
    jittered_latitude = base::jitter(latitude, amount = 0.05),
    jittered_longitude = base::jitter(longitude, amount = 0.05)
  )

#################################################################
# save the dataframe of 'keep' results and 'removed' results
#################################################################
# separate out dataframes for each round individually
# generate a new id column for each data subset

for (round_id in c(1:6, "6a", "6b", 
                   7, "7a", "7b", 
                   8, "8a", "8b", 
                   9, "9a", "9b",
                   10:13, "13a", "13b", 14:17)){
  print(round_id)
  
  mydata <- all_rounds %>%
    dplyr::filter(round == paste0("round_", round_id)) %>%
    mutate("id" = row_number())
  
  assign(paste0("round", round_id), mydata)
}

# round1 <- all_rounds %>%
#   dplyr::filter(round == "round_1") %>%
#   mutate("id" = row_number())
# 
# round2 <- all_rounds %>%
#   dplyr::filter(round == "round_2") %>%
#   mutate("id" = row_number())
# 
# round3 <- all_rounds %>%
#   dplyr::filter(round == "round_3") %>%
#   mutate("id" = row_number())
# 
# round4 <- all_rounds %>%
#   dplyr::filter(round == "round_4") %>%
#   mutate("id" = row_number())
# 
# round5 <- all_rounds %>%
#   dplyr::filter(round %in% c("round_5")) %>%
#   mutate("id" = row_number())
# 
# round6 <- all_rounds %>%
#   dplyr::filter(round %in% c("round_6")) %>%
#   mutate("id" = row_number())
# 
# round6a <- all_rounds %>%
#   dplyr::filter(round %in% c("round_6a")) %>%
#   mutate("id" = row_number())
# 
# round6b <- all_rounds %>%
#   dplyr::filter(round %in% c("round_6b")) %>%
#   mutate("id" = row_number())
# 
# round7 <- all_rounds %>%
#   dplyr::filter(round %in% c("round_7")) %>%
#   mutate("id" = row_number())
# 
# round7a <- all_rounds %>%
#   dplyr::filter(round %in% c("round_7a")) %>%
#   mutate("id" = row_number())
# 
# round7b <- all_rounds %>%
#   dplyr::filter(round %in% c("round_7b")) %>%
#   mutate("id" = row_number())
# 
# round8 <- all_rounds %>%
#   dplyr::filter(round %in% c("round_8")) %>%
#   mutate("id" = row_number())
# 
# round8a <- all_rounds %>%
#   dplyr::filter(round %in% c("round_8a")) %>%
#   mutate("id" = row_number())
# 
# round8b <- all_rounds %>%
#   dplyr::filter(round %in% c("round_8b")) %>%
#   mutate("id" = row_number())
# 
# round9 <- all_rounds %>%
#   dplyr::filter(round %in% c("round_9")) %>%
#   mutate("id" = row_number())
# 
# round9a <- all_rounds %>%
#   dplyr::filter(round %in% c("round_9a")) %>%
#   mutate("id" = row_number())
# 
# round9b <- all_rounds %>%
#   dplyr::filter(round %in% c("round_9b")) %>%
#   mutate("id" = row_number())
# 
# round10 <- all_rounds %>%
#   dplyr::filter(round %in% c("round_10")) %>%
#   mutate("id" = row_number())
# 
# round11 <- all_rounds %>%
#   dplyr::filter(round %in% c("round_11")) %>%
#   mutate("id" = row_number())
# 
# round12 <- all_rounds %>%
#   dplyr::filter(round %in% c("round_12")) %>%
#   mutate("id" = row_number())
# 
# round13 <- all_rounds %>%
#   dplyr::filter(round %in% c("round_13")) %>%
#   mutate("id" = row_number())
# 
# round13a <- all_rounds %>%
#   dplyr::filter(round %in% c("round_13a")) %>%
#   mutate("id" = row_number())
# 
# round13b <- all_rounds %>%
#   dplyr::filter(round %in% c("round_13b")) %>%
#   mutate("id" = row_number())
# 
# round14 <- all_rounds %>%
#   dplyr::filter(round %in% c("round_14")) %>%
#   mutate("id" = row_number())
# 
# round15 <- all_rounds %>%
#   dplyr::filter(round %in% c("round_15")) %>%
#   mutate("id" = row_number())
# 
# round16 <- all_rounds %>%
#   dplyr::filter(round %in% c("round_16")) %>%
#   mutate("id" = row_number())
# 
# round17 <- all_rounds %>%
#   dplyr::filter(round %in% c("round_17")) %>%
#   mutate("id" = row_number())



# save the removed data
# saveRDS(removed_data, "output_files/clean_data/dat_removed.rds")


# # save single round data
# saveRDS(round1, "output_files/clean_data/round_1_data.rds")
# saveRDS(round2, "output_files/clean_data/round_2_data.rds")
# saveRDS(round3, "output_files/clean_data/round_3_data.rds")
# saveRDS(round4, "output_files/clean_data/round_4_data.rds")
# saveRDS(round5, "output_files/clean_data/round_5_data.rds")
# saveRDS(round6, "output_files/clean_data/round_6_data.rds")
# saveRDS(round6a, "output_files/clean_data/round_6a_data.rds")
# saveRDS(round6b, "output_files/clean_data/round_6b_data.rds")
# saveRDS(round7, "output_files/clean_data/round_7_data.rds")
# saveRDS(round7a, "output_files/clean_data/round_7a_data.rds")
# saveRDS(round7b, "output_files/clean_data/round_7b_data.rds")
# saveRDS(round8, "output_files/clean_data/round_8_data.rds")
# saveRDS(round8a, "output_files/clean_data/round_8a_data.rds")
# saveRDS(round8b, "output_files/clean_data/round_8b_data.rds")
# saveRDS(round9, "output_files/clean_data/round_9_data.rds")
# saveRDS(round9a, "output_files/clean_data/round_9a_data.rds")
# saveRDS(round9b, "output_files/clean_data/round_9b_data.rds")
# saveRDS(round10, "output_files/clean_data/round_10_data.rds")
# saveRDS(round11, "output_files/clean_data/round_11_data.rds")
# saveRDS(round12, "output_files/clean_data/round_12_data.rds")
# saveRDS(round13, "output_files/clean_data/round_13_data.rds")
# saveRDS(round13a, "output_files/clean_data/round_13a_data.rds")
# saveRDS(round13b, "output_files/clean_data/round_13b_data.rds")
# saveRDS(round14, "output_files/clean_data/round_14_data.rds")
# saveRDS(round15, "output_files/clean_data/round_15_data.rds")
# saveRDS(round16, "output_files/clean_data/round_16_data.rds")
saveRDS(round17, "output_files/clean_data/round_17_data.rds")


