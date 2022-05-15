rm(list = ls())
path_to_scripts <- "E:/Group/report/round19/Scripts/"
setwd(path_to_scripts)

library(dplyr)
library(geosphere)
library(tidyr)
library(ggplot2)
library(sf)
library(snakecase)
library(readr)
library(janitor)
library(prevalence)
library(stringr)
library(ggspatial)

source(paste0(path_to_scripts, "functions/spatial_analysis_functions.R"))
source(paste0(path_to_scripts, "functions/overall_prev.R"))

recombinant_ids=c("MMRJQCBD", # 1st XE
                  "MMPDJF4G",
                  "MM943GHM",
                  #"LLVP9JWC",
                  #"LLWMYX3K",
                  #"LL6V3HKX",
                  "MMPTM9RV",
                  "MMR4BP9K",
                  "MMMBCHRJ",
                  "MMWX7N3Q",
                  "MM78KNJ3", # 1st XJ
                  "MMRNJLCY", # 1st XL
                  "MMWNP8LR")

mylineage="BA.2"
#mylineage="sequenced"

setwd("E:/Group/react1_spatial_analysis")

# select round:'round_1_2', 'round_1', 'round_2'
round <- "round_19"

# read in correct round data
dat <- readRDS(paste0("output_files/clean_data/", round, "_data.rds"))

# read in pairwise distance matrices
pd <- readRDS(paste0("output_files/", round, "_sample_all/pairwise_distances.rds"))

# read in ordered pairwise distance matrices
pd_ordered <- readRDS(paste0("output_files/", round, "_sample_all/ordered_pd_mat.rds"))

# read in sorted pairwise distance matrices
pd_sorted <- readRDS(paste0("output_files/", round, "_sample_all/sorted_pd_mat.rds"))

# set save directory
save_dir <- paste0("output_files/", round, "_sample_all/")

# Loading lineage data
lineage=readRDS("E:/Group/saved_objects/rep19_lineage.rds")
if (mylineage=="sequenced"){
  ids_lineage=lineage$u_passcode
} else {
  ids_lineage=lineage$u_passcode[which(substr(lineage$react_lineage, start=1, stop=4)==mylineage)]
}

# Re-defining positives and negatives
dat[which(!dat$u_passcode%in%ids_lineage), "estbinres"]=0
dat[which(dat$u_passcode%in%ids_lineage), "estbinres"]=1
dat[which(dat$u_passcode%in%recombinant_ids), "estbinres"]=0
print(table(dat$estbinres, useNA = "always"))


##################################################################
# find neighbours within X km
##################################################################

# set maximum distance
max_dist <- 30

# take the distance sorted matrix and set value greater than max distance to NA
nd_mat <- pd_sorted
nd_mat[nd_mat > max_dist] <- NA

# find number of nearest neighbours in each column
num_nn <- apply(nd_mat, 2, function(x) {
  length(which(!is.na(x)))
})

# find median nn
max_nn <- floor(median(num_nn))

dir.create(file.path(paste0(save_dir, max_dist, "km_", max_nn, "nn_", gsub("\\.", "", mylineage))))

##########################################################################
# find neighbourhood prevalence for points up to X neighbours
##########################################################################

sample_id <- pd$sample_ids

# calc nearest neighbour prevalence
nn <- vapply(1:length(pd$sample_ids),
             nn_func,
             numeric(1),
             results_data = dat,
             ordered_pairwise_distance_matrix = pd_ordered,
             pairwise_distance_samples = sample_id,
             num_nearest_neighbours = max_nn
)

nn_df <- tibble(
  id = sample_id,
  prevalence = nn
)

######################################################################
# save result
######################################################################

save_dat <- nn_df %>%
  left_join(dat, by = "id")

saveRDS(save_dat,
        file = paste0(save_dir, max_dist, "km_", max_nn, "nn_", gsub("\\.", "", mylineage), "/nn_prev_results.rds")
)