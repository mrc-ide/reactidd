library(dplyr)
library(geosphere)
library(tidyr)
library(ggplot2)
library(sf)
library(snakecase)
library(readr)
library(janitor)
library(prevalence)

setwd("E:/Group/react1_spatial_analysis")

# select round:'round_1_2', 'round_1', 'round_2'
round <- "round_17"

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

dir.create(file.path(paste0(save_dir, max_dist, "km_", max_nn, "nn")))

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
  file = paste0(save_dir, max_dist, "km_", max_nn, "nn/nn_prev_results.rds")
)















##################################################################
# find neighbours within X km
##################################################################

# set maximum distance
max_dist <- 20

# take the distance sorted matrix and set value greater than max distance to NA
nd_mat <- pd_sorted
nd_mat[nd_mat > max_dist] <- NA

# find number of nearest neighbours in each column
num_nn <- apply(nd_mat, 2, function(x) {
  length(which(!is.na(x)))
})

# find median nn
max_nn <- floor(median(num_nn))

dir.create(file.path(paste0(save_dir, max_dist, "km_", max_nn, "nn")))

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
  file = paste0(save_dir, max_dist, "km_", max_nn, "nn/nn_prev_results.rds")
)













##################################################################
# find neighbours within X km
##################################################################

# set maximum distance
max_dist <- 40

# take the distance sorted matrix and set value greater than max distance to NA
nd_mat <- pd_sorted
nd_mat[nd_mat > max_dist] <- NA

# find number of nearest neighbours in each column
num_nn <- apply(nd_mat, 2, function(x) {
  length(which(!is.na(x)))
})

# find median nn
max_nn <- floor(median(num_nn))

dir.create(file.path(paste0(save_dir, max_dist, "km_", max_nn, "nn")))

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
  file = paste0(save_dir, max_dist, "km_", max_nn, "nn/nn_prev_results.rds")
)
