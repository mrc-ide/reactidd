library(dplyr)
library(geosphere)
library(tidyr)
library(ggplot2)
library(sf)
library(snakecase)
library(readr)
library(janitor)

setwd("E:/Group/react1_spatial_analysis")

# select round:'round_1_2', 'round_1', 'round_2'
round <- "round_19"

# read in correct round data
dat <- readRDS(paste0("output_files/clean_data/", round, "_data.rds"))

# read in pairwise distance matrices
pd <- readRDS(paste0("output_files/", round, "_sample_all/pairwise_distances.rds"))

# set save directory
save_dir <- paste0("output_files/", round, "_sample_all/")

######################################################################
# order and sort the data
######################################################################

# column i of the pairwise distance matrix gives the vector of pairwise
#  distances between all sample data points and sample id at position i.

# apply the 'order' operation over each matrix column to find the sample id
#  closest (row 1) to furthest. Note that the order operation returns row
#  numbers, however here the row number corresponds to the sample ids in the data.

# apply the 'sort' operation over each matrix column to distances
#  closest (row 1) to furthest. Note that the ordered matrix will be required for
#  pairing up sample ids.


# extract the matrix
pd_mat <- pd$pairwise_distance_matrix

# order the results - matrix of positions
ordered_pd_mat <- apply(pd_mat, 2, order)

# sort the results - matrix of distances
sorted_pd_mat <- apply(pd_mat, 2, sort)

saveRDS(ordered_pd_mat,
  file = paste0(save_dir, "ordered_pd_mat.rds")
)

saveRDS(sorted_pd_mat,
  file = paste0(save_dir, "sorted_pd_mat.rds")
)
