library(dplyr)
library(geosphere)
library(tidyr)
library(ggplot2)
library(sf)
library(snakecase)

setwd("E:/Group/react1_spatial_analysis")

# select round:'round_1_2', 'round_1', 'round_2'
round <- "round_17"

# read in correct round data
dat <- readRDS(paste0("output_files/clean_data/", round, "_data.rds"))

# create a directory to save results for the results
dir.create(file.path(paste0("output_files/", round, "_sample_all")),
  showWarnings = FALSE
)

# set round save directory for the results
save_dir <- paste0("output_files/", round, "_sample_all/")

if (length(list.files(save_dir)) > 0) {
  tocopy <- list.files(save_dir)
  tocopy <- tocopy[!grepl("Previous", tocopy)]
  previous_dir <- paste0(save_dir, "Previous_", Sys.Date() - 1)
  dir.create(previous_dir,
    showWarnings = FALSE
  )
  file.copy(
    from = paste0(save_dir, tocopy),
    to = previous_dir, recursive = TRUE
  )
  unlink(paste0(save_dir, tocopy), recursive = TRUE)
}


######################################################################
# select sample ids to run pairwise distance calculations on
######################################################################

# want the same number of sample ids per ltla

samp_per_ltla <- 15

samp_dat <- dat %>%
  group_by(lacode) %>%
  sample_n(size = samp_per_ltla)

samps <- samp_dat %>% pull(id)


# ##################################################################
# # calculate and save pairwise distance matrix
# ##################################################################

pd <- pairwise_distances(
  df = dat,
  sample_ids = samps
)


saveRDS(pd, paste0(save_dir, "pairwise_distances.rds"))
