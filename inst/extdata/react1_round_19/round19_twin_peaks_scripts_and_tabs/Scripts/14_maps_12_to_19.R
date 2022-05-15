# Script to run spatial analyses and produce the maps
# Currently done for England and South West
# Outputs: 2 PDF documents (maps)
# and 1 csv with LTLAs ordered by decreasing prevalence
# stored in "E:/Group/react1_spatial_analysis/output_files/round_XXX_sample_all/30km_XXX"
# and automatically copied into the transfer folder

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

run_models <- FALSE

# Running spatial data wrangling and models
if (run_models) {
  print("Running step 1/4...")
  source(paste0(path_to_scripts, "spatial_analysis/1_data_wrangling_all_rounds.R"))
  
  print("Running step 2/4...")
  source(paste0(path_to_scripts, "spatial_analysis/2_pairwise_distance_sample_all.R"))
  
  print("Running step 3/4...")
  source(paste0(path_to_scripts, "spatial_analysis/3_order_sort_matrices_sample_all.R"))
  
  print("Running step 4/4...")
  source(paste0(path_to_scripts, "spatial_analysis/4_nearest_neighbour_distance_sample_all.R"))
}

for (round_id in 12:19){
  print(round_id)
  
  # Parameters
  direct_export <- TRUE
  # round_id <- 19
  dist <- "30"
  
  # Path to output folder
  output_folder <- paste0("E:/Group/react1_spatial_analysis/output_files/round_", round_id, "_sample_all/")
  myfiles <- list.files(output_folder)
  model_folder <- myfiles[grep(paste0(dist, "km_"), myfiles)]
  model_folder=model_folder[1]
  save_dir <- paste0(output_folder, model_folder, "/")
  
  # Extracting number of neighbours
  nn15 <- gsub("nn", "", gsub(".*_", "", model_folder))
  
  # read in correct round data
  dat <- readRDS(paste0(save_dir, "nn_prev_results.rds"))
  
  ##############################################################################
  # prevalence by ltla
  ##############################################################################
  
  # find number of positives and total number of samples per LTLA
  round_dat=readRDS(paste0("E:/Group/saved_objects/rep", round_id, ".rds"))
  # choose relevant columns
  names_to_keep <- c("id",
                     "estbinres",
                     "region",
                     "lacode")
  round_dat <- round_dat %>%
    dplyr::select(all_of(names_to_keep)) %>%
    dplyr::filter(!(is.na(estbinres)))
  dat_ltla <- round_dat %>%
    group_by(lacode, region) %>%
    summarise(positive = sum(estbinres),
              number_samples = length(estbinres)) %>%
    mutate(round = paste0("round ", round_id))
  
  prev <- propCI(x = dat_ltla$positive,
                 n = dat_ltla$number_samples,
                 method = "wilson"
  )
  
  dat_ltla <- dat_ltla %>%
    bind_cols(prev) %>%
    dplyr::select(-c(x, n, method, level)) %>%
    rename(prevalence = p)
  
  ############################################################################
  # load shapefiles for plotting
  ############################################################################
  
  # read in lad to region
  lad_to_region <- read_csv("E:/Group/aux_data/Local_Authority_District_to_Region_(April_2019)_Lookup_in_England.csv")
  lad_to_region <- clean_names(lad_to_region)
  
  # load england shapefile
  eng_shp <- st_read("E:/Group/aux_data/GBR/GBR_adm1.shp") %>% filter(NAME_1 == "England")
  
  # load lad shapefile
  lad_shp <- st_read("E:/Group/aux_data/Local_Authority_Districts__December_2019__Boundaries_UK_BGC-shp/Local_Authority_Districts__December_2019__Boundaries_UK_BGC.shp") %>%
    st_transform(crs = st_crs(4326)) %>%
    left_join(dplyr::select(lad_to_region, c(lad19cd, rgn19cd, rgn19nm)),
              by = "lad19cd"
    )
  
  # england regions convert to sf
  region_shp <- readRDS("E:/Group/aux_data/england_regions.rds") %>%
    st_as_sf(coords = c("longitude", "latitude"), remove = FALSE) %>%
    st_set_crs(4326)
  
  #########################################################################
  # nn numbers
  #########################################################################
  
  nn_nums <- tibble(
    distance = as.numeric(dist),
    round = c(paste0("round ", round_id)),
    nn = c(nn15)
  )
  
  write.csv(nn_nums,
            file = paste0(save_dir, "nn_size_", dist, "km.csv"),
            row.names = FALSE
  )
  
  ##############################################################################
  # joining up with shapefile for maps and LTLA names
  ##############################################################################
  
  
  dat_ltla_shp <- lad_shp %>%
    right_join(dat_ltla,
               by = c("lad19cd" = "lacode"))
  
  ##############################################################################
  # generate table
  ##############################################################################
  
  ordered_ltla_df <- dat_ltla_shp %>%
    st_drop_geometry() %>%
    dplyr::select(lad19nm, lad19cd, region, positive, number_samples, prevalence,
                  lower, upper) %>%
    rename(ltla_code = lad19cd,
           ltla = lad19nm) %>%
    arrange(desc(prevalence))
  
  
  write.csv(ordered_ltla_df,
            file = paste0(save_dir, "unwt_ordered_ltla_prev_", dist, "km_r",round_id,".csv"),
            row.names = FALSE)
  
  
  #########################################################################
  # ltla prevalence
  #########################################################################
  
  # average the neighbourhood prevalence over each ltla
  
  ltla_dat <- dat %>%
    group_by(round, lacode, laname) %>%
    summarise(ltla_prevalence = mean(prevalence)) %>%
    mutate(round = str_replace_all(round, pattern = "_", replacement = " ")) %>%
    left_join(lad_shp,
              by = c("lacode" = "lad19cd")
    ) %>%
    st_as_sf()
  
  
  ordered_ltla <- ltla_dat %>%
    st_drop_geometry() %>%
    dplyr::select(round, laname, rgn19nm, ltla_prevalence) %>%
    rename(
      region = rgn19nm,
      ltla = laname
    ) %>%
    group_by(round) %>%
    arrange(desc(ltla_prevalence)) %>%
    dplyr::select(ltla, region, ltla_prevalence, round)
  
  write.csv(ordered_ltla,
            file = paste0(save_dir, dist, "km_ordered_prev.csv")
  )
  
  
  # Map parameters
  ann_x <- -5.5
  ann_y <- 56
  region_ann_size <- 4
  annotate_size <- 6
  line_dat <- data.frame(
    x = c(0.6, 0),
    y = c(50.7, 51.4)
  )
  
  low_colour <- "#FEF0F0"
  high_colour <- "darkred"
  
  max_legend <- 0.09
  print(max_legend)
  
  ltla_dat$round <- rep("Round 15 partial", nrow(ltla_dat))
  
  # Preparing map of England
  ltla_map <- ggplot() +
    geom_sf(data = ltla_dat, aes(fill = ltla_prevalence), colour = NA) +
    geom_sf(data = region_shp, fill = NA) +
    scale_fill_gradient(low = low_colour, high = high_colour, limits = c(0, max_legend)) +
    theme_bw() +
    xlab("") +
    ylab("") +
    labs(fill = "prevalence") +
    # facet_wrap(~ round) +
    theme(
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      panel.grid.major = element_blank()
    ) +
    annotate("text", x = -1.9, y = 55.3, label = "NE", size = region_ann_size) +
    annotate("text", x = -3, y = 54.7, label = "NW", size = region_ann_size) +
    annotate("text", x = -0.8, y = 54.2, label = "YH", size = region_ann_size) +
    annotate("text", x = -2.5, y = 52.7, label = "WM", size = region_ann_size) +
    annotate("text", x = -0.2, y = 53.3, label = "EM", size = region_ann_size) +
    annotate("text", x = 1, y = 52.8, label = "EE", size = region_ann_size) +
    annotate("text", x = -0.5, y = 51.0, label = "SE", size = region_ann_size) +
    annotate("text", x = -3.7, y = 51.0, label = "SW", size = region_ann_size) +
    annotate("text", x = 0.75, y = 50.65, label = "L", size = region_ann_size) +
    geom_line(data = line_dat, aes(x = x, y = y)) +
    annotation_scale(location = "bl")
  
  ltla_map
  
  save_dir="T:/"
  ggsave(
    plot = ltla_map,
    filename = paste0(save_dir, dist, "km_nn_ltla_map_r", round_id, ".pdf"),
    device = "pdf",
    height = 10,
    width = 8
  )
  
  ggsave(
    plot = ltla_map,
    filename = paste0(save_dir, dist, "km_nn_ltla_map_r", round_id, ".png"),
    device = "png",
    height = 10,
    width = 8
  )
}
