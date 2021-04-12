# LOAD REQUIRED PACKAGES AND FUNCTIONS -----------------------------------------
if (!require("pacman")) install.packages("pacman")
pkgs = c("PrevMap", "raster", "sf", "dplyr") # package names
pacman::p_load(pkgs, character.only = T)

source("R/functions.R")

# LOAD DATA --------------------------------------------------------------------

# LTLA boundaries
ltla <- st_read("data/original/geodata/ltla.gpkg")

# Population raster
pop_raster <- raster("data/original/geodata/gbr_ppp_2020_1km_Aggregated.tif")

# Fitted spatio-temporal model
fit <- readRDS("output/models/geotime_binomial.rds")

# Compliance data
compliance <- readr::read_csv("data/processed/react_clean.csv") 

# DATA PROCESSING --------------------------------------------------------------

# Create 5km grid for predictions
pred_coords_sf <- st_make_grid(ltla, cellsize = 5000, what = "centers") %>% 
  st_as_sf() %>% 
  st_join(ltla, join = st_within, left = F) 

# If some small LTLA is left out, add pop_weighted centroid for that LTLA
if(length(unique(pred_coords_sf$lad19cd)) < length(unique(ltla$lad19cd))) {
  miss_ltlas <- setdiff(ltla$lad19cd, unique(pred_coords_sf$lad19cd))
  temp <- compliance[compliance$round == "round 1", ]
  pred_coords_miss <- temp[temp$lacode %in% miss_ltlas, c("utm_x", "utm_y")]
  names(pred_coords_miss) <- c("X", "Y")
  
  pred_coords <- pred_coords_sf %>% 
    st_coordinates() / 1000 
  
  pred_coords <- rbind(pred_coords, as.matrix(pred_coords_miss))
  
} else {
  
  pred_coords <- pred_coords_sf %>% 
    st_coordinates() / 1000 
  
}

# Create predictors:
# - extract LTLA for each point
# - assign round
# - extract compliance
admin <- c(as.character(pred_coords_sf$lad19cd), miss_ltlas)
predictors <- data.frame(lacode = rep(admin, times = 4),
                         round = rep(paste("round", 1:4), each = nrow(pred_coords)))

compliance$lacode[compliance$lacode == "E06000052"] <- "E06000052,E06000053"
compliance$lacode[compliance$lacode == "E09000033"] <- "E09000001,E09000033"

predictors <- predictors %>% 
  inner_join(compliance[c("lacode", "round", "elog_comprate")]) %>% 
  select(-lacode)

grid.pred <- rbind(pred_coords, pred_coords, pred_coords, pred_coords)

# OBTAIN PREDICTIONS -----------------------------------------------------------

cmcmc <- control.mcmc.MCML(n.sim = 11000, burnin = 1000, thin = 10)

predictions <- sptime_pred(object = fit, 
                           grid.pred = grid.pred, 
                           time.pred = rep(as.numeric(unique(fit$times)),
                                           each = nrow(pred_coords)), 
                           predictors = predictors, 
                           control.mcmc = cmcmc, 
                           type = "joint", 
                           save.sim = T, 
                           plot.correlogram = T, 
                           messages = T)
saveRDS(predictions, file = "output/predictions/sptime_pred.rds")

# SUMMARIES PREDICTIONS --------------------------------------------------------

times <- unique(predictions$times)

# Aggregate at LTLA level
for(i in 1:length(times)) {
  
  id <- which(predictions$times == times[i])
  
  # Transform prediction grid centroid to sf object
  grid_sf <- st_as_sf(as.data.frame(predictions$coords[id, ] * 1000), 
                      coords = c(1, 2), 
                      crs =  st_crs(ltla)) %>% 
    st_transform(crs = crs(pop_raster))
  
  # Extract population at prediction locations 
  pop <- raster::extract(pop_raster, grid_sf)
  
  # Calculate national prevalence
  pop_weigths <- pop / sum(pop)
  prev_samples <- plogis(predictions$samples[id, ])
  pmean <- rowMeans(prev_samples)
  assign(paste0("nationalp", i), sum(pmean * pop_weigths, na.rm = T))
  
  # Aggregate at LTLA level 
  samples <- data.frame(lad19cd = admin, pop, prev_samples) %>% 
    na.omit() %>% 
    joint_weighted_prev() 
  
  names(samples)[3] <- paste0("prev", i)
  samples$round <- paste("Round", i)
  
  # Store the samples
  assign(paste0("round", i, "_samples"), samples)
}

# Calculate mean prevalence at LTLA level

r1 <- round1_samples %>% 
  group_by(lad19cd, round) %>% 
  summarise(mean_prev = mean(prev1) * 100, pabove = mean(prev1 > nationalp1)) %>%
  ungroup() %>% 
  inner_join(ltla)

r2 <- round2_samples %>% 
  group_by(lad19cd, round) %>% 
  summarise(mean_prev = mean(prev2) * 100, pabove = mean(prev2 > nationalp2)) %>%
  ungroup() %>% 
  inner_join(ltla)

r3 <- round3_samples %>% 
  group_by(lad19cd, round) %>% 
  summarise(mean_prev = mean(prev3) * 100, pabove = mean(prev3 > nationalp3)) %>%
  ungroup() %>% 
  inner_join(ltla)

r4 <- round4_samples %>% 
  group_by(lad19cd, round) %>% 
  summarise(mean_prev = mean(prev4) * 100, pabove = mean(prev4 > nationalp4)) %>%
  ungroup() %>% 
  inner_join(ltla)

# Put everything together and save
rall <- rbind(r1, r2, r3, r4)
st_write(rall, "output/predictions/ltla_mean.gpkg", delete_dsn = T)


# Calculate Prob(P_t > P_{t - 1})
rchange <- round1_samples %>% 
  select(-round) %>% 
  inner_join(select(round2_samples, -round)) %>% 
  inner_join(select(round3_samples, -round)) %>% 
  inner_join(select(round4_samples, -round)) %>% 
  group_by(lad19cd) %>% 
  summarise(prob12 = mean(prev2 > prev1),
            prob23 = mean(prev3 > prev2),
            prob34 = mean(prev4 > prev3)) %>% 
  tidyr::gather(key = "round", value =  "pchange", -lad19cd) %>% 
  inner_join(ltla) %>% 
  mutate(round = factor(round, 
                        labels = c("Round 1 to 2", "Round 2 to 3", "Round 3 to 4"))) 
# Save
st_write(rchange, "output/predictions/ltla_change.gpkg", delete_dsn = T)

# Calculate E(P_t - P_{t - 1})
rdiff <- round1_samples %>% 
  select(-round) %>% 
  inner_join(select(round2_samples, -round)) %>% 
  inner_join(select(round3_samples, -round)) %>% 
  inner_join(select(round4_samples, -round)) %>% 
  group_by(lad19cd) %>% 
  summarise(diff12 = mean(prev2 - prev1) * 100,
            diff23 = mean(prev3 - prev2) * 100,
            diff34 = mean(prev4 - prev3) * 100) %>% 
  tidyr::gather(key = "round", value =  "diff", -lad19cd) %>% 
  inner_join(ltla) %>% 
  mutate(round = factor(round, 
                        labels = c("Round 1 to 2", "Round 2 to 3", "Round 3 to 4"))) 
# Save
st_write(rdiff, "output/predictions/ltla_diff.gpkg", delete_dsn = T)

# Probability that all LTLAs in Round 4 have a higher prevalence thant Round 3
samples3 <- round3_samples %>%
  arrange(lad19cd) %>% 
  pull(prev3) %>% 
  matrix(nrow = 1000, ncol = 315, byrow = F)

samples4 <- round4_samples %>%
  arrange(lad19cd) %>% 
  pull(prev4) %>% 
  matrix(nrow = 1000, ncol = 315, byrow = F)

mean(apply(samples4 - samples3, 1, function(x) all(x > 0)))

# Double check
select(round3_samples, -round) %>% 
  inner_join(select(round4_samples, -round)) %>% 
  mutate(diff = prev4 - prev3) %>% 
  group_by(sample_id) %>% 
  summarise(out = all(diff > 0)) %>% 
  summarise(mean(out))

