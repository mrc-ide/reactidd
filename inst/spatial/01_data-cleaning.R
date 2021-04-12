# LOAD REQUIRED PACKAGES AND FUNCTIONS -----------------------------------------
if (!require("pacman")) install.packages("pacman")
if (!require("pacman")) install.packages("readr")
pkgs = c("dplyr", "sf", "readr") # package names
pacman::p_load(pkgs, character.only = T)

source("R/functions.R")

# LOAD DATA --------------------------------------------------------------------

# Prevalence data from REACT1 rounds 1-4 at LTLA level
react <- readr::read_csv("data/original/r1_4_unwt_ltla_prev.csv")

# Population-weighted centroids at LTLA level
ltla_pw <- st_read("data/original/geodata/LTLA_popweightedcentroid.shp")

# Compliance data
compliance <- readr::read_csv("data/original/compliance_r1to4.csv")

# DATA PREPARATION -------------------------------------------------------------

# For round 1, aggregate:
# - Cornwall E06000052 and Isles of Sicily E06000053 
react$number_samples[react$lacode == "E07000009" & react$round == "round 1"] <- 404 + 1
react$number_samples[react$lacode == "E06000052" & react$round == "round 1"] <- 312 

react <- react[-which(react$lacode %in% c("E0700009", "E06000053") & react$round == "round 1"), ]

# Convert LTLA pop-weighted centroids from meters to KM
ltla_pw <- ltla_pw %>% 
  mutate(utm_x = X / 1000, utm_y = Y / 1000) %>% 
  st_drop_geometry() 

# Calculate empirical log-odds of compliance rate
compliance$elog_comprate <- elogit(compliance$swab, compliance$issued)
compliance$comp_prev <- compliance$swab / compliance$issued

# Merge centroids and compliance data with the main dataset
react_clean <- react %>% 
  inner_join(ltla_pw[c("lad19cd", "utm_x", "utm_y")], by = c("lacode" = "lad19cd")) %>% 
  inner_join(compliance) %>% 
  select(lacode, round, positive, tested = number_samples, utm_x, utm_y, 
         elog_comprate, comp_prev) %>% 
  mutate(elogit = elogit(num = positive, den = tested), prev = positive / tested)

# Save cleaned dataset
readr::write_csv(react_clean, "data/processed/react_clean.csv")
