# LOAD REQUIRED PACKAGES AND FUNCTIONS -----------------------------------------
if (!require("dplyr")) install.packages("pacman")
pkgs = c("") # package names
pacman::p_load(pkgs, character.only = T)

source("R/functions.R")

# LOAD DATA --------------------------------------------------------------------

# Fitted model
fit <- readRDS("output/models/geotime_binomial.rds")

# TABLES -----------------------------------------------------------------------
tab <- create_tab(fit)
knitr::kable(tab, booktabs = T, linesep = "", align = "l", escape = F,
             caption = "Monte Carlo maximum likelihood estimates and corresponding 95\\% confidence intervals.", 
             format = "latex",
             col.names = c("Parameter", rep(c("Estimate", "95\\% CI"), times = 3))) %>%
  kable_styling(latex_options = c("HOLD_position", "scale_down")) %>% 
  add_header_above(c(" " = 1, "Baseline" = 2, "Midterm" = 2,
                     "Endline" = 2), bold = T) %>% 
  row_spec(0, bold = T) 

