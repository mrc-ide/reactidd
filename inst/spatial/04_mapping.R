# LOAD REQUIRED PACKAGES AND FUNCTIONS -----------------------------------------
if (!require("pacman")) install.packages("pacman")
pkgs = c("sf", "tmap", "emojifont") # package names
pacman::p_load(pkgs, character.only = T)

source("R/functions.R")

# LOAD DATA --------------------------------------------------------------------

# Model predictions at LTLA level
rall <- st_read("output/predictions/ltla_mean.gpkg")
rchange <- st_read("output/predictions/ltla_change.gpkg")
rdiff <- st_read("output/predictions/ltla_diff.gpkg")

# LTLA boundaries
ltla <- st_read("data/original/geodata/ltla.gpkg")

# Region boundaries
regions <- st_read("data/original/geodata/Regions__December_2017__Boundaries.shp")
regions$labels <- c("NE", "NW", "YH", "EM", "WM", "EE", "L", "SE", "SW")

# MAPS -------------------------------------------------------------------------


# CHANGE OVER TIME PROBABILITIES -----------------------------------------------
palette = "-RdBu"
# 
# breaks <- c(0, 0.01, 0.05, 0.15, 0.30, 0.40, 0.60, 0.70, 0.85, 0.95, 0.99, 1)
# 
# # Generate labels for color legend
# labs <- c("< 0.01", "0.01 - 0.05", "0.05 - 0.15",
#           "0.15 - 0.30", "0.30 - 0.40", "0.40 - 0.60",
#           "0.60 - 0.70", "0.70 - 0.85",
#           "0.85 - 0.95", "0.95 - 0.99", "\U2265 0.99")

breaks <- c(0, 0.01, 0.05, 0.10, 0.90, 0.95, 0.99, 1)

# Generate labels for color legend
labs <- c("< 0.01", "0.01 - 0.05", "0.05 - 0.10", "0.10 - 0.90",
          "0.90 - 0.95", "0.95 - 0.99", "\U2265 0.99")

# Generate color palette
if (is.character(palette)) {
  pal <- tmaptools::get_brewer_pal(palette, n = length(labs), contrast = c(0, 1), plot = F)
} else {
  pal <- palette(length(labs))
}

ex <- tm_shape(rchange) +
  tm_fill(col = "pchange", palette = pal, style = "fixed",
          breaks = breaks, labels = labs, contrast = c(0, 1),
          title = bquote("Probability of an increase in prevalence")) +
  tm_facets(by = "round") +
  tm_shape(regions) +
  tm_borders(col = "black", lwd = 0.8) +
  tm_shape(st_centroid(regions)) +
  tm_text(text = "labels", size = 1, ymod = c(0, 2, 0, 0, 0, 0, 0, 0, 0), 
          xmod = c(0, 0, 0, 0, 0, 0, 0, -0.5, 0),
          col = "black") +
  tm_compass(size = 2, text.size = 1, position = c("right", "top")) +
  tm_scale_bar(position = c("right", "bottom"), text.size = .7) +
  tm_layout(outer.margins = 0, asp = 0, 
            legend.text.size = 0.9, legend.title.size = 1.35,
            panel.label.size = 1.4,
            legend.outside.size = .23, scale = 1.1)

tmap_save(ex, filename = paste0("figs/changeprev_LTLA.pdf"), 
          width = 5 * 3, height = 5)

# EXCEEDANCE PROBS. NATIONAL ---------------------------------------------------

palette = "-RdBu"

breaks <- c(0, 0.01, 0.05, 0.15, 0.30, 0.40, 0.60, 0.70, 0.85, 0.95, 0.99, 1)

# Generate labels for color legend
labs <- c("< 0.01", "0.01 - 0.05", "0.05 - 0.15",
          "0.15 - 0.30", "0.30 - 0.40", "0.40 - 0.60",
          "0.60 - 0.70", "0.70 - 0.85",
          "0.85 - 0.95", "0.95 - 0.99", "\U2265 0.99")

ex <- tm_shape(rall) +
  tm_fill(col = "pabove", palette = pal, style = "fixed",
          breaks = breaks, labels = labs, contrast = c(0, 1),
          title = "Probability of exceeding\nnational average") +
  tm_facets(by = "round") +
  tm_shape(regions) +
  tm_borders(col = "black", lwd = 0.8) +
  tm_shape(st_centroid(regions)) +
  tm_text(text = "labels", size = 1, ymod = c(0, 2, 0, 0, 0, 0, 0, 0, 0), 
          xmod = c(0, 0, 0, 0, 0, 0, 0, -0.5, 0),
          col = "black") +
  tm_compass(size = 2, text.size = 1, position = c("right", "top")) +
  tm_scale_bar(position = c("right", "bottom"), text.size = .7) +
  tm_layout(outer.margins = 0, asp = 0, 
            legend.text.size = 0.9, legend.title.size = 1.35,
            panel.label.size = 1.4, 
            legend.outside.size = .23, scale = 1.1)

tmap_save(ex, filename = paste0("figs/exnational_LTLA.pdf"), 
          width = 5 * 2, height = 5 * 2)

# MEAN PREVALENCE --------------------------------------------------------------

mp <- tm_shape(rall) +
  tm_fill(col = "mean_prev", palette = "Reds", style = "cont",
          title = "Predicted mean\nswab prevalence (%)") +
  tm_facets(by = "round") +
  tm_shape(regions) +
  tm_borders(col = "black", lwd = 0.8) +
  tm_shape(st_centroid(regions)) +
  tm_text(text = "labels", size = 1, ymod = c(0, 2, 0, 0, 0, 0, 0, 0, 0), 
          xmod = c(0, 0, 0, 0, 0, 0, 0, -0.5, 0),
          col = "black") +
  tm_compass(size = 2, text.size = 1, position = c("right", "top")) +
  tm_scale_bar(position = c("right", "bottom"), text.size = .7) +
  tm_layout(outer.margins = 0, asp = 0, 
            legend.text.size = 0.8, legend.title.size = 1.35,
            panel.label.size = 1.4, 
            legend.outside.size = .21, scale = 1.1)

tmap_save(mp, filename = paste0("figs/meanprev_LTLA.pdf"), 
          width = 5 * 2, height = 5 * 2)

# DIFFERENCE -------------------------------------------------------------------

df <- tm_shape(rdiff) +
  tm_fill(col = "diff", palette = "-RdBu", style = "cont",
          title = "Change in prevalence (%)") +
  tm_facets(by = "round") +
  tm_shape(regions) +
  tm_borders(col = "black", lwd = 0.8) +
  tm_shape(st_centroid(regions)) +
  tm_text(text = "labels", size = .7, ymod = c(0, 2, 0, 0, 0, 0, 0, 0, 0), 
          xmod = c(0, 0, 0, 0, 0, 0, 0, -0.5, 0),
          col = "black") +
  tm_compass(size = 2, text.size = 1, position = c("right", "top")) +
  tm_scale_bar(position = c("right", "bottom"), text.size = .7) +
  tm_layout(outer.margins = 0, asp = 0, 
            legend.text.size = 0.8, legend.title.size = 1.35,
            panel.label.size = 1.2,
            legend.outside.size = .21, scale = 1.1)

tmap_save(df, filename = paste0("figs/diffprev_LTLA.pdf"), 
          width = 5 * 3, height = 5)