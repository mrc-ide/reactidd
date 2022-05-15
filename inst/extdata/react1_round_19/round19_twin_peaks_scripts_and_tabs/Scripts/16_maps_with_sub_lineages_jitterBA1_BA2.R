rm(list=ls())

library("dplyr")
library("glmmfields")
library("ggplot2")
library("sf")
library("ggspatial")
library("viridis")
library("data.table")
source("E:/group/Ollies_workspace/bayesian_gam_functions.R")

round_ids=c(16,17,18,19)

for (round_id in round_ids) {
setwd(paste0("E:/Group/report/round_",round_id,"_lineage"))
df16 <- readRDS(paste0("E:/Group/saved_objects/rep",round_id,"_lineage.rds"))

df16 <- df16[c("u_passcode", "lineage_pan","uk_lineage","react_lineage","region","age_group_char","age","round","postcode","d_comb")]

ba2_ids=df16$u_passcode[which(df16$react_lineage%like%"BA.2")]
ba1_ids=df16$u_passcode[which(df16$react_lineage%like%"BA.1")]

df=readRDS(paste0("E:/Group/saved_objects/rep",round_id,".rds"))
df <- df[,c("u_passcode","region","age_group_char","age","postcode","d_comb")]
rownames(df)=df$u_passcode
df=df[c(ba2_ids, ba1_ids),]

df <- df %>%
  mutate(react_lineage = ifelse(u_passcode %in% ba2_ids, "BA.2", "BA.1"))

df <- mutate(df, del_omi = react_lineage)

tab_dates<-table( df$d_comb,df$del_omi, exclude = NULL)


df <- mutate(df, postcode = ifelse(nchar(postcode)==7, paste(substr(postcode,1,4),substr(postcode,5,7),sep=" "),
                                   ifelse(nchar(postcode)==6, paste(substr(postcode,1,3),substr(postcode,4,6),sep=" "),
                                          ifelse(nchar(postcode)==5, paste(substr(postcode,1,2),substr(postcode,3,5),sep=" "),NA))))

# Load and add lat lon data using postcodes
postcode_locs <- read.csv("E:/Group/aux_data/ukpostcodes.csv")
df <- inner_join(df, postcode_locs, by = "postcode")


# Load England shape file
eng_shp <- st_read("E:/Group/aux_data/GBR/GBR_adm1.shp") %>% filter(NAME_1 == "England")

# england regions convert to sf
rgn_shp <- readRDS("E:/Group/aux_data/england_regions.rds") %>%
  st_as_sf(coords = c("longitude","latitude"), remove = FALSE) %>%
  st_set_crs(4326)

ann_x <- -5.5
ann_y <- 56
region_ann_size <- 4
annotate_size <- 6

line_dat <- data.frame(x = c(0.6, 0),
                       y = c(50.7, 51.4))

df$Dates <-"Unknown"

# england regions convert to sf
region_shp <- readRDS("E:/Group/aux_data/england_regions.rds") %>%
  st_as_sf(coords = c("longitude", "latitude"), remove = FALSE) %>%
  st_set_crs(4326)

# Preparing map of England
ltla_map <- ggplot() +
  geom_sf(data = region_shp, fill = NA) +
  # scale_fill_gradient(low = low_colour, high = high_colour, limits = c(0, max_legend)) +
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

# ba1 background then switch ba2 background
if (round_id%in%c(17,18)) {
  raw_plot<-ltla_map+
    geom_jitter(data= df[df$del_omi=="BA.1"&is.na(df$Dates)==FALSE,], aes(x=longitude, y= latitude, color = as.factor(del_omi)),
                size=1.9,alpha=0.5, fill="#7873b4", stroke=0,
                position = position_jitter(width = 0.05, height =0.05, seed = 12)) +
    geom_jitter(data= df[df$del_omi=="BA.2"&is.na(df$Dates)==FALSE,], aes(x=longitude, y= latitude, color = as.factor(del_omi)),
                size=1.9,alpha=0.3, fill="#7874b5", stroke=0,
                position = position_jitter(width = 0.05, height =0.05, seed = 12))
} else {
  raw_plot<-ltla_map+
    geom_jitter(data= df[df$del_omi=="BA.2"&is.na(df$Dates)==FALSE,], aes(x=longitude, y= latitude, color = as.factor(del_omi)),
                size=1.9,alpha=0.3, fill="#7874b5", stroke=0,
                position = position_jitter(width = 0.05, height =0.05, seed = 12)) +
    geom_jitter(data= df[df$del_omi=="BA.1"&is.na(df$Dates)==FALSE,], aes(x=longitude, y= latitude, color = as.factor(del_omi)),
              size=1.9,alpha=0.5, fill="#7873b4", stroke=0,
              position = position_jitter(width = 0.05, height =0.05, seed = 12))
}

raw_plot

# if (round_id==17){
#   raw_plot=raw_plot+
#     scale_color_brewer(palette = "Set1")
# } else {
#   raw_plot=raw_plot+
#     scale_color_brewer(palette = "Dark2")
# }
#
# raw_plot

# save_dir=paste0("E:/Group/report/round",round_id,"/Figures/")
# ggsave(
#   plot = raw_plot,
#   filename = paste0(save_dir, "km_nn_ltla_map_other_r", round_id, ".pdf"),
#   device = "pdf",
#   height = 10,
#   width = 8
# )

ggsave(
  plot = raw_plot,
  filename = paste0("T:/km_nn_ltla_map_other_r", round_id, ".pdf"),
  device = "pdf",
  height = 10,
  width = 8
)

ggsave(
  plot = raw_plot,
  filename = paste0("T:/km_nn_ltla_map_other_r", round_id, ".png"),
  device = "png",
  height = 10,
  width = 8
)
}


raw_plot<-ltla_map+
  geom_jitter(data= df[df$del_omi=="BA.2"&is.na(df$Dates)==FALSE,], aes(x=longitude, y= latitude, color = as.factor(del_omi)),
              size=1.9,alpha=0.3, fill="#7874b5", stroke=0,
              position = position_jitter(width = 0.05, height =0.05, seed = 12)) +
  geom_jitter(data= df[df$del_omi=="BA.1"&is.na(df$Dates)==FALSE,], aes(x=longitude, y= latitude, color = as.factor(del_omi)),
              size=1.9,alpha=0.5, fill="#7873b4", stroke=0,
              position = position_jitter(width = 0.05, height =0.05, seed = 12))

raw_plot
