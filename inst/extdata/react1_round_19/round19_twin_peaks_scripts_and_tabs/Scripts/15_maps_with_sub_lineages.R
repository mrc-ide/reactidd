rm(list=ls())

library("dplyr")
library("glmmfields")
library("ggplot2")
library("sf")
library("ggspatial")
library("viridis")
source("E:/group/Ollies_workspace/bayesian_gam_functions.R")

round_id=19

#' Set the working directory and load up the current round 8 data
setwd(paste0("E:/Group/report/round_",round_id,"_lineage"))
df16 <- readRDS(paste0("E:/Group/saved_objects/rep",round_id,"_lineage.rds"))

df16 <- df16[c("u_passcode", "lineage_pan","uk_lineage","react_lineage","region","age_group_char","age","round","postcode","d_comb")]

if (round_id==16){
  ba2_ids=df16$u_passcode[which(df16$react_lineage%in%c("BA.1", "BA.1.1"))]
} else {
  ba2_ids=df16$u_passcode[which(df16$react_lineage=="BA.2")]
}

df=readRDS(paste0("E:/Group/saved_objects/rep",round_id,".rds"))
df <- df[,c("u_passcode","region","age_group_char","age","postcode","d_comb")]
rownames(df)=df$u_passcode
df=df[c(ba2_ids),]
df=cbind(df, react_lineage=c(rep("BA.2", length(ba2_ids))))

# df <- mutate(df, del_omi = ifelse(react_lineage=="BA.3","BA.3","Recombinant"))
df <- mutate(df, del_omi = react_lineage)

tab_dates<-table( df$d_comb,df$del_omi, exclude = NULL)


###################################################################
#df$round <- as.factor(rep("Round 17", nrow(df)))

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

raw_plot<-ltla_map+
  geom_jitter(data= df[df$del_omi=="BA.2"&is.na(df$Dates)==FALSE,], aes(x=longitude, y= latitude, color = as.factor(del_omi)),
              size=1.5,alpha=0.8,
              position = position_jitter(width = 0.05, height =0.05, seed = 12))

if (round_id==17){
  raw_plot=raw_plot+
    scale_color_brewer(palette = "Set1")
} else {
  raw_plot=raw_plot+
    scale_color_brewer(palette = "Dark2")
}

raw_plot

save_dir=paste0("E:/Group/report/round",round_id,"/Figures/")
ggsave(
  plot = raw_plot,
  filename = paste0(save_dir, "km_nn_ltla_map_other_r", round_id, ".pdf"),
  device = "pdf",
  height = 10,
  width = 8
)

ggsave(
  plot = raw_plot,
  filename = paste0("T:/km_nn_ltla_map_other_r", round_id, ".pdf"),
  device = "pdf",
  height = 10,
  width = 8
)

ggsave(
  plot = raw_plot,
  filename = paste0(save_dir, "km_nn_ltla_map_other_r", round_id, ".png"),
  device = "png",
  height = 10,
  width = 8
)
