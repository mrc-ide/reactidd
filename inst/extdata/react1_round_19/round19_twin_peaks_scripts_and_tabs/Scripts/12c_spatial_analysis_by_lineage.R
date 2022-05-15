rm(list=ls())

library("dplyr")
library("glmmfields")
library("ggplot2")
library("sf")
library("ggspatial")
library("viridis")
source("E:/group/Ollies_workspace/bayesian_gam_functions.R")

#' Set the working directory and load up the current round 8 data
setwd("E:/Group/report/round_19_lineage")
df16 <- readRDS("E:/Group/saved_objects/rep19_lineage.rds")

df16 <- df16[c("u_passcode", "lineage_pan","uk_lineage","react_lineage","region","age_group_char","age","round","postcode","d_comb")]

recombinant_ids=c("MMRJQCBD", # 1st XE
                  "MMPDJF4G",
                  "MM943GHM",
                  #"LLVP9JWC",
                  #"LLWMYX3K",
                  #"LL6V3HKX",
                  "MMPTM9RV",
                  "MMR4BP9K",
                  "MMMBCHRJ",
                  "MMWX7N3Q",
                  "MM78KNJ3", # 1st XJ
                  "MMRNJLCY", # 1st XL
                  "MMWNP8LR")

ba3_ids=df16$u_passcode[which(df16$react_lineage=="BA.3")]

df=readRDS("E:/Group/saved_objects/rep19.rds")
df <- df[,c("u_passcode","region","age_group_char","age","postcode","d_comb")]
rownames(df)=df$u_passcode
df=df[c(recombinant_ids, ba3_ids),]
df=cbind(df, react_lineage=c(rep("XE", 7),
                             rep("XJ", 1),
                             rep("XL", 2),
                             rep("BA.3", length(ba3_ids))))

# df <- mutate(df, del_omi = ifelse(react_lineage=="BA.3","BA.3","Recombinant"))
df <- mutate(df, del_omi = react_lineage)

tab_dates<-table( df$d_comb,df$del_omi, exclude = NULL)


###################################################################
df$round <- as.factor(rep("Round 19", nrow(df)))

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


#tab1<-table(df$react_lineage)
#lin_list <- names(tab1[tab1>9])
#df<-df[df$react_lineage%in%lin_list,]




df <- mutate(df, Dates = ifelse(d_comb < as.Date("2021-12-01"),"2021-11-23 to 2021-11-30",
                                ifelse(d_comb<as.Date("2021-12-07"),"2021-12-01 to 2021-12-06",
                                       ifelse(d_comb<as.Date("2021-12-13"),"2021-12-07 to 2021-12-12","2021-12-13 to 2021-12-17") )))

df$Dates <-"Unknown"


raw_plot<-ggplot()+
  geom_sf(data = rgn_shp)+
  theme_bw(base_size = 18)+
  geom_jitter(data= df[df$del_omi=="XE"&is.na(df$Dates)==FALSE,], aes(x=longitude, y= latitude, color = as.factor(del_omi)),
              size=1.5,alpha=0.8,
              position = position_jitter(width = 0.05, height =0.05, seed = 12))+
  geom_jitter(data= df[df$del_omi=="XJ"&is.na(df$Dates)==FALSE,], aes(x=longitude, y= latitude, color = as.factor(del_omi)),
              size=1.5,alpha=0.8,
              position = position_jitter(width = 0.05, height =0.05, seed = 12))+
  geom_jitter(data= df[df$del_omi=="XL"&is.na(df$Dates)==FALSE,], aes(x=longitude, y= latitude, color = as.factor(del_omi)),
              size=1.5,alpha=0.8,
              position = position_jitter(width = 0.05, height =0.05, seed = 12))+
  geom_jitter(data= df[df$del_omi=="BA.3"&is.na(df$Dates)==FALSE,], aes(x=longitude, y= latitude, color = as.factor(del_omi)),
              size=1.5,alpha=0.8,
              position = position_jitter(width = 0.05, height =0.05, seed = 12))+
  coord_sf(ylim=c(50,56))+
  scale_alpha(guide = "none")+
  labs(color = "Lineage")+
  #facet_wrap(.~Dates, ncol = 3)+
  ylim(c(50.2,55.5))+
  xlim(c(-5.5,1.8))+
  theme(legend.position = "bottom")+
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())+
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank()) +
  annotate("text", x = -1.9, y = 55.3, label = "NE", size = region_ann_size) +
  annotate("text", x = -3, y = 54.7, label = "NW", size = region_ann_size) +
  annotate("text", x = -0.8, y = 54.2, label = "YH", size = region_ann_size) +
  annotate("text", x = -2.5, y = 52.7, label = "WM", size = region_ann_size) +
  annotate("text", x = -0.2, y = 53.3, label = "EM", size = region_ann_size) +
  annotate("text", x = 1, y = 52.8, label = "EE", size = region_ann_size) +
  annotate("text", x = -0.5, y = 51.0, label = "SE", size = region_ann_size) +
  annotate("text", x = -3.7, y = 51.0, label = "SW", size = region_ann_size) +
  annotate("text", x = 0.75, y = 50.65, label = "L", size = region_ann_size) +
  geom_line(data = line_dat, aes(x = x, y = y))+
  scale_color_brewer(palette = "Dark2")


raw_plot

save_dir="E:/Group/report/round19/Figures/"
round_id=19
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

