rm(list=ls())
library(data.table)

round_id=17
data_file <- paste0("E:/Group/saved_objects/rep", round_id, ".rds")
df_round <- data.frame(readRDS(data_file))
bridge=fread("E:/Group/report/round16/Parameters/urban_vs_rural.csv", data.table=FALSE)

bridge$postcode=gsub(" ", "", bridge$postcode)
df_round$postcode=gsub(" ", "", df_round$postcode)

# table(!df_round$postcode%in%bridge$postcode)

mapping=bridge$URBAN
names(mapping)=bridge$postcode

df_round$urban=mapping[df_round$postcode]

saveRDS(df_round, paste0("E:/Group/saved_objects/rep", round_id, ".rds"))