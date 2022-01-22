# check vax status for omicron cases

library(tidyverse)
library(plyr)
library(openxlsx)

round_id <- 16

# Identifying delta and omicron cases
df_round <- readRDS(paste0("E:/dt20/linkedR", round_id, "datANG.rds"))
df_lineage=readRDS(paste0("E:/Group/saved_objects/rep",round_id,"_lineage.rds"))

ids_omicron=df_lineage$u_passcode[which(df_lineage$react_lineage=="BA.1")]
ids_delta=df_lineage$u_passcode[which(df_lineage$react_lineage!="BA.1")]

df_round$lineage=rep(NA,nrow(df_round))
df_round[ids_omicron, "lineage"]="omicron"
df_round[ids_delta, "lineage"]="delta"

df_omicron <- df_round %>%
  filter(u_passcode %in% ids_omicron)

df_delta <- df_round %>%
  filter(u_passcode %in% ids_delta)


table(df_omicron$link_vax1type)
table(df_omicron$link_vax_status_booster0days)
prop.table(table(df_omicron$link_vax_status_booster0days))

table(df_delta$link_vax1type)
table(df_delta$link_vax_status_booster0days)
prop.table(table(df_delta$link_vax_status_booster0days))*100

output_table <- rbind.fill(
  as.data.frame(rbind(table(df_delta$link_vax_status_booster0days))),
  as.data.frame(rbind(table(df_omicron$link_vax_status_booster0days))))

output_table[is.na(output_table)] <- 0

print("Chi-Sq test Delta / Omicron")
print(chisq.test(output_table))

chisq.test(output_table)[3]
chisq.test(output_table)$p.value

output_table <- cbind(output_table, chisq.test(output_table)$p.value)
rownames(output_table) <- c("delta", "omicron")
print(output_table)

output_filename <- "E:/Group/report/round16/Tables/Omicron_vs_delta_linked_vacc_chiSq"

write.xlsx(output_table,
  paste0(output_filename, "_r", paste(round_id, collapse=""), "_", Sys.Date(), ".xlsx"),
  rowNames = TRUE
)

file.copy(
  from = paste0(output_filename, "_r", paste(round_id, collapse=""), "_", Sys.Date(), ".xlsx"),
  to = "T:/", overwrite = TRUE
)

print(output_filename)
