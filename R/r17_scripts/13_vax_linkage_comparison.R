# 1. check link vaccine vs self declared vaccine status

# library(data.table)
library(tidyverse)
library(reaction)
library(progressr)


list_of_round_ids <- c(16,17)

for (round_id in list_of_round_ids) {

linked_df_path <- paste0("E:/dt20/linkedR", round_id, "datANG.rds")

linked_df   <- as.data.frame(readRDS(paste0("E:/dt20/linkedR", round_id, "datANG.rds")))
original_df <- as.data.frame(readRDS(paste0("E:/Group/saved_objects/rep", round_id, ".rds")))
cleaned_linked_df <- with_progress(clean_linked_df(linked_df_path, progressor))

output_path <- "E:/Group/report/round17/Tables"
#df17 <- as.data.table(readRDS("E:/dt20/linkedR17datANG.rds"))

# [586] "link_vax1date"                 "link_vax2date"                 "link_vax3date"
# [589] "link_vax1type"                 "link_vax2type"                 "link_boostertype"
# [592] "unlinked"                      "link_vax1toswab"               "link_vaxrecenttoswab"
# [595] "link_vax3toswab"               "diff_t1_t2"                    "diff_t2_t3"
# [598] "discard"                       "link_vax_status_booster0days"  "link_vax_status_booster7days"
# [601] "link_vax_status_booster14days" "link_v1time"                   "link_v1tov2"
# [604] "link_v2tov3"

calculate_difference <- function(df) {
  colnames(df) <- c("original", "linked")
  df <- df %>%
  mutate(difference = linked - original) %>%
  mutate(percentage_difference = difference / original * 100)
}
## first vaccination numbers -- tidyverse

no1 <- as.data.frame(cbind(table(original_df$vax_status),
                           table(linked_df$vax_status)))

rownames(no1) <- c("Unvaccinated", "Vaccinated", "Unknown")


no1 <- calculate_difference(no1)

no1

original_df$vax_status_noDate_v2 <- factor(
  recode(original_df$vax_status_noDate_v2,
         "One does" = "1 dose",
         "Two does" = "2 doses",
         "Three does" = "3 doses",
         "Not vaccinated" = "unvaccinated",
         "Unknown does" = "unknown"),
         c("1 dose", "2 doses", "3 doses",
           "unvaccinated", "unknown", "NA"))


# original_df <- original_df %>% fct_relevel(vax_status_noDate_v2,

no2 <- as.data.frame(t(bind_rows(table(original_df$vax_status_noDate_v2),
                           table(linked_df$link_vax_status_booster0days))))

no2 <- no2 %>%
  mutate(across(everything(), ~replace_na(.x, 0)))

no2 <- calculate_difference(no2)

# table(linked_df$link_vax_status_booster0days)
# table(linked_df$link_vax_status_booster14days)
# table(cleaned_linked_df$vax_status)

no2

## second vaccination brand

original_df$vax_type_char <- factor(
  recode(original_df$vax_type_char,
         "oxford" = "AZ",
         "pfizer" = "Pfizer",
         "moderna" = "Moderna",
         "NA" = "unknown"),
  c("AZ", "Pfizer", "Moderna",
    "unvaccinated", "unknown"))

no3 <-  as.data.frame(t(bind_rows(
  table(original_df$vax_type_char),
  table(cleaned_linked_df$vax_type))))

no3 <- calculate_difference(no3)

no3

# these do not have unvac and na numbers!
# table(linked_df$link_vax1type)
# table(linked_df$link_vax2type)
# table(linked_df$link_vax3type)

# TODO:: create output file ... in xlsx!


openxlsx::write.xlsx(
  rbind(no1,no2,no3),
  paste0(output_path, "/Vax_link_comparison_Round", round_id,
    Sys.Date(), ".xlsx"),
  rowNames = TRUE
)

}
