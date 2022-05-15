library(ggplot2)
library(tidyverse)
library(ggpubr)

l95wtPrev5 <- c(1.12, 0.56, 0.26, 0.07, 0.24, 0.74, 2.05, 4.16, 4.15, 7.10, 4.01, 7.89)
wtPrev5    <- c(1.37, 0.75, 0.41, 0.14, 0.38, 0.98, 2.44, 4.76, 4.74, 7.85, 4.69, 8.81)
u95wtPrev5 <- c(1.67, 1.02, 0.64, 0.28, 0.58, 1.29, 2.92, 5.44, 5.40, 8.69, 5.48, 9.82)

l95wtPrev12 <- c(1.98, 0.40, 0.13, 0.04, 0.09, 1.21, 1.98, 4.78, 1.91, 4.57, 3.00, 4.21)
wtPrev12    <- c(2.37, 0.57, 0.23, 0.11, 0.17, 1.50, 2.36, 5.35, 2.31, 5.20, 3.42, 4.71)
u95wtPrev12 <- c(2.84, 0.81, 0.40, 0.29, 0.33, 1.85, 2.81, 5.99, 2.80, 5.92, 3.91, 5.25)

# calculate unweighted percentage vaccinated
# percentageVac12 <-  NULL
# for (round_no in (13:19)){
#
#   df <- readRDS(paste0("E:/dt20/linkedR", round_no, "datANG.rds"))
#   df <-df %>%
#     filter(age_group_char3 == "12-17") %>%
#     mutate(vacc_res = ifelse(link_vax_status_booster0days %in% c("1 dose", "2 doses", "3 doses"), 1,
#                              ifelse(link_vax_status_booster0days %in% c("Not vaccinated", "unvaccinated", "Unknown does"), 0, NA))) %>%
#     filter(!is.na(vacc_res))%>% filter(!is.na(estbinres))
#
#   x <- sum(df$vacc_res) / length(df$vacc_res)
#   percentageVac12 <- c(percentageVac12, x)
#
# }


# get weighted percentage vaccinated
wtPercentageVac <- read.csv("E:/Group/report/round19/linked_vacc_prop_under18_r8_r19.csv")

wtPercentageVac12 <- wtPercentageVac %>%
  filter(Category == "12-17")
wtPercentageVac12$vacc_res


# setup dataframe for plotting
df <- as.data.frame(
  cbind(
    wtPrev5[6:12],
    l95wtPrev5[6:12],
    u95wtPrev5[6:12],
    wtPrev12[6:12],
    l95wtPrev12[6:12],
    u95wtPrev12[6:12]))

colnames(df) <- c(
  "weightedPrev5", "lower95wp5", "upper95wp5",
  "weightedPrev12", "lower95wp12", "upper95wp12")

df$ratio12to5 <- df$weightedPrev12 / df$weightedPrev5
df$ratio5to12 <- df$weightedPrev5 / df$weightedPrev12
df$roundVariance5 <- ((df$upper95wp5 - df$lower95wp5)/(2*1.96))^2
df$roundVariance12 <- ((df$upper95wp12 - df$lower95wp12)/(2*1.96))^2

df$lower95ratio <- df$ratio5to12 - 1.96 * sqrt(
  (df$roundVariance5 + df$ratio5to12^2 * df$roundVariance12)/df$weightedPrev12)
df$upper95ratio <- df$ratio5to12 + 1.96 * sqrt(
  (df$roundVariance5 + df$ratio5to12^2 * df$roundVariance12)/df$weightedPrev12)

df <- cbind(df, wtPercentageVac12$vacc_res*100)
names(df)[ncol(df)] <- "wtPercentageVac12"
df <- cbind(df, wtPercentageVac12$ci_l *100 )
df <- cbind(df, wtPercentageVac12$ci_u *100 )
names(df)[ncol(df)-1] <- "lowerwtPercentageVac12"
names(df)[ncol(df)] <- "upperwtPercentageVac12"

df$Round <- as.factor(c(13,14,15,16,17,18,19))
df

# REVISION: 13 May 2022 - only weighted prevalence
ggplot(df %>% filter(!Round %in% c(13,14,15)))+
  geom_bar(stat="identity", aes(x=Round, y=ratio5to12, color=Round, fill=Round), alpha=0.5, width=0.52) +
  geom_errorbar(aes(x=Round, y=ratio5to12, ymin=lower95ratio, ymax=upper95ratio, color=Round), width=0.11, size=0.75) +
  theme_bw() +
  theme(legend.position = c(0.94, 0.84),
        legend.background = element_rect(fill="white", color="black")) +
  ylab("Ratio of weighted prevalence in 5-11\nto weighted prevalence in 12-17 years") +
  xlab("Round Number") +
  theme(panel.grid.major.x = element_blank())

ggsave(paste0("T:/ratio_round_bar", Sys.Date(), ".pdf"), width=10, height=6, units="in")

# REVISION: 13 May 2022 - only weighted prevalence
ggplot(df %>% filter(!Round %in% c(13,14,15)))+
  geom_bar(stat="identity", aes(x=Round, y=ratio5to12, color=Round, fill=Round), alpha=0.5, width=1) +
  geom_point(aes(x=Round, y=ratio5to12, color=Round, fill=Round), alpha=1, size=2) +
  geom_errorbar(aes(x=Round, y=ratio5to12, ymin=lower95ratio, ymax=upper95ratio, color=Round), width=0.11, size=0.75) +
  theme_bw() +
  theme(legend.position = "none") +
  ylab("Ratio of weighted prevalence in 5-11\nto weighted prevalence in 12-17 years") +
  xlab("Round Number") +
  theme(panel.grid.major.x = element_blank())

ggsave(paste0("T:/ratio_round_bar_pointEstimates_ci", Sys.Date(), ".pdf"), width=5.5, height=6, units="in")


# dotplot
# ggplot(df %>% filter(!Round %in% c(13,14,15)))+
#   geom_point(size=1.3, aes(x=Round, y=ratio5to12, color=Round)) +
#   geom_errorbar(aes(x=Round, y=ratio5to12, ymin=lower95ratio, ymax=upper95ratio, color=Round), width=0.15) +
#   theme_bw() +
#   theme(legend.position = c(0.90, 0.78),
#         legend.background = element_rect(fill="white", color="black")) +
#   ylab("Ratio of weighted prevalence in 5-11\nto weighted prevalence in 12-17 years") +
#   xlab("Round Number")


# previous vaccination plot
# ggplot(df)+
#   geom_smooth(method="lm", formula=y~x, aes(x=wtPercentageVac12, y=ratio12to5), color="darkgray") +
#   geom_point(size=1.3, aes(x=wtPercentageVac12, y=ratio12to5, color=Round)) +
#   geom_errorbar(aes(x=wtPercentageVac12, y=ratio12to5, ymin=lower95ratio, ymax=upper95ratio, color=Round)) +
#   geom_errorbar(aes(x=wtPercentageVac12, y=ratio12to5, xmin=lowerwtPercentageVac12, xmax=upperwtPercentageVac12, color=Round)) +
#   theme_bw() +
#   theme(legend.position = c(0.90, 0.78),
#         legend.background = element_rect(fill="white", color="black")) +
#   ylab("Ratio of weighted prevalence in 12-17\nto weighted prevalence in 5-11 years") +
#   xlab("Weighted percentage vaccinated in 12-17 years (%)")
#
# ggsave("E:/Group/report/round19/Figures/weightedRatiovsweightedVaccKids.pdf",
#        width= 6.1, height=5.75, units="in")
#
# write.csv(df, "E:/Group/report/round19/kids_RatioWeightedPrevalence_WeightedVacc.csv", row.names = FALSE)
#
# lmKids <- lm(ratio12to5 ~ wtPercentageVac12, data=df)
# summary(lmKids)
# confint(lmKids)
#
# vcov(lmKids)
# cov2cor(vcov(lmKids))
# cov2cor(vcov(lmKids))[1,2]

# to deprecate
# noci <- ggplot(df)+
#   geom_point(size=2, aes(x=pctVac12, y=ratio12to5)) +
#   # geom_errorbar(aes(x=pctVac12, y=ratio12to5, ymin = l95, ymax = u95)) +
#   geom_smooth(method="lm", formula=y~x, aes(x=pctVac12, y=ratio12to5)) +
#   ylab("Ratio of weighted prevalence in 12-17\nto weighted prevalence in 5-11 years") +
#   xlab("(unweighted) Percentage vaccinated in 12-17 years")
#
# ggarrange(ci, noci, nrow=2)
#
# y <- ratio12to5
# x <- perVac12
#
# prev12 <- c(1, 5, 3, 8, 9 , 2)
#
# var()
#
# plot(x, y)
# abline(lm(y~x))
