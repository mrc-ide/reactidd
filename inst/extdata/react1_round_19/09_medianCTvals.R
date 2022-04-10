library(tidyverse)

df1 <- readRDS("E:/Group/saved_objects/rep18_lineage.rds")
df2 <- readRDS("E:/Group/saved_objects/rep17_lineage.rds")

df1$round <- 18
df2$round <- 17

df2$agprev5_new <- format(as.Date(df2$agprev5_new, format="%m/%d/%Y"), "%Y-%m-%d")
df1[setdiff(names(df2), names(df1))] <- NA
df2[setdiff(names(df1), names(df2))] <- NA


df3 <- rbind(df1,df2)
table(df3$react_lineage)

write.csv(x = df3, file = "E:/dt20/data/rep17_18_lineage.csv", row.names = FALSE)
t1 <- read.csv("E:/dt20/data/rep17_18_lineage.csv")

t2 <- read.csv("E:/Incoming _data/AG/Antigen18/AntigenR18_FINAL File_INTUSE_06-MAR-22.csv")

# calculate the median and iqr
# KRW test

getMedians=function(df){

  medianList <- NULL
  kruskalList <- NULL

  for (ctType in c("ct1","ct2")) {

    print(ctType)

    ctList <- NULL
    krct   <- NULL

    ba1 <- df %>%
      filter(react_lineage == "BA.1")
    ba11 <- df %>%
      filter(react_lineage == "BA.1.1")
    ba2 <- df %>%
      filter(react_lineage == "BA.2")

    ctList <-
      rbind(
        c(median(ba1[[ctType]]),
          quantile(ba1[[ctType]], 0.25),
          quantile(ba1[[ctType]], 0.75)),
        c(median(ba11[[ctType]]),
          quantile(ba11[[ctType]], 0.25),
          quantile(ba11[[ctType]], 0.75)),
        c(median(ba2[[ctType]]),
          quantile(ba2[[ctType]], 0.25),
          quantile(ba2[[ctType]], 0.75))
      )

    print(kruskal.test(
      c(ba1[[ctType]],ba11[[ctType]]),
      rep(c(1,2), c(length(ba1[[ctType]]), length(ba11[[ctType]]))))
    )
    print(kruskal.test(
      c(ba1[[ctType]],ba2[[ctType]]),
      rep(c(1,2), c(length(ba1[[ctType]]), length(ba2[[ctType]]))))
    )
    print(kruskal.test(
      c(ba11[[ctType]],ba2[[ctType]]),
      rep(c(1,2), c(length(ba11[[ctType]]), length(ba2[[ctType]]))))
    )

    medianList <- cbind(medianList, ctList)
  }

  return(medianList)
}


