# Author: Kylie Ainslie

# Update by Haowei Wang: enable specify method of calculating CIs
# Update by Barbara Bodinier: enable choice of outcome variable and application on factors

overall_prev <- function(dat, method, outcome="estbinres") {
  x=as.numeric(as.character(dat[,outcome]))
  overPrev <- propCI(
    x = sum(x, na.rm = TRUE),
    n = sum(x == 1, na.rm = TRUE) + sum(x == 0, na.rm = TRUE),
    method = method
  )

  rtn <- data.frame(
    "Positive" = overPrev$x,
    "Total" = overPrev$n,
    "Prevalence" = overPrev$p,
    "Lower" = overPrev$lower,
    "Upper" = overPrev$upper
  )
  return(rtn)
}
