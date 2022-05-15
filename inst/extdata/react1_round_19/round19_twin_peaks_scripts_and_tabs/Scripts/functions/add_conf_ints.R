#' Updated 08-04-2021 by Haowei Wang to add an argument "method" in add_conf_ints() so we can
#' specify the method for CIs

add_conf_ints <- function(tab, method, poscol = "Detected", negcol = "Not Detected",
                          spec = 1, sens = 1) {

  # browser()
  rtn <- tab
  tmp <- dim(tab)
  nrows <- tmp[1]
  ncols <- tmp[2]
  rowP <- vector(mode = "numeric", length = ncols)
  rowUB <- vector(mode = "numeric", length = ncols)
  rowLB <- vector(mode = "numeric", length = ncols)
  rowP_adj <- vector(mode = "numeric", length = ncols)
  rowUB_adj <- vector(mode = "numeric", length = ncols)
  rowLB_adj <- vector(mode = "numeric", length = ncols)
  rowAll <- vector(mode = "numeric", length = ncols)
  for (i in 1:ncols) {
    tmpbin <- propCI(
      x = tab[poscol, i],
      n = tab[poscol, i] + tab[negcol, i],
      method = method
    )
    rowP[i] <- as.numeric(tmpbin$p)
    rowUB[i] <- as.numeric(tmpbin$upper)
    rowLB[i] <- as.numeric(tmpbin$lower)

    # Peter Diggle's correction
    rowP_adj[i] <- max(0, min(1, (rowP[i] + spec - 1) / (sens + spec - 1)))
    rowUB_adj[i] <- max(0, min(1, (rowUB[i] - (1 - spec)) / (sens + spec - 1)))
    rowLB_adj[i] <- max(0, min(1, (rowLB[i] - (1 - spec)) / (sens + spec - 1)))

    rowAll[i] <- tab[poscol, i] + tab[negcol, i]
  }
  ## rtn[1:2,] <- round(rtn[1:2,])
  rtn <- rbind(rtn,
    all = rowAll, p = rowP, lb = rowLB, ub = rowUB,
    p_adj = rowP_adj, lb_adj = rowLB_adj, ub_adj = rowUB_adj
  )
  t(rtn)
}
