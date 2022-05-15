FormatCount <- function(x) {
  if (is.matrix(x)) {
    ncolumns <- ncol(x)
    x <- as.vector(x)
  } else {
    ncolumns <- NULL
  }
  output <- formatC(as.numeric(x), format = "f", digits = 0, big.mark = ",")
  if (!is.null(ncolumns)) {
    output <- matrix(output, ncol = ncolumns)
  }
  print(ncolumns)
  return(output)
}


FormatPrevalence <- function(x, digits = 2, add_symbol = TRUE) {
  if (is.matrix(x)) {
    ncolumns <- ncol(x)
    x <- as.vector(x)
  } else {
    ncolumns <- NULL
  }
  if (add_symbol) {
    output <- paste0(formatC(as.numeric(x) * 100, format = "f", digits = 2), "%")
  } else {
    output <- formatC(as.numeric(x) * 100, format = "f", digits = 2)
  }
  if (any(grepl("NA", output))) {
    output[grep("NA", output)] <- ""
  }
  if (!is.null(ncolumns)) {
    output <- matrix(output, ncol = ncolumns)
  }
  return(output)
}


FormatOR <- function(x, digits = 2) {
  if (is.matrix(x)) {
    ncolumns <- ncol(x)
    x <- as.vector(x)
  } else {
    ncolumns <- NULL
  }
  output <- formatC(as.numeric(x), format = "f", digits = 2)
  if (any(grepl("NA", output))) {
    output[grep("NA", output)] <- ""
  }
  if (!is.null(ncolumns)) {
    output <- matrix(output, ncol = ncolumns)
  }
  return(output)
}


FormatPvalue=function(x, digits=2){
  return(formatC(x, format="e", digits=digits))
}


FormatCI <- function(x, CI = c(" (", ", ", ")")) {
  if (!is.matrix(x)) {
    if (length(x) != 3) {
      stop("Please provide an argument x of length 3 with point estimate, lowerbound and upperbound.")
    }
    x <- matrix(x, ncol = 3)
  }
  output <- matrix(NA, nrow = nrow(x), ncol = 1)
  for (k in 1:nrow(x)) {
    output[k, 1] <- paste0(x[k, 1], CI[1], x[k, 2], CI[2], x[k, 3], CI[3])
  }
  if (any(x[, 1] == "")) {
    output[which(x[, 1] == "")] <- ""
  }
  return(output)
}


ExtendList <- function(x) {
  output <- rep(NA, length(x))
  for (i in 1:length(x)) {
    if (!is.na(x[i])) {
      tostore <- x[i]
    }
    output[i] <- tostore
  }
  return(output)
}


ExtractPrevalence <- function(df_round, covs, covs_names,
                              res_param, weight_params = NULL,
                              weighted = FALSE) {
  if (!weighted) {
    perc <- FALSE
    sig_figs <- 6
    
    # Computing unweighted prevalences
    prev_tables_r15 <- make_tables(
      dat = df_round, covariates = covs, sens = 1, spec = 1, method = "exact",
      result_var = res_param, suffix = "r15", sf = sig_figs, percent = perc
    ) %>%
      bind_rows(.) %>%
      rename(
        "Positive_r15" = "Positive", "Total_r15" = "Total", "Prevalence_r15" = "Prevalence",
        "Lower_r15" = "Lower", "Upper_r15" = "Upper"
      )
    
    # Adding covariate names
    prev_tables_r15[, 1] <- covs_names[prev_tables_r15[, 1]]
    rownames(prev_tables_r15) <- paste0(
      prev_tables_r15[, 1],
      "_",
      prev_tables_r15[, 2]
    )
    
    # Re-formatting table
    prev_tables_r15 <- as.matrix(prev_tables_r15)
    mytable <- cbind(
      prev_tables_r15[, 1:2, drop = FALSE],
      FormatCount(prev_tables_r15[, 3:4, drop = FALSE]),
      FormatCI(FormatPrevalence(prev_tables_r15[, 5:7]))
    )
  } else {
    dclus15g <- svydesign(
      id = as.formula(paste0("~", weight_params["id"])),
      strata = as.formula(paste0("~", weight_params["strata"])),
      weights = as.formula(paste0("~", weight_params["weights"])),
      data = df_round, nest = TRUE
    )
    
    prev_tables_r15 <- NULL
    for (covariate in covs) {
      print(covs_names[covariate])
      # if (all(table(df_round[,res_param], df_round[,covariate])>3)){
      prev_tab_g_r15 <- svyby(as.formula(paste0("~", res_param)),
                              by = as.formula(paste0("~", covariate)),
                              design = dclus15g,
                              FUN = svyciprop, vartype = "ci"
      )
      # } else {
      #   prev_tab_g_r15=matrix(NA, nrow=length(levels(df_round[,covariate])), ncol=4)
      #   prev_tab_g_r15[,1]=levels(df_round[,covariate])
      # }
      tmp <- cbind(rep(covariate, nrow(prev_tab_g_r15)), prev_tab_g_r15)
      colnames(tmp) <- c("Variable", "Category", "Estimate", "Lower", "Upper")
      prev_tables_r15 <- rbind(prev_tables_r15, tmp)
    }
    
    # Adding covariate names
    prev_tables_r15[, 1] <- covs_names[prev_tables_r15[, 1]]
    rownames(prev_tables_r15) <- paste0(
      prev_tables_r15[, 1],
      "_",
      prev_tables_r15[, 2]
    )
    
    # Re-formatting table
    prev_tables_r15 <- as.matrix(prev_tables_r15)
    mytable <- cbind(
      prev_tables_r15[, 1:2, drop = FALSE],
      FormatCI(FormatPrevalence(prev_tables_r15[, 3:5]))
    )
  }
  return(mytable)
}


DropLevels=function(x){
  if (!is.factor(x)){
    stop("Argument 'x' must be a factor.")
  }
  x=factor(x, levels=names(table(x))[table(x)!=0])
  return(x)
}


ProportionalJitter=function(x, nbreaks=40, alpha=0.3){
  x=x[!is.na(x)]
  mybreaks=seq(min(x)-1,max(x)+1, length.out=nbreaks)
  x_cat=cut(x, breaks = mybreaks, labels = (2:length(mybreaks))-1)
  counts=table(x_cat)
  
  jittered=rep(0, length(x))
  for (id_jitter in 1:length(x)){
    jittered[id_jitter]=jitter(0, amount=counts[as.character(x_cat[id_jitter])]/sum(counts))
  }
  
  jittered=jittered/max(abs(jittered))*alpha
  # jittered=scale(jittered)[,1]*alpha
  
  return(jittered)
}


ProportionalJitterContinuous <- function(x, alpha = 0.3) {
  x <- x[!is.na(x)]
  myd <- density(x)
  mydist <- ProportionalDistribution(x = x)
  
  jittered <- rep(0, length(x))
  for (id_jitter in 1:length(x)) {
    jittered[id_jitter] <- runif(n = 1, min = 0, max = mydist[id_jitter])
  }
  
  return(jittered)
}


ProportionalDistribution <- function(x, alpha = 0.3) {
  x <- x[!is.na(x)]
  myd <- density(x)
  
  jittered <- rep(0, length(x))
  for (id_jitter in 1:length(x)) {
    jittered[id_jitter] <- myd$y[which.min(abs(myd$x - x[id_jitter]))] / max(myd$y)
  }
  
  jittered <- jittered / max(abs(jittered)) * alpha
  
  return(jittered)
}

ViolinPlot <- function(mylist, mycolours, xlab = "", ylab = "", bty = "o", ylim = NULL,
                       tick_length = 0.1, segment_lwd = 2, density_lwd = 2, median_lwd = 5,
                       point_cex = 0.5, point_pch = 19) {
  if (is.null(ylim)) {
    ylim <- range(mylist)
  }
  
  plot(NA,
       xlim = c(0, length(mylist) / 2 + 1),
       ylim = ylim,
       xlab = xlab, ylab = ylab,
       las = 1, xaxt = "n",
       bty = bty
  )
  
  mysign <- rep(c(-1, 1), length(mylist) / 2)
  
  set.seed(1)
  for (i in 1:length(mylist)) {
    x <- sort(mylist[[i]])
    if (length(x) > 0) {
      if (i %% 2 != 0) {
        abline(v = ceiling(i / 2), lty = 3, col = "grey")
      }
      
      # Density distribution (bg)
      polygon(
        x = c(
          ceiling(i / 2) + mysign[i] * ProportionalDistribution(x),
          rep(ceiling(i / 2), length(x))
        ),
        y = c(x, rev(x)), border = NA,
        col = lighten(mycolours[i], amount = 0.9)
      )
      
      # Data points
      points(ceiling(i / 2) + mysign[i] * abs(ProportionalJitterContinuous(x)), x,
             pch = point_pch, cex = point_cex,
             col = lighten(mycolours[i], amount = 0.5)
      )
      
      # Density distribution (lines)
      lines(ceiling(i / 2) + mysign[i] * ProportionalDistribution(x), x,
            col = mycolours[i],
            lwd = density_lwd
      )
      
      # Range and median
      segments(
        x0 = ceiling(i / 2), x1 = ceiling(i / 2),
        y0 = min(x), y1 = max(x),
        col = mycolours[i],
        lend = 1,
        lwd = segment_lwd
      )
      segments(
        x0 = ceiling(i / 2),
        x1 = ceiling(i / 2) + mysign[i] * tick_length,
        y0 = min(x), y1 = min(x),
        col = mycolours[i],
        lend = 1,
        lwd = segment_lwd
      )
      segments(
        x0 = ceiling(i / 2),
        x1 = ceiling(i / 2) + mysign[i] * tick_length,
        y0 = max(x), y1 = max(x),
        col = mycolours[i],
        lend = 1,
        lwd = segment_lwd
      )
      segments(
        x0 = ceiling(i / 2),
        x1 = ceiling(i / 2) + mysign[i] * ProportionalDistribution(x)[which.min(abs(x - median(x)))],
        y0 = median(x), y1 = median(x),
        col = mycolours[i],
        lend = 1,
        lwd = median_lwd
      )
    }
  }
}


NonParamTest=function(mylist){
  pvalues=NULL
  for (i in 1:(length(mylist)/2)){
    mytest=kruskal.test(x=c(mylist[[2*i-1]], mylist[[2*i]]), 
                        g=c(rep(1, length(mylist[[2*i-1]])),
                            rep(2, length(mylist[[2*i]]))))
    pvalues=c(pvalues, mytest$p.value)
  }
  return(pvalues)
}


TTest=function(mylist){
  pvalues=NULL
  for (i in 1:(length(mylist)/2)){
    mytest=t.test(x=mylist[[2*i-1]], 
                  y=mylist[[2*i]])
    pvalues=c(pvalues, mytest$p.value)
  }
  return(pvalues)
}


MoodTest=function(mylist){
  pvalues=NULL
  for (i in 1:(length(mylist)/2)){
    mytest=mood.test(x=mylist[[2*i-1]], 
                     y=mylist[[2*i]])
    pvalues=c(pvalues, mytest$p.value)
  }
  return(pvalues)
}

