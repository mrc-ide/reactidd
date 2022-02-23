### Make prevalence tables
# Creator: Kylie Ainslie
# Last modified: 28/07/2020 11:24 AM

#' Updated 10-12-2020 by Matt Whitaker to avoid function failing when
#' there are zero pos_vals or zero neg_vals
#' Backup of function before this edit is saved in Group/Backup, just in case
#' Updated 15-2-2021 to include a wrapper function for stratified tables

#' Updated 08-04-2021 by Haowei Wang to add an argument "method" in make_tables() so we can
#' specify the method for CIs

## helper function
# quick function to make sure 2dp always shown
specifyDecimal <- function(x, k) trimws(format(round(x, k), nsmall = k))



#' Make prevalence tables
make_tables <- function(dat = dfRes,
                        result_var = "res",
                        pos_val = "1",
                        neg_val = "0",
                        covariates = covs, # character vector naming variables to create tables for
                        categories = NULL, # list of character vector of levels of each variable in covariates
                        sens = 0.844,
                        spec = 0.986,
                        weights = NULL,
                        suffix = NULL,
                        sf = 4,
                        method,
                        for_report = FALSE,
                        write_to_file = FALSE,
                        percent = TRUE) {

  #  browser()
  rtn <- list()

  for (i in 1:length(covariates)) {
    print(paste0("Now processing ", covariates[[i]]))
    # generate table of result_var by covariate[i]
    if (is.null(weights)) {
      tab <- table(pull(dat, result_var), pull(dat, covariates[i]))
    }
    else {
      tab <- round(wtd.table(
        x = pull(dat, result_var),
        y = pull(dat, covariates[i]),
        weights = pull(dat, weights),
        normwt = F,
        na.rm = T,
        na.show = F
      ), 0)
    }


    #' for instances where there are only zeros or ones in the results (ie all pos or all neg)
    if (dim(tab)[1] != 2) {
      nonmissing <- rownames(tab)
      missing <- setdiff(c("0", "1"), nonmissing)
      newtab <- (rep(0, ncol(tab)))
      tab <- rbind(tab, missing = newtab)
      rownames(tab) <- c(nonmissing, missing)
    }

    ### reorder
    tab <- as.table(tab[c("0", "1"), ])



    # calculate prevalence and confidence intervals
    tab_ci <- add_conf_ints(tab,
      poscol = pos_val,
      negcol = neg_val,
      sens = sens,
      spec = spec,
      method
    )

    # define category levels
    if (rownames(tab_ci)[1] == "0") {
      indx <- c(1, 2)
    } else {
      indx <- as.numeric(rownames(tab_ci))
    }
    if (is.null(categories)) {
      cats <- rownames(tab_ci)
    } else {
      cats <- categories[[i]][indx]
    }
    # } else {cats <- categories[[i]]}


    # make pretty data frame for output
    if (for_report) {
      df <- as.data.frame.matrix(tab_ci) %>%
        dplyr::mutate(
          Variable = covariates[i],
          Category = cats,
          p = paste0(specifyDecimal(p * 100, sf), " [", specifyDecimal(lb * 100, sf), "-", specifyDecimal(ub * 100, sf), "]"),
          p_adj = paste0(specifyDecimal(p_adj * 100, sf), " [", specifyDecimal(lb_adj * 100, sf), "-", specifyDecimal(ub_adj * 100, sf), "]")
        ) %>%
        dplyr::select(Variable, Category, `1`, all, p, p_adj) %>%
        dplyr::rename("Positive" = "1", "Total" = "all", "Prevalence" = p, "Prevalence_adjusted" = "p_adj", )
    } else {
      if (percent) {
        df <- as.data.frame.matrix(tab_ci) %>%
          dplyr::mutate(
            Variable = covariates[i],
            Category = cats,
            p_adj = round(p_adj, sf) * 100,
            lb_adj = round(lb_adj, sf) * 100,
            ub_adj = round(ub_adj, sf) * 100
          ) %>%
          dplyr::select(Variable, Category, `1`, all, p_adj, lb_adj, ub_adj) %>%
          dplyr::rename(
            "Positive" = "1", "Total" = "all", "Prevalence" = "p_adj",
            "Lower" = "lb_adj", "Upper" = "ub_adj"
          )
      } else {
        df <- as.data.frame.matrix(tab_ci) %>%
          dplyr::mutate(
            Variable = covariates[i],
            Category = cats,
            p_adj = round(p_adj, sf),
            lb_adj = round(lb_adj, sf),
            ub_adj = round(ub_adj, sf)
          ) %>%
          dplyr::select(Variable, Category, `1`, all, p_adj, lb_adj, ub_adj) %>%
          dplyr::rename(
            "Positive" = "1", "Total" = "all", "Prevalence" = "p_adj",
            "Lower" = "lb_adj", "Upper" = "ub_adj"
          )
      }
    }
    rtn[[i]] <- df


    # write df to csv file (for easy incorporation into rmarkdown report)
    if (write_to_file) {
      if (!is.null(suffix)) {
        write.csv(df, paste0(covariates[i], "_prev_", suffix, ".csv"), row.names = FALSE)
      } else {
        write.csv(df, paste0(covariates[i], "_prev.csv"), row.names = FALSE)
      }
    }
  }

  names(rtn) <- covariates
  return(rtn)
}




### Wrapper function to create tables as above but stratified by some covariate
makeStratifiedTables <- function(dat = dfRes,
                                 result_var = "res",
                                 strat_var = "sex",
                                 pos_val = "1",
                                 neg_val = "0",
                                 covariates = covs, # character vector naming variables to create tables for
                                 categories = NULL, # list of character vector of levels of each variable in covariates
                                 sens = 0.844,
                                 spec = 0.986,
                                 weights = NULL,
                                 suffix = NULL,
                                 sf = 4,
                                 for_report = FALSE,
                                 write_to_file = FALSE,
                                 percent = TRUE) {

  ### Get unique levels of a variable
  uniqueLevels <- unique(pull(dat, strat_var))
  uniqueLevels <- uniqueLevels[!is.na(uniqueLevels)]

  ### Loop round creating tables for each one
  res_list <- list()
  res_list_compressed <- list()

  for (i in 1:length(uniqueLevels)) {
    tempdat <- dat[dat[, strat_var] == uniqueLevels[[i]] &
      !is.na(dat[, strat_var]), ]
    tab <- make_tables(
      dat = tempdat,
      result_var = result_var,
      covariates = covariates,
      sens = sens,
      spec = spec,
      weights = weights,
      suffix = suffix,
      sf = sf,
      for_report = for_report,
      write_to_file = write_to_file,
      percent = percent
    )
    tabs_df <- do.call(rbind, tab)
    tabs_df$Level <- uniqueLevels[[i]]
    tabs_df_compressed <- data.frame(tabs_df[, 1:3],
      Total = tabs_df$Total,
      Prevalence = paste0(
        tabs_df$Prevalence, " [",
        tabs_df$Lower, "-", tabs_df$Upper,
        "]"
      ),
      level = tabs_df$Level
    )
    res_list[[i]] <- tabs_df
    res_list_compressed[[i]] <- tabs_df_compressed
  }
  ### Join data frames and do a bit of wrangling
  res_df <- do.call(rbind, res_list)
  res_df_compressed <- do.call(rbind, res_list_compressed)
  res_df_compressed <- res_df_compressed %>% pivot_wider(
    id_cols = c(Variable, Category),
    names_from = level, values_from = c(Positive, Total, Prevalence)
  )


  return(list(
    df_for_report = res_df_compressed,
    df_for_plotting = res_df
  ))
}
