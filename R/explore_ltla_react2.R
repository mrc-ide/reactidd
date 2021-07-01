if (false) {

    ## Clear out the memory
    rm(list=ls(all=TRUE))

    ## Load the data
    df_ab_prev <- read.csv("../inst/extdata/ltla_age_vax_ab_prev.csv")

    ## Make an initial histogram
    hist(df_ab_prev$Overall.n..sampled.in.LTLA)
    hist(df_ab_prev$Overall.n..antibody.positive,breaks = seq(0,2000,10))
    hist(df_ab_prev$Overall.weighted.and.adjusted.antibody.prevalence,
         breaks = seq(0,100,5),
         main = "High variance of prevalence of antibody positivity by LTLA",
         xlab = "Weighted adjusted prevalence",
         ylab = "Frequency of LTLAs")

    ## Load up LTLA case data from UK gov dashboard
    ## install_github("c97sr/idd")
    ## library(idd)
    x <- idd::load_uk_cov_data(at="ltla")

    head(x)
    head(df_ab_prev)

    ## need to merge on areaCode and LTLA.name once have made a
    ## growth rate calculation for the dashboard data
    
}
