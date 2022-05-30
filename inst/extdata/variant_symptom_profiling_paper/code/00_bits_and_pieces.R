covid_yesnos <- paste0("symptnowaw_",1:26)
covid_yesnos_firstsymp=paste0("symptfirst_",1:26)
covid_yesnos_month=c(paste0("sympt_any1_",1:4),
                     paste0("sympt_any2_",1:13),
                     paste0("sympt_any3_",1:9)
)

# add symptoms to cov_name_list
cov_name_list[covid_yesnos] <- react1_sympnames[1:26]
cov_name_list[covid_yesnos_firstsymp] <- react1_sympnames[1:26]
cov_name_list[covid_yesnos_month] <- react1_sympnames[1:26]
cov_name_list$one_of_four <- "Classic symptoms"
cov_name_list$estbinres <- "PCR positive at time of survey"
cov_name_list$estbinres_char <- "PCR positive at time of survey"



sympnames_type_df <- data.frame(symptom_code = paste0("symptnowaw_",1:26),
                                first_symptom_code = paste0("symptfirst_",1:26),
                                symptom = react1_sympnames[1:26],
                                symptom_type = c(rep("Smell/taste",2),
                                                 "Respiratory/cardiac",
                                                 "Influenza-like",
                                                 rep("Coryzal (cold-like)", 6),
                                                 "Influenza-like","Other",
                                                 rep("Gastrointestinal",4),
                                                 rep("Respiratory/cardiac", 3),
                                                 "Influenza-like",
                                                 rep("Fatigue", 3),
                                                 rep("Other",2),
                                                 "Influenza-like"
                                )
)


# sympnames_type_df <- data.frame(symptom_code = paste0("symptnowaw_",1:26), 
#                                 first_symptom_code = paste0("symptfirst_",1:26), 
#                                 symptom = react1_sympnames[1:26],
#                                 symptom_type = c(rep("Ear, nose and throat",3),
#                                                  rep("Systemic", 1),
#                                                  rep("Ear, nose and throat", 1),
#                                                  
#                                                  "Other","Other",
#                                                  rep("Gastrointestinal",4),
#                                                  rep("Respiratory/cardiac", 3),
#                                                  "Other",
#                                                  rep("Fatigue", 3),
#                                                  rep("Other",3)
#                                 )
# )


# add additional symptoms
cov_name_list$symptomatic_26 <- "Any of 26 symptoms"
cov_name_list$one_of_four <-  "Any of 4 classic symptoms"
sympnames_type_df <- rbind(sympnames_type_df,c("symptomatic_26","symptomatic_26","Any of 26 symptoms", "Overall"))
# sympnames_type_df <- rbind(sympnames_type_df,c("one_of_four","one_of_four","Any of 4 classic symptoms", "Overall"))


# Barbara's functions for conditional logistc LSS -------------------------


# Functions
ResamplePair <- function(data, pair_id, tau, ...) {
  # Ensuring that pair_id and data are aligned
  pair_id <- pair_id[rownames(data)]
  
  # Sampling pairs to include
  sampled_pairs <- sort(sample(unique(pair_id), size = tau * length(unique(pair_id))))
  
  # Identifying corresponding observation IDs
  s <- unique(which(pair_id %in% sampled_pairs))
  
  return(s)
}

LambdaGridCLogit <- function(xdata, ydata, Lambda_cardinal = 100) {
  # Data preparation as in clogitLasso
  y <- as.vector(ydata)
  x <- xdata[y == 1, ] - xdata[y == 0, ]
  
  # Using internal function from clogitLasso with default parameters
  Lambda <- clogitLasso:::frac(
    epsilon = 1e-04, log = TRUE,
    x = x, n = nrow(x), m = ncol(x),
    nbfraction = Lambda_cardinal
  )
  
  return(Lambda)
}

PenalisedCLogit <- function(xdata, ydata, pair_id, Lambda, ...) {
  # Storing extra arguments
  extra_args <- list(...)
  
  # Extracting relevant extra arguments
  tmp_extra_args <- MatchingArguments(extra_args = extra_args, FUN = clogitLasso::clogitLasso)
  tmp_extra_args <- tmp_extra_args[!names(tmp_extra_args) %in% c("X", "y", "strata", "fraction", "trace")]
  
  # Ensuring that pair ID is matched to correct observation
  pair_id <- pair_id[rownames(xdata)]
  
  # Re-ordering to have case first as required in clogitLasso
  ids <- do.call(c, lapply(split(ydata[, 1], f = pair_id), FUN = function(z) {
    names(sort(z, decreasing = TRUE))
  }))
  xdata <- xdata[ids, ]
  ydata <- ydata[ids, ]
  pair_id <- pair_id[ids]
  
  # Fitting conditional logistic LASSO
  fitLasso <- do.call(clogitLasso::clogitLasso, args = c(
    list(
      X = xdata, y = as.vector(ydata),
      strata = pair_id,
      fraction = Lambda,
      trace = FALSE
    ),
    tmp_extra_args
  ))
  
  # Extracting selection status and beta coefficients
  beta_full <- fitLasso$beta
  selected <- ifelse(beta_full != 0, yes = 1, no = 0)
  rownames(beta_full) <- rownames(selected) <- paste0("s", 0:(nrow(beta_full) - 1))
  colnames(beta_full) <- colnames(selected) <- colnames(xdata)
  
  return(list(selected = selected, beta_full = beta_full))
}




# round_dates -------------------------------------------------------------

round_dates=data.frame(
  stringsAsFactors = FALSE,
       check.names = FALSE,
               Round = c("1","2","3","4","5","6","7",
                       "8","9","10","11","12","13","14","15","16","17",
                       "18","19"),
  Month = c("May 2020","June 2020","July 2020",
                       "August 2020","September 2020","October 2020","November 2020",
                       "January 2021","February 2021","March 2021","April 2021",
                       "May 2021 *","June 2021 *","September 2021 *",
                       "October 2021 *","December 2021 *","January 2022 *",
                       "February 2022 *","March 2022 *"),
  Date_sample_drawn=as.Date(c("29/04/2020","29/05/2020",
                       "10/07/2020","06/08/2020","03/09/2020","21/09/2020",
                       "29/10/2020","07/12/2020","20/01/2021","19/02/2021",
                       "19/03/2021","23/04/2021","26/05/2021","20/08/2021",
                       "21/09/2021","07/11/2021","26/11/2021","17/01/2022",
                       "16/02/2022"),format = "%d/%m/%Y"),
  Date_invitation_letters_sent=as.Date(c("30/07/2020","11/06/2020","14/07/2020","11/08/2020",
                       "08/09/2020","06/10/2020","05/11/2020","10/12/2020",
                       "26/01/2021","02/03/2021","01/04/2021","04/05/2021",
                       "07/06/2021","23/08/2021","04/10/2021","11/11/2021",
                       "06/12/2021","24/01/2022","21/02/2022"),format = "%d/%m/%Y"),
  Fieldwork_dates=c("01/05/20 – 01/06/20",
                       "19/06/20 – 07/07/20","24/07/20 – 11/08/20",
                       "20/08/20 – 08/09/20","18/09/20 – 05/10/20","16/10/20 – 02/11/20",
                       "13/11/20 – 03/12/20","06/01/21 – 22/01/21",
                       "06/02/21 – 23/02/21","13/03/21 – 29/03/21","15/04/21 – 04/05/21",
                       "20/05/21 – 07/06/21","24/06/21 – 12/07/21",
                       "11/09/21 – 27/09/21","19/10/21 – 05/11/21",
                       "23/11/21 – 14/12/21","05/01/22 – 20/01/22","08/02/22 – 01/03/22",
                       "08/03/22 – 31/03/22"))

# Tidy up date table
round_dates=tidyr::separate(data = round_dates,col = "Fieldwork_dates",into = c("Fieldwork_start","Fieldwork_end"),sep = " – ")
round_dates$Fieldwork_start=as.Date(round_dates$Fieldwork_start,format = "%d/%m/%y")
round_dates$Fieldwork_end=as.Date(round_dates$Fieldwork_end,format = "%d/%m/%y")

# Add midpoint
round_dates$Fieldwork_midpoint=round_dates$Fieldwork_start+floor((round_dates$Fieldwork_end-round_dates$Fieldwork_start)/2)



