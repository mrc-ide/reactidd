# Save this file as `R/load_example_phe_data.R`

#' Load and format example data
#'
#' @export
#' @return PHE Pillar1&2 case numbers data for England only.

load_example_phe_data <- function(){
  phe <- read.csv(system.file("extdata", "phe_pillar1&2_20220204.csv", package = "reactidd"))
  phe$date <- as.Date(phe$date)
  names(phe)[names(phe) == "newCasesBySpecimenDate"] <- "n_cases"
  phe <- phe[c("date","n_cases")]

  return(phe)
}



