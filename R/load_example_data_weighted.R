# Save this file as `R/load_example_data.R`

#' Load and format example data
#'
#' @export
#' @return An object of class `list` with length 2. The first item is a dataframe of weighted number of positive swabs from the REACT data set by region and time. The second item is the same but for weighted total number of swab tests.
#'

load_example_data_weighted <- function(){
  pos <- read.csv(system.file("extdata", "positive_weighted.csv", package = "reactidd"))
  tot <- read.csv(system.file("extdata", "total_weighted.csv", package = "reactidd"))


  for(i in seq_len(nrow(pos))){
    pos$England[i] <- sum(pos[i,2:10])
    tot$England[i] <- sum(tot[i,2:10])
  }

  pos$X <- as.Date(pos$X)
  tot$X <- as.Date(tot$X)

  return(list(pos, tot))
}
