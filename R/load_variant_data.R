# Save this file as `R/load_example_data.R`

#' Load and format example data
#'
#' @export
#' @return An object of class `list` with length 2. The first item is a dataframe of number of positive swabs from the REACT data set by region and time. The second item is the same but for total number of swab tests.
#'

load_variant_data <- function(){
  pos <- read.csv(system.file("extdata", "positive.csv", package = "reactidd"))
  tot <- read.csv(system.file("extdata", "total.csv", package = "reactidd"))

  for(i in seq_len(length(unique(tot$X)))){
    if(unique(tot$X)[i] %in% pos$X == FALSE){
      new_row <- data.frame(unique(tot$X)[i], 0, 0, 0, 0, 0, 0, 0, 0, 0)
      colnames(new_row) <- colnames(pos)
      new_row
      pos<-rbind(pos, new_row)
    }
  }
  pos <- pos[order(as.Date(pos$X)),]
  for(i in seq_len(nrow(pos))){
    pos$England[i] <- sum(pos[i,2:10])
    tot$England[i] <- sum(tot[i,2:10])
  }

  pos$X <- as.Date(pos$X)
  tot$X <- as.Date(tot$X)

  return(list(pos, tot))
}
