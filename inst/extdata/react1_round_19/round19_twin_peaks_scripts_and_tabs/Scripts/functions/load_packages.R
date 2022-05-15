
### This function takes a vector of package names, installs any that aren't currently installed, then loads them all ###

load_packages <- function(package.list) {
  new.packages <- package.list[!(package.list %in% installed.packages()[, "Package"])]
  print(paste(length(new.packages), "packages require installation. Installing now"))
  if (length(new.packages)) install.packages(new.packages, repos = "http://se-r.sm.med.ic.ac.uk", dependencies = TRUE)
  print("Loading packages")
  lapply(package.list, require, character.only = TRUE)
}
