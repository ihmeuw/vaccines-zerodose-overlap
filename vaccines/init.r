#----HEADER-------------------------------------------------------------------------------------------------------------
# Author:  USERNAME
# Date:    DATE
# Purpose: Prep vaccination code, set repo and paths
# Run:     source("FILEPATH")
#***********************************************************************************************************************


#----PREP---------------------------------------------------------------------------------------------------------------
# Set OS flexibility
os       <- .Platform$OS.type
username <- 'USERNAME'
system   <- Sys.info()["sysname"]
if (system=="Linux") {
  j <- "FILEPATH"
  j_root <- "FILEPATH"
  h <- paste0("FILEPATH")
  h_root <- paste0("FILEPATH")
} else if (system=="Darwin") { 
  j <- "FILEPATH"
  j_root <- "FILEPATH"
  h <- "FILEPATH"
  h_root <- "FILEPATH"
} else { 
  j <- "FILEPATH"
  j_root <- "FILEPATH"
  h <- "FILEPATH"
  h_root <- "FILEPATH"
}

# Load packages
library(data.table)
library(dplyr)
library(magrittr)

# Path locals
ref_data_repo <- "FILEPATH"
setwd(ref_data_repo)
paths <- fread("paths.csv", header=TRUE)
paths <- paths[obj=="user_root", path5 := paste0("FILEPATH", username)]  
fwrite(paths, file=file.path(ref_data_repo, "paths.csv"))  # this line should save any username changes, so that when run path_loader() function it takes effect
source(paths[obj=="ubcov_tools", 2, with=FALSE] %>% gsub("FILEPATH", j, .) %>% unlist)
path_loader("paths.csv")
source(ubcov_tools)
if (username != "USERNAME") source(db_tools)  
if (os != "windows") source(file.path( "FILEPATH"))

### function to clean paths referenced in other scripts
path.clean <- function(path) {
  if (os=="windows") {
    path <- gsub("FILEPATH")
    path <- gsub(paste0("FILEPATH"), h, path)
  } else {
    path <- gsub("FILEPATH", j, path)
    path <- gsub("FILEPATH", h, path)
  }
  
  return(path)
}

# Set working directory to code_root as original 'FILEPATH' intended
setwd(code_root)

#***********************************************************************************************************************