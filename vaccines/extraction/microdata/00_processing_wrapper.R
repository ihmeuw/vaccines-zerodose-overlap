#----HEADER------------------------------------------------------------------------------------------------------------
#' Author:  USERNAME; USERNAME
#' Date:    DATE
#' Hub Explanation: <link>
#' Purpose: Process extracted survey microdata from UbCov
#'         - *00_processing_wrapper.R*: Parent script used to launch child processing scripts
#'         - 01_process.R: Dose determination - for each child and antigen of interest, 
#'                         determine number of doses received
#'         - 02_tabulate_gbd.R: Collapse data by GBD location ids 
#'         - 02_tabulate_lbd.R: Collapse by lat/long or most granular admin unit
#'         - 03_disaggregate_lbd.R: Disaggregate tabulated lbd data by antigen and 
#'                                  resample polygon data (if needed)
#'
#' Inputs:  nid  ->  NID from FILEPATH on which to launch processing
#' Qsub:    system("qsub -l m_mem_free=4G -P proj_geospatial -q all.q -l fthread=3 -l archive=TRUE -N PARENT_vax_processing FILEPATH -e s 
#'                  FILEPATH --action launch --project proj_geospatial --nid 'all' --lbd FALSE --gbd FALSE")
#***********************************************************************************************************************

#----SETUP-------------------------------------------------------------------------------------------------------------
# Set file paths
username        <- 'USERNAME'
core_repo       <- paste0("FILEPATH")
vaccines_repo   <- paste0("FILEPATH")
extraction_root <- "FILEPATH"

# Load packages available to the LBG Singularity image
source('FILEPATH')
packages_to_load_lbd <- c("proto", "findpython", "getopt", "argparse", "data.table", "magrittr", "survey", "parallel", "plyr", "dplyr",
                          "rgeos", "raster", "rgdal", "dismo", "gbm", "foreign", "doParallel", "grid", "gridExtra", "gtools", "ggplot2",
                          "assertthat", "INLA", "seegSDM", "seegMBG", "pacman", "glmnet", "RMySQL", "tictoc", "binom", "sf", "fasterize")
suppressMessages(mbg_setup(packages_to_load_lbd, repos=core_repo))

# Load functions from processing script for use by parent script (job_hold, check_nids)
source(paste0("FILEPATH"))

#**********************************************************************************************************************


#----LOAD ARGUMENTS----------------------------------------------------------------------------------------------------
# Get arguments from qsub or set defaults
parser <- ArgumentParser()
parser$add_argument("--action",   help="'launch' to launch a job to process each NID listed; 'process' to actually run processing", default="launch", type="character")
parser$add_argument("--nid",      help="NID(s) to process", type="character")
parser$add_argument("--project",  help="prod cluster project", default="proj_covariates", type="character")
parser$add_argument("--cores",    help="how many cores?", default=10, type="integer")
parser$add_argument("--date",     help="leave missing", default=as.character(Sys.Date()), type="character")
parser$add_argument("--lbd",      help="logical; process for LBD?", default="FALSE", type="character")
parser$add_argument("--gbd",      help="logical; process for GBD?", default="FALSE", type="character")
parser$add_argument("--resample", help="logical; resample polygon data for LBD?", default="FALSE", type="character")

# Save job arguments as variables in environment
args <- parser$parse_args()
list2env(args, environment())
lbd       <- as.logical(lbd)
gbd       <- as.logical(gbd)
resample  <- as.logical(resample)
if (nid=="all") {
  nids      <- gsub(".csv", "", gsub(".*_\\s*|_.*", "", list.files(file.path(extraction_root, "00_raw")))) %>% unique %>% as.numeric %>% sort
} else { 
  nids      <- as.numeric(strsplit(gsub(",", " ", nid), split=" +")[[1]])
}
message(""); print(as.data.table(args)); message(""); rm(args)

#**********************************************************************************************************************


#----LAUNCH------------------------------------------------------------------------------------------------------------
if (action == "launch") {
  
  # Launch one job per NID
  for (nid in nids) {
    
    # Set date
    date <- format(lubridate::with_tz(Sys.time(), tzone="America/Los_Angeles"), "%Y-%m-%d")
    
    # Launch
    job_name <- paste0("CHILD_", nid, "_vax_processing")
    shell    <- "FILEPATH"
    job      <- paste0("qsub -N ", job_name,
                       " -l m_mem_free=15G",
                       " -q long.q",
                       " -P ", project,
                       " -l fthread=", cores,
                       " -l archive=TRUE",
                       " -o FILEPATH", username, "FILEPATH",
                       " -e FILEPATH", username, "FILEPATH",
                       " ", shell, " -e s ",
                       vaccines_repo, "FILEPATH",
                       " --action process --nid ", nid, " --project ", project, 
                       " --date ", date, " --gbd ", gbd, " --lbd ", lbd, " --resample ", resample)
    system(job); print(job)
  } 
  
  # Check for success/failure
  source("FILEPATH")
  
  # Wait
  message("\nBeginning job_hold()")
  job_hold(paste0("CHILD_", nids, "_vax_processing"))
  
  # Check each job
  message("\nBeginning lapply()")
  invisible(lapply(nids, check_nids))
  
} else if (action == "process") {
  
  # Step 1: process NID (dose determination)
  process_successful <- process(nid)

  # Step 2: launch tabulation
  if (process_successful) {
    
    # GBD tabulation
    if (gbd) {
      source(paste0("FILEPATH"))
      tabulate_gbd(nid)
    }

    # LDB tabulation
    if (lbd) {
      
      # Tabulate
      source(paste0("FILEPATH"))
      tabulate_lbd(nid)
      
      # Disaggregate by antigen and resample
      source(paste0("FILEPATH"))
      disaggregate_and_resample(nid, resample_polys=resample, vaccine_stems='all')
      

    }
  }
}

#**********************************************************************************************************************
