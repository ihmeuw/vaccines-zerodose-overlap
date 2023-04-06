# HEADER -------------------------------------------------------------------------
# Author: USERNAME
# Date: DATE
# Project: Vaccines: Data Prep
# Purpose: Combine processed survey microdata with prepped literature extraction to produce MBG-ready inputs
# source("FILEPATH")
# Outputs: pointpoly-positioned csvs (in out_dir) ready for mbg
#********************************************************************************* 

# LOAD PACKAGES ---------------------------------------------------------------
# Set file paths
username       <- 'USERNAME'
core_repo      <- paste0("FILEPATH")
indic_repo     <- paste0("FILEPATH")
vaccines_repo  <- paste0("FILEPATH")
commondir      <- paste0("FILEPATH")

# Load packages
source("FILEPATH")
invisible(load_packages(c("proto", "findpython", "getopt", "argparse", "data.table", "magrittr", "survey", "binom", "parallel", "plyr", "dplyr", "haven",
                          "rgeos", "raster", "rgdal", "dismo", "gbm", "foreign", "doParallel", "grid", "gridExtra", "gtools", "ggplot2",
                          "sf", "fasterize", "assertthat")))

# Load functions
source(paste0("FILEPATH"))
source(paste0("FILEPATH"))
source(paste0("FILEPATH"))
source('FILEPATH')

# Avoid namespace conflicts
str_match <- stringr::str_match


# GLOBAL VARIABLES ------------------------------------------------------------
# Set folder date of literature extraction
lit_data_date <- "DATE"

# Get file locations
outdir                  <- "FILEPATH"
outlier_path            <- "FILEPATH"
lbd_regions_path        <- "FILEPATH"
processed_lit_data_path <- "FILEPATH"

# Load MICS decider results
mics_decision_path <- "FILEPATH"
mics_decisions     <- fread(mics_decision_path)


# HELPER FUNCTIONS ------------------------------------------------------------

# Shorthand "not in" function for readability
'%!in%' <- function(x,y)!('%in%'(x,y))

# Given vaccine prefix, load and rbind all data from resampled folder
load_processed_microdata <- function(vax_prefix, resample = FALSE) {
  # Get folder path 
  folder_path <- ifelse(resample == TRUE, 
                        paste0("FILEPATH"),
                        paste0("FILEPATH"))
  processed_files <- list.files(folder_path, full.names = TRUE)
  # Load processed microdata
  processed_microdata <- data.table()
  for (file in processed_files) {
    processed_microdata <- rbind(processed_microdata, 
                                 fread(file), 
                                 fill = TRUE)
  }
  # Return data
  return(processed_microdata)
}


# Given survey data and decider results, remove rows from survey data which will be replaced by the report data
merge_with_decisions <- function(data, mics_decisions, decider_vax, is_ratio){
  
  # Subset decider data to rows relevant to current vaccine in loop
  vaccines_of_interest        <- unique(mics_decisions$me_name)[grepl(paste(paste0("vacc_", decider_vax), collapse = "|"), 
                                                                      unique(mics_decisions$me_name))]
  mics_decisions_vax_specific <- mics_decisions[me_name %in% vaccines_of_interest, ]
  
  # Multi-dose vaccines must have same swap decision for all doses. For mcv and rcv, which aren't accounted for by the decider, 
  # set swap decision to decision from 1st dose. For ratios between vaccines, set the decision for each vaccine based on the 
  # first dose. If either vaccine is swapped, ratio should be swapped
  mics_decisions_vax_specific <- mics_decisions_vax_specific[order(nid, ihme_loc_id, me_name), ]
  if (!is_ratio) {
    
    # Check that swap/keep decision is the same for all doses of multi-dose vaccines (except mcv and rcv, which are dose-specific)
    multi_dose_decisions_standardized <-
      ifelse((mics_decisions_vax_specific[, .("N_decisions" = length(unique(decision))), by=.(ihme_loc_id, nid)][N_decisions > 1, .N] > 0), 
             FALSE, 
             TRUE)
    if (!multi_dose_decisions_standardized) {
      if (decider_vax %!in% c("mcv", "rcv")) {
        stop(paste0("Conflicting decider results for different doses of ", decider_vax, ". \nStop and fix (should be accounted for in decider code)\n "))
      } else {
        # Fix 1st/2nd dose disagreement for mcv/rcv, which aren't accounted for by the decider
        multi_decision_nids <- mics_decisions_vax_specific[, .("N_decisions" = length(unique(decision))), by=.(ihme_loc_id, nid)][N_decisions > 1, .(ihme_loc_id, nid)]
        for (i in 1:nrow(multi_decision_nids)) {
          location_id <- multi_decision_nids[i, ihme_loc_id]
          svy_id      <- multi_decision_nids[i, nid]
          first_dose_decision <- mics_decisions_vax_specific[ihme_loc_id == location_id & nid == svy_id & grepl("1", me_name), decision]
          mics_decisions_vax_specific[ihme_loc_id == location_id & nid == svy_id, decision := first_dose_decision]
        }
      }
    }
    
    # Reformat decider data to merge with vaccine data
    setnames(mics_decisions_vax_specific, c("ihme_loc_id", "nid"), c("country", "svy_id"))
    mics_decisions_vax_specific[, me_name := NULL]
    mics_decisions_vax_specific <- unique(mics_decisions_vax_specific)
    
    # Merge vaccine data with decider data
    data_with_decider_results <- merge(data, mics_decisions_vax_specific, by=c("country", "svy_id"), all.x=TRUE)
    return(data_with_decider_results)
  } else {
    stop("merge_with_decisions not configured to merge ratio decisions")
  }
}

# Add lit extractions except for nids that already exist
add_lit_extractions <- function(processed_microdata, literature_extractions, decider_vax){
  
  if ("ihme_loc_id" %in% names(literature_extractions)) {
    setnames(literature_extractions, old = "ihme_loc_id", new = "country")
  }
  
  literature_extractions <- subset(literature_extractions, select = names(literature_extractions) %in% names(processed_microdata))
  redundant_surveys      <- unique(literature_extractions$svy_id)[unique(literature_extractions$svy_id) %in% unique(processed_microdata$svy_id)]
  literature_extractions <- literature_extractions[!svy_id %in% redundant_surveys, ]
  processed_microdata    <- rbind(processed_microdata, literature_extractions, fill=TRUE)
  return(processed_microdata)
}

# Load resampled lit extraction based on vaccine and date provided. If no date is provided, use most recent
load_literature_extraction <- function(path, lit_data_date, vaccine_to_run, resample = FALSE) {
  
  # If date is provided, look for folder with appropriate date
  if(!resample) {
    folder_path <- file.path(path, lit_data_date, "not_resampled")
  } else {
    folder_path <- file.path(path, lit_data_date, "resampled")
  }
  if (!file.exists(folder_path)) {
    stop(paste0("Lit data folder does not exist: ", folder_path))
  }
  
  # Identify vaccine-specific resample file at file path
  file <- list.files(folder_path)[grepl(paste0(vaccine_to_run), list.files(folder_path))]
  file <- file[!grepl("pre_resample|ratio", file)]
  if (length(file) == 0) {
    stop("Vaccine not found in resampled file path\n", "Vaccine: ", vaccine_to_run, "\nPath: ", folder_path)
  } else if (length(file) > 1) { stop(paste0("Problem reading in resampled lit data at ", folder_path))
  } else {
    lit_extraction_path      <- file.path(folder_path, file)
    lit_extraction           <- fread(lit_extraction_path)
    lit_extraction$type      <- "report"
    message("||-------- ", lit_extraction_path)
    return(lit_extraction)
  }
}

# Apply outliers to data
apply_outliers <- function(data, outlier_path, outlier_names) {
  
  # HELPER FUNCTIONS-------------------------------------------------------------------------------
  
  # Given data table with unformatted year_id column, format year_id
  format_dates <- function(outlier) {
    
    # Helper function: given data table with unformatted year_id column and symbol, separate year_ids
    # based on symbol and extend data table for each year in range
    parse_date_range <- function(data, symbol) {
      dates_formatted <- data.table()
      for(i in 1:nrow(data)) {
        year_start_end        <- as.numeric(unlist(strsplit(data[i]$year_id, symbol)))
        year_range            <- year_start_end[1]:year_start_end[2]
        row_formatted         <- data[i][rep.int(1, length(year_range))]
        row_formatted$year_id <- year_range
        dates_formatted       <- rbind(dates_formatted, 
                                       row_formatted, 
                                       fill = TRUE)
      }
      return(dates_formatted)
    }
    
    # Helper function: given data table with unformatted year_id column and symbol, separate year_ids
    # based on symbol and extend data table for each year in year_id
    parse_date_specific <- function(data, symbol) {
      dates_formatted <- data.table()
      for(i in 1:nrow(data)) {
        years                 <- as.numeric(unlist(strsplit(data[i]$year_id, symbol)))
        row_formatted         <- data[i][rep.int(1, length(years))]
        row_formatted$year_id <- years
        dates_formatted       <- rbind(dates_formatted, 
                                       row_formatted, 
                                       fill = TRUE)
      }
      return(dates_formatted)
    }
    
    # Parse dates in to no date, formatted dates, and unformatted dates
    outliers_no_date          <- outlier[is.na(year_id) | year_id == ""]
    outliers_date_formatted   <- outlier[!is.na(as.numeric(year_id))]
    outliers_date_unformatted <- outlier[!is.na(year_id) & year_id != "" & is.na(as.numeric(year_id))]
    # Parse unformatted date in to date formatted by "-", "/" and ","
    hyphen_dates_unformatted <- outliers_date_unformatted[grep("-", year_id), ]
    slash_dates_unformatted  <- outliers_date_unformatted[grep("/", year_id), ]
    comma_dates_unformatted  <- outliers_date_unformatted[grep(",", year_id), ]
    # Format each type of unformatted date
    hyphen_dates_formatted <- parse_date_range(hyphen_dates_unformatted, "-")
    slash_dates_formatted  <- parse_date_specific(slash_dates_unformatted, "/")
    comma_dates_formatted  <- parse_date_specific(comma_dates_unformatted, ",")
    # Combine formatted dates
    outliers_date_unformatted_fixed <- rbind(hyphen_dates_formatted, 
                                             slash_dates_formatted, 
                                             comma_dates_formatted)
    # Combine newly formatted dates with dates that were originally formatted correctly 
    dates_formatted <- rbind(outliers_date_formatted, 
                             outliers_date_unformatted_fixed, 
                             outliers_no_date)
    # Return data
    return(dates_formatted)
  }
  
  # OUTLIERING-------------------------------------------------------------------------------------
  
  # Read in outliers
  outlier <- fread(outlier_path)
  
  # Format outlier dates
  outlier <- suppressWarnings(format_dates(outlier))
  
  # Currently, the outliering process isn't configured to outlier subnationally for LBD - only for GBD (feature may be added in the
  # future). Ignore subnationally outliered locations, except subnational locations that are batch outliered. In the GBD outliering
  # logic, batch outliers are removed for all locations associated with the nid, meaning that subnationally batch-outliered nids are 
  # interpreted the same as nationally outliered nids on the GBD side. To be consistent with GBD, remove nids that are
  # outliered subnationally and batch outliered completely.
  subnational_rows                   <- grepl(outlier$ihme_loc_id, pattern = "_")
  non_batch_outlier_rows             <- outlier$batch_outlier != 1
  non_batch_outlier_subnational_rows <- subnational_rows & non_batch_outlier_rows
  outlier                            <- outlier[!non_batch_outlier_subnational_rows, ]
  
  # Subset to relevant outliers
  outlier <- outlier[lbd == 1 & (me_name == "" | me_name %in% outlier_names)]
  
  # Get nids to be outliered for all years (batch outlier) 
  outlier_nids_years_all <- outlier[is.na(year_id) | year_id == "" | batch_outlier == 1, nid]
  
  # Get nids to be outliered year-specifically
  outlier_nids_years_specific <- outlier[!is.na(year_id) & year_id != "", paste0(nid, "_", year_id)]
  
  # Add outlier column to data
  data$outlier <- 0
  
  # Outlier data
  data[svy_id %in% outlier_nids_years_all, outlier := 1]
  data[paste0(svy_id, "_", year_id) %in% outlier_nids_years_specific, outlier := 1]
  
  # Return data
  return(data)
}

# Encapsulated in function due to file path complexity. Due to the logic of the processing pipeline, 
# microdata data drops are recorded as individual files (because surveys are processed individually), 
# whereas report data drops are recorded in master file (because all reports are processed together).
# As a result, saving the drops for a given file requires identifying whether the nid is from microdata
# vs report to identify the correct file path
save_data_drops_from_outlier <- function(vaccine_to_run, data_drops_from_outlier, microdata_drop_table, report_data_drop_table) {
  
  # Round data drops
  data_drops_from_outlier[, n_after_outlier := round(n_after_outlier, digits = 0)]
  
  # For the given antigen, grab either the microdata or report data drop record depending on which source is being used
  antigen_data_from_microdata <- data_drops_from_outlier[type == "microdata", ]
  antigen_data_from_report    <- data_drops_from_outlier[type == "report", ]
  
  # Merge microdata and report data drop records with data drops from outliering
  microdata_drop_table$n_after_outlier   <- NULL
  report_data_drop_table$n_after_outlier <- NULL
  
  antigen_microdata_drop_table <- merge(antigen_data_from_microdata, microdata_drop_table, all.x = TRUE)
  antigen_report_drop_table    <- merge(antigen_data_from_report, report_data_drop_table, all.x = TRUE)
  
  # Combine microdata and report data drop tables to get the complete drop table for the given antigen
  final_data_drop_table <- rbind(antigen_microdata_drop_table, antigen_report_drop_table)

  # Format
  setcolorder(final_data_drop_table, 
              c("nid", "total", "n_with_age_data", "n_0_59_months", "antigen", "n", 
                "n_after_geomatch", "resampled", "n_after_resample", "n_after_outlier", "type"))  
  
  # Save
  folder_path <- file.path("FILEPATH", gsub("-", "_", Sys.Date()))
  file_path   <- file.path(folder_path, paste0(vaccine_to_run, ".csv"))
  message("||-------- ", file_path)
  if(!dir.exists(folder_path)) dir.create(folder_path)
  write.csv(final_data_drop_table, file = file_path, row.names = FALSE)

}

#***********************************************************************************************************************


#----GENERATE CSVS----------------------------------------------------------------------------------------------------

generate_csvs <- function(resample = FALSE) {
  
  # LOAD DATA DROP DATA--------------------------------------------------------------
  # Data drops were recorded during data processing and lit extraction. Load records
  # of data drops and incorporate data dropped from outliering at the end of the function
  
  # Get report data drops
  report_data_drop_path  <- file.path("FILEPATH")
  report_data_drop_table <- fread(report_data_drop_path)
  
  # Get microdata data drops and merge in to one table (recorded survey-specifically due to parallelized processing logic)
  microdata_data_drop_path  <- file.path("FILEPATH", ifelse(resample, "resampled", "not_resampled"))
  microdata_drop_table      <- lapply(list.files(microdata_data_drop_path, full.names = TRUE), fread) %>% rbindlist
  
  # BEGIN LOOP ------------------------------------------------------------------

  ratios <- c("hib3_dpt3_ratio", "mcv2_mcv1_ratio", "pcv3_dpt3_ratio", "rotac_dpt3_ratio") # Not currently producing ratios
  
  # Begin vaccine loop
  message("||-- Begin vaccine loop:")
  for (vaccine_to_run in c("hib", "bcg", "hepb", "pcv", "polio", "yfv", "rota", "rcv", "mcv", "dpt")) {

    # Identify whether or not running ratio
    is_ratio <- grepl("ratio", vaccine_to_run)
    
    # Make table to record drops for outliering, to be saved in the nid-specific data drop files
    data_drops_from_outlier <- data.table(svy_id             = NA_integer_,
                                          antigen            = NA_character_,
                                          n_after_outlier    = NA_integer_)
    
    # Set up the vaccine depending on the "vaccine_to_run" parameter
    message(paste0("||\n||---- ", vaccine_to_run))
    
    vaccine_list 	    	  <- NULL
    graph_vaccines    	  <- NULL
    final_model_vaccines  <- NULL
    graph_vaccine_titles  <- NULL
    
    if ("hepb3_dpt3_ratio" %in% vaccine_to_run) {
      outlier_names = c("vacc_hepb3", "vacc_dpt3")
      decider_vax   = c("hepb", "dpt")
    }
    
    if ("hib3_dpt3_ratio" %in% vaccine_to_run) {
      outlier_names = c("vacc_hib3", "vacc_dpt3")
      decider_vax   = c("hib", "dpt")
    }
    
    if ("pcv3_dpt3_ratio" %in% vaccine_to_run) {
      outlier_names = c("vacc_pcv3", "vacc_dpt3")
      decider_vax   = c("pcv", "dpt")
    }
    
    if ("rotac_dpt3_ratio" %in% vaccine_to_run) {
      outlier_names = c("vacc_rotac", "vacc_dpt3")
      decider_vax   = c("rota", "dpt")
    }
    
    if ("mcv2_mcv1_ratio" %in% vaccine_to_run) {
      outlier_names = c("vacc_mcv1", "vacc_mcv2")
      decider_vax   = c("mcv")
    }
    
    if("dpt3_timeliness_ratio" %in% vaccine_to_run) {
      outlier_names = c("vacc_dpt3")
      decider_vax   = c("dpt")
    }
    
    if ("polio" %in% vaccine_to_run) {
      if(is.null(vaccine_list)) vaccine_list <- list()
      vaccine_list$polio <- add_vaccine(prefix  = "polio", 
                                        title   = "Polio",
                                        doses   = 3,
                                        age_min = 0,
                                        age_max = 59)$polio
      
      vaccines <- names(vaccine_list)
      if(is.null(graph_vaccines))       graph_vaccines        <- character()
      if(is.null(final_model_vaccines)) final_model_vaccines  <- character()
      if(is.null(graph_vaccine_titles)) graph_vaccine_titles  <- character()
      graph_vaccines <- c(graph_vaccines, "polio3_cov")
      final_model_vaccines <- c(final_model_vaccines, "polio3_cov", "polio1_cov", "polio1_3_abs_dropout", "polio1_3_rel_dropout")
      graph_vaccine_titles <- c(graph_vaccine_titles, "Polio3")
      outlier_names <- c("vacc_polio1", "vacc_polio2", "vacc_polio3")
    }
    
    if ("bcg" %in% vaccine_to_run) {
      
      vaccine_list <- add_vaccine(prefix = "bcg", 
                                  title = "BCG",
                                  doses = 1,
                                  age_min = 0,
                                  age_max = 59)
      
      vaccines <- names(vaccine_list)
      
      graph_vaccines <- c("bcg_cov")
      final_model_vaccines <- c("bcg_cov")
      graph_vaccine_titles <- c("BCG")
      outlier_names <- c("vacc_bcg")
      
    }
    
    if ("rcv" %in% vaccine_to_run) {
      
      vaccine_list <- add_vaccine(prefix = "rcv", 
                                  title = "RCV",
                                  doses = 1,
                                  age_min = 0,
                                  age_max = 59)
      
      vaccines <- names(vaccine_list)
      
      graph_vaccines <- c("rcv_cov")
      final_model_vaccines <- c("rcv_cov")
      graph_vaccine_titles <- c("RCV")
      outlier_names <- c("vacc_rcv")
      
    }
    if ("yfv" %in% vaccine_to_run) {
      
      vaccine_list <- add_vaccine(prefix = "yfv", 
                                  title = "YFV",
                                  doses = 1,
                                  age_min = 0,
                                  age_max = 59)
      
      vaccines <- names(vaccine_list)
      
      graph_vaccines <- c("yfv_cov")
      final_model_vaccines <- c("yfv_cov")
      graph_vaccine_titles <- c("YFV")
      outlier_names <- c("vacc_yfv")
      
    }
    if ("mcv" %in% vaccine_to_run) {
      
      vaccine_list <- add_vaccine(prefix = "mcv", 
                                  title = "MCV",
                                  doses = 2,
                                  age_min = 0,
                                  age_max = 59)
      
      vaccines <- names(vaccine_list)
      
      graph_vaccines       <- c("mcv1_cov", "mcv2_cov")
      final_model_vaccines <- c("mcv12_cond", "mcv1_cov", "mcv2_cov")
      graph_vaccine_titles <- c("MCV1", "MCV2")
      outlier_names <- c("vacc_mcv1") 
      vaccine_list$mcv$all_doses <- c("1", "2")
      
    }  
    if ("dpt" %in% vaccine_to_run) {
      if(is.null(vaccine_list)) vaccine_list <- list()
      vaccine_list$dpt <- add_vaccine(prefix = "dpt", 
                                      title = "DPT",
                                      doses = 3,
                                      age_min = 0,
                                      age_max = 59)$dpt
      
      vaccine_list$dpt["cond_doses"] <- c(12)
      vaccine_list$dpt["cond_vaccines"] <- "dpt12_cond"
      
      vaccines <- names(vaccine_list)
      if(is.null(graph_vaccines))       graph_vaccines        <- character()
      if(is.null(final_model_vaccines)) final_model_vaccines  <- character()
      if(is.null(graph_vaccine_titles)) graph_vaccine_titles  <- character()
      graph_vaccines <- c(graph_vaccines, "dpt3_cov")
      final_model_vaccines <- c(final_model_vaccines, "dpt12_cond", "dpt3_cov", "dpt1_cov", "dpt1_3_abs_dropout", "dpt1_3_rel_dropout")
      graph_vaccine_titles <- c(graph_vaccine_titles, "DPT3")
      outlier_names <- c("vacc_dpt1", "vacc_dpt2", "vacc_dpt3")
    }
    if ("pcv" %in% vaccine_to_run) {
      if(is.null(vaccine_list)) vaccine_list <- list()
      vaccine_list$pcv <- add_vaccine(prefix = "pcv", 
                                      title = "PCV",
                                      doses = 3,
                                      age_min = 0,
                                      age_max = 59)$pcv
      
      vaccines <- names(vaccine_list)
      if(is.null(graph_vaccines))       graph_vaccines        <- character()
      if(is.null(final_model_vaccines)) final_model_vaccines  <- character()
      if(is.null(graph_vaccine_titles)) graph_vaccine_titles  <- character()
      graph_vaccines <- c(graph_vaccines, "pcv3_cov")
      final_model_vaccines <- c(final_model_vaccines, "pcv3_cov", "pcv1_cov", "pcv1_3_abs_dropout", "pcv1_3_rel_dropout")
      graph_vaccine_titles <- c(graph_vaccine_titles, "PCV3")
    }
    if ("hib" %in% vaccine_to_run) {
      if(is.null(vaccine_list)) vaccine_list <- list()
      vaccine_list$hib <- add_vaccine(prefix = "hib", 
                                      title = "hib",
                                      doses = 3,
                                      age_min = 0,
                                      age_max = 59)$hib
      
      vaccines <- names(vaccine_list)
      if(is.null(graph_vaccines))       graph_vaccines        <- character()
      if(is.null(final_model_vaccines)) final_model_vaccines  <- character()
      if(is.null(graph_vaccine_titles)) graph_vaccine_titles  <- character()
      graph_vaccines <- c(graph_vaccines, "hib3_cov")
      final_model_vaccines <- c(final_model_vaccines, "hib3_cov", "hib1_cov", "hib1_3_abs_dropout", "hib1_3_rel_dropout")
      graph_vaccine_titles <- c(graph_vaccine_titles, "hib3")
      outlier_names <- c("vacc_hib", "vacc_hib1", "vacc_hib2", "vacc_hib3")
      
    }
    if ("hepb" %in% vaccine_to_run) {
      if(is.null(vaccine_list)) vaccine_list <- list()
      vaccine_list$hepb <- add_vaccine(prefix = "hepb", 
                                       title = "hepb",
                                       doses = 3,
                                       age_min = 0,
                                       age_max = 59)$hepb
      
      vaccines <- names(vaccine_list)
      if(is.null(graph_vaccines))       graph_vaccines        <- character()
      if(is.null(final_model_vaccines)) final_model_vaccines  <- character()
      if(is.null(graph_vaccine_titles)) graph_vaccine_titles  <- character()
      graph_vaccines <- c(graph_vaccines, "hepb3_cov")
      final_model_vaccines <- c(final_model_vaccines, "hepb3_cov", "hepb1_cov", "hepb1_3_abs_dropout", "hepb1_3_rel_dropout")
      graph_vaccine_titles <- c(graph_vaccine_titles, "hepb3")
    }
    if ("rota" %in% vaccine_to_run) {
      if(is.null(vaccine_list)) vaccine_list <- list()
      vaccine_list$rota <- add_vaccine(prefix = "rota", 
                                       title = "ROTA",
                                       doses = 1,
                                       age_min = 0,
                                       age_max = 59)$rota
      vaccine_list$rota$all_titles <- c("rota0", "rotac")
      vaccine_list$rota$cond_vaccines <- c("rotac_cond", "rota0_cond")
      vaccine_list$rota$cond_doses <- c("c", "0")
      vaccine_list$rota$all_doses <- c("0", "c")
      
      vaccines <- names(vaccine_list)
      if(is.null(graph_vaccines))       graph_vaccines        <- character()
      if(is.null(final_model_vaccines)) final_model_vaccines  <- character()
      if(is.null(graph_vaccine_titles)) graph_vaccine_titles  <- character()
      graph_vaccines <- c(graph_vaccines, "rotac_cov")
      final_model_vaccines <- c(final_model_vaccines, "rotac_cov")
      graph_vaccine_titles <- c(graph_vaccine_titles, "rotac")
      outlier_names <- c("vacc_rota", "vacc_rotac", "vacc_rota1", "vacc_rota2", "vacc_rota3")
    }
    
    # Define objects
    if (!grepl("ratio", vaccine_to_run)) {
      vax_prefix  <- vaccine_list[[vaccine_to_run]][["prefix"]]
      vax_doses   <- vaccine_list[[vaccine_to_run]][["all_doses"]]
      cond_vax    <- vaccine_list[[vaccine_to_run]][["cond_vaccines"]]
      decider_vax <- vaccine_to_run
    } else {
      vax_prefix <- vaccine_to_run
    }
    
    # LOAD DATA -----------------------------------------------------------------------
    
    # Load resampled microdata
    message("||------ Load processed survey data")
    processed_microdata_raw <- load_processed_microdata(vax_prefix, resample)
    
    # Merge microdata with decider results (pass copy to preserve original)
    message("||------ Incorporate Decider results")
    processed_microdata_with_decisions <- merge_with_decisions(copy(processed_microdata_raw), mics_decisions, decider_vax, is_ratio)
    
    # Drop rows which will be swapped out for resampled literature extractions. Distinguish between data from reports vs microdata
    removal_decisions            <- c("use_report_extractions", "drop")
    processed_microdata          <- processed_microdata_with_decisions[!decision %in% removal_decisions, ]  
    processed_microdata$decision <- NULL
    processed_microdata$type     <- "microdata"
    
    # Merge microdata with literature extractions
    if (!is_ratio) {
      message("||------ Load processed literature extractions:")
      literature_extractions <- load_literature_extraction(processed_lit_data_path, lit_data_date, vaccine_to_run, resample)
      data                   <- add_lit_extractions(processed_microdata, literature_extractions, vaccine_to_run)
    } else {
      data <- processed_microdata
    }
    
    # APPLY OUTLIERING ----------------------------------------------------------------
    # Remove outliers
    message("||------ Identify outliers")
    data <- apply_outliers(data, outlier_path, outlier_names)
    
    
    # ADD MBG REGIONS -----------------------------------------------------------------
    # Add country regions (e.g. "cssa", "crbn", "essa")
    lbd_regions <- fread(lbd_regions_path)
    lbd_regions <- lbd_regions[, .(country, region)]
    data        <- merge(data, lbd_regions, all.x = TRUE, by = "country")
    
    
    # CALCULATE 0 DOSE CHILDREN--------------------------------------------------------
    
    # Add a rowID common between all derivative data sets (for use in creation of holdouts)
    data$row_id <- seq.int(nrow(data))
    
    # Define list of vaccines to save as separate csv's. 
    # Need _cov for last dose, _cond for intermediate doses, and _cov for
    # any doses that are going to be the ultimate models (for validation)
    if(!grepl("ratio", vaccine_to_run)){
      # Message User
      message("||------ Calculate 0-dose children")
      # Get vaccines to save
      last_dose <- max(vax_doses)
      csv_vaccines <- unique(c(cond_vax, 
                               paste0(vax_prefix, last_dose, "_cov"), 
                               final_model_vaccines))
      if(vaccine_to_run == "rota") { csv_vaccines <- "rotac_cov" }
      
      # Calculate 0-dose children
      if(vax_prefix=="rota")  {data[, rota_dose_0  := N - (rota_dose_1  + rota_dose_2  + rota_dose_3)]}
      if(vax_prefix=="polio") {data[, polio_dose_0 := N - (polio_dose_1 + polio_dose_2 + polio_dose_3)]}
      if(vax_prefix=="dpt")   {data[, dpt_dose_0   := N - (dpt_dose_1   + dpt_dose_2   + dpt_dose_3)]}
      if(vax_prefix=="dpt")   {data[, dpt_dose_0   := N - (dpt_dose_1   + dpt_dose_2   + dpt_dose_3)]}
      if(vax_prefix=="pcv")   {data[, pcv_dose_0   := N - (pcv_dose_1   + pcv_dose_2   + pcv_dose_3)]}
      if(vax_prefix=="hepb")  {data[, hepb_dose_0  := N - (hepb_dose_1  + hepb_dose_2  + hepb_dose_3)]}
      if(vax_prefix=="hib")   {data[, hib_dose_0   := N - (hib_dose_1   + hib_dose_2   + hib_dose_3)]}
      if(vax_prefix=="bcg")   {data[, bcg_dose_0   := N - (bcg_dose_1)]}
      if(vax_prefix=="yfv")   {data[, yfv_dose_0   := N - (yfv_dose_1)]}
      if(vax_prefix=="rcv")   {data[, rcv_dose_0   := N - (rcv_dose_1)]}
      if(vax_prefix=="mcv") {
        data[, mcv_dose_0 := N - mcv_dose_1]
        data[!is.na(mcv_dose_2), mcv_dose_0 := N - (mcv_dose_1 + mcv_dose_2)]
      }
    }
    
    
    # SAVE CSVS ----------------------------------------------------------------------
    
    # Create output directory
    if(!dir.exists(outdir)) {
      dir.create(outdir)
    }
    
    message("||------ Save CSV's")
    setnames(data, c("survey_name", "year_id"), c("source", "year"))
    
    if (!is_ratio){
      
      # Catch if single dose
      if (last_dose == 1) csv_vaccines <- paste0(vax_prefix, last_dose, "_cov")
      
      for (ind in csv_vaccines) {
        
        message(paste0("||-------- ", ind))
        
        # Create csv for the last dose -----------------------------------------
        
        df_temp <- copy(data)
        
        # Extract dose number
        if(ind == "rotac_cov") {
          dose <- 1
          last_dose <- 1
        } else {
          dose <- str_match(ind, "^[a-z]*([0-9]*)_.*")[2] %>% as.numeric 
        }
        
        if (grepl("_cov", ind)) {
          
          # Calculate coverage for dose of interest as sum of dose-specific doses <= dose of interest
          if (vaccine_to_run == "rota"){
            drop_vars <- ""
            names(df_temp)[names(df_temp) == "rota_dose_c"] <- "outcome"
          } else if (vax_prefix == "mcv" & dose == 1) {
            df_temp[, outcome := mcv_dose_1]
            df_temp[!is.na(mcv_dose_2), outcome := rowSums(.SD),
                    .SDcols = paste0(vaccine_to_run, "_dose_", dose:last_dose)] 
          } else {
            df_temp[, outcome := rowSums(.SD),
                    .SDcols = paste0(vaccine_to_run, "_dose_", dose:last_dose)]
          }
          
          # Drop extra columns
          drop_vars <- names(df_temp)[grepl(paste0(vax_prefix, "_dose"), names(df_temp))]
          df_temp   <- subset(df_temp, select = !names(df_temp) %in% drop_vars)
          
        } else if (grepl("_cond", ind)) {
          
          # Change to account for "dpt12_cond". We don't want to create
          # dpt cond for doses 1-12. Consider changing naming convention
          # to something like dpt_1_2_cond to avoid confusion
          if (dose == 12){
            # DPT1/2 conditional is the proportion who received 1 or 2 doses out of the people who didn't receive 3 doses
            if (ind == "dpt12_cond") dose <- 2
            # MCV 1/2 conditional is the proportion of people in regions where mcv2 has been introduced who receive
            # MCV1 given that they didn't receive MCV2. Subset to survey with MCV2 data
            if (ind == "mcv12_cond") {
              df_temp <- df_temp[!is.na(mcv_dose_2), ]
              dose <- 1
            }
            
            # Identify doses up to conditional dose and rename
            setnames(df_temp, paste0(vax_prefix, "_dose_", 0:dose), paste0("d", 0:dose))
            
            # Remove doses greater than conditional dose
            drop_vars <- names(df_temp)[grepl(paste0(vax_prefix, "_dose"), names(df_temp))]
            df_temp   <- subset(df_temp, select = !(names(df_temp) %in% drop_vars))
            
            # Overwrite N with sum of remaining variables
            df_temp[, N := rowSums(.SD), .SDcols = paste0("d", 0:dose)]
            
            # Add dose-specific 1 to dose-specific 2 for 1-2 conditional
            df_temp[, paste0("d", dose) := rowSums(.SD), .SDcols = paste0("d", 1:dose)]
            
          } else {
            # Identify doses up to conditional dose and rename
            setnames(df_temp, paste0(vax_prefix, "_dose_", 0:dose), paste0("d", 0:dose))
            
            # Drop other dose columns (ignoring: conditional)
            drop_vars <- names(df_temp)[grepl(paste0(vax_prefix, "_dose"), names(df_temp))]
            
            # Subset to just the columns of interest
            df_temp <- subset(df_temp, select = !(names(df_temp) %in% drop_vars))
            
            # Overwrite N with sum of remaining variables
            df_temp[, N := rowSums(.SD), .SDcols = paste0("d", 0:dose)]
          }
          
          # Drop extra columns
          df_temp <- subset(df_temp, select = !(names(df_temp) %in% paste0("d", 0:(dose-1))))
          
          # Rename to "outcome" for standardization
          setnames(df_temp, paste0("d", dose), "outcome")
          
          # Drop if N = 0 (none in that row have [dose] or fewer doses)
          df_temp <- subset(df_temp, N > 0)
          
        } else if (grepl("1_3_abs_dropout", ind)) {
          
          setnames(df_temp, paste0(vax_prefix, "_dose_", 0:last_dose), paste0("d", 0:last_dose))
          df_temp[, outcome := d1 + d2]  
          
          # Drop extra columns
          df_temp <- subset(df_temp, select = !(names(df_temp) %in% paste0("d", 0:last_dose)))
          
        } else if (grepl("1_3_rel_dropout", ind)) {
          
          # Get the denominator (P doses > 1)
          df_temp[, denom := rowSums(.SD),
                  .SDcols = paste0(vaccine_to_run, "_dose_", 1:last_dose)]
          
          setnames(df_temp, paste0(vax_prefix, "_dose_", 0:last_dose), paste0("d", 0:last_dose))
          
          # Calculate relative dropout
          df_temp[, outcome := (d1 + d2) / denom]  
          
          # Drop extra columns
          df_temp <- subset(df_temp, select = !(names(df_temp) %in% paste0("d", 0:last_dose)))
          df_temp[, denom := NULL]
          
        }
        
        # Deal with 100% plus coverage
        df_over <- df_temp[as.numeric((outcome/N) > 1.02),]
        if(nrow(df_over) > 0){
          greater_than_100_coverage_nids <- unique(df_over$svy_id)
          message(paste0("||---------- ", "WARNING for ", ind, ": nids have observations with > 100% coverage and a difference of more than 2%."))
          message(paste0("||---------- NIDS: ", paste(greater_than_100_coverage_nids, collapse = ", ")))
        } 
        df_temp[(outcome - N) > .01, outcome := N]
        
        # Drop negatives and missing values
        df_temp <- df_temp[!is.na(outcome), ]
        df_temp <- df_temp[outcome >= 0, ]
        
        # Rename outcome to indicator
        setnames(df_temp, "outcome", ind)
        
        # Write
        if(ind %in% c("mcv1_cov", "bcg_cov", "polio3_cov", "dpt1_cov", "dpt3_cov", "hib3_cov", "hepb_cov", "pcv_cov", "rotac_cov", "bcg1_cov", "hepb3_cov", "pcv3_cov")){
          write.csv(df_temp, paste0(outdir, ind, "_outliers_inc.csv"))
        }
        
        # Drop outliers
        df_temp <- df_temp[outlier == 0, ]
        if (df_temp[is.na(point), .N] > 0) message("data is being saved for which point == NA. Consider removing")
        
        # Record outlier data drops to be saved in data drop table
        data_drops_from_outlier <- rbind(data_drops_from_outlier, 
                                         df_temp[, .("antigen" = ind, "n_after_outlier" = sum(weight * N)), by = c("svy_id", "type")], 
                                         fill = TRUE)
        
        message(paste0("||---------- ", paste0(outdir, ind, ".csv")))
        fwrite(df_temp, paste0(outdir, ind, ".csv"))
      }
    } else {
      # For Ratios, drop outliers, "infinite" ratios (!0/0 or 0/0), and NaN ratios
      message(paste0("||-------- ", vaccine_to_run))
      df_vax <- df_vax[outlier == 0, ]
      df_vax <- df_vax[get(vaccine_to_run) != Inf & get(vaccine_to_run) != -Inf & !is.nan(get(vaccine_to_run)), ]
      write.csv(df_vax, paste0(outdir, vaccine_to_run, ".csv"))
    }
    
    # Save data drops from outliering to the nid-specific data-drop files
    message("||------ Save data drops")
    data_drops_from_outlier <- data_drops_from_outlier[grepl("cov", antigen), ] # Subset to only proper antigens (no conditionals, dropouts, or ratios)
    data_drops_from_outlier[, antigen := gsub("_cov", "", antigen)]             # Change antigen naming convention to match data drop table
    setnames(data_drops_from_outlier, "svy_id", "nid")
    save_data_drops_from_outlier(vaccine_to_run, data_drops_from_outlier, microdata_drop_table, report_data_drop_table)
    
  }
  # Message user
  message("******************************************************\n")
}

#***********************************************************************************************************************

generate_csvs(resample = TRUE)
