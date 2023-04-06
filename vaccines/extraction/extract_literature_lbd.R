#----HEADER------------------------------------------------------------------------------------------------------------
# Author:  USERNAME
# Date:    DATE
# Purpose: This script is responsible for interpretting and preparing all literature and report extractions for LBD. The data
#          is read from the Google sheet and can optionaly be resampled. Later, the prepared literature data is incorporated 
#          into the larger mbg-ready dataset by generate_csvs.R. Report data is incorporated if there is no survey microdata, 
#          or if the survey microdata differs from the report data by >10% (as decided by the Decider).
# Run:     source("FILEPATH")
# Qsub:    qsub -l m_mem_free=6G -P proj_geospatial -q all.q -l fthread=8 -l archive=TRUE -N Lit_Extraction 
#          FILEPATH -e s FILEPATH
#***********************************************************************************************************************

#----ENVIRONMENT-------------------------------------------------------------------------------------------------------
# Set file paths
username        <- "USERNAME"
core_repo       <- paste0("FILEPATH")
vaccine_repo    <- paste0("FILEPATH")
extraction_root <- "FILEPATH"

# Load packages available to the LBG Singularity image
source('FILEPATH')
mbg_setup(c("proto", "findpython", "getopt", "argparse", "data.table", "magrittr", "survey", "parallel", "plyr", "dplyr", "haven",
            "rgeos", "raster", "rgdal", "dismo", "gbm", "foreign", "doParallel", "grid", "gridExtra", "gtools", "ggplot2",
            "sf", "assertthat", "INLA", "seegSDM", "seegMBG", "pacman", "glmnet", "RMySQL", "tictoc"), repos=core_repo)
message("Successfully loaded packages")

# Get needed functions from processing scripts
source(paste0("FILEPATH"))
source(paste0("FILEPATH"))

#----GLOBAL VARIABLES--------------------------------------------------------------------------------------------------
# Set Global Variables
lit_extraction_key         <- "FILEPATH"
lit_extraction_sheet_name  <- FILEPATH
outlier_filepath           <- "FILEPATH"
vax_info_table_path        <- "FILEPATH"
lit_extraction_output_root <- "FILEPATH"

# Make vaccine rename table to distinguish stem from vaccine from dose-specific vaccine
vax_info     <- fread(vax_info_table_path)
rename_table <- vax_info[, .("stem"          = vaccine_stem, 
                             "vaccine"       = paste0("vacc_", vaccine), 
                             "dose_specific" = paste0(vaccine_stem, "_dose_", dose))]

#----HELPER FUNCTIONS -------------------------------------------------------------------------------------------------

# Shorthand "not in" function for readability
'%!in%' <- function(x,y)!('%in%'(x,y))

# Get literature data
gs_download <- function(key, sheet_name=NULL) {
  if (!is.null(sheet_name)) {
    sheet <- paste0("&gid=", sheet_name)
  } else {
    sheet <- ""
  }
  google_sheet_url <- paste0("URLPATH", sheet)
  google_sheet     <- fread(google_sheet_url)
  return(google_sheet)
}

# Rename columns if they are present
check_and_rename_col <- function(data, old_names, new_names) {
  # Confirm argument lengths are same
  if (length(old_names) != length(new_names)) {
    stop("Length of 'old_names' argument must equal length of 'new_names' argument")
  } 
  # Rename column names where present
  for (i in 1:length(old_names)) {
    if (old_names[i] %in% names(data)) {
      setnames(data, old_names[i], new_names[i])
    }
  }
}

# Subset to lbd reports in correct year range with proper geodata, distribute multi-vacc, save data drops for publication flowcharts 
prepare_report_data <- function(report_data){
  
  # FORMAT DATA----------------------------------------------------------------
  # Message user
  message("||-- Format Data")
  
  # Make copy to preserve original (data tables are passed by reference)
  data <- copy(report_data)
  
  # Rename columns
  column_names <- c("sample_size"     = "N_obs", 
                    "nid"             = "svy_id", 
                    "parent_location" = "ihme_loc_id")
  check_and_rename_col(data      = data, 
                       old_names = names(column_names), 
                       new_names = unname(column_names))
  
  # Subset columns
  keep_cols       <- c("svy_id", "file_path", "ihme_loc_id", "cv_admin", "location_code", "location_name_short_ihme_loc_id", 
                       "shapefile", "smaller_site_unit", "site_memo", "year_start", "year_end", "age_start", "age_end", "N_obs")
  vaccine_columns <- names(data)[grepl("vacc", names(data))]
  data            <- subset(data, select = c(keep_cols, vaccine_columns))
  
  # Set missing values to NA
  data[location_code  == "", location_code := NA]
  data[shapefile      == "", shapefile := NA]
  data[N_obs          == "", N_obs := NA]
  data[svy_id         == "", svy_id := NA]
  data[site_memo      == "", site_memo := NA]
  data[location_name_short_ihme_loc_id == "", location_name_short_ihme_loc_id := NA]
  
  # Standardize data
  data$N_obs       <- suppressWarnings(as.numeric(as.character(data$N_obs)))
  data$survey_name <- "Custom"
  data$point       <- 0
  data[, year_id := floor((year_start + year_end)/2) - (floor(age_end/12))]

  # Get vaccine columns and set data type to numeric
  for (column in vaccine_columns) {
    suppressWarnings(data[, (column) := as.numeric(get(column))])
  }
  
  # Fix any sample sizes separated with commas
  data[grepl(",", N_obs), N_obs := as.numeric(gsub(",", "", N_obs))]
  
  # Drop administrative data
  data <- data[!which(cv_admin == 1), ]
  
  # Drop data that start before year 2000
  data <- data[!is.na(year_end) & year_end > 1999, ]
  
  # Drop data with no observations (sample size can't be drawn if missing)
  data <- data[!is.na(N_obs), ]
  
  # DISTRIBUTE MULTI-DOSE------------------------------------------------------
  # Message user
  message("||-- Distribute Multi-Vaccines")
          
  # If have higher combination vaccine coverage, replace component vaccine doses with dose from combination vaccine
  duples <- list(c("vacc_dpt1",  "vacc_tetra1"),
                 c("vacc_dpt1",  "vacc_pent1"),
                 c("vacc_dpt2",  "vacc_pent2"),
                 c("vacc_dpt2",  "vacc_tetra2"),
                 c("vacc_dpt3",  "vacc_tetra3"),
                 c("vacc_dpt3",  "vacc_pent3"),
                 c("vacc_hib1",  "vacc_pent1"), 
                 c("vacc_hib1",  "vacc_tetra1"),
                 c("vacc_hib2",  "vacc_pent2"), 
                 c("vacc_hib2",  "vacc_tetra2"),
                 c("vacc_hib3",  "vacc_pent3"),
                 c("vacc_hib3",  "vacc_tetra3"),
                 c("vacc_hepb1", "vacc_pent1"),
                 c("vacc_hepb2", "vacc_pent2"),
                 c("vacc_hepb3", "vacc_pent3"),
                 c("vacc_mcv1",  "vacc_mmr1"),
                 c("vacc_mcv2",  "vacc_mmr2"),
                 c("vacc_rcv1",  "vacc_mmr1"),
                 c("vacc_rcv2",  "vacc_mmr2"))
  for (i in duples) {
    if(i[2] %in% names(data)){
      data <- data[, (i[1]) := ifelse((!is.na(get(i[2])) & get(i[2]) > get(i[1])) | (!is.na(get(i[2])) & is.na(get(i[1]))), get(i[2]), get(i[1]))]
    }
  }
  
  # RECORD DATA DROPS----------------------------------------------------------
  #' Preserve number of children, by survey, before and after each data drop, to be used in manuscript
  #' flowcharts. Processed microdata data drops are preserved in the `microdata` folder, while report
  #' data drops are preserved in the `report` folder. Depending on which version of a given survey
  #' is used (based on Decider results), `generate_csvs.R` will move the proper data drop record in to
  #' the `final_folder`.
  #' Explanation: Many report extractions have overlapping age bins (ie 12-23 months and 12-59 months). 
  #'              The sample size of the reports with overlapping age bins relies on the largest age bin.
  #'              Sample sizes from surveys without overlapping age bins is taken as the sum of sample sizes
                
  # Message user
  message("||-- Record Data Drops")
  
  # Identify admin1 vs admin0-lvl extractions
  subnational_ihme_locations <- lapply(data$location_name_short_ihme_loc_id, function(x) {
    if(!is.na(x)) {
      split_name <- strsplit(x, "\\|")[[1]][2]
      if(grepl("[0-9]", split_name)) {
        return(x)
      } else {
        return(NULL)
      }
    } else {
      return(NULL)
    }
  }) %>% unlist
  data[location_name_short_ihme_loc_id %in% subnational_ihme_locations | !is.na(location_code) | !is.na(shapefile) | !is.na(site_memo), admin_level := 1]
  data[is.na(admin_level) | admin_level != 1, admin_level := 0]
  
  # Drop admin0-level data
  data <- data[admin_level == 1, ]

  # Identify surveys with overlapping age bins to avoid double counting kids in sample size
  temp       <- unique(data[, .(svy_id, location_code, shapefile, site_memo, year_start, year_end, age_start, age_end)])
  age_range  <- data.table()
  row_length <- nrow(temp)
  for(i in seq.int(row_length)) {
    # Message user
    if(i %% 100 == 0) message(paste0("||---- ", round(i / row_length, digits = 2) * 100, "%"))
    nid                             <- temp[i, svy_id]
    location_code                   <- temp[i, location_code]
    shapefile                       <- temp[i, shapefile]
    site_memo                       <- temp[i, site_memo]
    year_start                      <- temp[i, year_start]
    year_end                        <- temp[i, year_end]
    age_start                       <- temp[i, age_start]
    age_end                         <- temp[i, age_end]
    if(!is.na(age_start) & !is.na(age_end)) {
      ages       <- seq.int(age_start, age_end)
    } else {
      ages       <- NA
    }
    age_range <- rbind(age_range, 
                       data.table(svy_id        = nid, 
                                  location_code = location_code, 
                                  shapefile     = shapefile, 
                                  site_memo = site_memo,
                                  year_start    = year_start, 
                                  year_end      = year_end, 
                                  age_start     = age_start, 
                                  age_end       = age_end, 
                                  age_range     = ages), 
                       fill = TRUE)
  }
  age_overlap <- age_range[, .("age_overlap" = any(duplicated(age_range))), by = .(svy_id, location_code, shapefile, site_memo, year_start, year_end)]
  svy_ids_with_age_overlap    <- age_overlap[age_overlap == TRUE, unique(svy_id), ]
  svy_ids_with_no_age_overlap <- age_overlap[svy_id %!in% svy_ids_with_age_overlap, unique(svy_id)]
  data <- merge(data, age_overlap, all.x = TRUE)
  
  # For survey/location/years with overlapping age bins, preserve only the largest age bin to calculate sample size 
  data_svys_with_age_overlap <- data[svy_id %in% svy_ids_with_age_overlap, ]
  data_svys_with_age_overlap[, age_span := age_end - age_start]
  largest_age_span <- data_svys_with_age_overlap[, .("largest_age_span" = max(age_span)), by = c("svy_id", "shapefile", "location_code", "site_memo", 
                                                                                                 "age_start", "age_end", "year_start", "year_end")]
  data_age_overlap <- merge(data_svys_with_age_overlap, largest_age_span, all.x = TRUE, by = c("svy_id", "shapefile", "location_code", "site_memo", 
                                                                                               "age_start", "age_end", "year_start", "year_end"))
  data_age_overlap_largest_only <- data_age_overlap[age_overlap == FALSE | (age_overlap == TRUE & age_span == largest_age_span), ]
  data_age_overlap_largest_only$age_span <- NULL
  data_age_overlap_largest_only$largest_age_span <- NULL
  
  # After dropping the smaller age bins from survey/location/years with overlapping age bins, combine with data from surveys without overlapping age bins
  data_no_age_overlap        <- data[svy_id %in% svy_ids_with_no_age_overlap, ]
  age_overlap_corrected_data <- rbind(data_age_overlap_largest_only, data_no_age_overlap)
  
  # Record sample sizes at each stage of data drops
  total           <- age_overlap_corrected_data[, .("total" = sum(N_obs)), by = svy_id] 
  n_with_age_data <- age_overlap_corrected_data[!is.na(age_start) & !is.na(age_end), .("n_with_age_data" = sum(N_obs)), by = svy_id] 
  n_age_0_59_data <- age_overlap_corrected_data[!is.na(age_start) & !is.na(age_end) & age_end <= 60, .("n_0_59_months" = sum(N_obs)), by = svy_id] 
  
  # Reshape wide to long to get antigen-specific data drops
  id_variables       <- c("svy_id", "file_path", "ihme_loc_id", "location_code", "location_name_short_ihme_loc_id", "shapefile", "site_memo", "year_start", "year_end", "age_start", "age_end", "year_id", "N_obs")
  data_long          <- melt.data.table(age_overlap_corrected_data, id.vars=id_variables, measure.vars = vaccine_columns, variable.name="me_name", variable.factor = F, value.name = "value")
  data_long          <- data_long[!is.na(value), ]
  n_antigen_specific <- data_long[!is.na(age_start) & !is.na(age_end) & age_end <= 60, .("n" = sum(N_obs, na.rm = TRUE)), by = c("svy_id", "me_name")]
  n_after_geomatch   <- data_long[!is.na(age_start) & !is.na(age_end) & age_end <= 60 & !is.na(shapefile) & !is.na(location_code), .("n_after_geomatch" = sum(N_obs, na.rm = TRUE)), by = c("svy_id", "me_name")] 

  # Merge to get drop table
  data_drop_table <- merge(merge(merge(total, n_with_age_data, all = TRUE), 
                                 merge(n_age_0_59_data, n_antigen_specific, all = TRUE), all = TRUE), 
                           n_after_geomatch, by = c("svy_id", "me_name"), all = TRUE)
  
  # Subset to only antigens that are present 
  data_drop_table <- data_drop_table[!is.na(n), ]
  
  # Set missings to 0's
  data_drop_table[!is.na(total) & is.na(n_with_age_data), n_with_age_data := 0]
  data_drop_table[!is.na(n_with_age_data) & is.na(n_0_59_months), n_0_59_months := 0]
  data_drop_table[!is.na(n_0_59_months) & is.na(n), n := 0]
  data_drop_table[!is.na(n) & is.na(n_after_geomatch), n_after_geomatch := 0]

  # Add columns for recording resample data drops
  setcolorder(data_drop_table, c("svy_id", "total", "n_with_age_data", "n_0_59_months", "me_name", "n", "n_after_geomatch"))
  data_drop_table <- cbind(data_drop_table, 
                           data.table(resampled            = FALSE,
                                      n_after_resample     = NA_integer_,
                                      n_after_outlier      = NA))
  
  # Formatting changes to match format of data drop table from microdata
  data_drop_table[, me_name := gsub("vacc_", "", me_name)]
  setnames(data_drop_table, c("svy_id","me_name"),  c("nid", "antigen"))
  
  data_drop_table_archive_path <- paste0("FILEPATH")
  data_drop_table_regular_path <- paste0("FILEPATH")
  
  message(paste0("||---- Saving Drops: ", data_drop_table_archive_path))
  write.csv(data_drop_table, file = data_drop_table_archive_path, row.names = FALSE)
  write.csv(data_drop_table, file = data_drop_table_regular_path, row.names = FALSE)
  
  
  # RETURN DATA----------------------------------------------------------------
  return(data)
}

# Add vaccine-dpt3 ratios
add_ratios <- function(data, ratio_antigen_1, ratio_antigen_2){
  if (!all(c(ratio_antigen_1, ratio_antigen_2) %in% names(data))) {
    stop("One of the provided ratio_antigens in add_ratios is missing from the extracted report data")
  }
  ratio_name <- gsub("vacc_", "", paste0(ratio_antigen_1, "_", ratio_antigen_2, "_ratio"))
  message(paste0("||---- ", ratio_name))
  data[!is.na(get(ratio_antigen_1)) & !is.na(get(ratio_antigen_2)), (ratio_name) := get(ratio_antigen_1) / get(ratio_antigen_2)]
  data[get(ratio_name) %in% c(Inf, -Inf), (ratio_name) := NA]
  data[is.nan(get(ratio_name)), (ratio_name) := NA]
  return(data)
}

#' Report data with age bins that don't map neatly on to the regular age bins (15-26 months and 18-29 months), 
#' are matched to the age bins they fall within (inclusive), signified by the `age_bin_agg_id` variable. 
#' Later in the modeling process, the microdata age distribution is relied upon to distribute the data in to
#' individual age_bins. This makes the assumption that the microdata age distribution roughly matches that of
#' the report data. This function checks that the raw microdata sample sizes are reasonably close to the report 
#' sample size
check_mismatched_age_bins <- function(data) {
  
  # Subset to mismatched data
  mismatch <- data[age_bin_agg_id == TRUE, ]
  # If no mismatched data, end function
  if(nrow(mismatch) == 0) return(NULL)
  # Collapse data 
  mismatch <- mismatch[, .("N_obs" = sum(N_obs, na.rm = TRUE)), by = c("svy_id", "ihme_loc_id", "year_start", "year_end", "year_id", "age_start", "age_end")]
  # Drop countries that aren't being modeled
  drop_countries <- c("KAZ", "MNE", "WSM", "SLB", "TUV")
  mismatch <- mismatch[!ihme_loc_id %in% drop_countries, ]
  # Drop surveys with no sample size
  mismatch <- mismatch[N_obs != 0, ]
  # For each mismatched survey, get the sample size of the mismatched age cohort from the raw microdata
  raw_files       <- list.files("FILEPATH", full.names = TRUE)
  nids            <- unique(mismatch$svy_id)
  raw_survey_size <- data.table()
  for(nid in nids) {
    # Get raw microdata
    file_path <- raw_files[grepl(paste0(nid, ".csv"), raw_files)]
    if (length(file_path) > 1) {
      file_path <- file_path[!grepl("_CH_", file_path)]
    } else if (length(file_path) == 0) {
      mismatched_age_range <- paste(unlist(mismatch[svy_id == nid, .(age_start, age_end)]), collapse = "-")
      message(paste0("\nUnable to compare sample size of mismatched age cohort ", mismatched_age_range, " from nid ", nid, " against microdata - no raw microdata available\n"))
      next
    }
    temp <- fread(file_path)
    # Get microdata sample size for ages in mismatched report age range
    age_start_age_end <- unlist(mismatch[svy_id == nid, .(age_start, age_end)])
    age_start         <- age_start_age_end[[1]]
    age_end           <- age_start_age_end[[2]]
    N                 <- temp[age_month >= age_start & age_month <= age_end, .N]
    raw_survey_size   <- rbind(raw_survey_size, 
                               data.table(svy_id = nid, 
                                          N_micro = N))
  }
  
  # Combine data
  mismatch <- merge(mismatch, raw_survey_size, all = TRUE, by = "svy_id")
  
  # Get difference
  mismatch[, diff := round(abs(N_obs - N_micro) / N_obs, digits = 2)]
  mismatch <- mismatch[order(-diff), ]
  setnames(mismatch, "N_obs", "N_report")
  
  # Message user
  message("\nLit data with age cohorts that don't map neatly on to the regular age cohorts are matched to overlapping age_bins. Mismatched age cohorts are signified by age_bin_agg_id == TRUE.
Later in the modeling process, the data from the mismatched age cohorts are distributed according to the distribution of the microdata. If the sample size of the microdata differ
too greatly from the sample size of the report data, this method might not be valid. Inspect mismatched data below for differences of greater than 10% between report and microdata:\n")
  print(mismatch)
  
  # Return data
  return(mismatch)
}


#----BEGIN FUNCTION ---------------------------------------------------------------------------------------------------
extract_report_data <- function(resample = FALSE){
  
  #----GET DATA--------------------------------------------------------------------------------------------------------
  # Message user
  message(paste0("======================================================\n||   EXTRACT REPORT DATA\n||-- Load data"))
  
  # Download data from Google Sheet
  raw_report_data      <- gs_download(lit_extraction_key, lit_extraction_sheet_name)
  report_data          <- prepare_report_data(raw_report_data)
  report_data          <- add_age_bins(report_data)
  mismatch_comparison  <- check_mismatched_age_bins(report_data)
  
  # Get outlier data (used in vaccine loop)
  outlier <- fread(outlier_filepath)
  
  #----PREP FOR RESAMPLE-----------------------------------------------------------------------------------------------
  # Get necessary resample functions from related scripts if resampling
  if (resample == TRUE) {
    source(paste0("FILEPATH"))
    source(paste0("FILEPATH"))
    source(paste0("FILEPATH"))
  }
  

  #----ADD RATIOS------------------------------------------------------------------------------------------------------
  # Message user
  message("||-- Add Ratios")
  report_data    <- add_ratios(report_data, "vacc_hib3", "vacc_dpt3")
  report_data    <- add_ratios(report_data, "vacc_hepb3", "vacc_dpt3")
  report_data    <- add_ratios(report_data, "vacc_pcv3",  "vacc_dpt3")
  report_data    <- add_ratios(report_data, "vacc_mcv2", "vacc_mcv1")
  
  
  #----VACCINE LOOP----------------------------------------------------------------------------------------------------
  # Message User
  message("||-- Vaccine Loop")
  vaccines_to_run <- c("bcg", "dpt", "hepb", "hib", "mcv", "pcv", "polio", "rcv", "rota", "yfv",
                       "hib3_dpt3_ratio", "hepb3_dpt3_ratio", "pcv3_dpt3_ratio", "mcv2_mcv1_ratio")
  
  for (vaccine_to_run in vaccines_to_run){
    
    message(paste0("||---- ", vaccine_to_run))
    
    # ----SUBSET DATA--------------------------------------------------------------------------------------------------
    
    is_ratio <- grepl("ratio", vaccine_to_run)
    
    # Get vaccine names for given vaccine stem
    if(!is_ratio) {
      dose_vars <- rename_table[stem == vaccine_to_run, unique(vaccine)]
    } else {
      dose_vars <- vaccine_to_run
    }
    
    # Keep only columns specific to vaccine in vaccine loop 
    keep_cols <- c("svy_id", "file_path", "point", "ihme_loc_id", "location_code", "shapefile", 
                   "year_id", "year_start", "year_end", "age_start", "age_end", 
                   "age_bin", "age_bin_agg", "age_bin_agg_id", "N_obs",
                   dose_vars)
    df_lbd <- subset(report_data, select = keep_cols[keep_cols %in% names(report_data)])
    
    # For multi-dose vaccines, drop rows missing data for any of the doses
    variable_dose_vaccines <- c("mcv", "rcv", "rota")
    if (vaccine_to_run %!in% variable_dose_vaccines){
      for(var in dose_vars) {
        df_lbd <- df_lbd[!is.na(get(var)), ]
      }
    } else {
      # For variable dose vaccines, only drop if 1st dose is missing
      df_lbd <- df_lbd[!is.na(get(paste0("vacc_", vaccine_to_run, "1"))), ]
    }
    
    
    # ----MAKE ROTAC---------------------------------------------------------------------------------------------------
    
    # Make RotaC based on location-specific dose schedule
    if (vaccine_to_run == "rota"){
      
      # Message user
      message("||------ Make RotaC")
      
      # Get rota dose schedule and intro years
      doses <- readRDS("FILEPATH")
      intro <- readRDS("FILEPATH")
      
      # Extract for each country in data
      df_rotac  <- data.table()
      countries <- sort(unique(df_lbd$ihme_loc_id))
      for (country in countries){
        data                   <- df_lbd[ihme_loc_id == country, ]
        country_specific_doses <- doses[me_name=="vacc_rotac" & ihme_loc_id==unique(data$ihme_loc_id), .(doses)]
        country_specific_intro <- intro[me_name=="vacc_rotac" & ihme_loc_id==unique(data$ihme_loc_id), .(cv_intro)] %>% unique
        
        # Cap doses at 3; if doses == 0, set to 2
        if (nrow(country_specific_doses) == 0) country_specific_doses <- 0
        country_specific_doses <- ifelse(country_specific_doses >= 3, 3, 2)
        rota_dose <- paste0("vacc_rota", country_specific_doses)
        
        # Set rota c doses if rotavirus has been introduced
        if (nrow(intro) == 0) stop(paste0("BREAK | Missing introduction year for ", unique(data$ihme_loc_id), "; need to prep introduction frame for this geography before continuing"))
        if (country_specific_intro < 9999) data$rotac <- data[, rota_dose, with=FALSE]
        df_rotac <- rbind(df_rotac, data, fill = TRUE)
      }
      df_lbd <- df_rotac
    }
    
    # ----APPLY DESIGN EFFECT------------------------------------------------------------------------------------------
    
    # Message user
    message("||------ Apply Design Effect")
    
    # Remove outliered nids from design effect
    if (is_ratio) {
      outlier_names <- gsub("_ratio", "", vaccine_to_run)
      outlier_names <- unlist(strsplit(outlier_names, split = "_"))
      outlier_names <- unique(outlier$me_name)[grep(paste(outlier_names, collapse = "|"), unique(outlier$me_name))]
      outlier_names <- outlier_names[!grepl("ratio",  outlier_names)]
      outlier_lbd   <- outlier[lbd == 1 & me_name %in% outlier_names,] 
    } else {
      outlier_names  <- unique(outlier$me_name)[grep(vaccine_to_run, unique(outlier$me_name))]
      outlier_names  <- outlier_names[!grepl("ratio",  outlier_names)]
      outlier_names  <- c(outlier_names, "")
      outliered_nids <- outlier[lbd == 1 & me_name %in% outlier_names, unique(nid)]
    }
    
    df_de         <- fread("FILEPATH")
    df_de         <- df_de[vacc == vaccine_to_run, ]
    df_de         <- df_de[nid %!in% outliered_nids, ]
    
    # Get vaccine-specific design effect as the median design effect
    de <- median(df_de$de, na.rm=T)

    # Apply design effect
    df_lbd[ , N := N_obs * de]
    
    # Convert from percentages to decimals
    for (var in dose_vars) {
      df_lbd[ ,(var) := (get(var)/100) * N]
    }
     
    # Set weight
    df_lbd[ , weight := 1]
    
    # ----CALCULATE DOSE-SPECIFIC COVERAGE-----------------------------------------------------------------------------
    
    # Message user
    message("||------ Calculate Dose-Specific")
    
    # Rename dose variables to indicate switch from dose aggregate to dose-specific
    setnames(df_lbd, rename_table$vaccine, rename_table$dose_specific, skip_absent = TRUE)
    
    # Convert from aggregate to dose-specific
    if(vaccine_to_run %in% c("dpt", "pcv", "hib", "hepb", "polio")){
      vars <- paste0(vaccine_to_run, c("_dose_1", "_dose_2", "_dose_3"))
      df_lbd[ , (vars[1]) := get(vars[1]) - get(vars[2])]
      df_lbd[ , (vars[2]) := get(vars[2]) - get(vars[3])]
      
      # Drop data from multi-dose vaccines which report non-first doses higher than prior doses
      df_lbd <- df_lbd[get(vars[1]) > 0, ]
      df_lbd <- df_lbd[get(vars[2]) > 0, ]
    }
    
    # Calculate dose-specific MCV values for locations where mcv2 is present
    if(vaccine_to_run == "mcv" & "mcv_dose_2" %in% names(df_lbd)) {
      vars <- paste0(vaccine_to_run, c("_dose_1", "_dose_2"))
      df_lbd[!is.na(get(vars[2])) , (vars[1]) := get(vars[1]) - get(vars[2])]
    }

    
    # ----RESAMPLE POLYGONS (LONG)-------------------------------------------------------------------------------------
    
    # Set output path
    datestamp <- gsub("-", "_", Sys.Date())
    outdir    <- ifelse(resample == TRUE, 
                        paste0(lit_extraction_output_root, datestamp, "/resampled/"), 
                        paste0(lit_extraction_output_root, datestamp, "/not_resampled/"))
    if (!dir.exists(outdir)){
      dir.create(outdir, recursive = TRUE)
    }
    
    # Resample
    if (resample == TRUE) {
      
      # Message user
      message("||------ Resampling")
      
      # Save pre-resampled data
      write.csv(df_lbd, paste0(outdir, vaccine_to_run, "_pre_resample.csv"), row.names=FALSE)

      # Resample polygon data
      if (nrow(df_lbd[point == 0]) > 0) {
        # Get population raster for resampling. If can't be read in from file location (quicker), 
        # create the raster and resave to file location for faster future loading
        if (!exists("popraster")) {
          if (file.exists("FILEPATH")) {
            popraster <- readRDS("FILEPATH")
          } else {
            popraster <- load_pop_raster()
          }
        }
        
        # Resample
        df_lbd <- resample_polygons(data = df_lbd,
                                    cores = 15,
                                    indic = vaccine_to_run,
                                    density = 0.001,
                                    gaul_list = get_adm0_codes('all', shapefile_version = 'current'),
                                    pull_poly_method = "fast")
        
        # Record number of children after resampling for data drop table
        data_drops_from_resample <- data.table()
        dose_specific_variables  <- rename_table[stem == vaccine_to_run, dose_specific]
        for(dose_specific_variable in dose_specific_variables) {
          n_after_resample <- df_lbd[, .("dose_specific"    = dose_specific_variable,
                                         "resampled"        = TRUE,
                                         "n_after_resample" = sum(weight * N, na.rm = TRUE)), by = svy_id]
          data_drops_from_resample <- rbind(data_drops_from_resample, n_after_resample, fill = TRUE)
        }
        data_drops_from_resample <- merge(rename_table[stem == vaccine_to_run, ], data_drops_from_resample, 
                                          all.x = TRUE, by = "dose_specific")
        
        # Format for merge with drop table
        data_drops_from_resample[, vaccine := gsub("vacc_", "", vaccine)]

        # Merge with data drop table
        data_drop_table_path <- paste0("FILEPATH")
        
        data_drop_table      <- fread(data_drop_table_path)
        data_drop_table      <- merge(data_drop_table, 
                                      data_drops_from_resample[, .("nid"                       = svy_id, 
                                                                   "antigen"                   = vaccine, 
                                                                   "n_after_resample_recorded" = n_after_resample, 
                                                                   "resample_recorded"         = resampled)], 
                                      all.x = TRUE, by = c("nid", "antigen"), sort = FALSE)
        
        # Format and drop redundant columns
        data_drop_table$n_after_resample <- as.numeric(data_drop_table$n_after_resample)
        data_drop_table[is.na(n_after_resample) & !is.na(n_after_resample_recorded), n_after_resample := n_after_resample_recorded]
        data_drop_table[resampled == FALSE & resample_recorded == TRUE, resampled := resample_recorded]
        data_drop_table$n_after_resample_recorded <- NULL
        data_drop_table$resample_recorded         <- NULL
        setcolorder(data_drop_table, c("nid", "total", "n_with_age_data", "n_0_59_months", "antigen"))
        
        # Save data drop table (saved in each loop to preserve data in case of error)
        resampled_data_drop_table_path         <- paste0("FILEPATH")
        resampled_data_drop_table_archive_path <- paste0("FILEPATH")
        
        write.csv(data_drop_table, file = resampled_data_drop_table_path, row.names = FALSE)
        write.csv(data_drop_table, file = resampled_data_drop_table_archive_path, row.names = FALSE)
      }
    }

    
    # ----SAVE DATA----------------------------------------------------------------------------------------------------

    # Message user
    message("||------ Saving")
    
    # Set column order
    final_vax_variables <- names(df_lbd)[grepl(vaccine_to_run, names(df_lbd))]
    setcolorder(df_lbd, c(names(df_lbd)[!names(df_lbd) %in% final_vax_variables], final_vax_variables))
    
    # Save data
    file_name <- ifelse(resample == TRUE, 
                        paste0(vaccine_to_run, "_resampled.csv"),
                        paste0(vaccine_to_run, ".csv"))
    write.csv(df_lbd, paste0(outdir, file_name), row.names=FALSE)
    message(paste0("||-------- ", outdir, file_name))
  }
  message("||-- Lit Extraction Complete: ")
  message("******************************************************\n")
}

#***********************************************************************************************************************

# Run
extract_report_data(resample = FALSE)




