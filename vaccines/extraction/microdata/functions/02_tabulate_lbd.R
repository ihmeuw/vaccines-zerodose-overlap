#----HEADER------------------------------------------------------------------------------------------------------------
#' Author:  USERNAME
#' Date:    DATE
#' Hub Explanation: <link>
#' Purpose: Process extracted survey microdata from UbCov
#'         - 00_processing_wrapper.R: Parent script used to launch child processing scripts
#'         - 01_process.R: Dose determination - for each child and antigen of interest, 
#'                         determine number of doses received
#'         - 02_tabulate_gbd.R: Collapse data by GBD location ids 
#'         - *02_tabulate_lbd.R*: Collapse by lat/long or most granular admin unit
#'         - 03_disaggregate_and_resample.R: For LBD - reshape data from long to wide, 
#'                                           resample polygons to points (if resample == TRUE), and save by vaccine.
#'                                           
#' Run:     source("FILEPATH")
#' Inputs:  nid  ->  NID from dataset in 'FILEPATH' to tabulate
#**********************************************************************************************************************

#----ENVIRONMENT-------------------------------------------------------------------------------------------------------
# Source scripts
source("FILEPATH")
source("FILEPATH")

# Load additional packages
mbg_packages_path <- "FILEPATH"
package_list      <- fread(mbg_packages_path, header = FALSE)[[1]]
suppressMessages(mbg_setup(package_list=package_list, repos=core_repo))

# Specific function loads (avoid namespace conflicts)
str_match <- stringr::str_match

#----HELPER FUNCTIONS--------------------------------------------------------------------------------------------------
# Given pair of antigens, return dataset with ratio of antigens created at cluster-level.
# Ratios not generated where denominator (ratio_antigen_2) is 0.
add_ratios <- function(data, ratio_antigen_1, ratio_antigen_2) {
  vaccines_in_dataset    <- unique(data$me_name)
  ratio_antigen_1_exists <- ratio_antigen_1 %in% vaccines_in_dataset
  ratio_antigen_2_exists <- ratio_antigen_2 %in% vaccines_in_dataset
  
  if (ratio_antigen_1_exists & ratio_antigen_2_exists) {
    
    ### Reshape long to wide in order to add ratios, then wide back to long
    # value   me_name   (Long to Wide)   DPT3   HepB3  HepB3_DPT3_Ratio  (Wide to Long)   value   me_name   
    #   2      DPT3          --->          2      1          .5              --->          2      DPT        
    #   1      HepB3         --->                                            --->          1      HepB       
    #                                                                                     .5      HepB3_DPT3_Ratio 
    
    # Reshape long to wide
    lhs_names <- names(data)[!names(data) %in% c("me_name", "value")]
    lhs       <- paste(lhs_names, collapse = " + ")
    rhs       <- "me_name"
    lhs_rhs   <- paste(c(lhs, rhs), collapse = " ~ ")
    
    data <- data.table(data.table::dcast(data, eval(parse(text = lhs_rhs))))
    
    # Add ratios
    ratio_column_name <- gsub("vacc_", "", paste0(ratio_antigen_1, "_", ratio_antigen_2, "_ratio"))
    data[, (ratio_column_name) := numeric() ]
    data[!is.na(get(ratio_antigen_2)) & get(ratio_antigen_2) != 0 & !is.na(get(ratio_antigen_1)), (ratio_column_name) := get(ratio_antigen_1) / get(ratio_antigen_2)]
    
    # Reshape wide back to long
    vax_cols     <- names(data)[grep("vacc|ratio", names(data))]
    id_variables <- names(data)[!names(data) %in% vax_cols]
    data <- melt.data.table(data, id.vars=id_variables, measure.vars = vax_cols, variable.name="me_name", variable.factor = F, value.name = "value")
    
    # Drop NAs resulting from reshaping process  
    data <- data[!is.na(value), ]
  }
  return(data)
}


#----TABULATE----------------------------------------------------------------------------------------------------------

tabulate_lbd <- function(nid, collapse_method="kish") {
  
  #----GET DATA--------------------------------------------------------------------------------------------------------
  # Message user
  message(paste0("======================================================\n||   TABULATE (LBD): ", nid))
  
  # Get geodata
  geo <- suppressMessages(get_geocodebooks(nids = nid))
  
  # Get survey data, check for required variables, and assign year_id and age_bins
  message("||-- Load and Prep Data")
  dataset <- prep_for_tabulation(nid, team="lbd", vaccines.=c(vaccines, "rotac", "dpt3_timeliness_ratio")) 
  
  # Data from prep for tabulation is returned as a list. The first element is T/F indicating whether the data is usable,
  # and the 2nd element is the data
  if (dataset[[1]] != FALSE){ 
    data            <- dataset[[2]]
    data_drop_table <- dataset[[3]] 
    
    # If we are missing these variables, add them as NAs to make it through reshapes/tabulations. Only drop if these are missing for polygons later
    for (variable in c("pweight", "psu", "strata")){   
      if (!variable %in% names(data)){
        data[ ,(variable):=NA]
      }
    }
    
    if(max(data$year_end) >= 2000){  # Include surveys from 2000 onwards
      
      #----GEOMATCH DATA---------------------------------------------------------------------------------------------
      # Message User
      message("||-- Geomatch Data")
      
      # Prep survey data and geodata for merge
      data$geospatial_id <- as.character(data$geospatial_id)
      data$ihme_loc_id   <- as.character(data$ihme_loc_id)
      data$ihme_loc_id   <- substr(data$ihme_loc_id, 0, 3)
      setnames(geo, "iso3", "ihme_loc_id")
      geo$ihme_loc_id    <- substr(geo$ihme_loc_id, 0, 3)
      geo$geospatial_id  <- as.character(geo$geospatial_id)
      
      # Merge survey data with geodata
      data <- merge(geo, data, 
                    by=c("nid", "geospatial_id", "ihme_loc_id"), 
                    all.x=FALSE, all.y=TRUE, allow.cartesian=TRUE)  
      
      # Rename columns
      rename_table <- data.table(rbind(c("ihme_loc_id", "country"),
                                       c("nid", "svy_id"),
                                       c("geospatial_id", "geo_id"), 
                                       c("lat", "latitude"),
                                       c("long", "longitude"), 
                                       c("loc_name1", "admin1_name"),
                                       c("loc_code1", "loc_code_admin1"),
                                       c("admin_level1", "admin_level1"),
                                       c("point", "point")))
      names(rename_table) <- c("old", "new")
      invisible(lapply(1:nrow(rename_table), function(x) check_and_rename_col(data, rename_table$old[x], rename_table$new[x])))
      

      
      #----PREP DATA FOR POINT/POLY COLLAPSE-----------------------------------------------------------------------------
      
      # Ensure column formats correct
      if (class(data$latitude)  != "numeric") data[, latitude := as.numeric(as.character(latitude))]
      if (class(data$longitude) != "numeric") data[, longitude := as.numeric(as.character(longitude))]
      if (class(data$strata)    != "numeric") data[, strata := as.numeric(as.character(strata))]
      if (class(data$pweight)   != "numeric") data[, pweight := as.numeric(as.character(pweight))]
      if ("cluster_id" %in% names(data))      data[, cluster_id := NULL]
      
      # Fix country iso3 codes
      data[, country := tstrsplit(country, "_")[1]]
      
      # Set empty shapefiles to missing
      data[shapefile == "", shapefile := NA]
      
      # Get shapefiles into character format (in case factor)
      data$shapefile <- as.character(data$shapefile)
      
      # Identify missing shapefiles and exclude a priori (otherwise will break resample_polygons)
      shapefile_list       <- unique(data[!is.na(shapefile), shapefile])
      shapefile_dir        <- "FILEPATH"
      shapefiles_available <- list.files(shapefile_dir, ".shp$") %>% gsub(".shp", "", .)
      no_shapefile         <- shapefile_list[!(shapefile_list %in% shapefiles_available)]
      
      if (length(no_shapefile) > 0) {
        message("The following shapefiles are missing from the shapefile directory. Associated data will be dropped:")
        message(paste(paste0("  ", no_shapefile, ".shp"), collapse = " \n"))
        data <- data[!(shapefile %in% no_shapefile), ]
      }
      
      # Geo missingness log file
      geo_missingness <- nrow(data[is.na(latitude) & is.na(shapefile)]) / nrow(data)
      vet_log         <- fread(file.path(extraction_root, "log/details", paste0(nid, ".csv")))
      vet_log[, missingness_geomatch := geo_missingness]
      write.csv(vet_log, file.path("FILEPATH"), row.names=FALSE)
      
      # Confirm that points are appropriately categorized
      data[is.na(point) & !is.na(latitude) & !is.na(longitude), point := 1]
      data[is.na(point) & !is.na(shapefile) & !is.na(pweight),  point := 0]
      
      # Drop data with geomissingness
      data <- data[!(is.na(latitude) & is.na(shapefile))] # geo_drop
      data <- data[!(is.na(latitude) & is.na(pweight))]   # pweight_drop if polygon data does not have pweight
      data <- data[!is.na(year_id), ]                     # age_drop
      data <- data[!is.na(value), ]                       # missing vaccination_data
      
      # Record geomissingness in data drop table
      geo_missingness <- data[, .N, by = me_name]
      geo_missingness <- geo_missingness[me_name %in% paste0("vacc_", data_drop_table$antigen), ]
      geo_missingness[, antigen := gsub("vacc_", "", me_name)]
      data_drop_table <- merge(data_drop_table, geo_missingness[, .(antigen, N)], all.x = TRUE, all.y = FALSE)
      data_drop_table[, n_after_geomatch := N]
      data_drop_table$N <- NULL
      setcolorder(data_drop_table, c("nid", "total", "n_with_age_data", "n_0_59_months", "antigen"))
      write.csv(data_drop_table, file = paste0("FILEPATH"), row.names = FALSE)
      
      #----COLLAPSE POINT/POLYGON DATA-----------------------------------------------------------------------------------
      # Message User
      message("||-- Collapse Data")
      
      # Placeholder value to sum over
      data[, N := 1]
      
      # Split in to point and poly data
      df_point <- data[point == 1]
      df_poly  <- data[point == 0]
      
      # Collapse point data
      if (nrow(df_point) > 0) {
        keep_vars <- c("svy_id", "survey_name", "country", "point", "svy_year", "year_start", "year_end", 
                       "location_code", "shapefile", 
                       "year_id", "age_year", "age_bin", "age_bin_agg", "age_bin_agg_id",
                       "psu", "latitude", "longitude", "outlier", "me_name", "value", "N")
        df_point <- subset(df_point, select = names(df_point) %in% keep_vars)
        df_point[, N_obs := 1]
        
        # Collapse to clusters
        df_point <- df_point[, lapply(.SD, sum), 
                             by = .(svy_id, survey_name, country, point, 
                                    location_code, shapefile, 
                                    year_start, year_end, year_id, age_year, age_bin, age_bin_agg, age_bin_agg_id, 
                                    psu, latitude, longitude, me_name),
                             .SDcols = c("value", "N", "N_obs")]
      }
      
      # Collapse the polygon data
      if (nrow(df_poly) > 0) {
        if (collapse_method == "kish") {
          # Collapse using the Kish approximation to generate an effective sample size (N)
          # Keeps the original sample size in N_obs
          df_poly <- df_poly[, c(list(N_obs = sum(N),
                                      SSW = sum(pweight),
                                      N = sum(pweight) ^ 2 / sum(pweight ^ 2)),
                                 lapply(.SD, FUN = weighted.mean, w = pweight)), 
                             by = list(svy_id, survey_name, country, point, 
                                       location_code, shapefile, 
                                       year_start, year_end, year_id, age_year, age_bin, age_bin_agg, age_bin_agg_id, me_name),
                             .SDcols = "value"]
          
          # Convert these back to counts
          df_poly[, ("value") := lapply(.SD, '*', N), .SDcols = "value"]                            
          
        } else {
          # Throw error if not Kish method
          stop(paste0("Polygon collapse is only configured with the 'kish' method. '", collapse_method, "' method isn't configured"))
        }
      }
      
      # Combine point data with collapsed poly data
      df_pointpoly <- rbind(df_point, df_poly, fill=TRUE)
      if ("age_month" %in% names(df_pointpoly)) df_pointpoly$age_month <- NULL
      
      # Add a weights column to be resampled in resample_polygons
      df_pointpoly$weight <- 1
      
      # ADD RATIOS -------------------------------------------------------------
      # Message User
      message("||-- Add Vaccine Ratios")
      
      df_pointpoly[, value := as.double(value)]
      
      df_pointpoly <- add_ratios(df_pointpoly, "vacc_mcv2",  "vacc_mcv1")
      df_pointpoly <- add_ratios(df_pointpoly, "vacc_pcv3",  "vacc_dpt3")
      df_pointpoly <- add_ratios(df_pointpoly, "vacc_rotac", "vacc_dpt3")
      df_pointpoly <- add_ratios(df_pointpoly, "vacc_hib3",  "vacc_dpt3")
      df_pointpoly <- add_ratios(df_pointpoly, "vacc_hepb3", "vacc_dpt3")
      
      
      # SAVE DATA --------------------------------------------------------------
      
      # Message user
      message("||-- Save Data")
      write.csv(df_pointpoly, paste0("FILEPATH"), row.names=FALSE)
      message(paste0("||-- Tabulate LBD Complete: ", nid))
      message("******************************************************\n")      
      
    } else {
      
      message(paste0(nid, ": All data is pre-2000, skipping tabulation"))
      if (paste0(nid, ".txt") %in% list.files(file.path("FILEPATH"))) {
        cat("\nDataPre2000\n", file=file.path("FILEPATH"), append=TRUE)
      }
    } 
  } else {
    message(paste0(nid, ": Tabulation prep failed"))
  }
}

#***********************************************************************************************************************
