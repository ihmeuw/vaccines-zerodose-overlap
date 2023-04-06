#----HEADER------------------------------------------------------------------------------------------------------------
#' Author:  USERNAME
#' Date:    DATE
#' Hub Explanation: <link>
#' Purpose: Process extracted survey microdata from UbCov
#'         - 00_processing_wrapper.R: Parent script used to launch child processing scripts
#'         - 01_process.R: Dose determination - for each child and antigen of interest, 
#'                         determine number of doses received
#'         - 02_tabulate_gbd.R: Collapse data by GBD location ids 
#'         - 02_tabulate_lbd.R: Collapse by lat/long or most granular admin unit
#'         - *03_disaggregate_and_resample.R*: For LBD - reshape data from long to wide, 
#'                                             resample polygons to points (if resample == TRUE), and save by vaccine.
#'                                             
# Run:     source("FILEPATH")
# Inputs:  nid  ->  NID from dataset in 'FILEPATH' to tabulate; 
#***********************************************************************************************************************

#----ENVIRONMENT--------------------------------------------------------------------------------------------------------
# Load required packages (should be redundant)
require(fasterize, lib.loc = "FILEPATH") # Fasterize is used for a modification in the resampling process.
require(sf, lib.loc = "FILEPATH" )       # Pairs with fasterize

#----HELPER FUNCTIONS------------------------------------------------------------------------------------------------
# Function to build simple population raster
build_simple_raster_pop_fast <- function(subset_shape, field="ADM0_CODE") {
  
  # Get raster brick
  master_pop <- brick('FILEPATH')
  
  # Get rasters of interest 
  rasters_of_interest <- c(paste0("worldpop_total_1y_", c(2000, 2005, 2010, 2015), "_00_00.tif"))
  raster_paths_all <- list.files(path = "FILEPATH", 
                                 pattern = ".tif$", full.names = TRUE)
  raster_paths     <- raster_paths_all[grepl(paste(rasters_of_interest, collapse = "|"), raster_paths_all)]
  
  # Combine population rasters in to raster brick 
  master_pop <- do.call(brick, lapply(raster_paths, raster))
  
  # Subset raster brick to regions of interest
  cropped_pop <- crop(master_pop, extent(subset_shape), snap="out")
  
  # Rasterize polygons
  initial_raster <- fasterize(st_as_sf(subset_shape), cropped_pop[[1]], field = field)
  
  if(length(subset(subset_shape, subset = !(get(field) %in% unique(initial_raster)))) != 0) {
    rasterized_shape <- merge(fasterize(st_as_sf(subset(subset_shape, !(get(field) %in% unique(initial_raster)))), cropped_pop[[1]], field = field), initial_raster)
  }
  if(length(subset(subset_shape, !(get(field) %in% unique(initial_raster)))) == 0) {
    rasterized_shape <- initial_raster
  }
  
  masked_pop <- raster::mask(x=cropped_pop, mask=rasterized_shape)
  
  raster_list <- list()
  raster_list[['simple_raster']] <- rasterized_shape
  raster_list[['pop_raster']]    <- masked_pop
  
  return(raster_list)
}

# Load population raster manually and save to file location for quicker future reading
load_pop_raster <- function(shapefile_version = 'current') {
  
  gaul_list         <- get_adm0_codes('all', shapefile_version = shapefile_version)
  use_1k_popraster  <-  TRUE
  
  # Loads popraster. Takes a long time, so saves to file location which can be read from in future runs.
  # File location (might) be intermittantly cleaned by mbg_setup, which empties the tmp files of any
  # file older than a week old, impacting the saved popraster. As of now, will need to intermittantly
  # load and save.
  if(!exists("popraster")){
    
    # Load polygon
    simple_polygon <- suppressMessages(suppressWarnings(load_simple_polygon(gaul_list = gaul_list, buffer=0.4, subset_only=TRUE,
                                                                            shapefile_version = shapefile_version)))
    # Merge polygons with population raster
    raster_list <- suppressMessages(suppressWarnings(build_simple_raster_pop_fast(simple_polygon$subset_shape)))
    
    # If using 1k popraster, create new, higher-resolution raster layer, disaggregating from 5km density to 1 km (.001) density
    if (use_1k_popraster){
      popraster <- disaggregate(raster_list[['pop_raster']], 5) 
    } else {
      popraster <-  raster_list[['pop_raster']]
    }
  }
  
  # Save popraster for quicker future loading
  saveRDS(popraster, "FILEPATH")
  
  # Return popraster 
  return(popraster)
}


#----DISAGGREGATE-------------------------------------------------------------------------------------------------------

disaggregate_and_resample <- function(nid, resample_polys = FALSE, vaccine_stems){
  
  #----GET DATA------------------------------------------------------------------------------------------------------
  # Message user
  if (resample_polys == TRUE){
    message(paste0("======================================================\n||   DISAGGREGATE & RESAMPLE: ", nid, "\n||-- Load data"))
  } else {
    message(paste0("======================================================\n||   DISAGGREGATE: ", nid, "\n||-- Load data"))
  }
  
  # Load data from tabulated folder
  df_tabulated <- load_nid(nid, folder="tabulated")
  
  # Load table for recording data drops from resampling
  if(resample_polys == TRUE){
    data_drop_table <- fread(file.path("FILEPATH"))
    data_drop_table[, n_after_resample := as.integer(n_after_resample)]
  }
  
  # If no "vaccine_stems" argument provided, get vaccine stems (and ratios) from data to be disaggregated
  if (vaccine_stems == "all") {
    vax_info              <- fread("FILEPATH")
    indicators_in_data    <- unique(df_tabulated$me_name)
    vaccines_in_data      <- indicators_in_data[!grepl("ratio", indicators_in_data)]
    ratios_in_data        <- indicators_in_data[grepl("ratio", indicators_in_data)]
    vaccine_stems_in_data <- vax_info[vaccine %in% gsub("vacc_", "", vaccines_in_data) & multi_vaccine == FALSE, unique(vaccine_stem)]
    ratios_stems_in_data  <- gsub("vacc_", "", ratios_in_data)
    vaccine_stems         <- c(vaccine_stems_in_data, ratios_stems_in_data)
  }
  
  # Clean up data types and survey names
  df_tabulated[, latitude := as.numeric(latitude)]
  df_tabulated[, longitude := as.numeric(longitude)]
  df_tabulated[, survey_name := gsub("/", "_", survey_name)] 
  
  # Standardize columns in data to allow for reshaping when disaggregating
  for (n in c("SSW", "N_obs")){  
    if (!n %in% names(df_tabulated)) {
      df_tabulated[ ,(n) := NA]
    }
  }
  
  #----DISAGGREGATE ANTIGENS-----------------------------------------------------------------------------------------
  
  if (max(df_tabulated$year_end) >= 2000){ 
    
    # Message user
    if (resample_polys == TRUE){
      message("||-- Disaggregate and Resample")
    } else {
      message("||-- Disaggregate")
    }
    
    # Vaccine loop for disaggregating antigens and ratios (and resampling if the resample argument is TRUE)
    if("mcv" %in% vaccine_stems) {
      vaccine_stems <- c("mcv", vaccine_stems[which(vaccine_stems != "mcv")])
    }
    
    for (stem in vaccine_stems){  
      
      # Message vaccine in loop
      message(paste0("||---- ", stem))
      
      # Get copy of data to subset
      df_input <- copy(df_tabulated)
      
      # Figure out which doses are present in the data and grab max
      dose_vars <- unique(df_input$me_name)[grepl(stem, unique(df_input$me_name))] 
      
      # If vaccine in loop isn't a ratio antigen, remove related ratio antigens from dose_vars
      if (!grepl("ratio", stem)) {
        dose_vars <- dose_vars[!grepl("ratio", dose_vars)]
      }
      
      if (length(dose_vars) > 0) { 
        
        # Subset to vaccine of interest
        df_input <- df_input[me_name %in% dose_vars, ]
        
        # Reshape wide
        df_input <- data.table(data.table::dcast(df_input, svy_id + country + point + latitude + longitude + location_code +
                                                   shapefile + survey_name + year_start + year_end + age_year + age_bin + age_bin_agg + age_bin_agg_id + 
                                                   year_id + psu + weight ~ me_name, sum, value.var = c("value", "N", "N_obs", "SSW"), fill = NA))
        
        # Rename dose columns from "value_vax_<vaccine>" to "vacc_<vaccine>"
        val_vars <- paste0("value_", dose_vars)
        setnames(df_input, val_vars, dose_vars)
        
        
        #----SWITCH FROM DOSE-AGGREGATE TO DOSE-SPECIFIC-------------------------------------------------------------------
        
        # Rename columns from reshape and calculate dose-specific doses from dose-aggregate
        if(stem == "dpt"){
          
          # Rename columns
          names(df_input)[names(df_input) == "vacc_dpt1"]       <- "dpt_dose_1"
          names(df_input)[names(df_input) == "vacc_dpt2"]       <- "dpt_dose_2"
          names(df_input)[names(df_input) == "vacc_dpt3"]       <- "dpt_dose_3"
          names(df_input)[names(df_input) == "N_vacc_dpt1"]     <- "N"
          names(df_input)[names(df_input) == "N_obs_vacc_dpt1"] <- "N_obs"
          names(df_input)[names(df_input) == "SSW_vacc_dpt1"]   <- "SSW"
          
          # Switch from dose-aggregate to dose-specific
          df_input[ , dpt_dose_1 := dpt_dose_1 - dpt_dose_2]
          df_input[ , dpt_dose_2 := dpt_dose_2 - dpt_dose_3]
        }
        
        if(stem == "pcv"){
          
          # Rename columns
          names(df_input)[names(df_input) == "vacc_pcv1"]       <- "pcv_dose_1"
          names(df_input)[names(df_input) == "vacc_pcv2"]       <- "pcv_dose_2"
          names(df_input)[names(df_input) == "vacc_pcv3"]       <- "pcv_dose_3"
          names(df_input)[names(df_input) == "N_vacc_pcv1"]     <- "N"
          names(df_input)[names(df_input) == "N_obs_vacc_pcv1"] <- "N_obs"
          names(df_input)[names(df_input) == "SSW_vacc_pcv1"]   <- "SSW"
          
          # Switch from dose-aggregate to dose-specific
          df_input[ , pcv_dose_1 := pcv_dose_1 - pcv_dose_2]
          df_input[ , pcv_dose_2 := pcv_dose_2 - pcv_dose_3]
        }
        
        if(stem == "polio"){
          
          # Rename columns
          names(df_input)[names(df_input) == "vacc_polio1"]       <- "polio_dose_1"
          names(df_input)[names(df_input) == "vacc_polio2"]       <- "polio_dose_2"
          names(df_input)[names(df_input) == "vacc_polio3"]       <- "polio_dose_3"
          names(df_input)[names(df_input) == "N_vacc_polio1"]     <- "N"
          names(df_input)[names(df_input) == "N_obs_vacc_polio1"] <- "N_obs"
          names(df_input)[names(df_input) == "SSW_vacc_polio1"]   <- "SSW"
          
          # Switch from dose-aggregate to dose-specific
          df_input[ , polio_dose_1 := polio_dose_1 - polio_dose_2]
          df_input[ , polio_dose_2 := polio_dose_2 - polio_dose_3]
        }
        
        if(stem == "hepb"){
          
          # Rename columns
          names(df_input)[names(df_input) == "vacc_hepb1"]       <- "hepb_dose_1"
          names(df_input)[names(df_input) == "vacc_hepb2"]       <- "hepb_dose_2"
          names(df_input)[names(df_input) == "vacc_hepb3"]       <- "hepb_dose_3"
          names(df_input)[names(df_input) == "N_vacc_hepb1"]     <- "N"
          names(df_input)[names(df_input) == "N_obs_vacc_hepb1"] <- "N_obs"
          names(df_input)[names(df_input) == "SSW_vacc_hepb1"]   <- "SSW"
          
          # Switch from dose-aggregate to dose-specific
          df_input[ , hepb_dose_1 := hepb_dose_1 - hepb_dose_2]
          df_input[ , hepb_dose_2 := hepb_dose_2 - hepb_dose_3]
        }
        
        if(stem == "hib"){
          
          # Rename columns
          names(df_input)[names(df_input) == "vacc_hib1"]       <- "hib_dose_1"
          names(df_input)[names(df_input) == "vacc_hib2"]       <- "hib_dose_2"
          names(df_input)[names(df_input) == "vacc_hib3"]       <- "hib_dose_3"
          names(df_input)[names(df_input) == "N_vacc_hib1"]     <- "N"
          names(df_input)[names(df_input) == "N_obs_vacc_hib1"] <- "N_obs"
          names(df_input)[names(df_input) == "SSW_vacc_hib1"]   <- "SSW"
          
          # Switch from dose-aggregate to dose-specific
          df_input[ , hib_dose_1 := hib_dose_1 - hib_dose_2]
          df_input[ , hib_dose_2 := hib_dose_2 - hib_dose_3]
        }
        
        if(stem == "bcg"){
          
          # Rename columns
          names(df_input)[names(df_input) == "vacc_bcg"]       <- "bcg_dose_1"
          names(df_input)[names(df_input) == "N_vacc_bcg"]     <- "N"
          names(df_input)[names(df_input) == "N_obs_vacc_bcg"] <- "N_obs"
          names(df_input)[names(df_input) == "SSW_vacc_bcg"]   <- "SSW"
        }
        
        if(stem == "yfv"){
          
          # Rename columns
          names(df_input)[names(df_input) == "vacc_yfv"]       <- "yfv_dose_1"
          names(df_input)[names(df_input) == "N_vacc_yfv"]     <- "N"
          names(df_input)[names(df_input) == "N_obs_vacc_yfv"] <- "N_obs"
          names(df_input)[names(df_input) == "SSW_vacc_yfv"]   <- "SSW"
        }
        
        if(stem == "mcv"){
          
          # Rename columns
          names(df_input)[names(df_input) == "vacc_mcv1"]       <- "mcv_dose_1"
          names(df_input)[names(df_input) == "vacc_mcv2"]       <- "mcv_dose_2"
          names(df_input)[names(df_input) == "N_vacc_mcv1"]     <- "N"
          names(df_input)[names(df_input) == "N_obs_vacc_mcv1"] <- "N_obs"
          names(df_input)[names(df_input) == "SSW_vacc_mcv1"]   <- "SSW"
          
          # Switch from dose-aggregate to dose-specific
          if ("mcv_dose_2" %in% names(df_input)) {
            df_input[!is.na(mcv_dose_2), mcv_dose_1 := mcv_dose_1 - mcv_dose_2]
          }
        }
        
        if(stem == "rota"){
          
          # Rename columns (surveys before intro year do not have rota complete)
          names(df_input)[names(df_input) == "vacc_rota1"]       <- "rota_dose_1"
          names(df_input)[names(df_input) == "vacc_rota2"]       <- "rota_dose_2"
          names(df_input)[names(df_input) == "vacc_rota3"]       <- "rota_dose_3"
          names(df_input)[names(df_input) == "vacc_rotac"]       <- "rota_dose_c"
          
          if("rota_dose_c" %in% names (df_input)) {  
            names(df_input)[names(df_input) == "N_vacc_rotac"]     <- "N"
            names(df_input)[names(df_input) == "N_obs_vacc_rotac"] <- "N_obs"
            names(df_input)[names(df_input) == "SSW_vacc_rotac"]   <- "SSW"
          } else {
            names(df_input)[names(df_input) == "N_vacc_rota1"]     <- "N"
            names(df_input)[names(df_input) == "N_obs_vacc_rota1"] <- "N_obs"
            names(df_input)[names(df_input) == "SSW_vacc_rota1"]   <- "SSW"
          }  
          
          # Switch from dose-aggregate to dose-specific
          if ("rota_dose_2" %in% names(df_input)) {
            df_input[ , rota_dose_1 := rota_dose_1 - rota_dose_2]
          }
          if ("rota_dose_3" %in% names(df_input)) {
            df_input[ , rota_dose_2 := rota_dose_2 - rota_dose_3]
          }
        }
        
        if(stem == "rcv"){
          
          # Rename columns
          names(df_input)[names(df_input) == "vacc_rcv1"]       <- "rcv_dose_1"
          names(df_input)[names(df_input) == "vacc_rcv2"]       <- "rcv_dose_2"
          names(df_input)[names(df_input) == "N_vacc_rcv1"]     <- "N"
          names(df_input)[names(df_input) == "N_obs_vacc_rcv1"] <- "N_obs"
          names(df_input)[names(df_input) == "SSW_vacc_rcv1"]   <- "SSW"
          
          # Switch from dose-aggregate to dose-specific
          if ("rcv_dose_2" %in% names(df_input)) {
            df_input[!is.na(rcv_dose_2), rcv_dose_1 := rcv_dose_1 - rcv_dose_2]
          }
        }
        
        if(stem == "hib3_dpt3_ratio"){
          names(df_input)[names(df_input) == "N_hib3_dpt3_ratio"]      <- "N"
          names(df_input)[names(df_input) == "N_obs_hib3_dpt3_ratio"]  <- "N_obs"
          names(df_input)[names(df_input) == "SSW_hib3_dpt3_ratio"]    <- "SSW"
        }
        
        if(stem == "hepb3_dpt3_ratio"){
          names(df_input)[names(df_input) == "N_hepb3_dpt3_ratio"]      <- "N"
          names(df_input)[names(df_input) == "N_obs_hepb3_dpt3_ratio"]  <- "N_obs"
          names(df_input)[names(df_input) == "SSW_hepb3_dpt3_ratio"]    <- "SSW"
        }
        
        if(stem == "mcv2_mcv1_ratio"){
          names(df_input)[names(df_input) == "N_mcv2_mcv1_ratio"]      <- "N"
          names(df_input)[names(df_input) == "N_obs_mcv2_mcv1_ratio"]  <- "N_obs"
          names(df_input)[names(df_input) == "SSW_mcv2_mcv1_ratio"]    <- "SSW"
        }
        
        if(stem == "pcv3_dpt3_ratio"){
          names(df_input)[names(df_input) == "N_pcv3_dpt3_ratio"]      <- "N"
          names(df_input)[names(df_input) == "N_obs_pcv3_dpt3_ratio"]  <- "N_obs"
          names(df_input)[names(df_input) == "SSW_pcv3_dpt3_ratio"]    <- "SSW"
        }
        
        if(stem == "rotac_dpt3_ratio"){
          names(df_input)[names(df_input) == "N_rotac_dpt3_ratio"]      <- "N"
          names(df_input)[names(df_input) == "N_obs_rotac_dpt3_ratio"]  <- "N_obs"
          names(df_input)[names(df_input) == "SSW_rotac_dpt3_ratio"]    <- "SSW"
        }
        
        if(stem == "dpt3_timeliness_ratio"){
          names(df_input)[names(df_input) == "vacc_dpt3_timeliness_ratio"]       <- "dpt3_timeliness_ratio"
          names(df_input)[names(df_input) == "N_vacc_dpt3_timeliness_ratio"]     <- "N"
          names(df_input)[names(df_input) == "N_obs_vacc_dpt3_timeliness_ratio"] <- "N_obs"
          names(df_input)[names(df_input) == "SSW_vacc_dpt3_timeliness_ratio"]   <- "SSW"
        }
        
        
        # Remove other vaccine variables
        if (grepl("ratio", stem)) {
          df_input <- subset(df_input, select = !(names(df_input) %in% c(paste0("N_", dose_vars),
                                                                         paste0("N_obs_", dose_vars),
                                                                         paste0("SSW_", dose_vars))))
        } else {
          df_input <- subset(df_input, select = !(names(df_input) %in% c(dose_vars, 
                                                                         paste0("N_", dose_vars),
                                                                         paste0("N_obs_", dose_vars),
                                                                         paste0("SSW_", dose_vars))))
        }
        
        
        #----SAVE DESIGN EFFECT-----------------------------------------------------------------------------------------
        
        # Save design effect information (weighted N)/(unweighted N) of polygon data
        if (nrow(df_input[point==0, ]) > 0){  
          design_effect_file_path <- "FILEPATH"
          de_data   <- fread(design_effect_file_path)    # Existing DE data to append to/modify
          svy_id    <- nid                               # To avoid namespace confusion
          N_sum     <- sum(df_input[point == 0, N])      # Sum of weighted N in polygon data
          N_obs_sum <- sum(df_input[point == 0, N_obs])  # Sum of unweighted N in polygon data
          de_calc   <- N_sum/N_obs_sum                   # Survey-level design effect is (weighted N)/(unweighted N)
          
          # If nid/vaccine exists in sheet, modify row, otherwise append
          if (paste0(nid, stem) %in% unique(paste0(de_data$nid, de_data$vacc))) { 
            de_data[nid == svy_id  & vacc == stem, N := N_sum]
            de_data[nid == svy_id  & vacc == stem, N_obs := N_obs_sum]
            de_data[nid == svy_id  & vacc == stem, de := de_calc]
          } else {  
            de_data <- rbind(de_data, 
                             data.table("nid"   = svy_id, 
                                        "vacc"  = stem, 
                                        "N"     = N_sum, 
                                        "N_obs" = N_obs_sum, 
                                        "de"    = de_calc)) 
          }
          
          # Save updated design effect sheet
          write.csv(de_data, file = design_effect_file_path, row.names=FALSE)         
        }
        
        #----RESAMPLE POLYGON DATA--------------------------------------------------------------------------------------
        
        # Resample polygon data
        if (resample_polys == TRUE) {
          
          if (nrow(df_input[point == 0]) > 0) {
            # Get population raster for resampling. If can't be read in from file location (quicker), 
            # create the raster and resave to file location for faster future loading
            if (!exists("popraster")) {
              if (file.exists("FILEPATH")) {
                popraster <- readRDS("FILEPATH")
              } else {
                popraster <- load_pop_raster()
              }
            }
            
            # Get necessary functions from related scripts and resample
            source(paste0("FILEPATH"))
            source(paste0("FILEPATH"))
            source(paste0("FILEPATH"))
            
            df_pointpoly <- resample_polygons(data = df_input,
                                              cores = 10,
                                              indic = stem,
                                              density = 0.001,
                                              gaul_list = get_adm0_codes('all', shapefile_version = 'current'),
                                              pull_poly_method = "fast")
            
          } else {
            # If no polys, don't resample!
            df_pointpoly <- df_input
          }
          
          # Record any data dropped due to resampling
          data_drop_table[grepl(stem, antigen), resampled := TRUE]
          data_drop_table[grepl(stem, antigen), n_after_resample := df_pointpoly[, sum(N * weight)]]
          write.csv(data_drop_table, file.path("FILEPATH"), row.names = FALSE)

        } else {
          df_pointpoly <- df_input
        }
        
        #----ADD VACCINE SCHEDULE---------------------------------------------------------------------------------------
        # Add vaccine cohort schedule:
        # - Regular coverage modeling: children from age of recommended vaccination will be dropped due to incomplete coverage
        # - Age-cohort modeling: all data is preserved
        
        # Get vaccine cohort schedule for survey country
        vax_targets  <- set_target_ages(df_pointpoly, vaccines.=vaccines)
        setnames(vax_targets, "age_cohort", "me_cohort_schedule")   
        
        # Identify antigen in loop and get vaccine cohort schedule for 1st dose
        vaccine_name       <- paste0("vacc_", vax_info[vaccine_stem == stem & dose == 1, vaccine])
        me_cohort_schedule <- unique(vax_targets[me_name == vaccine_name, me_cohort_schedule])
        
        # Add to dataset
        df_pointpoly$me_cohort_schedule <- me_cohort_schedule

        #----SAVE DISAGGREGATED DATA------------------------------------------------------------------------------------
        
        # Save disaggregated (vax-specific) data
        outdir <- ifelse(resample_polys == TRUE,
                         file.path(extraction_root, "03_disaggregated", "resampled", stem),
                         file.path(extraction_root, "03_disaggregated", "not_resampled", stem))
        if (!dir.exists(outdir)) dir.create(outdir)
        write.csv(df_pointpoly, file.path(outdir, paste0(nid, ".csv")), row.names=FALSE)

      } else {
        message(paste0("||------ No dose variables, skipping.."))
      }
    }
    
    # Disaggregation complete, message user
    if (resample == TRUE){
      message(paste0("||-- Disaggregate & resample complete: ", nid))
    } else {
      message(paste0("||-- Disaggregate complete: ", nid))
    } 
    message("******************************************************\n")
  } else {
    message(paste0("||-- All pre-2000, skipping disaggregate"))
    message("******************************************************\n")
  }
} 

