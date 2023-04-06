##########################################################################
### Function to read in, clean, and save .rds of all GBD 2019 extracted MICS tabulations from Fix The MICS for the Decider
##########################################################################

# based off of preparation of literature.rds
mics_read_save <- function(...) {
  #----SETUP--------------------------------------------------------------------------------------------------------------
  ### load packages
  username <- Sys.info()[["user"]]
  pacman::p_load(data.table, dplyr, parallel, readxl)  # make sure these already loaded at top of full script maybe or just add new ones
  
  ### source functions
  source(paste0("FILEPATH"))
  source(db_tools)
  source("FILEPATH")
  locations <- get_location_metadata(location_set_id=location_set_id, gbd_round_id=gbd_round)[level >= 3]  
  
  df <- fread("FILEPATH")
  #***********************************************************************************************************************
  
  #----PREP---------------------------------------------------------------------------------------------------------------
  ### if have higher combination vaccine coverage, swap out component vaccines
  duples <- list(c("vacc_dpt1",  "vacc_tetra1"),
                 c("vacc_dpt1",  "vacc_pent1"),
                 c("vacc_dpt2",  "vacc_pent2"),
                 c("vacc_dpt2",  "vacc_tetra2"),
                 c("vacc_dpt3",  "vacc_tetra3"),
                 c("vacc_dpt3",  "vacc_pent3"),
                 c("vacc_hib3",  "vacc_pent3"),
                 c("vacc_hib3",  "vacc_tetra3"),
                 c("vacc_hepb3", "vacc_pent3"),
                 c("vacc_mcv1",  "vacc_mmr1"),
                 c("vacc_mcv2",  "vacc_mmr2"),
                 c("vacc_rcv1",  "vacc_mmr1"),
                 c("vacc_rcv2",  "vacc_mmr2"))
  for (i in duples) {
    ## e.g. if !is.na(penta) & penta > dpt | !is.na(penta) & is.na(dpt), replace dpt with penta
    df <- df[, (i[1]) := ifelse((!is.na(get(i[2])) & get(i[2]) > get(i[1])) | (!is.na(get(i[2])) & is.na(get(i[1]))), get(i[2]), get(i[1]))]
  }
  
  ### location
  df[, ihme_loc_id := tstrsplit(location_name_short_ihme_loc_id, "[|]")[[2]] ]
  df[is.na(ihme_loc_id), ihme_loc_id := parent_location]

  ### reshape
  id <- c("nid", "file_path", "ihme_loc_id", "year_start", "year_end", "age_start", "age_end", "sample_size")
  vacc <-  grep("vacc", names(df), value=TRUE)
  id <- c("nid", "file_path", "ihme_loc_id", "year_start", "year_end", "age_start", "age_end", "sample_size")
  df <- melt(df, id=id, measure=patterns("^vacc|^maternal"), value.name="data", variable.name="me_name")
  
  ## going to have to match column classes, etc.
  #***********************************************************************************************************************
  
  #----CLEANUP------------------------------------------------------------------------------------------------------------
  ### drop rows without antigen-specific data
  df <- df[!is.na(data)]
  
  ### divide coverage by 100
  df[, data := as.numeric(data) / 100]
  # cap coverage at 1
  df[data > 1, data := 1]
  
  ### drop anything without mapped ihme_loc_id
  df <- df[!is.na(ihme_loc_id)]
  # remove sample size if > 5000, usually refers to target population
  df[, sample_size := gsub(",", "", sample_size) %>% as.numeric]
  df[sample_size > 5000, sample_size := NA]
  # ages --> Round to single year ages, assume age_year = 1 otherwise
  df[, `:=` (age_start = round(age_start / 12), 
             age_end   = round(age_end / 12))]
  df[, age_length := age_end - age_start]
  df[age_length %in% c(0, 1), age_year := age_start]
  df[is.na(age_length), age_year := 1]
  df[age_length > 1, age_year := 1]
  
  df[, year_id := year_start-age_year] 
  df[, c("year_start", "year_end") := NULL]
  
  parsed <- df[, tstrsplit(file_path, "_", fixed=TRUE)]
  parsed <- parsed[, survey_name := ifelse(nchar(V3)==3, paste0(V3, "/", V4), V3)]
  df <- cbind(df, parsed$survey_name)
  setnames(df, "V2", "survey_name")
  df <- df[survey_name=="MICS", survey_name := "UNICEF_MICS"]
  #***********************************************************************************************************************
  
  #----SAVE---------------------------------------------------------------------------------------------------------------
  ### save
  saveRDS(df, paste0("FILEPATH"))
  print("MICS report tabulations prepped and saved for use in GBD pipeline")
  #***********************************************************************************************************************
  
  return(df)
}


### run
saved_tabs <- mics_read_save()




##########################################################################
### Function to clean MICS report data for MBG December Deliverables DPT3 timeliness ratio
##########################################################################

mics_timeliness_reports <- function(...) {
  #----SETUP--------------------------------------------------------------------------------------------------------------
  ### load packages
  pacman::p_load(data.table, dplyr, parallel, readxl)  

  ### source functions
  source(paste0(j, "FILEPATH"))
  source(db_tools)
  source("FILEPATH")
  locations <- get_location_metadata(location_set_id=location_set_id, gbd_round_id=gbd_round)[level >= 3]  
  
  df <- fread("FILEPATH")
  #***********************************************************************************************************************
    
  #----CLEANUP------------------------------------------------------------------------------------------------------------
  ### standardize col names and row identifiers
  setnames(df, c("parent_location.x", "timeliness_dpt3"), c("ihme_loc_id", "data"))
  df[, me_name := "vacc_dpt3_timeliness_ratio"]
  ### drop rows without antigen-specific data or logical data
  df <- df[!is.na(data)]
  df <- df[data < 1]
  
  ### drop anything without mapped ihme_loc_id
  df <- df[!is.na(ihme_loc_id)]
  
  # # remove sample size if > 5000, usually refers to target population
  df[, sample_size := gsub(",", "", sample_size) %>% as.numeric]

  df[, `:=` (age_start = round(age_start / 12),
             age_end   = round(age_end / 12))]
  df[, age_length := age_end - age_start]
  df[age_length %in% c(0, 1), age_year := age_start]
  df[is.na(age_length), age_year := 1]
  df[age_length > 1, age_year := 1]

  ### Matching year_id for MICS tabulations to how done in lit (from year_start to year_start - cohort)
  df[, year_id := year_start-age_year]
  df[, c("year_start", "year_end") := NULL]
  
  df <- df[, survey_name := "UNICEF_MICS"]
  #***********************************************************************************************************************
  
  #----SAVE---------------------------------------------------------------------------------------------------------------
  ### save
  saveRDS(df, paste0("FILEPATH"))
  print("MICS timeliness report tabulations prepped and saved for use in GBD pipeline")
  #***********************************************************************************************************************
  
  return(df)
}

### run
saved_timeliness <- mics_timeliness_reports()

