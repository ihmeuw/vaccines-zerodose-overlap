#----HEADER------------------------------------------------------------------------------------------------------------
# Author:  USERNAME
# Date:    DATE
# Purpose: Prep literature data on vaccine coverage
# Run:     source("FILEPATH") 
#**********************************************************************************************************************

#----ENVIRONMENT-------------------------------------------------------------------------------------------------------
# Clear workspace
rm(list=ls())

# Load packages
library(data.table)
library(dplyr)
library(parallel)
library(readxl)

# Source functions
username <- 'USERNAME'
source(paste0("FILEPATH"))
source(db_tools)
source("FILEPATH")

#----FUNCTIONS---------------------------------------------------------------------------------------------------------
# Classic `nin` function
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

#----GLOBAL VARIABLES--------------------------------------------------------------------------------------------------
# Get literature data
key         <- "URLPATH" 
sheet_name  <- URLPATH   
df          <- gs_download(key, sheet_name)

# Get locations
decomp_step <- "iterative"
locations   <- get_location_metadata(location_set_id = location_set_id, 
                                     gbd_round_id    = gbd_round, 
                                     decomp_step     = decomp_step)[level >= 3]

# Get vaccine schedule
cohorts <- readRDS(file.path("FILEPATH"))

# Get country-specific RotaC doses
rotac_schedule <- readRDS(file.path("FILEPATH"))
#**********************************************************************************************************************

#----1. INDICATOR ASSIGNMENT AND CLEANING------------------------------------------------------------------------------

# If have higher combination vaccine coverage, swap out component vaccines
message("Swapping out coverage from component vaccine where combination vaccine is higher")
duples <- list(c("vacc_dpt1",  "vacc_tetra1"),
               c("vacc_dpt1",  "vacc_pent1"),
               c("vacc_dpt2",  "vacc_pent2"),
               c("vacc_dpt2",  "vacc_tetra2"),
               c("vacc_dpt3",  "vacc_tetra3"),
               c("vacc_dpt3",  "vacc_pent3"),
               c("dpt3_timeliness",  "pent3_timeliness"),  
               c("vacc_hib3",  "vacc_pent3"),
               c("vacc_hib3",  "vacc_tetra3"),
               c("vacc_hepb3", "vacc_pent3"),
               c("vacc_mcv1",  "vacc_mmr1"),
               c("vacc_mcv2",  "vacc_mmr2"),
               c("vacc_rcv1",  "vacc_mmr1"),
               c("vacc_rcv2",  "vacc_mmr2"))
duples_alt_1 <- list(c("vacc_dpt1",  "vacc_tetra1"),
                     c("vacc_dpt1",  "vacc_pent1"),
                     c("vacc_dpt2",  "vacc_pent2"),
                     c("vacc_dpt2",  "vacc_tetra2"),
                     c("vacc_dpt3",  "vacc_tetra3"),
                     c("vacc_dpt3",  "vacc_pent3"),
                     c("dpt3_timeliness",  "pent3_timeliness"),  
                     c("vacc_hib3",  "vacc_pent3"),
                     c("vacc_hib3",  "vacc_tetra3"),
                     c("vacc_polio3", "vacc_pent3"),  # polio3 instead of HEPB3 in pentavalent
                     c("vacc_mcv1",  "vacc_mmr1"),
                     c("vacc_mcv2",  "vacc_mmr2"),
                     c("vacc_rcv1",  "vacc_mmr1"),
                     c("vacc_rcv2",  "vacc_mmr2"))
duples_alt_2 <- list(c("vacc_dpt1",  "vacc_tetra1"),
                     c("vacc_dpt1",  "vacc_pent1"),
                     c("vacc_dpt2",  "vacc_pent2"),
                     c("vacc_dpt2",  "vacc_tetra2"),
                     c("vacc_dpt3",  "vacc_tetra3"),
                     c("vacc_dpt3",  "vacc_pent3"),
                     c("dpt3_timeliness",  "pent3_timeliness"),  
                     c("vacc_polio3",  "vacc_pent3"),  # polio3 instead of HIB3 in pentavalent
                     c("vacc_hib3",  "vacc_tetra3"),
                     c("vacc_hepb3", "vacc_pent3"),
                     c("vacc_mcv1",  "vacc_mmr1"),
                     c("vacc_mcv2",  "vacc_mmr2"),
                     c("vacc_rcv1",  "vacc_mmr1"),
                     c("vacc_rcv2",  "vacc_mmr2"))

# Account for alternative pentavalent in select countries/nids!
message("-- Account for unconentional pentavalent formulation (country-specific)")
unconventional_penta_filepath <- file.path(ref_data_repo, "penta_unconventional_component.csv")
alt_penta                     <- fread(unconventional_penta_filepath)
# What countries/nids should be processed with duples_alt_1?
alt_penta_1 <- alt_penta[DTAP==1 & HepB==0 & HiB==1 & IPV==1]
# What countries/nids should be processed with duples_alt_2?
alt_penta_2 <- alt_penta[DTAP==1 & HepB==1 & HiB==0 & IPV==1]

for (i in duples_alt_1) {  
  # e.g. if !is.na(penta) & penta > dpt | !is.na(penta) & is.na(dpt), replace dpt with penta
  df <- df[nid %in% unique(alt_penta_1$survey_id), 
           (i[1]) := ifelse((!is.na(get(i[2])) & get(i[2]) > get(i[1])) | (!is.na(get(i[2])) & is.na(get(i[1]))), get(i[2]), get(i[1]))]
}

for (i in duples_alt_2) {
  # e.g. if !is.na(penta) & penta > dpt | !is.na(penta) & is.na(dpt), replace dpt with penta
  df <- df[nid %in% unique(alt_penta_2$survey_id), 
           (i[1]) := ifelse((!is.na(get(i[2])) & get(i[2]) > get(i[1])) | (!is.na(get(i[2])) & is.na(get(i[1]))), get(i[2]), get(i[1]))]
}

# original
for (i in duples) {
  # e.g. if !is.na(penta) & penta > dpt | !is.na(penta) & is.na(dpt), replace dpt with penta
  df <- df[!nid %in% unique(alt_penta$survey_id), 
           (i[1]) := ifelse((!is.na(get(i[2])) & get(i[2]) > get(i[1])) | (!is.na(get(i[2])) & is.na(get(i[1]))), get(i[2]), get(i[1]))]
}

# Get ihme_loc_id from location_name
df[, ihme_loc_id := tstrsplit(location_name_short_ihme_loc_id, "[|]")[[2]]]

# Change to correct ratio indicator name for dpt3 timeliness used in GBD modeling
setnames(df, "dpt3_timeliness", "vacc_dpt3_timeliness_ratio")  

# Reshape wide to long
id   <- c("nid", "file_path", "ihme_loc_id", "year_start", "year_end", "age_start", "age_end", "sample_size", "cv_admin", "cv_do_not_shift") 
vacc <-  grep("vacc", names(df), value=TRUE)
df   <- melt(df, id=id, measure=patterns("^vacc|^maternal"), value.name="data", variable.name="me_name")  

#----2. CLEANING DOSE AND SAMPLE SIZE ASSIGNMENT-----------------------------------------------------------------------

# Message user
message("Cleaning data")

# Drop rows without antigen-specific data
df <- df[!is.na(data)]

# Divide coverage by 100 (except for timeliness ratio, which is already in proportion space)
df[me_name != "vacc_dpt3_timeliness_ratio", data := (as.numeric(data) / 100)] 

# Cap coverage at 1
df[data > 1, data := 1]

# Again, drop rows without antigen-specific data that were blank and are now NA
df <- df[!is.na(data)]

# Drop anything without mapped ihme_loc_id
df <- df[!is.na(ihme_loc_id)]

# Convert 'data' col class to numeric 
df$data <- as.numeric(df$data)

# Calculate country schedule-specific rotac rows and add to primary dataset 
setnames(rotac_schedule, "me_name", "rota_complete")
rotac_schedule[, me_name := paste0("vacc_rota", doses)]
rotac <- merge(df, rotac_schedule, by=c("ihme_loc_id", "me_name"), 
               all.x = TRUE, all.y = TRUE)[!is.na(rota_complete) & !is.na(data)]
rotac[, `:=` (me_name = rota_complete, location_id = NULL, doses = NULL, rota_complete = NULL)]
df    <- rbind(df, rotac)

# Remove sample size if > 5000, usually refers to target population 
# (will be arbitrarily assigned a sample size of 100 in prep_exp.r) 
df[, sample_size := gsub(",", "", sample_size) %>% as.numeric]
df[sample_size > 5000, sample_size := NA]  

# Save the df at this point -- still with age_start and age_end in months -- for age-specific modeling prep
df_months <- copy(df)

# If age_start and age_end not listed, assuming 12-23mo as done below after round cols (save here for age-specific data prep)
df_months <- df_months[is.na(age_start) & is.na(age_end) & is.na(cv_admin), `:=` (age_start=12, age_end=23)]  
saveRDS(df_months, file.path(ref_data_repo, "vaccine_lit_months.rds"))

# Preserve original month data
df[, age_start_months  := age_start]
df[, age_end_months    := age_end]
df[, age_length_months := age_end_months - age_start_months]

# Identify special age extractions (ie age-ranges that don't fit neatly into single years)
df$special <- TRUE
df[paste(age_start_months, age_end_months, sep = "-") %in% c("12-23", "24-35", "36-47", "48-59", "60-71"), special := FALSE]

# Calculate age in years:
# - If age-range is 1 year or less, use age_start as age_year
# - If age-range is greater than 1 year, take rounded mean of age_start and age_end
# - If age_start or age_end is missing, assume data is of 12-23 month olds
df[, `:=` (age_start = round(age_start / 12),
           age_end   = round(age_end / 12))]
df[, age_length := age_end - age_start]  

df[age_length %in% c(0, 1), age_year := age_start]
df[is.na(age_length), age_year := 1]
df[age_length > 1, age_year := round((age_start + age_end)/2)]

# Get year_id from difference between survey year and age year. Do not backshift for admin data
df[cv_admin==1 | cv_do_not_shift==1, year_id := year_start]
df[is.na(cv_admin) & is.na(cv_do_not_shift), year_id := year_start-age_year]

# Drop non-admin data outside of 1-5 age range
df <- df[(age_year > 0 & age_year < 5) |
           cv_admin==1 |
           cv_do_not_shift==1]

# Account for multiple data points for given location-antigen-year using rule-based approach. Multiple data points can be due to:
# - Duplicate extractions in lit extraction sheet
# - Result of age processing logic (e.g. 12-35 and 0-59 months both evaluate to age_year 3)
df[, row_id := 1:nrow(df)]
df[, keep := TRUE]
overlap <- df[is.na(cv_admin), .N, by = c("nid", "ihme_loc_id", "me_name", "year_id")][N > 1, ]

message("Resolving duplicate extractions")

for(i in 1:nrow(overlap)) {
  
  # Message user
  if(i %% 10 == 0) {
    percentage_complete <- round(i / nrow(overlap), digits = 2) * 100
    message(paste0("-- ", percentage_complete, "%"))
  }
  
  # Set needed variables
  overlap_fixed <- FALSE
  
  i_nid         <- overlap[i, nid]
  i_ihme_loc_id <- overlap[i, ihme_loc_id]
  i_me_name     <- overlap[i, me_name]
  i_year_id     <- overlap[i, year_id]
  
  overlap_row_ids <- df[nid == i_nid & ihme_loc_id == i_ihme_loc_id & me_name == i_me_name & year_id == i_year_id, row_id]
  
  # 1. If duplicate extraction from different file path (eg. national and subnational MICS), drop duplicate
  duplicates <- df[row_id %in% overlap_row_ids, setdiff(names(df), c("file_path", "row_id")), with = FALSE] %>% duplicated
  if(any(duplicates)) {
    duplicated_row_ids <- overlap_row_ids[duplicates]
    df[row_id %in% duplicated_row_ids, keep := FALSE] 
    overlap_row_ids <- df[row_id %in% overlap_row_ids & keep == TRUE, row_id]
    if(length(overlap_row_ids) == 1) {
      overlap_fixed <- TRUE
    }
  } 
  
  # 2. If regular and special tabulation exist for the overlapping nid/location/antigen/year, use regular
  if(!overlap_fixed){
    regular_cohort_exists <- df[row_id %in% overlap_row_ids, any(!special)]
    if(regular_cohort_exists) {
      df[row_id %in% overlap_row_ids & special == TRUE, keep := FALSE]
      overlap_row_ids <- df[row_id %in% overlap_row_ids & keep == TRUE, row_id]
      if(length(overlap_row_ids) == 1) {
        overlap_fixed <- TRUE
      }
    }
  } 
  
  # 3. If overlapping tabulations are both special, use tabulation with age length closest to 1 year
  #    (e.g. 12-35 and 0-59 both evaluate to age_year 3, use 12-35)
  if(!overlap_fixed){
    special_cohorts_only <- df[row_id %in% overlap_row_ids, all(special)]
    if(special_cohorts_only) {
      different_age_lengths <- length(df[row_id %in% overlap_row_ids, unique(age_length_months)]) > 1
      if(different_age_lengths) {
        incorrect_age_length_row_ids <- df[row_id %in% overlap_row_ids, .SD[which(abs(age_length_months - 11) != min(abs(age_length_months - 11))), row_id], 
                                           .SDcols = c("nid", "row_id", "age_length_months")]
        df[row_id %in% incorrect_age_length_row_ids, keep := FALSE] 
        overlap_row_ids <- df[row_id %in% overlap_row_ids & keep == TRUE, row_id]
      }
      if(length(overlap_row_ids) == 1) {
        overlap_fixed <- TRUE
      }
    }
  }
  
  # 4. If overlapping tabulations are both special and have same age-length (e.g. 18-29 and 30-41), 
  #    use tabulation with larger sample size
  
  if(!overlap_fixed){
    special_cohorts_only <- df[row_id %in% overlap_row_ids, all(special)]
    if(special_cohorts_only) {
      same_age_lengths <- length(df[row_id %in% overlap_row_ids, unique(age_length_months)]) == 1
      if(same_age_lengths) {
        overlap_row_id_largest_sample_size <- df[row_id %in% overlap_row_ids, .SD[sample_size == max(sample_size), row_id]]
        df[row_id %in% overlap_row_ids & row_id %!in% overlap_row_id_largest_sample_size, keep := FALSE]
        overlap_row_ids <- df[row_id %in% overlap_row_ids & keep == TRUE, row_id]
        if(length(overlap_row_ids) == 1) {
          overlap_fixed <- TRUE
        }
      }
    }
  }
  
  # 5. If duplicative lit extraction still unresolved, throw error and flag user
  if(!overlap_fixed){
    error_message <- paste0("Multiple data points for ", i_year_id, " ", i_ihme_loc_id, " ", i_me_name, " from lit extraction sheet. Review lit extraction age-processing logic.")
    stop(error_message)
  }
}

message("Duplicate extractions resolved")

# Drop duplicate extractions and remove unneeded indicators
df <- df[keep == TRUE, ]
df[, c("age_start_months", "age_end_months", "year_start", "year_end", "row_id", "keep") := NULL]

#----4. SCHEDULE CHECK-------------------------------------------------------------------------------------------------

# Drop rows/cohorts of MCV/RCV survey data from younger cohorts pre-schedule
df <- merge(df, cohorts, by=c("ihme_loc_id", "me_name"), all.x=T)
# if !is.na(age_cohort) & age_start < age_cohort, mark 1 in "drop" column  # <<<<<<<< but have after admin data !!
df[!is.na(age_cohort) & ((age_end==1 & age_start != age_end) | age_end < age_cohort) & is.na(cv_admin), drop := "drop"]  # use age_end instead of age_start or age_year since aggregate age groups includes older, eligible kids
# keep all rows where drop != "drop"
df <- df[is.na(drop)]
# get rid of arbitrary columns age_cohort, drop
df[, c("age_cohort", "drop") := NULL]


#----5. SURVEY NAME ASSIGN---------------------------------------------------------------------------------------------

# Create survey name
# This needs to me improved. If survey doesn't have UNICEF_MICS in file path it isn't used in prep_exp.R...
parsed <- df[, tstrsplit(file_path, "/", fixed=TRUE)]
parsed <- parsed[, survey_name := ifelse(nchar(V3)==3, paste0(V3, "/", V4), V3)]  # <<<<<< Re-write this line so that who_csv/xxxxxx survey names still populate!
df     <- cbind(df, parsed$survey_name)
setnames(df, "V2", "survey_name")

#----6. CUSTOM ADD-ONS-------------------------------------------------------------------------------------------------

# Prepare USA NIS Kindergarten report 2009-10 through 2016-17
prep.nisreport <- function() {
  
  # NIS coverage data
  length <- year.est.end - 2009
  path   <- paste0(data_root, "FILEPATH")  # MMR data only
  nis    <- read_excel(path, sheet="SVV Coverage Trend 2016-17 Data", skip=2, 
                       col_types=c("text", rep(c("text", "skip", "text", "text", "text", "text"), length))) %>% data.table
  
  # Remove unnecessary rows and columns
  nis <- nis[!Names %in% c("HP 2020 Target", "Median")]
  nis <- nis[, colnames(nis)[grep("SURVEY TYPE|TARGET|TOTAL KINDERGARTEN POPULATION|PERCENT SURVEYED", colnames(nis))] := NULL]
  
  # Add years to colnames
  yrs <- 2009:(year.est.end - 1)
  names <- c("name", paste0("x_", yrs))
  colnames(nis) <- names
  
  # Reshape wide to long
  nis[, (paste0("coverage_", yrs)) := lapply(paste0("x_", yrs), function(cols) as.numeric(get(cols)))]
  nis <- melt(nis[, c("name", paste0("coverage_", yrs)), with=FALSE], value.name="data")
  
  # Add in year of fall of academic year
  nis[, year_start := NA_integer_]
  nis[, year_end := NA_integer_]
  nis[, year_start := substring(variable, 10, 13) %>% as.integer]
  nis[, year_end := year_start + 1]
  nis[, year_id := floor((year_start + year_end) / 2)]
  nis[, variable := NULL]
  setnames(nis, "name", "location_name")
  nis <- merge(nis[!grep("Median", location_name)], locations[parent_id==102, .(ihme_loc_id, location_name)], by="location_name", all.x=TRUE) %>%  
    .[, location_name := NULL]
  
  # Calculate population included
  nis[, data := data / 100]
  
  # Add vars
  nis[, nid := 334866]
  nis[, file_path := path]
  nis[, variance := data * (1 - data) / 50]
  nis[, cv_survey := 0]
  nis[, age_group_id := 22]
  nis[, sex_id := 3]
  nis[, survey_name := "School Vaccination Assessment Program"]
  nis <- nis[!is.na(data), ]
  nis_mmr <- copy(nis)[, me_name := "vacc_mmr2"]
  nis_mcv <- copy(nis)[, me_name := "vacc_mcv2"]
  nis <- rbind(nis_mmr, nis_mcv)
  
  # Retur data
  return(nis)
}

year.est.end <- 2017  # need to hard-code year because NID/dataset does not include after 2017
nis <- prep.nisreport()

df  <- rbind(df, nis, fill=TRUE)

#----7. SAVE------------------------------------------------------------------------------------------------------------
# Save prepped literature extractions as RDS
message(paste0("Literature extractions prepped and saved: ", file.path(ref_data_repo, "vaccine_lit.rds")))
saveRDS(df, file.path(ref_data_repo, "vaccine_lit.rds"))

#***********************************************************************************************************************
