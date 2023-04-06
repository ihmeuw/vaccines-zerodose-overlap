#----HEADER------------------------------------------------------------------------------------------------------------
# Author:  USERNAME, USERNAME, USERNAME
# Date:    DATE
# Purpose: Compare estimates produced by processing wrapper with literature reports for MICS surveys. If internal
#          estimates differ too greatly from report coverage, use literature extraction. Record decisions in flat
#          file to be used in LBD and GBD pipeline.
# 
# Section 1: Match reports with processed tabulations 
#   1) Load and prepare literature extractions
#   2) Get processed estimates for surveys with special age groups
#   3) Get processed estimates for surveys with regular age groups
# Section 2: Decider 
#   1) 
#**********************************************************************************************************************
library(data.table)
library(ggplot2)
library(magrittr)

#----GLOBAL VARIABLES--------------------------------------------------------------------------------------------------
# File paths
raw_file_path                 <- "FILEPATH"
ref_data_repo                 <- "FILEPATH"
unconventional_penta_filepath <- "FILEPATH"

# Google sheet for lit extraction download
lit_extraction_key         <- "FILEPATH"
lit_extraction_sheet_name  <- FILEPATH

# Get date for date stamped file paths
date <- Sys.Date()

#----FUNCTIONS---------------------------------------------------------------------------------------------------------
# Classic 'nin' function
'%!in%' <- function(x,y)!('%in%'(x,y))

# Modified min/max function to return NA if all NA
min_na_rm <- function(x,...) {
  if (all(is.na(x))) {
    return(as.numeric(NA))
  } else {
    return(min(x, na.rm = T))
  }
}

max_na_rm <- function(x,...) {
  if (all(is.na(x))) {
    return(as.numeric(NA))
  } else {
    return(max(x, na.rm = T))
  }
}

# Get literature data
google_sheet_download <- function(key, sheet_name) {
  if (!is.null(sheet_name)) sheet <- paste0("&gid=", sheet_name) else sheet <- ""
  return(fread(paste0("FILEPATH")))
}

# Return integer vector of extracted MICS nids
get_mics_nids <- function(raw_file_path) {
  
  # Given raw file string, return nid as integer
  get_nid_from_filename <- function(raw_file) {
    string_split_file <- unlist(strsplit(raw_file, "_"))
    nid_csv           <- string_split_file[length(string_split_file)]
    nid               <- as.integer(gsub(".csv", "", nid_csv))
    return(nid)
  }
  
  # Get list of MICS files
  all_raw_files  <- list.files(raw_file_path)
  mics_raw_files <- all_raw_files[grepl("UNICEF_MICS", all_raw_files)]
  
  # Extract nids from file names
  mics_nids <- unlist(lapply(mics_raw_files, get_nid_from_filename))
  return(mics_nids)
}

#**********************************************************************************************************************

# ===================================================
# === Section 1: Match Report with Processed Data ===
# ===================================================


# -- 1. Load and Prepare Literature Extractions -----------------------------------------------------------------------

# Get lit extraction sheet and subset to the NIDs with microdata
lit_extraction  <- google_sheet_download(lit_extraction_key, lit_extraction_sheet_name)
mics_nids       <- get_mics_nids(raw_file_path)
lit_extraction  <- lit_extraction[nid %in% mics_nids, ]

# Get location id 
lit_extraction[, ihme_loc_id := tstrsplit(location_name_short_ihme_loc_id, "[|]")[[2]]]

# -- 1.1 Distribute Multi-Vaccines ------------------------------------------------------------------------------------

# Identify surveys with alternate pentavalent formulation (IPV instead of Hib or Hepb)
alt_penta   <- fread(unconventional_penta_filepath)
alt_penta_1 <- alt_penta[DTAP==1 & HepB==0 & HiB==1 & IPV==1]
alt_penta_2 <- alt_penta[DTAP==1 & HepB==1 & HiB==0 & IPV==1]

# Set multi-dose combinations
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
duples_alt_1 <- list(c("vacc_dpt1",  "vacc_tetra1"),
                     c("vacc_dpt1",  "vacc_pent1"),
                     c("vacc_dpt2",  "vacc_pent2"),
                     c("vacc_dpt2",  "vacc_tetra2"),
                     c("vacc_dpt3",  "vacc_tetra3"),
                     c("vacc_dpt3",  "vacc_pent3"),
                     c("vacc_hib3",  "vacc_pent3"),
                     c("vacc_hib3",  "vacc_tetra3"),
                     c("vacc_polio3", "vacc_pent3"),  # polio3 instead of hepb3 in pentavalent
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
                     c("vacc_polio3", "vacc_pent3"),  # polio3 instead of hib3 in pentavalent
                     c("vacc_hib3",  "vacc_tetra3"),
                     c("vacc_hepb3", "vacc_pent3"),
                     c("vacc_mcv1",  "vacc_mmr1"),
                     c("vacc_mcv2",  "vacc_mmr2"),
                     c("vacc_rcv1",  "vacc_mmr1"),
                     c("vacc_rcv2",  "vacc_mmr2"))


# Distribute multi-vax: if multi-vax is greater than component vaccine or vaccine is missing, replace with multi-vax dose
for (i in duples) {
  lit_extraction <- lit_extraction[!nid %in% unique(alt_penta$survey_id), (i[1]) := ifelse((!is.na(get(i[2])) & get(i[2]) > get(i[1])) | (!is.na(get(i[2])) & is.na(get(i[1]))), 
                                                                                           get(i[2]), get(i[1]))]
}
for (i in duples_alt_1) {
  lit_extraction <- lit_extraction[nid %in% unique(alt_penta_1$survey_id), (i[1]) := ifelse((!is.na(get(i[2])) & get(i[2]) > get(i[1])) | (!is.na(get(i[2])) & is.na(get(i[1]))), 
                                                                                            get(i[2]), get(i[1]))]
}
for (i in duples_alt_2) {
  lit_extraction <- lit_extraction[nid %in% unique(alt_penta_2$survey_id), (i[1]) := ifelse((!is.na(get(i[2])) & get(i[2]) > get(i[1])) | (!is.na(get(i[2])) & is.na(get(i[1]))), 
                                                                                            get(i[2]), get(i[1]))]
}

# -- 1.2 Convert Wide to Long ------------------------------------------------------------------------------------------

# Preserve lit extraction for later use. Use "report" variable to wide to long
report <- suppressWarnings(melt(lit_extraction, 
                                id.vars = c("age_start", "age_end",  "ihme_loc_id",  "location_code", "nid", "parent_location", "sample_size","year_end", "year_start"),
                                measure.vars = c("vacc_bcg",    
                                                 "vacc_yfv", 
                                                 "vacc_meng", 
                                                 "vacc_mcv1",   "vacc_mcv2",
                                                 "vacc_mmr1",   "vacc_mmr2", 
                                                 "vacc_rcv1",   "vacc_rcv2", 
                                                 "vacc_dpt1",   "vacc_dpt2",   "vacc_dpt3",                       
                                                 "vacc_hepb1",  "vacc_hepb2",  "vacc_hepb3",                      
                                                 "vacc_hib1",   "vacc_hib2",   "vacc_hib3",        
                                                 "vacc_pcv1",   "vacc_pcv2",   "vacc_pcv3",
                                                 "vacc_pent1",  "vacc_pent2",  "vacc_pent3",
                                                 "vacc_polio1", "vacc_polio2", "vacc_polio3",                    
                                                 "vacc_rota1",  "vacc_rota2",  "vacc_rota3",
                                                 "vacc_tetra1", "vacc_tetra2", "vacc_tetra3")))


# Drop missing data
report <- report[!is.na(value) & value != "", ]   
# Also drop data that isn't tagged to an ihme_loc_id (i.e. lit-extracted admin2 data)
report <- report[!is.na(ihme_loc_id), ]
# Standardize columns to merge with processed data
setnames(report, "year_start", "survey_yr")
# report[, survey_yr := year_start]

# -- 1.3 Calculate RotaC ----------------------------------------------------------------------------------------------

# Calculate rota complete and add to report data
rotac_schedule <- readRDS(file.path(ref_data_repo, "vaccine_schedule.rds"))
setnames(rotac_schedule, "me_name", "rota_complete")
rotac_schedule[, variable := paste0("vacc_rota", doses)]
rotac <- merge(report, rotac_schedule, by=c("ihme_loc_id", "variable"), all.x=T)[!is.na(rota_complete) & !is.na(value)]
rotac[, `:=` (variable=rota_complete, location_id=NULL, doses=NULL, rota_complete=NULL)]
report <- rbind(report, rotac)

# Get copy of report to preserve original
for_now <- copy(report)  


# -- 2. Get Processed Estimates for Special Age Groups ----------------------------------------------------------------

# *EXPLANATION* 
# To compare processed estimates against the report estimates, it's necessary that the processed estimates match the 
# same age group as the reports. Identify special age groups in report (ie 18-27 months) and process survey microdata 
# to tabulate to the same age group (as opposed to the default processing to year-based age groups)

# Get surveys with special age groups
special_age_lit_extractions <- lit_extraction[age_start > 12, ]
special_age_lit_extractions <- unique(special_age_lit_extractions[, .(nid, age_start, age_end)])

# -- 2.1 Prepare for Processing ---------------------------------------------------------------------------------------

# Setup
username        <- 'USERNAME'
core_repo       <- paste0("FILEPATH")
vaccines_repo   <- paste0("FILEPATH")
extraction_root <- "FILEPATH"

# Load packages
source('FILEPATH')
lapply(c("proto", "findpython", "argparse", "data.table", "magrittr", "survey", "parallel", "plyr", "dplyr", "haven",
         "rgeos", "raster", "rgdal", "dismo", "gbm", "foreign", "doParallel", "grid", "gridExtra", "gtools", "ggplot2",
         "sf", "assertthat", "INLA", "seegSDM", "seegMBG", "pacman", "glmnet", "RMySQL"), library, character.only = TRUE)

# Load packages directly from file path for packages not yet available to LBD Singularity image
library("fasterize", lib.loc="FILEPATH", character.only=TRUE)
library("binom", lib.loc="FILEPATH", character.only=TRUE)

# Load custom functions for special age group processing
source(file.path("FILEPATH"))

#--- 2.2 Process Special Age Groups -----------------------------------------------------------------------------------

# Settings
topic       <- "vaccines"
config.path <- "FILEPATH"
parallel    <- FALSE


# Process, tabulate, and rbind all MICS surveys with unconventional age brackets (long)
special_ages <- data.table()
for (i in 1:nrow(special_age_lit_extractions)){
  
  # Process
  message("\nNow working on ", i)
  nid       <- special_age_lit_extractions[i, nid]
  age_start <- special_age_lit_extractions[i, age_start]
  age_end   <- special_age_lit_extractions[i, age_end]
  test      <- process_mics(nid)
  
  # Tabulate
  if (length(vaccines[vaccines %in% names(test)]) > 0) {
    dataset     <- prep_for_tabulation_mics_special_age(nid, team="gbd", vaccines.=c(vaccines, "rotac"), age_start=age_start, age_end=age_end)
    # launch collapse
    dataset  <- dataset[[2]] %>% as.data.table
    config   <- load.settings(config.path, topic)
    tab_data <- collapse.run.mics.special.age(dataset, config=config)
  } else {
    tab_data <- data.table()
  }
  
  # Combine tabulated data
  special_ages <- rbind(special_ages, tab_data, fill=T)
}

# Explictly check to make sure all NIDs got processed to avoid unnoticed permissions issues
if (length(unique(special_ages$nid)) != length(unique(special_age_lit_extractions$nid))){
  stop("Not all NIDs got processed. Look at your error message. Probably a permissions issue?")
}

#--- 2.3 Combine and Save Special Age Groups --------------------------------------------------------------------------

# Combine report data with specially processed microdata for special age groups
special_ages <- merge(special_ages, special_age_lit_extractions, by="nid")
setnames(special_ages, c("mean"), c("data"))

# Save processed estimates of special ages
dir.create(paste0("FILEPATH"))
fwrite(special_ages, paste0("FILEPATH"))  


#--- 3. Get Processed Estimates for Regular Age Groups ----------------------------------------------------------------
# *EXPLANATION* 
# Determine processed estimates for regular age groups. 

# Setup
username    <- Sys.info()[["user"]]
decomp_step <- "iterative"
source("FILEPATH")  
source(paste0("FILEPATH"))

# Get location metadata
locs <- get_location_metadata(location_set_id=location_set_id, release_id = release_id)[level >= 3, ]  
# locs <- get_location_metadata(location_set_id=location_set_id, gbd_round_id=gbd_round, decomp_step=decomp_step)[level >= 3, ]  

#--- 3.1 Get Data -----------------------------------------------------------------------------------------------------
# Load processed data for regular age groups
regular_ages <- data.table()
for (mics_nid in mics_nids) {
  tabulated_survey_path <- paste0("FILEPATH")
  if (file.exists(tabulated_survey_path)) {
    temp         <- fread(tabulated_survey_path)
    regular_ages <- rbind(regular_ages, temp, fill = TRUE)  
  }
}
setnames(regular_ages, "mean", "data")

#--- 3.2 Prep Data ----------------------------------------------------------------------------------------------------

regular_ages[, variance := standard_error ^ 2]  # Takes sample size into account here
regular_ages[, age_group_id := 22]
regular_ages[, cv_survey := 1]

# Drop messed up locations
drop.locs <- regular_ages[!(ihme_loc_id %in% locs$ihme_loc_id), ihme_loc_id] %>% unique
if (length(drop.locs) > 0) {
  print(paste0("UNMAPPED LOCATIONS (DROPPING): ", toString(drop.locs)))
  regular_ages <- regular_ages[!(ihme_loc_id %in% drop.locs)]
}

# If mean==0/1 use Wilson Interval Method 
if (nrow(regular_ages[data %in% c(0, 1)]) > 0) {
  regular_ages.w <- regular_ages[data %in% c(0, 1)]
  sample_size <- regular_ages.w$sample_size
  n  <- ifelse(regular_ages.w$data==0, 0, sample_size)
  ci <- binom.confint(n, sample_size, conf.level=0.95, methods="wilson")
  se <- (ci$upper - ci$lower) / 3.92
  variance <- se ^ 2 * 2.25 ## Inflate by design effect according to DHS official/unofficial rule of thumb 
  regular_ages[data %in% c(0, 1)]$variance <- variance
}

# If design effect is unreasonably small, assume 2.25 and readjust
regular_ages <- regular_ages[design_effect < 1, variance := variance * 2.25 / design_effect]
regular_ages <- regular_ages[design_effect < 1, design_effect := 2.25]

regular_ages <- regular_ages[age_year == 1, ]
regular_ages[, `:=` (age_bin=NULL, age_bin_agg=NULL, age_bin_agg_id=NULL)]
regular_ages$age_start <- 12
regular_ages$age_end   <- 23
regular_ages <- regular_ages[survey_name == "UNICEF_MICS", ]

# Which surveys have only weird tabulations that we want to be dropping?
both                <- unique(for_now[, .(nid, age_start)])
two_ages            <- both$nid[which(duplicated(both$nid))]
single_ages         <- both[nid %!in% two_ages, ]
single_ages_to_drop <- single_ages[age_start != 12, nid]

# Remove single, non-conventional age groups since those are accounted for in 'special_ages' 
regular_ages <- regular_ages[nid %!in% single_ages_to_drop, ] 


#--- 4. Combine Report Data with Processed Data -----------------------------------------------------------------------
# Combine processed data for regular and special age groups
new <- rbind(regular_ages, special_ages, fill=T)
new <- new[me_name %in% unique(for_now$variable), ]
new[, parent_location := gsub("_*[0-9]", "", ihme_loc_id)]
new <- unique(new) 

for_now <- unique(for_now)

test <- merge(for_now, new,  # merge the cleaned lit data with the combined special + conventional microdata tabulations
              by.x=c("parent_location","ihme_loc_id", "survey_yr", "variable", "nid", "age_start", "age_end"), 
              by.y=c("parent_location","ihme_loc_id", "year_start", "me_name", "nid", "age_start", "age_end"), all.x=T, all.y=T)

keep <- test[, .(parent_location, ihme_loc_id, year_start, nid, variable, age_start, age_end, sample_size.x, value, data)]
colnames(keep) <- c("parent_location", "ihme_loc_id", "year_start", "nid", "me_name","age_start", "age_end", "sample_size", "report", "tabulation")


keep[, report := (as.numeric(report) / 100)] 
keep[, difference := (report - tabulation)]

hist(keep$difference, breaks=100)

keep <- keep[!is.na(tabulation) | !is.na(report), ]


# Save merged special and regular estimates merged with report estimates in date-versioned folder and reference files repo
fwrite(keep, paste0("FILEPATH"))
fwrite(keep, file.path("FILEPATH"))

#**********************************************************************************************************************

# ===================================================
# === Section 2: Decider                          ===
# ===================================================

#--- 1. Decider Set-Up ------------------------------------------------------------------------------------------------

# Settings
threshold        <- 0.1 # 10% threshold for abs difference for now
subnat_threshold <- 0.9 # spearman rank corr coeff of 0.9 or higher

# Get prepared data
mics_compare_df <- fread(paste0("FILEPATH"))

# Folder for diagnostics
diagnostic_output_dir <- paste0("FILEPATH")
ifelse(!dir.exists(diagnostic_output_dir), dir.create(diagnostic_output_dir), FALSE)

# Folder for final outputs
final_output_dir <- paste0("FILEPATH")
ifelse(!dir.exists(final_output_dir), dir.create(final_output_dir), FALSE)

#--- 2. Archive & Preprocess ------------------------------------------------------------------------------------------

# Identify subnational ids that are for a single location and don't include
# national to rbind to the national set - should be treated the same

# Archive mics compare df
fwrite(mics_compare_df, paste0("FILEPATH"))

# Split into national and subnational
mics_nat_df    <- subset(mics_compare_df, !grepl("_", ihme_loc_id))
mics_subnat_df <- subset(mics_compare_df, grepl("_", ihme_loc_id))

# Treat single-location, subnational-only MICS NIDS (e.g. 104042 and 104043 from IDN with NO national IDN) the same as national MICS to get independent decisions
mics_nat_df_nids <- unique(mics_nat_df$nid)

# Identify the subnational-only NIDs but remove those with multiple subnationals
add_to_nat_df <- mics_subnat_df[!nid %in% mics_nat_df_nids]
actually_still_subnat <- data.table()
for (n in unique(add_to_nat_df$nid)) {
  test <- add_to_nat_df[nid==n]
  test <- test[, count := length(unique(test$ihme_loc_id))]
  if(unique(test$count) > 1){
    actually_still_subnat <- rbind(actually_still_subnat, test)
  }
}
actually_still_subnat[, count := NULL]

# NIDs that should be considered subnational (all those with a national tabulation and those representing multiple subnational units without national)
mics_nat_df_nids <- c(mics_nat_df_nids, unique(actually_still_subnat$nid))

# Remove the subnational-only NIDs from mics_subnat_df
add_to_nat_df  <- add_to_nat_df[!nid %in% mics_nat_df_nids]
mics_subnat_df <- mics_subnat_df[nid %in% mics_nat_df_nids]

# Rbind the subnational-only NIDs to mics_nat_df
mics_nat_df <- rbind(mics_nat_df, add_to_nat_df)

if (nrow(mics_nat_df) + nrow(mics_subnat_df) != nrow(mics_compare_df)) {
  stop("Something is going wrong in subnational assignments!")
}

# For ease of processing, remove any extraneous tabulations -- i.e. if 12-23 month tabs
# are in the comparison DF, but there is only a later age cohort with report data,
# we don't need the 12-23 month tabs for the national comparisons

mics_nat_df[is.na(report) & !is.na(tabulation), tab_only := TRUE]
mics_nat_df[, rows_per_me_nid := .N, by = c("nid", "me_name")]
mics_nat_df[, any_comparisons := sum(!is.na(difference)), by = c("nid", "me_name")]
mics_nat_df[tab_only == TRUE & rows_per_me_nid > 1 & any_comparisons > 0, remove_extraneous := TRUE]

mics_nat_df <- subset(mics_nat_df, is.na(remove_extraneous))
mics_nat_df[, tab_only := NULL]
mics_nat_df[, rows_per_me_nid := NULL]
mics_nat_df[, remove_extraneous := NULL]


#--- 3. Calculate Subnational Correlations ----------------------------------------------------------------------------

# Note that the subnational data set is really just for these sub-analyses
# Everything will be flagged at the national level, however, then collapsed to
# one decision per NID / me_name pair for final outputs

# Identify subnational NIDS
subnat_nids <- unique(mics_subnat_df$nid)

# Flag these in the national data set (which will become our final data set)
mics_nat_df[nid %in% subnat_nids, subnational := TRUE]
mics_nat_df[subnational != TRUE, subnational := FALSE]

# Double check this to ensure no additional surveys from the countries
# that have subnational surveys are present
subnat_loc_ids <- unique(mics_subnat_df$ihme_loc_id)
if (nrow(mics_nat_df[ihme_loc_id %in% subnat_loc_ids & subnational == FALSE])) {
  warning("There are some surveys from GBD subnational countries with only national report extractions - check CSV")
  fwrite(mics_nat_df[ihme_loc_id %in% subnat_loc_ids & subnational == FALSE],
         file = paste0("FILEPATH"))
}

# For each subnational NID, make a plot of report vs tabulation
dir.create(paste0("FILEPATH"))

df_cors <- lapply(unique(mics_subnat_df$nid), function(a_nid) {
  
  df_subset <- subset(mics_subnat_df, nid == a_nid)
  parent_loc <- unique(df_subset$parent_location)
  df_subset[, age_cohort := paste0(age_start, "_", age_end)]
  
  age_cors <- lapply(unique(df_subset$age_cohort), function(ac) {
    
    df_age <- df_subset[age_cohort == ac]
    
    # cor() doesn't handle all-NA vectors well,
    # so here's an improved version
    cor_na_rm <- function(x, y, ...) {
      if (all(is.na(x) + is.na(y))) {
        return(as.numeric(NA))
      } else {
        return(cor(x, y, use="complete.obs", ...))
      }
    }
    
    df_age[, subnat_sp_cor := cor_na_rm(x=report, y=tabulation,
                                        method = 'spearman'),
           by = "me_name"]
    
    # Table of correlations for plotting
    # Plotting here is somewhat inelegant but we have the data already prepped...
    cors <- unique(subset(df_age, select = c("me_name", "subnat_sp_cor",
                                             "age_start", "age_end",
                                             "nid")))
    cors[, a_label := paste0("sp_cor: ", round(subnat_sp_cor, 2))]
    
    if (nrow(df_age[!is.na(subnat_sp_cor)]) > 0) {
      
      gg_age <- ggplot(df_age,
                       aes(x = report, y = tabulation)) +
        geom_point(alpha = 0.5) +
        geom_abline() +
        geom_text(data = cors, x = 1, y = 0, aes(label = a_label), hjust=1) +
        theme_bw() +
        geom_smooth(method = "lm") +
        theme_bw() +
        xlim(0,1) +
        ylim(0,1) +
        facet_wrap(~me_name) +
        coord_equal() +
        labs(title = "Report vs tabulation",
             subtitle = paste(parent_loc , " | ", a_nid, " | ages: ", ac),
             x = "Report",
             y = "Tabulation")
      
      png(file = paste0("FILEPATH"),
          height = 10,
          width = 18,
          units = "in",
          res = 300)
      print(gg_age)
      dev.off()
      
    }
    
    # Drop label
    cors[, a_label := NULL]
    
    # Return cors object with all correlations for the age cohort / nid / me_name combos
    return(cors)
    
  })
  
  age_cors <- rbindlist(age_cors)
  
  return(age_cors)
  
})

# Rbind these rows together to create a data frame of correlations by survey & nid
df_cors <- rbindlist(df_cors)

# Now merge the correlation coefficients back on to the subnational surveys
mics_nat_df <- merge(mics_nat_df, df_cors,
                     by = c("me_name", "age_start", "age_end", "nid"),
                     all.x = TRUE)

# 3. Apply thresholds (national and subnational) by row

# 3a. National
mics_nat_df[abs(difference) > threshold & !is.na(difference), row_national_ok := FALSE]
mics_nat_df[abs(difference) <= threshold & !is.na(difference), row_national_ok := TRUE]

# 3b. Subnational
mics_nat_df[subnational == TRUE & subnat_sp_cor >= subnat_threshold & !is.na(subnat_sp_cor),
            row_subnational_ok := TRUE]

mics_nat_df[subnational == TRUE & subnat_sp_cor < subnat_threshold & !is.na(subnat_sp_cor),
            row_subnational_ok := FALSE]

#--- 4. Create Helper Metrics -----------------------------------------------------------------------------------------
# Create some related metrics to help with decision-making where limited data available

# 4a. Working with nid/me combos with more than one row, i.e.  more than one age cohort for a given nid and me_name
mics_nat_df[, n_rows_by_me_nid := .N, by = c("nid", "me_name")]

# Calculate maximum difference (national level) by me and NID
mics_nat_df[, max_abs_diff_by_me_nid := max_na_rm(abs(difference), na.rm = T), by = c("nid", "me_name")]

# Apply threshold to this nid/me_name level variable
mics_nat_df[abs(max_abs_diff_by_me_nid) > threshold & !is.na(max_abs_diff_by_me_nid), me_nid_national_ok := FALSE]
mics_nat_df[abs(max_abs_diff_by_me_nid) <= threshold & !is.na(max_abs_diff_by_me_nid), me_nid_national_ok := TRUE]

# 4b. Similarly, do this for subnationals & correlations

# Calculate the lowest spearman correlation for each nid/me combo (ie by age level)
mics_nat_df[, min_subnat_sp_cor_by_me_nid := min_na_rm(subnat_sp_cor), by = c("nid", "me_name")]

# Apply threshold to the lowest spearman correlation
mics_nat_df[subnational == TRUE & min_subnat_sp_cor_by_me_nid >= subnat_threshold & !is.na(min_subnat_sp_cor_by_me_nid),
            me_nid_subnational_ok := TRUE]
mics_nat_df[subnational == TRUE & min_subnat_sp_cor_by_me_nid < subnat_threshold & !is.na(min_subnat_sp_cor_by_me_nid),
            me_nid_subnational_ok := FALSE]


# 4c. Calculate how many MEs are classified as OK by the above procedures out of all MEs per NID
national_by_me_nid <- unique(subset(mics_nat_df, select = c("nid", "me_name", "me_nid_national_ok")))
national_by_me_nid[, n_mes_in_nid := .N, by = nid]
national_by_me_nid[, n_mes_ok_in_nid := sum(me_nid_national_ok == TRUE, na.rm = T), by = nid]
national_by_me_nid[, n_mes_not_ok_in_nid := sum(me_nid_national_ok == FALSE, na.rm = T), by = nid]

# Merge this back on to our master data set
mics_nat_df <- merge(mics_nat_df, national_by_me_nid, all.x = T, by = c("nid", "me_name", "me_nid_national_ok"))

# Calculate percent of mes that are OK out of those with comparisons available, by NID
mics_nat_df[, pct_ok_mes := n_mes_ok_in_nid / (n_mes_ok_in_nid + n_mes_not_ok_in_nid)]

#--- 5. Make Decisions ------------------------------------------------------------------------------------------------
# Make classifications for a given me/nid combo

# 5a. If OK by national threshold by nid/me group and not a subnational survey, use microdata
mics_nat_df[me_nid_national_ok == TRUE & (subnational %in% c(FALSE, NA)), decision := "use_microdata"]

# 5b. If not OK by national threshold, use report extractions regardless of subnational status
mics_nat_df[me_nid_national_ok == FALSE, decision := "use_report_extractions"]

# 5c. If OK by national threshold and OK by subnational threshold, use microdata
mics_nat_df[me_nid_national_ok == TRUE & me_nid_subnational_ok == TRUE & subnational == TRUE, decision := "use_microdata"]

# 5d. If OK by national threshold but fails subnational threshold, use report extraction
mics_nat_df[me_nid_national_ok == TRUE & me_nid_subnational_ok == FALSE & subnational == TRUE, decision := "use_report_extractions"]

# 5e. If OK for national threshold but NA for subnational threshold, use microdata <- CHALLENGE THIS?
mics_nat_df[me_nid_national_ok == TRUE & is.na(me_nid_subnational_ok) & subnational == TRUE, decision := "use_microdata"]

# 5f(1). If there is a tabulation but no report data available, OK to use microdata
# only if > 90% of other mes meet thresholds for that NID <- CHALLENGE THIS?
mics_nat_df[is.na(me_nid_national_ok) & is.na(report) & !is.na(tabulation) & pct_ok_mes >= 0.9 & !is.nan(pct_ok_mes),
            decision := "use_microdata"]

# 5f(2). In the above case, if the percent agreement is < 90% for other mes in this NID, drop this data <- CHALLENGE THIS?
mics_nat_df[is.na(me_nid_national_ok) & is.na(report) & !is.na(tabulation) & pct_ok_mes < 0.9 & !is.nan(pct_ok_mes),
            decision := "drop"]

# 5f(3). In the same case, if there are no comparisons available, preserve this data <- CHALLENGE THIS?
mics_nat_df[is.na(me_nid_national_ok) & is.na(report) & !is.na(tabulation) & is.nan(pct_ok_mes),
            decision := "use_microdata"]

# 5g. If we have report and not tabulation data and no other comparisons available, first check to see if there are
# any tabulations -- use these and drop report data if available tabulations & passing checks, but otherwise
# use the report data

# If there is tabulation data for other age cohorts, use tabulation if passing checks and drop report. Otherwise, 
# If the other tabulations aren't passing or if there are no tabulations, use report

if (nrow(mics_nat_df[is.na(tabulation) & !is.na(report)]) > 0) {
  
  fwrite(mics_nat_df[is.na(tabulation) & !is.na(report)],
         file = paste0(diagnostic_output_dir, "surveys_with_reports_but_no_tabs.csv"))
  
  # Label temporarily for convenience
  mics_nat_df[is.na(tabulation) & !is.na(report), report_no_tab := TRUE]
  
  # If there are other ages for the same me, use those decisions for the report-no-tab rows
  mics_nat_df[n_rows_by_me_nid > 1,
              any_valid_microdata := any(decision == "use_microdata", na.rm = T),
              by = c("me_name", "nid")]
  
  # Drop rows where there are valid microdata for that nid/me in other age cohorts
  mics_nat_df[n_rows_by_me_nid > 1 & report_no_tab == TRUE &  is.na(decision) & any_valid_microdata == TRUE,
              decision := "drop"]
  
  # Otherwise, use the report extractions (if there are no other data for that me,
  # or if the decision was made not to use microdata for other rows for that me/nid combo
  # where comparisons were possible)
  mics_nat_df[report_no_tab == TRUE & (any_valid_microdata == FALSE | n_rows_by_me_nid == 1),
              decision := "use_report_extractions"]
  
  # Clean up
  mics_nat_df[, report_no_tab := NULL]
  
}

# Check to ensure everything classified at the national level
if (nrow(mics_nat_df[is.na(decision)]) > 0) {
  warning("There are unclassified surveys! See csv files.")
  fwrite(mics_nat_df[is.na(decision)],
         file = paste0(diagnostic_output_dir, "unclassified_surveys.csv"))
}

#--- 6. Final Decisions -----------------------------------------------------------------------------------------------
# Final decisions, add on subnationals 

# Write  me/nid files with all data from the national and subnational level processing
fwrite(mics_nat_df, file = paste0(diagnostic_output_dir, "mics_compare_national_processed.csv"))
fwrite(mics_subnat_df, file = paste0(diagnostic_output_dir, "mics_compare_subnational_processed.csv"))

# Generate the final output with decisions
decision_df <-  subset(mics_nat_df,
                       select = c("me_name", "ihme_loc_id",
                                  "nid", "decision")) %>% unique()

# 6a. Check for NID - antigen combos not found in the original data set,
# particularly subnationals, and make decisions about those as needed.
# This could be integrated earlier in the script, probably
final_nids <- unique(subset(mics_nat_df, select = c("nid", "me_name")))
final_nids[, me_nid := paste0(me_name, "_", nid)]
final_me_nids <- unique(final_nids$me_nid)

original_nids <- unique(subset(mics_compare_df, select = c("nid", "me_name")))
original_nids[, me_nid := paste0(me_name, "_", nid)]
orig_me_nids <- unique(original_nids$me_nid)

subnational_nids <- unique(subset(mics_subnat_df, select = c("nid", "me_name")))
subnational_nids[, me_nid := paste0(me_name, "_", nid)]
subnat_me_nids <- unique(subnational_nids$me_nid)

# Missing subnational me_nids
missing_subnat_me_nids <- subset(subnational_nids, me_nid %in% subnat_me_nids[!(subnat_me_nids %in% final_me_nids)])

# Run a function to make a decision about each subnational data point that isn't
# in the national level data. The logic is this:
# 1) If there is good correlation or no ability to calculate correlation,
#    then use microdata as long as you have at least one tabulation row (indicating that
#    there is some microdata)
# 2) If there is poor correlation, or if there are no tabulations (no microdata), then
#    use the report data instead

subnat_decisions <- lapply(1:nrow(missing_subnat_me_nids), function(i) {
  
  a_me <- missing_subnat_me_nids[i, me_name]
  a_nid <- missing_subnat_me_nids[i, nid]
  
  df_subset <- mics_subnat_df[me_name == a_me & nid == a_nid]
  cor_subset <- df_cors[me_name == a_me & nid == a_nid]
  
  # Determine the minimal non-na spearman correlation coefficient
  if (any(!is.na(cor_subset$spearman))) {
    a_min_cor <- min(cor_subset$spearman, na.rm = T)
  } else {
    a_min_cor <- as.numeric(NA)
  }
  
  # Determine if these are mostly report or tabulation data
  tab_rows <- nrow(df_subset[!is.na(tabulation)])
  
  # If there is tabulation data, try to use it -- unless the minimum spearman
  # correlation is below our threshold (CHALLENGE THIS - may want instead to try to
  # use as much data as possible, ie if 2 tab rows as 20 report rows?)
  
  if ((a_min_cor >= subnat_threshold | is.na(a_min_cor)) & tab_rows >= 1) {
    decision <- "use_microdata"
  } else if ((!is.na(a_min_cor) & a_min_cor < subnat_threshold) | (is.na(a_min_cor) & tab_rows == 0)) {
    decision <- "use_report_extractions"
  }
  return(data.table(me_name = a_me,
                    nid = a_nid,
                    decision = decision,
                    ihme_loc_id = unique(df_subset$parent_location)))
})

subnat_decisions <- rbindlist(subnat_decisions)

decision_df <- rbind(decision_df, subnat_decisions, fill = T, use.names = T)
decision_df <- subset(decision_df, nid != 156268) 
decision_df <- decision_df[which(decision_df$ihme_loc_id != ""),]
# 5b. Process each NID to make decisions for ratios
decision_df <- subset(decision_df, nid %!in% c(27022, 141914)) 
ratio_decisions <- function(n) {
  message(n)
  # A. Pull decisions for each element of a ratio
  df_subset <- decision_df[nid == n]
  
  ratio_table <- data.table(rbind(c("hepb3", "dpt3"),
                                  c("dpt1", "dpt3"),
                                  c("hib3", "dpt3"),
                                  c("pcv3", "dpt3"),
                                  c("rotac", "dpt3"),
                                  c("mcv2", "mcv1"),
                                  c("hepb3", "dpt3"),
                                  c("rcv1", "mcv1")))
  names(ratio_table) <- c("numerator", "denominator")
  
  ratio_decisions <- lapply(1:nrow(ratio_table), function(i) {
    num <- ratio_table[i, numerator]
    den <- ratio_table[i, denominator]
    
    num_decision <- df_subset[me_name == paste0("vacc_", num), decision]
    den_decision <- df_subset[me_name == paste0("vacc_", den), decision]
    
    # NAs if both elements of ratios don't exist
    if ((length(num_decision) == 0) | (length(den_decision) == 0)) {
      num_decision <- NA
      den_decision <- NA
    }
    
    return(data.table(num = num,
                      den = den,
                      num_decision = num_decision,
                      den_decision = den_decision))
  }) %>% rbindlist()
  
  ratio_decisions[, me_name := paste0("ratio_", num, "_", den)]
  
  # B.Classify the ratios into a single decision
  
  # Agreement: report
  ratio_decisions[num_decision == "use_report_extractions" & den_decision == "use_report_extractions", decision := "use_report_extractions"]
  
  # Agreement: microdata
  ratio_decisions[num_decision == "use_microdata" & den_decision == "use_microdata", decision := "use_microdata"]
  
  # Numerator: report; denominator: microdata
  ratio_decisions[num_decision == "use_report_extractions" & den_decision == "use_microdata", decision := "use_report_extractions"]
  
  # Numerator: microdata; denominator: report
  ratio_decisions[num_decision == "use_microdata" & den_decision == "use_report_extractions", decision := "use_report_extractions"]
  
  # Format & return decisions for each ratio
  ratio_decisions[, nid := n]
  ratio_decisions[, ihme_loc_id := unique(df_subset$ihme_loc_id)]
  
  ratio_decisions <- subset(ratio_decisions,
                            !is.na(decision),
                            select = c("me_name", "nid", "ihme_loc_id", "decision"))
  
  return(ratio_decisions)
  
}

ratio_decision_df <- lapply(unique(decision_df$nid), ratio_decisions) %>% rbindlist



# Add these ratio decisions to the overall list of decisions
decision_df <- rbind(decision_df, ratio_decision_df, fill = T, use.names = T)
decisions_pre_multi_dose_standardizing <- copy(decision_df)

##### Standardize decisions across multi-dose vaccines
# DPT: If either 1st or last dose is 'report', set all doses to report
# Non-DPT multi-dose vaccines (polio, hib, hepb, pcv, rota): set all doses to decision from final dose
# MCV & RCV: dose-independent decisions
# BCG/YFV: N/A

# DPT
dpt_report_nids    <- decision_df[me_name %in% c("vacc_dpt1", "vacc_dpt3") & decision == "use_report_extractions", unique(nid)]
dpt_microdata_nids <- decision_df[me_name %in% c("vacc_dpt1", "vacc_dpt3") & !(nid %in% dpt_report_nids), unique(nid)]

decision_df[nid %in% dpt_report_nids & me_name %in% c("vacc_dpt1", "vacc_dpt2", "vacc_dpt3"), decision := "use_report_extractions"]
decision_df[nid %in% dpt_microdata_nids & me_name %in% c("vacc_dpt1", "vacc_dpt2", "vacc_dpt3"), decision := "use_microdata"]

# Hib, HepB, PCV, Rota
stems <- c("hib", "hepb", "pcv", "polio")
for (stem in stems) {
  
  # Get decision from last dose
  final_dose                    <- paste0("vacc_", stem, "3")
  final_dose_decision_microdata <- decision_df[me_name == final_dose & decision == "use_microdata", unique(nid)]
  final_dose_decision_report    <- decision_df[me_name == final_dose & decision == "use_report_extractions", unique(nid)]
  
  # Set decision for all doses to decision from last dose
  decision_df[grepl(stem, me_name) & !grepl("ratio", me_name) & nid %in% final_dose_decision_microdata, decision := "use_microdata"]
  decision_df[grepl(stem, me_name) & !grepl("ratio", me_name) & nid %in% final_dose_decision_report, decision := "use_report_extractions"]
  
  # If the third dose is missing, set decision of 1st and 2nd dose to decision from 1st dose (edge case)
  first_dose               <- paste0("vacc_", stem, "1")
  nids                     <- decision_df[grepl(stem, me_name) & !grepl("ratio", me_name), unique(nid)]
  nids_with_final_dose     <- decision_df[me_name == final_dose, unique(nid)]
  nids_missing_final_dose  <- nids[!nids %in% nids_with_final_dose]
  
  first_dose_decision_microdata <- decision_df[grepl(stem, me_name) & !grepl("ratio", me_name) & nid %in% nids_missing_final_dose & decision == "use_microdata", unique(nid)]
  first_dose_decision_report    <- decision_df[grepl(stem, me_name) & !grepl("ratio", me_name) & nid %in% nids_missing_final_dose & decision == "use_report_extractions", unique(nid)]
  
  decision_df[grepl(stem, me_name) & !grepl("ratio", me_name) & nid %in% first_dose_decision_microdata, decision := "use_microdata"]
  decision_df[grepl(stem, me_name) & !grepl("ratio", me_name) & nid %in% first_dose_decision_report, decision := "use_report_extractions"]
}

# Rota
# Identify last dose for each nid (in some cases it's "c", in others it's 3 or 2)
rota_last_dose_c_nids <- decision_df[me_name == "vacc_rotac", unique(nid)]
rota_last_dose_3_nids <- decision_df[me_name == "vacc_rota3" & !(nid %in% rota_last_dose_c_nids), unique(nid)]
rota_last_dose_2_nids <- decision_df[me_name == "vacc_rota2" & !(nid %in% c(rota_last_dose_c_nids, rota_last_dose_3_nids)), unique(nid)]
rota_last_dose_1_nids <- decision_df[me_name == "vacc_rota1" & !(nid %in% c(rota_last_dose_c_nids, rota_last_dose_3_nids, rota_last_dose_2_nids)), unique(nid)]

# Identify nids where decision for last dose is report vs. microdata
rota_last_dose_c_decision_microdata <- 
  decision_df[grepl("rota", me_name) & !grepl("ratio", me_name) & nid %in% rota_last_dose_c_nids & decision == "use_microdata", unique(nid)]
rota_last_dose_c_decision_report <- 
  decision_df[grepl("rota", me_name) & !grepl("ratio", me_name) & nid %in% rota_last_dose_c_nids & decision == "use_report_extractions", unique(nid)]

rota_last_dose_3_decision_microdata <- 
  decision_df[grepl("rota", me_name) & !grepl("ratio", me_name) & nid %in% rota_last_dose_3_nids & decision == "use_microdata", unique(nid)]
rota_last_dose_3_decision_report <- 
  decision_df[grepl("rota", me_name) & !grepl("ratio", me_name) & nid %in% rota_last_dose_3_nids & decision == "use_report_extractions", unique(nid)]

rota_last_dose_2_decision_microdata <- 
  decision_df[grepl("rota", me_name) & !grepl("ratio", me_name) & nid %in% rota_last_dose_2_nids & decision == "use_microdata", unique(nid)]
rota_last_dose_2_decision_report <- 
  decision_df[grepl("rota", me_name) & !grepl("ratio", me_name) & nid %in% rota_last_dose_2_nids & decision == "use_report_extractions", unique(nid)]

rota_last_dose_1_decision_microdata <- 
  decision_df[grepl("rota", me_name) & !grepl("ratio", me_name) & nid %in% rota_last_dose_1_nids & decision == "use_microdata", unique(nid)]
rota_last_dose_1_decision_report <- 
  decision_df[grepl("rota", me_name) & !grepl("ratio", me_name) & nid %in% rota_last_dose_1_nids & decision == "use_report_extractions", unique(nid)]

# Set all doses to decision from last dose
decision_df[grepl("rota", me_name) & !grepl("ratio", me_name) & nid %in% rota_last_dose_c_decision_microdata, decision := "use_microdata"]
decision_df[grepl("rota", me_name) & !grepl("ratio", me_name) & nid %in% rota_last_dose_c_decision_report, decision := "use_report_extractions"]

decision_df[grepl("rota", me_name) & !grepl("ratio", me_name) & nid %in% rota_last_dose_3_decision_microdata, decision := "use_microdata"]
decision_df[grepl("rota", me_name) & !grepl("ratio", me_name) & nid %in% rota_last_dose_3_decision_report, decision := "use_report_extractions"]

decision_df[grepl("rota", me_name) & !grepl("ratio", me_name) & nid %in% rota_last_dose_2_decision_microdata, decision := "use_microdata"]
decision_df[grepl("rota", me_name) & !grepl("ratio", me_name) & nid %in% rota_last_dose_2_decision_report, decision := "use_report_extractions"]

decision_df[grepl("rota", me_name) & !grepl("ratio", me_name) & nid %in% rota_last_dose_1_decision_microdata, decision := "use_microdata"]
decision_df[grepl("rota", me_name) & !grepl("ratio", me_name) & nid %in% rota_last_dose_1_decision_report, decision := "use_report_extractions"]


# MCV/RCV: Nothing to be done (already dose-independent)

### make new country schedule-specific vacc_rotac rows
rotac_schedule <- readRDS(file.path("FILEPATH"))
setnames(rotac_schedule, "me_name", "rota_complete")
rotac_schedule[, me_name := paste0("vacc_rota", doses)]
rotac <- merge(decision_df, rotac_schedule, by=c("ihme_loc_id", "me_name"), all.x=T)[!is.na(rota_complete) & !is.na(decision)]
rotac[, `:=` (me_name=rota_complete, location_id=NULL, doses=NULL, rota_complete=NULL)]
decision_df <- rbind(decision_df, rotac)

### Make comparison of decisions against previous decision set	
previous_decisions  <- fread("FILEPATH")	
decision_comparison <- merge(decision_df, previous_decisions, by = c("nid", "ihme_loc_id", "me_name"), all.x = TRUE)	
setnames(decision_comparison, c("decision.x", "decision.y"), c("new", "old"))	
changed_decisions   <- decision_comparison[!is.na(old) & old != "NA" & !is.na(new) & new != old, ]	
# Save changed decisions to date stamped folder	
fwrite(changed_decisions, paste0("FILEPATH"))
# Save final decisions to date stamped folder
fwrite(decision_df, paste0("FILEPATH"))
# Save non-date stamped version as "current" version to reference repo
fwrite(decision_df, file.path("FILEPATH"))  

#--- 7. Diagnostic Plots ----------------------------------------------------------------------------------------------
# 7a. Summarize completeness of comparisons at the national level
summarize_me_match <- function(me) {
  df_subset <- subset(mics_nat_df, me_name == me)
  
  # Create a summary table
  df_summary <- data.table(me_name = me,
                           total_rows = nrow(df_subset),
                           both_complete = nrow(subset(df_subset, !is.na(tabulation) & !is.na(report))),
                           both_missing = nrow(subset(df_subset, is.na(tabulation) & is.na(report))),
                           missing_tabulation_only = nrow(subset(df_subset, is.na(tabulation) & !is.na(report))),
                           missing_report_only = nrow(subset(df_subset, !is.na(tabulation) & is.na(report))))
  return(df_summary)
}

me_match_summaries <- rbindlist(lapply(unique(mics_nat_df$me_name), summarize_me_match))

fwrite(me_match_summaries, paste0("FILEPATH"))

# 7b. Examine raw threshold effects by me

summarize_me_threshold <- function(me) {
  df_subset <- subset(mics_nat_df, me_name == me)
  
  # Create a summary table
  df_summary <- data.table(me_name = me,
                           threshold = threshold,
                           total_rows = nrow(df_subset),
                           rows_with_difference = nrow(subset(df_subset, !is.na(difference))),
                           rows_dropped_by_threshold = nrow(subset(df_subset, !is.na(difference) & abs(difference) > threshold)),
                           rows_with_difference_post_2000 = nrow(subset(df_subset, !is.na(difference) & year_start >= 2000)),
                           rows_dropped_by_threshold_post_2000 = nrow(subset(df_subset, !is.na(difference) & abs(difference) > threshold & year_start >= 2000)))
  df_summary[, pct_dropped_by_threshold := round(rows_dropped_by_threshold / rows_with_difference, 2)]
  df_summary[, pct_dropped_by_threshold_post_2000 := round(rows_dropped_by_threshold_post_2000 / rows_with_difference_post_2000, 2)]
  
  # Create a histogram
  return(df_summary)
}

me_threshold_summaries <- rbindlist(lapply(unique(mics_nat_df$me_name), summarize_me_threshold))

fwrite(me_threshold_summaries, paste0("FILEPATH"))

mics_nat_df[abs(difference) > threshold, above_threshold := TRUE]
mics_nat_df[abs(difference) <= threshold, above_threshold := FALSE]

# 7c. Plot report vs tabs by antigen and show outliers
gg_threshold <- ggplot(mics_nat_df[!is.na(difference)],
                       aes(x = report, y = tabulation, group = me_name)) +
  geom_point(alpha = 0.5, aes(color = above_threshold)) +
  geom_abline() +
  geom_smooth(method = "lm") +
  theme_bw() +
  xlim(0,1) +
  ylim(0,1) +
  facet_wrap(~me_name, nrow = 4) +
  coord_equal() +
  labs(title = "Report vs tabulation",
       subtitle = paste0("Threshold: ", threshold),
       x = "Report",
       y = "Tabulation",
       color = "Above threshold?") +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red"))

png(file = paste0("FILEPATH"),
    height = 10,
    width = 16,
    units = "in",
    res = 300)
print(gg_threshold)
dev.off()

# 7d. Plot histograms for each me_id with threshold
for (me in unique(mics_nat_df$me_name)) {
  gg_histogram <- ggplot(mics_nat_df[!is.na(difference) & me_name == me]) +
    geom_histogram(aes(x = difference), binwidth = 0.01) +
    geom_vline(xintercept = threshold, color = "red") +
    geom_vline(xintercept = -threshold, color = "red") +
    theme_bw() +
    xlim(-1,1) +
    labs(title = paste0("Histogram of differences: ", me),
         subtitle = paste0("Threshold: ", threshold),
         x = "Difference",
         y = "Count (Rows)")
  
  png(file = paste0("FILEPATH"),
      height = 7,
      width = 12,
      units = "in",
      res = 300)
  print(gg_histogram)
  dev.off()
}

#**********************************************************************************************************************


