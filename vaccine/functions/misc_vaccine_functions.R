# HEADER ------------------------------------------------------------------
# Author: <AUTHOR>
# Date: <DATE>
# Project: MBG
# Purpose: various R functions for use in vaccine data processing

# FUNCTIONS ---------------------------------------------------------------

# Timer Functions ---------------------------------------------------------
#
## Functions to time within the child stacking regions
##
## General usage:
##    require(tictoc)
## 
##    tic("Step 1")
##    **your code here**
##    toc(log = T)
##
##    ticlog <- tic.log(format = F)
##    generate_time_log(ticlog)   
##
##  Returns: data table with two columns
##     "step": names of events (e.g. "Step 1")
##     "time": time elapsed (as text: Xh Xm Xs)
##
##  Note: can nest tic/toc pairs

generate_time_log <- function(ticlog) {
  
  # Set up packages
  require(magrittr)
  require(data.table)
  
  # Functions in functions
  strip_time <- function(x) {
    sec <- as.numeric(x$toc - x$tic)
    time <- format_time(sec)
    name <- x$msg
    
    df <- c(name, time) %>%
            t %>%
            as.data.table
    
    names(df) <- c("step", "time")
    
    return(df)
  }
  
  format_time <- function(run_time) {
    run_time <- round(as.numeric(run_time),0)
    
    hours <- run_time %/% 3600
    remainder <- run_time %% 3600
    minutes <- remainder %/% 60 
    seconds <- remainder %% 60
    
    run_time <- paste0(hours, "h ", minutes, "m ", seconds, "s")
    return(run_time)
  }
  
  df_out <- lapply(ticlog, strip_time) %>% rbindlist
  
  return(df_out)

}



# pushover_notify() -------------------------------------------------------

# Function to notify via pushover when job done on cluster

pushover_notify <- function(message, title = NULL, image = NULL) {
  
  ## Function to enable custom push notifications from within R on cluster
  ## Uses command line curl call 
  ## (the httr package would be a good alternative, but breaks on our cluster)
  
  ## Inputs: message (text you want to send) 
  #          title (optional)
  #          image (path to image to include)
  
  ## Outputs: command line curl call
  
  ## Setup: ---------------------------------------------------------------
  
  ## 1) Download the pushover app on your phone and sign up for an account
  ##    https://pushover.net
  
  ## 2) Find your user token and enter in the "keyfile" 
  ##    in your home directory (<FILEPATH>) 
  ##    (plain text, token only, no spaces)
  
  ## 3) Register an application (https://pushover.net/apps/build) with pushover
  ##    and find your app token, then enter that app token in the 'appfile'
  ##    in your home directory (<FILEPATH>)
  ##    (again, plain text, no spaces)
  
  ## 4) Enjoy custom push notifications from your R code with the syntax
  ##    `pushover_notify("my_message", "my_title")`.  You can even send images
  ##    if you point towards a saved .jpg or .png on the cluster!
  
  # -----------------------------------------------------------------------
  

    h_root <- "FILEPATH"

  
  url <- "<ADDRESS>"

  keyfile <- paste0(h_root, "/lib/pushover_key.txt")
  
  # Check if file exists to continue
  if (!file.exists(keyfile)) {
    return("Pushover not enabled")
  }
  
  key <- readChar(keyfile, file.info(keyfile)$size)
  
  appfile <- paste0(h_root, "/lib/pushover_app.txt")
  app <- readChar(appfile, file.info(appfile)$size)
  

  
  if (is.null(image)) {
    img_addin <- "" 
  } else if (!is.null(image)) {
    img_addin <- paste0("-F \"attachment=@", image, "\" \\\n")
  }

  curl_cmd <- paste0("curl -s ",
                      "--form-string \"token=", app, "\" \\\n", 
                      "--form-string \"user=", key, "\" \\\n",
                      "--form-string \"message=", message, " \" \\\n",
                     "--form-string \"title=", title, " \" \\\n",
                     img_addin,
                      "<ADDRESS>")
  
  system(curl_cmd, ignore.stdout = T)

}




# Create raster of GBD estimates

make_gbd_rasterbrick <- function(gbd_data, 
                                 indicator, 
                                 gaul_list,
                                 years = c(2000, 2005, 2010, 2015)) {

  message(paste0("Creating GBD RasterBrick for ", indicator))

  simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list,
                                        buffer = 0.4, shapefile_version=modeling_shapefile_version)

  subset_shape   <- simple_polygon_list[[1]]
  raster_list    <- build_simple_raster_pop(subset_shape)
  simple_raster  <- raster_list[['simple_raster']]

  # create template raster
  ext <- extent(simple_raster)
  ncol <- ncol(simple_raster)
  nrow <- nrow(simple_raster)
  r <- raster(ext, ncol, nrow)

  gbd_brick_list <- lapply(years, function(yr) {

    my_spdf <- merge(subset_shape, gbd_data[year == yr], by.x = "GAUL_CODE", by.y = "name", all.x = T, all.y = F)
    ras <- rasterize(my_spdf, r, "mean")

  })

  gbd_brick <- brick(gbd_brick_list)

  return(gbd_brick) 
  
}






clean_model_results <- function(rd   = run_date,
                                regs = Regions,
                                ages = 0,
                                nm   = '',
                                indicator = indicator, 
                                indicator_group = indicator_group){

  require(magrittr)

  sharedir <- paste0('<FILEPATH>/', indicator_group, '/', indicator, '/output/', rd, '/')

  # make loopvars
  lv <- expand.grid(regs,ages)

  # grab model fit objects
  mods <- model_fit_table(lv=lv,rd=rd,nullmodel=nm, indicator = indicator, indicator_group = indicator_group)

  # one table for each region
  tables=list()
  for(r in regs){

    # get parameter names
    tempbig <- lapply(ages, function(age) row.names(mods[[paste0(r, '_', age)]])) %>%
                  unlist %>%
                  unique %>%
                  data.table(parameter = .)

     tempbig = tempbig[!parameter%in%c('Theta1 for space','Theta2 for space','GroupRho for space'),]
     assign(r,tempbig)
      for(a in ages){
          tmp=data.table(cbind(parameter=row.names(mods[[paste0(r,'_',a)]]),
                                         mods[[paste0(r,'_',a)]]))[,c(1,2,4,6)]
          setnames(tmp, c("parameter",paste0("a",a,"_mean"),paste0("a",a,"_lb"),paste0("a",a,"_ub") ))

          # tranform spatial covariates to be more readbale
          tmp$parameter=as.character(tmp$parameter)
          tmp<-tmp[!parameter%in%c('Theta1 for space','Theta2 for space'),]
          tmp$parameter[tmp$parameter=="GroupRho for space"]="GPRandom Rho for time"
          tmp <- rbind(tmp,transform_spatial_effects(params.mat=mods[[paste0(r,'_',a)]],age=a))

          # round
          cols <- names(tmp)[2:4]
          tmp[,(cols) := round(.SD,4), .SDcols=cols]

          assign(r,merge(get(r),tmp,by='parameter',all=TRUE)) #,all.x=TRUE))

        #  get(r)$parameter[get(r)$parameter=="GroupRho for space"]
      }
      tables[[r]]<-get(r)

  }

  for(r in regs){
    write.csv(tables[[r]],sprintf('%s/model_results_table_%s.csv',sharedir,r))
  }
  return(tables)

}




# Load alternative raking targets

load_newest_gbd_vax  <- function(vaccine,
                                 gaul_list = NULL,
                                 return_cis = F,
                                 return_mode = "national", # national, subnational, or both
                                 gbd_year = 2020,
                                 years = c(2000:2021),
                                 gbd_date = "best", 
                                 return_field = "location_id",
                                 shapefile_version = modeling_shapefile_version,
                                 age_specific = FALSE,
                                 covid_version='') {

  source('<FILEPATH>/get_location_metadata.R')

  str_match <- stringr::str_match

  gaul_to_loc_id <- get_location_metadata(location_set_id = 22, gbd_round_id=6)
  names(gaul_to_loc_id)[names(gaul_to_loc_id)=="loc_id"] <- "location_id"
  
  if(gbd_date=='best') gbd_date <-'2023-02-12'
  
  if(age_specific==TRUE) age_folder = 'age_specific/' else age_folder = ''
  
  if(covid_version !='') covid_version <- paste0('_',covid_version)

  gbd_file <- paste0('<FILEPATH>/gbd', gbd_year, '/',age_folder, gbd_date,  '/vacc_',vaccine,covid_version,'.rds')

  gbd <- readRDS(gbd_file)
  gbd[,ihme_loc_id:=NULL]  #In case this is already included, the merge will lead to incorrect names
  gbd <- merge(gaul_to_loc_id, gbd, by="location_id",  all.y = T)  

  if (return_mode == "national") {
    gbd <- subset(gbd, !grepl('_', ihme_loc_id))
  }

  if (return_field == "GAUL_CODE") {
    loc_map <- get_location_code_mapping(shapefile_version=shapefile_version)
    gbd <- merge(gbd, subset(loc_map, select = c("loc_id", "GAUL_CODE")), by.x = "location_id", by.y = "loc_id", all.x = T)
  }

  # subset to year & gaul lists
  gbd <- subset(gbd, year_id %in% years)
  if (!is.null(gaul_list)) {
    gbd <- subset(gbd, get(return_field) %in% gaul_list)
  }
  
  setnames(gbd, 
           c(return_field, "year_id", "gpr_mean", "gpr_lower", "gpr_upper"),
           c("name", "year", "value", "lower", "upper"))
  if (return_cis == F)  return_vars <- c("name", "year", "value")
  if (return_cis == T)  return_vars <- c("name", "year", "value", "lower", "upper")
  if (age_specific ==T) return_vars <- c(return_vars, 'age_group_id')

  if (return_mode != "national") {
    gbd[, parent_ihme_loc_id := str_match(ihme_loc_id,"(.*)_")[,2]]

    if (return_mode == "subnational") {
      # Drop national estimates if we have subnational estimates
      ihme_loc_ids_to_drop <- c( NA , "CHN", "GBR", "ZAF", "USA", "MEX", "IDN", "BRA", "IND", "NOR", "SWE", "JPN",
                                 "KEN" ,"NZL", "ETH", "IRN", "RUS", "PAK", "NGA","PHL")
      ihme_loc_ids_to_drop <- ihme_loc_ids_to_drop[!is.na(ihme_loc_ids_to_drop)]
      gbd <- subset(gbd, !(ihme_loc_id %in% ihme_loc_ids_to_drop))  
    }
  }

  gbd <- subset(gbd, select = return_vars)

  return(gbd)
}


## check_input_data ################################################

#' Checks to see if input data has been saved into sharedir
#' If not, saves an input data object for the given region
#' This is useful if post-estimating indicators that are calculated
#' from other modeled indicators (e.g. if using continuation-ratio
#" approach for ordinal indicators)
#'
#' @param indicator Indicator
#' @param indicator_group Indicator group
#' @param age Age group
#' @param reg Region
#' @param run_date Run date
#' @param use_share Should the file look in share for master input_data csv?
#'
#' @return Does not return anything; saves input_data files by region
#'         into standard directory for this run
#'
#' @examples
#' 

check_input_data <- function(indicator,
                             indicator_group,
                             age,
                             reg,
                             holdout = 0,
                             run_date,
                             use_share = F,
                             ylist = year_list,
                             shapefile_version=modeling_shapefile_version,
                             sunmational_raking=subnational_raking) {

  message(paste0("Checking input data for indicator: ", indicator,
                 " | region: ", reg, 
                 " | holdout: ", holdout))
  # make a pathaddin
  pathaddin <- paste0('_bin',age,'_',reg,'_',holdout)
  input_file <- paste0("<FILEPATH>/", indicator_group, "/", 
                       indicator, "/output/", run_date, "/input_data", pathaddin, ".csv")
  
  if (!file.exists(input_file)) {
    ## Load simple polygon template to model over
    gaul_list           <- get_adm0_codes(reg, shapefile_version=modeling_shapefile_version, subnational_raking=subnational_raking)
    simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list, buffer = 1, tolerance = 0.4, use_premade = T, shapefile_version=modeling_shapefile_version)
    subset_shape        <- simple_polygon_list[[1]]
    simple_polygon      <- simple_polygon_list[[2]]

    message('Loading in full dataset using load_input_data() and saving...')

    df <- load_input_data(indicator   = gsub(paste0('_age',age),'',indicator),
                          simple      = simple_polygon,
                          agebin      = age,
                          removeyemen = TRUE,
                          pathaddin   = pathaddin,
                          years       = yearload,
                          withtag     = as.logical(withtag),
                          datatag     = datatag,
                          use_share   = as.logical(use_share),
                          yl          = ylist)
  }
}

## distribute_config ################################################

#' Distribute config file to directories for a set of indicators with the same run date
#'
#' After running this, use load_config() with post_est_only == T to pull this config
#' Allows the config to be loaded once only from the master script, e.g. in case
#' you want to change it and launch another model shortly afterwards
#'
#' @param cfg config object (from `load_config()`)
#' @param indicators vector of indicator names to distribute config to
#' @param indicator_group indicator group
#' @param run_date shared run date for all indicators
#'
#' @return does not return any objects; writes CSV with config to each of the 
#'         indicator directories on `<FILEPATH>`
#' @examples
#' 
# Distribute config file to directories for a set of indicators with the same run date.
# After running this, use `load_config()` with `post_est_only == T` to pull this config.
# Allows the config to be loaded once only from the master script, e.g. in case
#   you want to change it and launch another model shortly afterwards

distribute_config <- function(cfg = config,
                              indicators,
                              ig = indicator_group,
                              rd = run_date) {
  for (ind in indicators) {
    ind_dir <- paste0("<FILEPATH>/",
                      ig, "/", ind,
                      "/output/", rd, "/")
    dir.create(ind_dir, recursive = T, showWarnings = F)
    write.csv(cfg, paste0(ind_dir,'/config.csv'), row.names = FALSE)
  }
}



# Function to set up indicators for multiple-dose runs

set_up_indicators <- function(stem,
                              doses,
                              single_dose = T,
                              ir = indic_repo,
                              save_doses = F,
                              save_2_cov = F) {
  # Define objects
  vaccine            <- stem

  # Indicators to model in INLA
  model_indicators <- paste0(stem, doses, "_cov")
  if (single_dose == F & vaccine != "dpt") model_indicators <- c(paste0(stem, 1:(doses-1), "_cond"), model_indicators)
  if (single_dose == F & vaccine == "dpt") model_indicators <- c(paste0(stem, "12_cond"), model_indicators)

  # Dose specific indicators (i.e. p(exactly that # of doses))
  dose_indicators <- paste0(stem, doses, "_cov")
  if (single_dose == F & save_doses == T) dose_indicators <- c(paste0(stem, 0:(doses-1), "_dose"), dose_indicators)                          

  # Indicators to rake
  # This will need to be changed once we have full dose-specific information from GBD
  rake_indicators    <- c(unique(paste0(stem, doses, "_cov")))
  if (single_dose == F & vaccine == "dpt") rake_indicators <- c(rake_indicators, 
                                             paste0(stem, "12_cond") )

  postest_indicators <- rake_indicators  
  
  all_indicators     <- unique(c(model_indicators, rake_indicators, dose_indicators, postest_indicators))

  # Message outputs
  message("\nGenerating Indicator Lists:")
  message(paste0(" model_indicators:   ", paste(sort(model_indicators), collapse = ", ")))
  message(paste0(" dose_indicators:    ", paste(sort(dose_indicators), collapse = ", ")))
  message(paste0(" rake_indicators:    ", paste(sort(rake_indicators), collapse = ", ")))
  message(paste0(" postest_indicators: ", paste(sort(postest_indicators), collapse = ", ")))
  message(paste0(" all_indicators:     ", paste(sort(all_indicators), collapse = ", ")))

  # Check to ensure all files present
  message("\nChecking to ensure input data files present")

  config <- fread(paste0(ir, "config_", vaccine, ".csv"), header = F)
  load_dir  <- '<FILEPATH>'

  # For padding in output
  char_space <- max(nchar(all_indicators)) + 2

  # Loop over indicators for check
  for (ind in sort(unique(c(model_indicators)))) {
    check_exists <- file.exists(paste0(load_dir, ind, ".csv"))
    message(paste0("  ", ind, ":", 
                   paste0(rep(" ", char_space - nchar(ind)), collapse = ""), 
                   ifelse(check_exists, "SUCCESS", "FAILURE")))
    if (!check_exists) stop(paste0("CSV file for ", ind, " does not exist!"))
  }

  # Assign everything to global environment
  message("\nAssigning to global environment...")
  for (param in c("vaccine", "doses", "save_2_cov", "save_doses", 
                  "model_indicators", "dose_indicators", 
                  "rake_indicators", "postest_indicators", 
                  "all_indicators")) {
     assign(param, get(param), envir=globalenv())
  }

  message("  Completed successfully!")
}

# Function to load profiling information for a given vaccine, region, and draw combo

load_profiling <- function(vax,
                           draws,
                           reg,
                           ir = indic_repo,
                           default_mem = 300,
                           default_slots = 4) {

  profile_csv <- paste0(indic_repo, "profiling_", vax, ".csv")
  if (file.exists(profile_csv)) {
    profile <- fread(profile_csv)
    mem <- profile[region == reg & samples == draws, mem]
    slots <- profile[region == reg & samples == draws, slots]
  } else {
    mem <- NULL
    slots <- NULL
  }

  if (length(mem) == 0) {
    mem <- default_mem
    message(paste0("No profiling csv entry found for memory -- using default value of ", default_mem))
  }

  if (length(slots) == 0) {
    slots <- default_slots
    message(paste0("No profiling csv entry found for slots -- using default value of ", default_slots))
  }

  return(c(mem = mem, slots = slots))

}
