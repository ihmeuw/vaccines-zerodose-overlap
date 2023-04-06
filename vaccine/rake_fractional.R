# Parallel script for post-estimation and aggregation

## SETUP ################################################################################

## clear environment
rm(list=ls())

## Set repo location and indicator group
user               <- Sys.info()['user']
core_repo          <- '<FILEPATH>/lbd_core/'
indic_repo         <- '<FILEPATH>/vaccine/'


## sort some directory stuff
sharedir <- sprintf("<FILEPATH>/%s/%s", indicator_group, indicator)
commondir      <- sprintf('<FILEPATH>/lbd_core/mbg_central/share_scripts/common_inputs')
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))# source in functions
repo <- paste0(core_repo, '/mbg_central')
# load the functions from the repo
load_mbg_functions(repo)


lapply(package_list, require, character.only = TRUE)

# Custom load indicator-specific functions
source(paste0(indic_repo,'functions/misc_vaccine_functions.R'))

## Script-specific stuff begins here ##########################################
indicator_group <- 'vaccine'

# Load objects from qsub
load_from_parallelize() # run_date, indicator, vaccine, indicator_group, reg

interval_mo <- 12
metric_space <- "rates"
indicator <- rake_indicators


# Define directories
main_dir <- paste0("<FILEPATH>/", indicator_group, "/", indicator, "/output/", run_date, "/")
temp_dir <- paste0(main_dir, "temp_post_est/")


## Read config file and save all parameters in memory
config <- load_config(repo            = indic_repo,
                      indicator_group = indicator_group,
                      indicator       = indicator,
                      post_est_only   = TRUE,   
                      run_tests       = FALSE,        
                      run_date        = run_date)


# Get the necessary variables out from the config object into global env
rake_countries <- eval(parse(text = config[V1 == "rake_countries", V2]))
rake_subnational <- eval(parse(text = config[V1 == "subnational_raking", V2]))
modeling_shapefile_version <- config[V1 == "modeling_shapefile_version", V2]
raking_shapefile_version <- config[V1 == "raking_shapefile_version", V2]
countries_not_to_rake <- config[V1 == "countries_not_to_rake", V2]
countries_not_to_subnat_rake <- config[V1 == "countries_not_to_subnat_rake", V2]
year_list <- eval(parse(text = config[V1 == "year_list", V2]))
countries_not_to_rake <- c("ESF+GUF+VIR+CHN+BRA") 
countries_not_to_subnat_rake <- c("ESF+GUF+VIR+CHN+BRA") 

reg <- region


# Print some settings to console
message(indicator)
message(indicator_group)
message(run_date)
message(reg)
message(pop_measure)
message(gbd_date)
message(paste0("Summary stats: ", paste0(summstats, collapse = ", ")))

# Print raking info
print(paste0("Metric Space                       : ", metric_space))
print(paste0("Subnational raking                 : ", rake_subnational))
print(paste0("Countries not to rake at all       : ", countries_not_to_rake))
print(paste0("Countries not to rake subnationally: ", countries_not_to_subnat_rake))



# For now just assign stratum to reg (will need to modify the below for strata beyond reg)
age <- 0
holdout <- 0



# PULL gbd values
message("\n################################################################################")
message(paste0("Raking ", vaccine, doses, "_cov | Region: ", reg))
message("################################################################################")

if(indicator =='dpt12_cond') doses <-'12_cond' #Change to set up bcg to match gbd indicator name

# Load GBD estimates
gbd_max_cov <- load_newest_gbd_vax(vaccine = paste0(vaccine, doses),
                               gaul_list = NULL,
                               years = year_list,
                               gbd_date = gbd_date,
                               gbd_year = 2020,
                               return_mode = "subnational",
                               return_field = "location_id",
                               covid_version = '')

rake_to <- subset(gbd_max_cov, select = c("name", "year", "value"))

#Grab the national values also
gbd_max_cov_nat <- load_newest_gbd_vax(vaccine = paste0(vaccine, doses),
                                   gaul_list = NULL,
                                   years = year_list,
                                   gbd_date = gbd_date,
                                   gbd_year = 2020,
                                   return_mode = "national",
                                   return_field = "location_id",
                                   covid_version = '')
gbd_max_cov_nat <- subset(gbd_max_cov_nat, select = c("name", "year", "value"))
gbd_max_cov_nat <- gbd_max_cov_nat[!(name %in% rake_to$name)]

rake_to <- rbind(rake_to,gbd_max_cov_nat)


# For now, treat everything < 0.001 (that is, 0.1%) as zero
# This helps logit raking converge
rake_to[value < 0.001, value := 0] 

# Conform to rake target defaults
setnames(rake_to, "value", "mean")


gbd <- rake_to


calculate_fractional_rfs <- function(ndraws = ndraws,
                                     cell_pred = cell_pred,
                                     gbd = gbd,
                                     sharedir = sharedir,
                                     run_date = run_date,
                                     indicator = indicator,
                                     shapefile_version = shapefile_version,
                                     stratum = stratum,
                                     countries_not_to_rake = NULL,
                                     custom_output_folder = NULL,
                                     countries_not_to_subnat_rake = NULL,
                                     MaxJump = 10,
                                     MaxIter = 80, 
                                     FunTol = 1e-5,
                                     approx_0_1 = T,
                                     zero_heuristic = T,
                                     iterate = T,
                                     rake_method = "linear") {
  message("converting from prevalence to counts")
  # set the variables to aggregate
  overs <- paste0("V", 1:ndraws)
  
  # convert to counts
  cell_pred <- cell_pred[, (overs) := lapply(overs, function(x) get(x) * pop_raked)]
  
  # do calculations!
  rake_geo <- cell_pred[, lapply(c(overs, "pop_raked"), function(x) sum(get(x), na.rm = T)), by = c("year", "location_id")]
  setnames(rake_geo, grep("V[0-9]", names(rake_geo), value = T), c(overs, "pop_raked"))
  
  # convert back to rates/prevalence
  rake_geo <- rake_geo[, (overs) := lapply(overs, function(x) get(x) / pop_raked) ]
  
  # merge to admin estimates
  rake_geo <- merge(rake_geo, gbd, by.x = c("location_id", "year"), by.y = c("name", "year"), all.x = TRUE)
  
  # finding the mean of the draws at the raking geographies
  rake_geo$gbd_prev <- rake_geo$mean
  rake_geo[, mbg_rate:= rowMeans(.SD), .SDcols = paste0('V', c(1:ndraws))]
  rake_geo$rf <- rake_geo$gbd_prev / rake_geo$mbg_rate
  
  rake_geo <- subset(rake_geo, location_id != -1)   


  # Clear Out
  message("creating fractional raking factors table")
  fractional_rf <- rake_geo
  fractional_rf$mbg_prev <- fractional_rf$mbg_rate
  fractional_rf <- fractional_rf[, c("location_id", "year", "mbg_prev", "gbd_prev", "rf")]
  
  
  if(rake_method == "logit"){
    cell_pred <- cell_pred[, (overs) := lapply(overs, function(x) get(x) / pop_raked)]
    #for each country year, find the logit raking factor
    #redefine cys
    cys <- unique(rake_geo[,.(location_id, year)])
    rf <- lapply(seq(nrow(cys)), function(x) {
      
      #year and location
      theloc = cys[x,location_id]
      theyear = cys[x,year]
      
      message("now on ", theloc, " in ", theyear)

      if (nrow(rake_geo[location_id == theloc & year == theyear]) == 0) {
        return(0)
      } else if (rake_geo[location_id == theloc & year == theyear, .(mean)] == 0 & zero_heuristic == T) {
        # catch true zeroes (i.e. pre-introduction of an intervention) and return -9999. This will speed things up & will replace later with 0
        return(-9999)        
      } else {
        ret <- try(LogitFindK(gbdval     = rake_geo[location_id == theloc & year == theyear,mean],
                              pixelval   = as.matrix(cell_pred[location_id == theloc & year == theyear & !is.na(V1), paste0('V', 1:ndraws), with = F]), #pass the cell pred rows that corrospond to this country year
                              weightval  = as.numeric(cell_pred[location_id == theloc & year == theyear & !is.na(V1), pop_raked]),
                              MaxJump    = MaxJump,
                              MaxIter    = MaxIter,
                              FunTol     = FunTol,
                              approx_0_1 = approx_0_1))
        
        
        if(iterate) {
          # Iterate over larger values of MaxJump and NumIter if needed
          if (is.null(ret) | "try-error" %in% class(ret)) {
            message(paste0("Your GBD and MBG estimates are quite far apart for location ", theloc, " | year ", theyear))
            message("Increasing MaxJump and NumIter by 10-fold, but you should investigate this...")
            
            ret <- try(LogitFindK(gbdval     = rake_geo[location_id == theloc & year == theyear,mean],
                                  pixelval   = as.matrix(cell_pred[location_id == theloc & year == theyear & !is.na(V1), paste0('V', 1:ndraws), with = F]), #pass the cell pred rows that corrospond to this country year
                                  weightval  = as.numeric(cell_pred[location_id == theloc & year == theyear & !is.na(V1), pop_raked]),
                                  MaxJump    = MaxJump*10,
                                  MaxIter    = MaxIter*10,
                                  FunTol     = FunTol,
                                  approx_0_1 = approx_0_1))
          }
          
          # If we get this far, the estimates are generally *really* far apart
          if (is.null(ret) | "try-error" %in% class(ret)) {
            message(paste0("Your GBD and MBG estimates are REALLY far apart for location ", theloc, " | year ", theyear))
            message("Increasing MaxJump and NumIter by 100-fold, but you should investigate this...")
            
            ret <- try(LogitFindK(gbdval     = rake_geo[location_id == theloc & year == theyear,mean],
                                  pixelval   = as.matrix(cell_pred[location_id == theloc & year == theyear & !is.na(V1), paste0('V', 1:ndraws), with = F]), #pass the cell pred rows that corrospond to this country year
                                  weightval  = as.numeric(cell_pred[location_id == theloc & year == theyear & !is.na(V1), pop_raked]),
                                  MaxJump    = MaxJump*100,
                                  MaxIter    = MaxIter*100,
                                  FunTol     = FunTol,
                                  approx_0_1 = approx_0_1))
          }
          
          # Throw error if previous two didn't work
          if (is.null(ret) | "try-error" %in% class(ret)) {
            message(paste0("Your GBD and MBG estimates are WAY TOO far apart for location ", theloc, " | year ", theyear, " - stopping."))   #!# Core code has a 'stop' here
          }
        }
        
        return(ret)
      }
    })
    rf <- unlist(rf)
    fractional_rf$rf <- rf
  }
  
  # Don't rake if countries_not_to_rake is provided
  if (!is.null(countries_not_to_rake)) {

    # Get the GBD location IDs from the ISO codes to set raking
    # factors to 1
    nonrake_table <- get_gbd_locs(
      reg = countries_not_to_rake,
      rake_subnational = FALSE, shapefile_version = raking_shapefile_version
    )
    if (nrow(nonrake_table) > 0) {
      fractional_rf[location_id %in% nonrake_table$location_id, rf := 1]
    }
    rm(nonrake_table)
  }
  
  # Don't subnational rake if countries_not_to_subnat_rake is provided
  if (!is.null(countries_not_to_subnat_rake)) {
    
    # Get the GBD location IDs from the ISO codes to set raking
    # factors to 1. 
    # The main difference here from the above case is that we are pulling
    # in the subnational location IDs if they exist, and also filtering for the
    # case where rak_level (raking level) is >=1 (anything more detailed than ADM0)
    nonrake_table <- get_gbd_locs(
      reg = countries_not_to_subnat_rake,
      rake_subnational = TRUE,
      shapefile_version = raking_shapefile_version
    )[rak_level >= 1]
    if (nrow(nonrake_table) > 0) {
      fractional_rf[location_id %in% get_gbd_locs(
        reg = countries_not_to_subnat_rake,
        rake_subnational = TRUE, shapefile_version = raking_shapefile_version
      )$location_id, rf := 1]
    }
    rm(nonrake_table)
  }
  
  # Fill NAs with 1 BUT WITH A SOLEMN WARNING
  if(nrow(fractional_rf[is.na(rf)]) >= 1 ) {
    warning("You have NAs in your output raking factor. Please check your GBD outputs to make sure you didn't miss any locations.")
    warning("Forcing those missing raking factors to 1")
    fractional_rf[is.na(rf), rf:= 1]
  }
  
  # saving the raking factors
  if (!is.null(custom_output_folder)) {
    write.csv(fractional_rf, file = paste0(custom_output_folder, "/", indicator, "_", stratum, "_rf.csv"))
  } else {
    write.csv(fractional_rf, file = paste0(sharedir, "/output/", run_date, "/", indicator, "_", stratum, "_rf.csv"))
  }
  
  return(fractional_rf)
}



## PREPARE RASTERS, ETC. ################################################################

# Load cell draws
message("Loading Data...")
load(paste0(main_dir, indicator, "_cell_draws_eb_bin0_", reg, "_0.RData"))

# Check if load was successful; stop if not
if (!exists("cell_pred")) {
  message(filename_rds)
  stop("Unable to load cell_pred object! Check to make sure that the relevant object exists.")
}

# Rake estimates
if (rake_countries) {
  if (!exists("gbd")) {
    stop("rake_countries was specified as T in config, gbd raking targets must be provided.")
  }

  ## determine if a crosswalk is needed
  if (modeling_shapefile_version == raking_shapefile_version) crosswalk <- F else crosswalk <- T

  # Assume linear raking unless specified as logit
  if (rake_transform == "logit") rake_method <- "logit" else rake_method <- "linear"

'%!in%' <- function(x,y)!('%in%'(x,y))

  
  ## Get the simple and new_simple rasters prepped up for us
  print("Getting simple and prepped rasters")
  raster_outputs <- prep_shapes_for_raking(
    reg = reg,
    modeling_shapefile_version = modeling_shapefile_version,
    raking_shapefile_version = raking_shapefile_version,
    field = "loc_id"
  )

  ## Take out the objects from the list that actually matters to us:
  simple_raster <- raster_outputs[["simple_raster"]]
  new_simple_raster <- raster_outputs[["new_simple_raster"]]

  simple_polygon <- raster_outputs[["simple_polygon"]]
  new_simple_polygon <- raster_outputs[["new_simple_polygon"]]

  pixel_id <- raster_outputs[["pixel_id"]]

  ##### Using fractional raking #####
    print("Get GBD populations")
    gbd_pops <- prep_gbd_pops_for_fraxrake(pop_measure = pop_measure, reg = reg_gbd, year_list = year_list, gbd_round_id=7, decomp_step='step3')

    print("Using the rates raking and aggregation functions:")


####### now lets fractionally rake!
    ## First, create all the fractional rake factors
    fractional_rake_rates(
      cell_pred = cell_pred,
      simple_raster = simple_raster,
      simple_polygon = simple_polygon,
      pixel_id = pixel_id,
      shapefile_version = raking_shapefile_version,
      reg = reg,
      pop_measure = pop_measure,
      year_list = year_list,
      interval_mo = interval_mo,
      rake_subnational = rake_subnational,
      age_group = age_group,
      sex_id = sex_id,
      sharedir = sharedir,
      run_date = run_date,
      indicator = indicator,
      gbd = gbd,
      rake_method = "logit",
      gbd_pops = gbd_pops,
      countries_not_to_rake = countries_not_to_rake,
      countries_not_to_subnat_rake = countries_not_to_subnat_rake
    )

    ## Now, create the raked cell pred files!
    outputs <- fractional_agg_rates(
      cell_pred = cell_pred,
      simple_raster = simple_raster,
      simple_polygon = simple_polygon,
      pixel_id = pixel_id,
      shapefile_version = raking_shapefile_version,
      reg = reg,
      pop_measure = pop_measure,
      year_list = year_list,
      interval_mo = interval_mo,
      rake_subnational = rake_subnational,
      sharedir = sharedir,
      run_date = run_date,
      indicator = indicator,
      main_dir = main_dir,
      rake_method = "logit",
      age = age,
      holdout = holdout,
      countries_not_to_subnat_rake = countries_not_to_subnat_rake,
      return_objects = TRUE
    )

    ## Get the necessary outputs and rename the columns
    rf <- data.table(outputs[["rf"]])[, .(loc = location_id, year, start_point = mbg_prev, target = gbd_prev, raking_factor = rf)]
    raked_cell_pred <- outputs[["raked_cell_pred"]]


    ## Raked simple raster has been made above
    raked_simple_raster <- simple_raster

  
  gaul_list           <- get_adm0_codes(reg, shapefile_version = modeling_shapefile_version)
  simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list, buffer = 0.4, 
                                               shapefile_version = modeling_shapefile_version)
}


## SAVE THE RESULTS #####################################################################
message("Saving results...")

## save RF
save_post_est(rf, "csv", paste0(reg, "_rf"))

## save raked cell preds
save(raked_cell_pred, file = paste0(
  sharedir, "/output/", run_date, "/",
  indicator, "_raked_cell_draws_eb_bin0_", reg, "_0.RData"
))

# make and save summaries


# ## SAVE THE RESULTS #####################################################################
# make and save summaries

save_cell_pred_summary <- function(summstat, raked, ...) {
  message(paste0("Making summmary raster for: ", summstat, " (", raked, ")"))
  
  if (raked == "unraked") {
    cpred <- "cell_pred"
    mask_raster <- "simple_raster"
  }
  if (raked == "raked") {
    cpred <- "raked_cell_pred"
    mask_raster <- "raked_simple_raster"
  }
  if (raked == "raked_c") {
    cpred <- "raked_cell_pred_c"
    load(paste0(sharedir, "/output/", run_date, "/", indicator, "_raked_c_cell_draws_eb_bin0_", reg, "_0.RData" ))
    mask_raster <- "raked_simple_raster"
  }
  ras <- make_cell_pred_summary(
    draw_level_cell_pred = get(cpred),
    mask = get(mask_raster),
    return_as_raster = TRUE,
    summary_stat = summstat,
    ...
  )
  save_post_est(ras,'raster',paste0(reg, ifelse(raked == "raked", "_raked", ifelse(raked == 'raked_c', '_raked_c', '')), '_', summstat, '_raster'))
}

# Do this as lapply to not fill up memory in global env with big obs
if (is.null(gbd)) {
    rake_list <- c("unraked")
} else {
    rake_list <- c("unraked", "raked")
}


summ_list <- expand.grid(summstats, rake_list)


summstats <- c('mean', 'cirange', 'upper', 'lower', 'cfb')

for (r in rake_list) {
  if ("p_below" %in% summstats) {
    save_cell_pred_summary(
      summstat = "p_below",
      raked = r,
      value = 0.8,
      equal_to = F
    )
  }
  if ("mean" %in% summstats) {
    save_cell_pred_summary(
      summstat = "mean",
      raked = r
    )
  }
  if ("cirange" %in% summstats) {
    save_cell_pred_summary(
      summstat = "cirange",
      raked = r
    )
  }
  if ("upper" %in% summstats) {
    save_cell_pred_summary(
      summstat = "upper",
      raked = r
    )
  }
  if ("lower" %in% summstats) {
    save_cell_pred_summary(
      summstat = "lower",
      raked = r
    )
  }
    if ("cfb" %in% summstats) {
    save_cell_pred_summary(
      summstat = "cfb",
      raked = r
    )
  }
}


# All done
message(paste0("Done with post-estimation and aggregation for ", reg))