###############################################################################
###############################################################################
## MBG Master Launch Script
##
## Author: <AUTHOR>
## Vaccine: DPT
## Date: <DATE>
##
## This is the main launch script.  Will run individual models for each needed
## vaccine or conditional vaccine data set for ordinal regression using a
## continuation-ratio model, then combine these arithmetically in order to
## the desired vaccine coverage metrics. 
##
###############################################################################
###############################################################################

###############################################################################
## SETUP
###############################################################################

## clear environment
rm(list=ls())

## Set repo location and indicator group
user               <- Sys.info()['user']
core_repo          <- '<FILEPATH>/lbd_core/'
indic_repo         <- '<FILEPATH>/vaccine/'


## sort some directory stuff
commondir      <- sprintf('<FILEPATH>/lbd_core/mbg_central/share_scripts/common_inputs')
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))



# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
library(raster)
load_mbg_functions(core_repo)
lapply(package_list, require, character.only = TRUE)

# Custom load indicator-specific functions
source(paste0(indic_repo,'functions/misc_vaccine_functions.R'))
source(paste0(indic_repo,'functions/parallelize_hack.R'))
'%!in%' <- function(x,y)!('%in%'(x,y))

## Script-specific stuff begins here ##########################################
indicator_group <- "vaccine"
coverage <- TRUE 
vacc <- 'dpt'

if(coverage == TRUE){
  vaccine <- vacc
  set_up_indicators(stem = vaccine, 
                    doses = 3, 
                    single_dose = F, 
                    save_doses = F,
                    save_2_cov = F)
}


slots              <- 4

## Create run date in correct format
run_date <- paste0(make_time_stamp(TRUE))


# Define a log directory 
log_dir <- paste0("<FILEPATH>/logs/", run_date, "/")
dir.create(log_dir, recursive = T)

if(vaccine %in% c("dpt")){
  config <- load_config(repo            = indic_repo,
                        indicator_group = "",
                        indicator       = NULL,
                        config_name     = paste0("config_", vaccine), 
                        covs_name       = paste0(vaccine, "_covs"),
                        run_test        = FALSE) ### TO DO: fix this!
}



## Ensure you have defined all necessary settings in your config
check_config()

# distribute config to all indicator directories
distribute_config(cfg = config, indicators = all_indicators)

# parse formatting
if (class(Regions) == "character" & length(Regions) == 1) Regions <- eval(parse(text=Regions))
if (class(summstats) == "character" & length(summstats) == 1) summstats <- eval(parse(text=summstats))
if (class(year_list) == "character") year_list <- eval(parse(text=year_list))

###############################################################################
## Launch parallel models
###############################################################################

launch_lv = list(indicator = model_indicators,
                 use_gn = TRUE)

launch_para_output <- parallelize_hack(lv = launch_lv,
                                  save_objs = c("run_date", "indicator_group", "vaccine", "core_repo", "shapefile_version"),
                                  script = "launch.R",
                                  prefix='launch',
                                  script_dir = indic_repo,
                                  output_log_dir = paste0('<FILEPATH>/output/%x.o%j'),
                                  error_log_dir = paste0('<FILEPATH>/errors/%x.e%j'),
                                  project              = '<PROJECT>',
                                  queue             = '<QUEUE>',
                                  memory = 10,
                                  cores = 1,
                                  time = '16-00:00:00',
                                  sing_image = 'default')


###############################################################################
## Rake indicators
###############################################################################
pushover_notify("Master script: starting raking",
                title = "Raking")

multiple_raking_indicators <-rake_indicators  #Dealing with multiple indicators w DPT


for(r in multiple_raking_indicators){
    rake_indicators <- r


  rake_lv <- list(region = Regions,
                  holdout = 0)


  rake_para_output <- parallelize_hack(lv = rake_lv,
                                         save_objs = c("run_date", "indicator_group", "vaccine", "doses", "shapefile_version","rake_indicators"),
                                         script = "rake_fractional.R",
                                         prefix='rake',
                                         script_dir = indic_repo,
                                         output_log_dir = paste0('<FILEPATH>/output/%x.o%j'),
                                         error_log_dir = paste0('<FILEPATH>/errors/%x.e%j'),
                                         project              = '<PROJECT>',
                                         queue             = '<QUEUE>',
                                         memory = 400,
                                         cores = 4,
                                         time = '2-00:00:00',
                                         sing_image = 'default',
                                         raking=TRUE)

}

####################################################################################
#Unraked aggregation
####################################################################################

#Run in parallel #######################

for(r in multiple_raking_indicators){
  indicator <- r
  expand_lv <- list(region = Regions)


  expand_para_output <- parallelize_hack(lv = expand_lv,
                                         save_objs = c("Regions", "indicator_group",
                                                       "coverage", "vacc", "doses", "indicator",
                                                       "run_date","summstats","shapefile_version","year_list",
                                                       "pop_measure","interval_mo","pop_release","countries_not_to_subnat_rake",'use_IND_states'),
                                         script = "unraked_aggregate_cell_draws.R",
                                         prefix='unraked_agg',
                                         script_dir = indic_repo,
                                         output_log_dir = paste0('<FILEPATH>/output/%x.o%j'),
                                         error_log_dir = paste0('<FILEPATH>/errors/%x.e%j'),
                                         project              = '<PROJECT>',
                                         queue             = '<QUEUE>',
                                         memory = 200,
                                         cores = 4,
                                         time = '3-00:00:00',
                                         sing_image = 'default')

}

####################################################################################
#processing of conditional DTP indicators (DTP1 & DTP3<12)
####################################################################################

## Run in parallel #######################
 expand_lv <- list(region = Regions)
if(vaccine == 'dpt'){
  run_dropout=F
  dir.create(paste0("<FILEPATH>/vaccine/dpt1_cov/output/", run_date,'/output/'),recursive=T)
  dir.create(paste0("<FILEPATH>/vaccine/dpt1_cov/output/", run_date,'/errors/'),recursive=T)
  expand_para_output <- parallelize_hack(lv = expand_lv,
                                         save_objs = c("Regions", "indicator_group",
                                                       "coverage", "vacc", "doses",
                                                       "run_date", "shapefile_version", "summstats","run_dropout"),
                                         script = "post_hoc_combine_dpt1.R",
                                         prefix='combine_dpt1',
                                         script_dir = indic_repo,
                                         output_log_dir = paste0('<FILEPATH>/output/%x.o%j'),
                                         error_log_dir = paste0('<FILEPATH>/errors/%x.e%j'),
                                         project              = '<PROJECT>',
                                         queue             = '<QUEUE>',
                                         memory = 100,
                                         cores = 4,
                                         time = '3-00:00:00',
                                         sing_image = 'default')

}


###############################################################################
## Post predict (summarize, combine, save files)
###############################################################################
## Convenience--make sure dpt1 is added to the mix
all_indicators <- c('dpt3_cov', 'dpt1_cov')
postest_indicators <- c('dpt3_cov', 'dpt1_cov')

for(indicator in all_indicators){

  sharedir <- sprintf('<FILEPATH>/%s/%s',indicator_group, indicator)
  strata <- Regions
  #######################################################################################
  ## Merge!!
  #######################################################################################

  if(indicator %!in% c('dpt1_cov')) rr <- c("raked", "unraked") else rr <- 'raked'
  if(indicator %!in% c('dpt1_cov')) rr_tf <- c(T,F) else rr_tf <- T
  if(indicator %!in% c('dpt1_cov')) rf_tab <- T else rf_tab <- F
  run_summary <- T

  post_load_combine_save(regions    = Regions,
                         summstats  = summstats,
                         raked      = rr,
                         rf_table   = FALSE,
                         run_summ   = FALSE,
                         indic      = indicator,
                         ig         = indicator_group,
                         sdir       = sharedir)

  # Clean up / delete unnecessary files
  clean_after_postest(indicator             = indicator,
                      indicator_group       = indicator_group,
                      run_date              = run_date,
                      strata                = strata,
                      delete_region_rasters = F)


  #######################################################################################
  ## Combine aggreagations
  #######################################################################################
  holdouts <- 0
  combine_aggregation(rd = run_date, indic = indicator, ig = indicator_group,
                      ages = 0,
                      regions = strata,
                      holdouts = holdouts,
                      raked = rr_tf,
                      delete_region_files = F)


  #######################################################################################
  ## Summarize
  #######################################################################################


  summarize_admins(summstats = c("mean", "lower", "upper", "cirange", "cfb"),
                   ad_levels = c(0,1,2),
                   raked = rr_tf)

}

########################################################################################
#compile input data
for(indicator in all_indicators[1]){
  #indicator <- all_indicators[1]

  outputdir <- paste0("<FILEPATH>/", indicator_group, "/",indicator,"/output/", run_date, "/")

  csvs <- list.files(outputdir, pattern = "input_data_(.*).csv", full.names = T)
  csv_master <- rbindlist(lapply(csvs, fread),fill=T)
  csv_master[, V1 := NULL]
  csv_master<-unique(csv_master)
  write.csv(csv_master, file=paste0(outputdir, '/input_data.csv'))
}


########################################################################################
message('Done with master script!')
## END OF FILE
###############################################################################