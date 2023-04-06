###############################################################################
###############################################################################
## MBG Launch Script
##
## Author: AUTHOR
## Indicator: can be any within the DPT vaccine group
## Date: DATE
##
## Source:
##   source("<FILEPATH>")
###############################################################################
###############################################################################

###############################################################################
## SETUP
###############################################################################

## clear environment
rm(list=ls())

## Set repo location and indicator group
user               <- Sys.info()['user']
core_repo          <- '<FILEPATH>/lbd_core'
indic_repo         <- '<FILEPATH>/vaccine'

## sort some directory stuff
setwd(core_repo)
commondir      <- '<FILEPATH>/lbd_core/mbg_central/share_scripts/common_inputs/'
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
load_mbg_functions(core_repo)

# Custom load indicator-specific functions
source(paste0(indic_repo,'functions/misc_vaccine_functions.R'))

library(stringr)

## Script-specific stuff begins here ##########################################

# Load from qsub
load_from_parallelize()

# Definitions
indicator_group <- 'vaccine'
sharedir <- sprintf('<FILEPATH>/%s/%s',indicator_group, indicator)

## Read config file (from sharedir) and save all parameters in memory
config <- load_config(repo            = indic_repo,
                      indicator_group = indicator_group,
                      indicator       = indicator,
                      post_est_only   = TRUE,      
                      run_tests       = FALSE,     
                      run_date        = run_date)

## Ensure you have defined all necessary settings in your config
check_config()

## Create a few objects from options above
if (class(Regions) == "character" & length(Regions) == 1) try(Regions <- eval(parse(text=Regions)))
if (class(summstats) == "character" & length(summstats) == 1) summstats <- eval(parse(text=summstats))
gaul_list <- get_adm0_codes(Regions, shapefile_version=modeling_shapefile_version, subnational_raking=subnational_raking)
if (class(year_list) == "character") year_list <- eval(parse(text=year_list))
test <- as.logical(test)


###############################################################################
## Launch Parallel Script
###############################################################################

## Make loopvars aka strata grid (format = regions, ages, holdouts)
if(as.logical(makeholdouts)) loopvars <- expand.grid(Regions, 0, 0:n_ho_folds) else loopvars <- expand.grid(Regions, 0, 0)


## loop over them, save images and submit qsubs
## first check to make sure that there are models to run (for resum_broken)
if (nrow(loopvars) > 0) { 

  for(i in 1:nrow(loopvars)){

    message(paste(loopvars[i,2],as.character(loopvars[i,1]),loopvars[i,3]))

    # Try to load profiling information from csv

    mem_cores <- load_profiling(vax = vaccine,
                                draws = samples,
                                reg = as.character(loopvars[i, 1]))


    dir.create(paste0('<FILEPATH>',indicator_group,'/',indicator,'/output/',run_date,'/output/'))
    dir.create(paste0('<FILEPATH>',indicator_group,'/',indicator,'/output/',run_date,'/errors/'))
    
    output_log_dir = paste0('<FILEPATH>',indicator_group,'/',indicator,'/output/',run_date,'/output/%x.o%j')
    error_log_dir  = paste0('<FILEPATH>',indicator_group,'/',indicator,'/output/',run_date,'/errors/%x.e%j')
    project        = '<PROJECT>'
    queue          = '<QUEUE>'
    memory         = as.numeric(mem_cores["mem"])
    cores          = as.numeric(mem_cores["slots"])
    time           = '15-00:00:00'
    sing_image     = '<FILEPATH>/<IMAGE>'
    script_file    = paste0('<FILEPATH>/lbd_core/mbg_central/share_scripts/parallel_model.R')
    age            = loopvars[i,2]
    reg            = loopvars[i,1]
    holdout        = loopvars[i,3]
    test           = as.numeric(as.logical(test))
    indicator      = indicator
    job_name       = paste0("job__", reg, "_", age, "_", holdout)
    
    save.image(sprintf("<FILEPATH>/%s/%s/model_image_history/pre_run_tempimage_%s_bin%s_%s_%s.RData",
                       indicator_group,indicator, run_date, age, reg, holdout))


    sbatch <- paste0('<ADDRESS> -s ',script_file,
                 ' --reg ',reg,' --age ',age,' --run_date ',run_date, ' --test ',test,' --holdout ',holdout,
                 ' --indicator ',indicator,' --indicator_group ',indicator_group)
    
    system(sbatch)

  }
  ## check to make sure models are done before continuing
  waitformodelstofinish(lv = cbind(as.character(loopvars[,1]),loopvars[,3]),sleeptime=60)
}



##############################################################################
## Indicate that this model is done
##############################################################################

message(paste0("Done with launch script for ", indicator))

## END OF FILE
###############################################################################