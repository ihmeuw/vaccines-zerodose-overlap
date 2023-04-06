
############################################################################################################################
##### post-hoc make dpt1 coverage

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
load_mbg_functions(core_repo)
lapply(package_list, require, character.only = TRUE)

# Custom load indicator-specific functions
source(paste0(indic_repo,'functions/misc_vaccine_functions.R'))
source(paste0(indic_repo,'functions/parallelize_hack.R')) 
## Script-specific stuff begins here ##########################################
indicator_group <- "vaccine"
vaccine <- 'dpt'
load_from_parallelize()


# Define directories
dpt3_dir <- paste0("<FILEPATH>/", indicator_group, "/dpt3_cov/output/", run_date, "/")
dpt12_dir <- paste0("<FILEPATH>/", indicator_group, "/dpt12_cond/output/", run_date, "/")
dpt1_dir <- paste0("<FILEPATH>/", indicator_group, "/dpt1_cov/output/", run_date, "/")

dir.create(dpt1_dir)

## LOAD ADMIN DRAWS ##############


slots              <- 4

# Define a log directory and clean out any files that haven't been touched in the last week
log_dir <- paste0( "<FILEPATH>/logs/", run_date, "/")

## Read config file and save all parameters in memory
config <- load_config(repo            = indic_repo,
                      indicator_group = "",
                      indicator       = NULL,
                      config_name     = paste0("config_", vaccine), 
                      covs_name       = paste0(vaccine, "_covs"),
                      run_test=FALSE)
# distribute config to all indicator directories
distribute_config(cfg = config, indicators = 'dpt1_cov')


# parse formatting
#if (class(Regions) == "character" & length(Regions) == 1) Regions <- eval(parse(text=Regions))
if (class(summstats) == "character" & length(summstats) == 1) summstats <- eval(parse(text=summstats))
if (class(year_list) == "character") year_list <- eval(parse(text=year_list))
samples <-as.numeric(samples)


load_from_parallelize() #Again, to grab regions from before the config

## LAUNCH MBG MODEL #########################################################################################################

  combine_lv <- list(region = region,
                  holdout = 0)
log_dir <- "<FILEPATH>/dpt_transforms/combine_dpt1_cells_rasters/"
dir.create(log_dir, recursive = T, showWarnings = F)
ratio_combine_para_output <- parallelize_hack(lv = combine_lv,
                                       save_objs = c("indicator_group", "run_date", "vaccine", "summstats", "shapefile_version"),
                                       script = "calc_other_vax/combine_dpt1_cells_rasters.R",
                                       prefix='combine_dpt1_cells_rasters',
                                       script_dir = indic_repo,
                                       output_log_dir = paste0('<FILEPATH>/output/%x.o%j'),
                                       error_log_dir = paste0('<FILEPATH>/errors/%x.e%j'),
                                       project              = '<PROJECT>',
                                       queue             = '<QUEUE>',
                                       memory = 450,
                                       cores = 3,
                                       time = '5-00:00:00',
                                       sing_image = 'default')


########################################################################
Regions=region
for(reg in Regions){
# First, load in the cell_preds for the base vaccine
message("Loading Data...")
load(paste0(dpt3_dir, "dpt3_cov_raked_admin_draws_eb_bin0_", reg, "_0.RData"))

dpt3_ad0 <- copy(admin_0)
dpt3_ad1 <- copy(admin_1)
dpt3_ad2 <- copy(admin_2)


# Now, load in the cell_preds for the ratios
message("Loading Data...")
load(paste0(dpt12_dir, "dpt12_cond_raked_admin_draws_eb_bin0_", reg, "_0.RData"))

dpt12_ad0 <- copy(admin_0)
dpt12_ad1 <- copy(admin_1)
dpt12_ad2 <- copy(admin_2)

overs<-paste0("V", 1:samples)

### getting in same order?
dpt3_ad0 <- dpt3_ad0[order(year),] 
dpt3_ad0 <- dpt3_ad0[order(ADM0_CODE),] 
ad0_keeps <- dpt3_ad0[,c("year","ADM0_CODE","pop")]  #s depend on number of draws
dpt3_ad0 <- dpt3_ad0[,lapply(overs, function(x) get(x))]

dpt3_ad1 <- dpt3_ad1[order(year),] 
dpt3_ad1 <- dpt3_ad1[order(ADM1_CODE),] 
ad1_keeps <- dpt3_ad1[,c("year","ADM1_CODE","pop")]
dpt3_ad1 <- dpt3_ad1[,lapply(overs, function(x) get(x))]

dpt3_ad2 <- dpt3_ad2[order(year),] 
dpt3_ad2 <- dpt3_ad2[order(ADM2_CODE),] 
ad2_keeps <- dpt3_ad2[,c("year","ADM2_CODE","pop")]
dpt3_ad2 <- dpt3_ad2[,lapply(overs, function(x) get(x))]

dpt12_ad0 <- dpt12_ad0[order(year),] 
dpt12_ad0 <- dpt12_ad0[order(ADM0_CODE),] 
dpt12_ad0 <- dpt12_ad0[,lapply(overs, function(x) get(x))]

dpt12_ad1 <- dpt12_ad1[order(year),] 
dpt12_ad1 <- dpt12_ad1[order(ADM1_CODE),] 
dpt12_ad1 <- dpt12_ad1[,lapply(overs, function(x) get(x))]

dpt12_ad2 <- dpt12_ad2[order(year),] 
dpt12_ad2 <- dpt12_ad2[order(ADM2_CODE),] 
dpt12_ad2 <- dpt12_ad2[,lapply(overs, function(x) get(x))]



## COMBINE CELL PREDS ######################################################################################

cell_pred_ad0 <- ( dpt12_ad0 * (1 - dpt3_ad0) ) + dpt3_ad0
cell_pred_ad1 <- ( dpt12_ad1 * (1 - dpt3_ad1) ) + dpt3_ad1
cell_pred_ad2 <- ( dpt12_ad2 * (1 - dpt3_ad2) ) + dpt3_ad2

admin_0 <- cbind(ad0_keeps[,c(1,2)], cell_pred_ad0, ad0_keeps[,c(3)])
admin_1 <- cbind(ad1_keeps[,c(1,2)], cell_pred_ad1, ad1_keeps[,c(3)])
admin_2 <- cbind(ad2_keeps[,c(1,2)], cell_pred_ad2, ad2_keeps[,c(3)])



## SAVE CELL PRED IN NEW RUN_DATE FOLDER ######################################################################################

## save raked cell preds
save(admin_0, admin_1, admin_2, sp_hierarchy_list, file = paste0(
  dpt1_dir, "dpt1_cov_raked_admin_draws_eb_bin0_", reg, "_0.RData"
))

}

 
 message('done!')