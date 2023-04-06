

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
lapply(package_list, require, character.only = TRUE)
load_mbg_functions(core_repo)

# Custom load indicator-specific functions
source(paste0(indic_repo,'functions/misc_vaccine_functions.R'))

## Script-specific stuff begins here ##########################################

load_from_parallelize()


if(coverage == TRUE){
  vaccine <- vacc
  set_up_indicators(stem = vaccine, 
                    doses = doses, 
                    single_dose = T,
                    save_doses = F,
                    save_2_cov = F)
}


# parse formatting
if (class(summstats) == "character" & length(summstats) == 1) summstats <- eval(parse(text=summstats))
if (class(year_list) == "character") year_list <- eval(parse(text=year_list))

raking_shapefile_version <- modeling_shapefile_version <-shapefile_version


# Define a log directory and clean out any files that haven't been touched in the last week
log_dir <- paste0("<FILEPATH>/logs/", run_date, "/")
dir.create(log_dir, recursive = T, showWarnings = F)



load(paste0('<FILEPATH>/vaccine/', indicator, '/output/', run_date, '/', indicator,'_cell_draws_eb_bin0_',region,'_0.RData'))

main_dir <- paste0('<FILEPATH>/vaccine/', indicator, '/output/', run_date, '/')



  # setting a reference for the number of draws
  ndraws <- ncol(cell_pred)
  overs <- paste0("V", 1:ndraws)
  
  # setup the output folder
  output_dir <- main_dir
  
  # setup values to test against
  cell_pred_dims <- dim(cell_pred)

  reg <- region
  
  
  '%!in%' <- function(x,y)!('%in%'(x,y))
  
  ##### Prep input data into raking:

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
  
  #####################################################################
  # collect and load the population data from the WorldPop rasters
  covdt <- load_populations_cov(reg, pop_measure, measure = 'count', simple_polygon, simple_raster, year_list, interval_mo, pixel_id = pixel_id)
  
  #####################################################################
  #load the cell id to admin units link
  link_table <- get_link_table(simple_raster, shapefile_version = modeling_shapefile_version)
  
  
  
  #####################################################################
  # Prepping the cell_pred and link table to be linked by making sure they have the appropriate identifiers.  Also performs a 
  # zippering at the region boundary where cells that have a portion of their area outside of the modeling region are reintegrated
  # as a whole cell and considered to be only in one region.  This works becasue a given cell is only modeled in one region.
  link <- prep_link_table(link_table = link_table,
                          simple_raster = simple_raster,
                          pixel_id = pixel_id)
  
  cell_ids <- link_table[[2]]
  
  # getting the connector for sub-national or national raking, This connector gets the IHME location_code for our
  # gbd targets and connects that to the ADM0_CODE or ADM1_CODE as nessecary 
  
  connector <- get_gbd_locs(rake_subnational = TRUE,
                            reg = reg,
                            shapefile_version = modeling_shapefile_version)
  
  # getting the connector for sub-national raking - used to implement countries_not_to_subnat_rake
  nat_connector <- get_gbd_locs(rake_subnational = F,
                                reg = reg,
                                shapefile_version = modeling_shapefile_version)
  
  # merge the connectors on to the link table
  link <- sub_nat_link_merge(rake_subnational=T,
                             link,
                             connector,
                             nat_connector,
                             countries_not_to_subnat_rake)
  
  # set cell pred as a data table, and rename things
  cell_pred <- prep_cell_pred(cell_pred = cell_pred,
                              cell_ids  = cell_ids,
                              pixel_id  = pixel_id,
                              covdt     = covdt)
  
  # merge cell_pred on the link
  cell_pred = merge(link, cell_pred, by.x = 'ID', by.y = 'cell_id',allow.cartesian = TRUE)
  
  
  ## Raking Population ###################################################################
  # This is done to ensure that the total pop in each raking geography is the same as GBD
  message("raking population")
  
  #convert to fractional population 
  cell_pred = cell_pred[,pop := pop * area_fraction] 
  
  #NA out population where the pixel value is NA (to prevent weirdness with denominators)
  cell_pred = cell_pred[is.na(V1), pop := NA]
  
  
  ##### Using fractional raking #####
  scalars <- fread(paste0(main_dir, indicator, "_", reg, "_pop_rf.csv"))
  
  cell_pred <- merge(cell_pred, scalars, by = c("location_id", "year"))
  
  # rake fractional populations
  cell_pred$pop_raked <- 0
  cell_pred = cell_pred[,pop := pop * pop_scalar] 
  
  ##Unraked Counts#########################################
  # creating unraked counts aggregations 
  message("creating a unraked counts aggregations")
  cell_pred$V1 <- cell_pred$V1.x
  #convert cell_pred to counts
  cell_pred <- cell_pred[, (overs) := lapply(.SD, function(x) x * pop), .SDcols=overs]
  
  ## Create unraked counts agregation
  admin_2 <- cell_pred[, lapply(.SD, sum, na.rm=T), .SDcols=c(overs, 'pop'), by = .(year, ADM2_CODE)]
  admin_1 <- cell_pred[, lapply(.SD, sum, na.rm=T), .SDcols=c(overs, 'pop'), by = .(year, ADM1_CODE)]
  admin_0 <- cell_pred[, lapply(.SD, sum, na.rm=T), .SDcols=c(overs, 'pop'), by = .(year, ADM0_CODE)]
  
  sp_hierarchy_list <- 
    link[ADM0_CODE %in% unique(admin_0[, ADM0_CODE]), 
         .(ADM0_CODE, ADM1_CODE, ADM2_CODE, ADM0_NAME, ADM1_NAME, ADM2_NAME)] %>% 
    unique %>% 
    .[, region := reg]
  
  ## save unraked counts aggregations
    paste0(output_dir, "/", indicator, "_unraked_", "_c_admin_draws_eb_bin", 
           "0_", reg, "_0.RData") %>%  
      save(admin_0, admin_1, admin_2, sp_hierarchy_list, file = .)

  
  ##Raked Rates Aggregation#########################################
  # creating unraked rates aggregations
  message("creating a unraked rates aggregations")
  
  # convert back to rates
  # TODO test equivalence to using weighted.mean on pop for this
  admin_0 <- admin_0[, (overs) := lapply(.SD, function(x) x / pop), .SDcols=overs]
  admin_1 <- admin_1[, (overs) := lapply(.SD, function(x) x / pop), .SDcols=overs]
  admin_2 <- admin_2[, (overs) := lapply(.SD, function(x) x / pop), .SDcols=overs]
  
  ## save unraked rates aggregations
    paste0(output_dir, "/", indicator, "_unraked_admin_draws_eb_bin0", 
           "_", reg, "_0.RData") %>%
      save(admin_0, admin_1, admin_2, sp_hierarchy_list, file = .)

  ######## END