#' Directory where central computation keeps R environments for the shared functions
CC_ENV_DIR <- '<FILEPATH>'
#################################################################################
### Pull GBD estimates from database and return as a data table
## Inputs:
##    gbd_type = "covariate" or "output", depending which database you need to pull from
##    gbd_name = GBD cov_name_short if "covariate" and GBD cause_id if "output"
##    gaul_list = list of GAUL codes you want to pull
##    measure_id = if "output", which measure you want to pull (defaults to incidence)
##    age_group_id = if "output", which age group you want to pull (defaults to Under 5)
##    metric_id = if "output", which metric you want to pull (defaults to rate)
##    named_location_field: string specifying which name of numerical location identifier that should be returned in 'name' column. could be 'GAUL_CODE' or 'location_id'
## Outputs:
##    Returns 3-column data.table where "name" = (usually) GAUL code, "year" = year, and
##      "mean" = value. Consistent with expected input for calc_raking_factors.
#################################################################################
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param gbd_type PARAM_DESCRIPTION
#' @param gbd_name PARAM_DESCRIPTION
#' @param gaul_list PARAM_DESCRIPTION
#' @param measure_id PARAM_DESCRIPTION, Default: 6
#' @param age_group_id PARAM_DESCRIPTION, Default: 1
#' @param metric_id PARAM_DESCRIPTION, Default: 3
#' @param year_ids PARAM_DESCRIPTION, Default: c(2000:2017)
#' @param return_by_age_sex PARAM_DESCRIPTION, Default: 'no'
#' @param collapse_age_sex PARAM_DESCRIPTION, Default: FALSE
#' @param shapefile_version PARAM_DESCRIPTION, Default: 'current'
#' @param gbd_round_id PARAM_DESCRIPTION, Default: 5
#' @param named_location_field PARAM_DESCRIPTION, Default: 'GAUL_CODE'
#' @param ... any other fields to be passed to all of the get_* functions from GBD CC
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @rdname load_gbd_data
#' @export
load_gbd_data     <- function(gbd_type,
                              gbd_name,
                              gaul_list,
                              measure_id = 6,
                              age_group_id = 1,
                              metric_id = 3,
                              year_ids = c(2000:2017),
                              return_by_age_sex = "no",
                              collapse_age_sex = FALSE,
                              shapefile_version = 'current',
                              gbd_round_id = 6,
                              decomp_step = step1,
                              named_location_field = "GAUL_CODE") {


  # get GAUL to location_id mapping
  gaul_to_loc_id <- get_location_code_mapping(shapefile_version = shapefile_version)
  gaul_to_loc_id <- gaul_to_loc_id[, list(location_id = loc_id, GAUL_CODE)]

  loc_ids <- gaul_to_loc_id[GAUL_CODE %in% gaul_list, location_id]

  # load covariate data
  if (gbd_type == "covariate") {

    # get covariate metadata from covariate_name_short
    metadata <- get_covariate_metadata()
    covariate_id <- metadata[covariate_name_short == tolower(gbd_name), covariate_id]
    covariate_by_age <- metadata[covariate_name_short == tolower(gbd_name), by_age]
    covariate_by_sex <- metadata[covariate_name_short == tolower(gbd_name), by_sex]

    # get covariate data
    source(path_join(CC_ENV_DIR, 'get_covariate_estimates.R'))
    if (covariate_by_age) {
        gbd_estimates <- get_covariate_estimates(covariate_id = covariate_id, location_id = loc_ids, year_id = year_ids, age_group_id = age_group_id, gbd_round_id = gbd_round_id, decomp_step = decomp_step)
    } else {
        gbd_estimates <- get_covariate_estimates(covariate_id = covariate_id, location_id = loc_ids, year_id = year_ids, gbd_round_id = gbd_round_id, decomp_step = decomp_step)
    }

    # collapse to all age, both sexes (if specified, and if there are multiple age and/or sex groups)
    if (collapse_age_sex & gbd_estimates[, uniqueN(age_group_name) > 1 | uniqueN(sex_id) > 1]) {

      # get population data
      source(path_join(CC_ENV_DIR, 'get_population.R'))
      gbd_pops <- get_population(age_group_id = gbd_estimates[, unique(age_group_id)],
                                 location_id = gbd_estimates[, unique(location_id)],
                                 year_id = gbd_estimates[, unique(year_id)],
                                 sex_id = gbd_estimates[, unique(sex_id)], decomp_step = decomp_step)

      # population-weight the covariate data
      gbd_estimates <- merge(gbd_estimates, gbd_pops, by=c('location_id', 'sex_id', 'age_group_id', 'year_id'))
      gbd_estimates <- gbd_estimates[, list(mean_value = weighted.mean(mean_value, population, na.rm=T)), by='location_id,year_id']
    }

    # format and return
    gbd_estimates <- merge(gbd_estimates, gaul_to_loc_id, by="location_id")
    setnames(gbd_estimates, c(named_location_field, "year_id", "mean_value"), c("name", "year", "mean"))

    if(return_by_age_sex=='no') gbd_estimates <- gbd_estimates[, list(name, year, mean)]
    if(return_by_age_sex=='yes') gbd_estimates <- gbd_estimates[, list(name, year, mean, sex_id, age_group_id)]

    return(gbd_estimates)
  }

  # load cause data
  if (gbd_type == "output") {

    # get cause metadata
    source(path_join(CC_ENV_DIR, 'get_cause_metadata.R'))
    metadata <- get_cause_metadata(cause_set_id = 2, gbd_round_id = gbd_round_id, decomp_step = decomp_step)
    cause_id <- suppressWarnings(as.numeric(gbd_name))
    if (is.na(cause_id)) cause_id <- metadata[acause == gbd_name, cause_id]

    # get cause data
    source(path_join(CC_ENV_DIR, 'get_outputs.R'))
    gbd_estimates <- get_outputs(topic = "cause",
                                 version = "best",
                                 gbd_round_id = gbd_round_id,
                                 cause_id = cause_id,
                                 measure_id = measure_id,
                                 metric_id = metric_id,
                                 age_group_id = age_group_id,
                                 location_id = loc_ids,
                                 year_id = year_ids, 
                                 decomp_step = decomp_step)

    all_data <- merge(gaul_to_loc_id, gbd_estimates, by="location_id")
    setnames(all_data, c(named_location_field, "year_id", "val"), c("name", "year", "mean"))
    all_data <- all_data[, list(name, year, mean)]

    return(all_data)

  }

}


#' @title Get GBD locations
#' @description Loads in the shapefiles to get a complete set of codes for GBD
#' locations. this is particularly useful if you want to get
#' subnational GBD location codes which are a combination of different
#' administrative levels
#'
#' @param reg region you want to pull codes for. e.g. 'africa'
#' @param rake_subnational logical - do you want to pull location codes for subnational raking or just national raking
#' @param shapefile_version string specifying version date of shapefile to use
#'
#' @return returns a dataframe containing, at a minimum, an ADM_CODE column and an loc_id (ihme loc id) column. If rake_subnational=T, it will also contain ADM0*, ADM1* and rak_level (admin level for raking) columns

get_gbd_locs <- function(reg,
                         rake_subnational = T,
                         shapefile_version = raking_shapefile_version){
  require(sf)
  if(rake_subnational == T) {
    connector <- 
      st_read(get_admin_shapefile(admin_level = 0, raking = T, version = shapefile_version), quiet = T) %>% 
      st_set_geometry(NULL) %>% 
      mutate(ADM0_CODE = as.numeric(as.character(ADM0_CODE))) %>% 
      mutate(ADM1_CODE = as.numeric(as.character(ADM1_CODE))) %>% 
      filter(ADM0_CODE %in% get_adm0_codes(reg, shapefile_version = shapefile_version, subnational_raking = TRUE))
    
    ## get the lowest raking level and add it onto connector
    if('ad_level' %in% colnames(connector)){
      connector$rak_level = connector$ad_level ## gadm (ad_level) and gaul (rak_level) use different names for this!
    }
    
    connector <- connector %>%
      dplyr::select(ADM0_CODE, ADM1_CODE, loc_id, rak_level) %>% 
      mutate(rak_level = as.character(rak_level)) %>% 
      dplyr::rename(location_id = loc_id) %>% 
      mutate(location_id = as.numeric(as.character(location_id)))

    ADM_CODE <- connector$ADM0_CODE
    for(rl in unique(connector$rak_level)){
      ADM_CODE[which(rl == connector$rak_level)] <- connector[[paste0('ADM', rl, '_CODE')]][which(rl == connector$rak_level)]
    }

    connector <- cbind(connector, ADM_CODE)
    
  } else {
    connector <-
    get_location_code_mapping(shapefile_version = shapefile_version)[ADM_CODE %in% get_adm0_codes(reg, shapefile_version = shapefile_version),
                                                                     list(location_id = loc_id, ADM_CODE)]
  }
  return(data.table(connector))
}


#' @title Get GBD estimates
#' @description Loads national or subnational gbd estimates to be used for raking.
#'
#' @param gbd_name = GBD cov_name_short if "covariate" and GBD cause_id if "output"
#' @param region = name of region you want to pull results for. e.g. 'africa'
#' @param measure_id = if "output", which measure you want to pull (defaults to incidence)
#' @param age_group_id = if "output", which age group you want to pull (defaults to Under 5)
#' @param metric_id = if "output", which metric you want to pull (defaults to rate)
#' @param year_ids = numeric vector of years to pull
#' @param shapefile_version string of dated shapefile version to use when generating list of ihme loc ids and ADM code
#' @param rake_subnational Logical. do you want subnational estimates or just national ones
#' @param gbd_round_id numeric gbd round id to pull from
#' @param decomp_step character decomp_step to pull from, only relevant for gbd_round_id 5 and above
#' 
#' @return Returns 3-column data.table where "name" = ihme loc id, "year" = year, and
#'      "mean" = value. Consistent with expected input for rake_cell_pred and calculate_raking_factors. 

get_gbd_estimates <- function(gbd_type,
                              gbd_name,
                              region,
                              measure_id = 6,
                              age_group_id = 1,
                              metric_id = 3,
                              year_ids = c(2000:2017),
                              shapefile_version = 'current',
                              rake_subnational = TRUE, 
                              gbd_round_id = 6,
                              decomp_step = NULL) {

  ## get GAUL to location_id mapping
  gaul_to_loc_id <- get_gbd_locs(
    reg = region,
    rake_subnational = rake_subnational,
    shapefile_version = shapefile_version
  )
  
  ## If we had subnational raking on, then we additionally pull in the national estimates
  ## because fractional rates raking needs it
  if(rake_subnational) {
    gaul_to_loc_id_nats <- get_gbd_locs(
      reg = region,
      rake_subnational = FALSE,
      shapefile_version = shapefile_version
    )
    
    gaul_to_loc_id <- rbindlist(list(as.data.table(gaul_to_loc_id)[, list(location_id, ADM_CODE)],
                                     as.data.table(gaul_to_loc_id_nats)), use.names = TRUE)
    
  } else {
    gaul_to_loc_id <- as.data.table(gaul_to_loc_id)[, list(location_id, ADM_CODE)]
  }
  
  gaul_to_loc_id <- unique(gaul_to_loc_id)
  loc_ids <- gaul_to_loc_id[, location_id]
  
  ## get cause metadata
  source(path_join(CC_ENV_DIR, 'get_cause_metadata.R'))
  metadata <- get_cause_metadata(cause_set_id = 2, gbd_round_id = gbd_round_id)
  cause_id <- suppressWarnings(as.numeric(gbd_name))
  if (is.na(cause_id)) cause_id <- metadata[acause == gbd_name, cause_id]

  ## get cause data
  source(path_join(CC_ENV_DIR, 'get_outputs.R'))
  gbd_estimates <- get_outputs(topic = "cause",
                               version = "latest",
                               gbd_round_id = gbd_round_id,
                               cause_id = cause_id,
                               measure_id = measure_id,
                               metric_id = metric_id,
                               age_group_id = age_group_id,
                               location_id = loc_ids,
                               year_id = year_ids,
                               decomp_step = decomp_step)

  all_data <- merge(gaul_to_loc_id, gbd_estimates, by="location_id")
  setnames(all_data, c("location_id", "year_id", "val"), c("name", "year", "mean"))
  all_data <- all_data[, list(name, year, mean)]

  return(all_data)

}


#################################################################################
### Pull population data and return as a raster
## Inputs:
# simple_raster:
## Outputs:
## Notes:
#################################################################################
get_population_data     <- function(simple_raster){

  # load population raster
  # RB 27SEPT: Note this is a temporary location, and is only Africa so updates will be necessary
  root <- '<FILEPATH>'
  pop<- brick(sprintf('<FILEPATH>/pop_stack.tif',root))

  # make sure population is cropped and extented to same as simple_raster
  # this is important, otherwise IDX wont work.
  if(!is.null(simple_raster)){
    pop <- mask(crop(pop,simple_raster),simple_raster)
    extent(pop)=extent(simple_raster)
  }
  return(pop)
}



#################################################################################
### load admin raster
## Inputs:
# admin_level: 0,1, or 2 are accepted.
# disag_fact: factor to increase resolution (needed whith small districts), 50 makes it 100m
# simple_raster: keeps resolution and extent
# shapefile_version: string indicating which version of shapefile to pull
## Outputs: returns raster layer with gaul codes for values
#################################################################################
load_admin_raster  <- function(admin_level, simple_raster, disag_fact=NULL,
                               shapefile_version = 'current'){

  if(!admin_level %in% c(0,1,2)) stop("admin_level must be either 0, 1, or 2")

  # load admin raster
  if(!is.null(disag_fact)){
    sr = disaggregate(simple_raster,fact=disag_fact)
  } else {
    sr = simple_raster
  }

  # UPDATED: master gaul admin shapefiles
  shapes <- shapefile(get_admin_shapefile(admin_level, version = shapefile_version))

  # The variable we rasterize on must be numeric.
  shapes@data[[paste0('ADM', admin_level,'_CODE')]] <- as.numeric(as.character(shapes@data[[paste0('ADM', admin_level,'_CODE')]]))

  # crop
  cropped_shapes <- crop(shapes, extent(sr), snap="out")

  ## Fix rasterize
  initial_raster <- rasterize_check_coverage(cropped_shapes, sr, field = paste0('ADM', admin_level,'_CODE'))
  if(length(cropped_shapes[!cropped_shapes[[paste0('ADM', admin_level,'_CODE')]]%in%unique(initial_raster),])!=0) {
    rasterized_shape <- merge(rasterize_check_coverage(cropped_shapes[!cropped_shapes[[paste0('ADM', admin_level,'_CODE')]]%in%unique(initial_raster),], sr, field = paste0('ADM', admin_level,'_CODE')), initial_raster)
  }
  if(length(cropped_shapes[!cropped_shapes[[paste0('ADM', admin_level,'_CODE')]]%in%unique(initial_raster),])==0) {
    rasterized_shape <- initial_raster
  }
  #rasterized_shape <- rasterize(cropped_shapes, sr, field=paste0('ADM', admin_level,'_CODE'))
  masked_shapes <- mask(x=rasterized_shape, mask=sr)

  return(masked_shapes)

}


make_population_weights  <- function(admin_level,
                                     simple_raster,
                                     pop_raster,
                                     gaul_list,
                                     custom_admin_raster = NULL) {

  #use custom admin raster for raking
  if(!is.null(custom_admin_raster)) {
    adm <- custom_admin_raster
  } else {
    # load admin raster, crop down to simple_raster
    adm <- load_admin_raster(admin_level, simple_raster) #had to add this here, not working in call otherwise ers
  }

  adm <- crop(adm,simple_raster)
  adm <- mask(adm,simple_raster)

  # load population data
  pop <- pop_raster

  # do using cell index (vectorized)
  cell_idx <- seegSDM:::notMissingIdx(simple_raster)
  adm_cell <- raster::extract(adm, cell_idx)

  pop_layer_names <- names(pop_raster)

  make_weights_vector <- function(x) {

    ad2_code <- raster::extract(adm$layer, cell_idx)
    ad2_code[is.na(ad2_code)] <- 0
    pop_cell <- raster::extract(pop[[x]], cell_idx)
    pop_cell[is.na(pop_cell)] <- 0

    message(paste0(length(unique(ad2_code)), " unique admin codes in this level."))

    pop_totals_adm <- aggregate(pop_cell~ad2_code, FUN=sum)
    pop_totals_adm$pop_cell[pop_totals_adm$pop_cell==0] <- 0.0001

    ## replicate totals for all cells in area
    pop_totals_adm_cell <- as.vector(pop_totals_adm)[match(adm_cell,
                                                           pop_totals_adm$ad2_code),]
    pop_totals_adm_cell$ad2_code <- NULL

    ## get  population weights for each cell
    pop_wt_adm <- pop_cell / pop_totals_adm_cell
    pop_wt_adm[is.na(pop_wt_adm)] <- 0
    pop_wt_adm_vector <- pop_wt_adm$pop_cell

    wt_sum_ad2 <- tapply(pop_wt_adm_vector, ad2_code, sum)
    stopifnot(all.equal(wt_sum_ad2, round(wt_sum_ad2)))

    names(pop_wt_adm_vector) <- x
    return(pop_wt_adm_vector)

  }

  pop_weights_by_layer <- lapply(pop_layer_names, make_weights_vector)
  pop_wt_matrix <- do.call(cbind, pop_weights_by_layer)
  rownames(pop_wt_matrix) <- cell_idx
  colnames(pop_wt_matrix) <- pop_layer_names

  # return
  return(list('pop_wt'     =as.matrix(pop_wt_matrix),
              'admin_level'=admin_level,
              'adm_code'   =adm_cell))

}






#################################################################################
### Split country/years (in format 'country_year') out of rownames of a condSim df
## Inputs:
# condSim_object: data matrix that make_condSim() spits out
## Outputs: returns matrix with name and year in two columns of same length
#################################################################################
split_geo_names <- function(condSim_object){
  str_match <- stringr::str_match
  splits <- strsplit(gsub("C\xf4te d'Ivoire","Cote dIvoire" ,rownames(condSim_object)), split = '_')
  ctry   <- sapply(splits, '[', 1)
  year   <- as.numeric(sapply(splits, '[', 2))
  ret    <- cbind('name' = ctry,
                  'year' = as.numeric(year))
  rownames(ret) <- NULL
  return (ret)
}





#################################################################################
### Custom summary functions
#################################################################################

cirange = function(x){
  z=quantile(x,probs=c(.025,.975),na.rm=T)
  return(z[2]-z[1])
}




p_below <- function(x, value, equal_to = F) {

  # probability <  (or <= if equal_to = T) target value
  value <- as.numeric(value)
  if (equal_to == T) output <- sum(x <= value)
  if (equal_to == F) output <- sum(x < value)

  output <- output / length(x)
  return(output)

}

p_above <- function(x, value, equal_to = F) {

  # probability > (or >= if equal_to = T) target value
  value <- as.numeric(value)
  if (equal_to == T) output <- sum(x >= value)
  if (equal_to == F) output <- sum(x > value)

  output <- output / length(x)
  return(output)

}

lower <- function(x) {
  # Simply get and return a percentile
  output <- quantile(x, 2.5 / 100, na.rm = T)
  return(output)
}

upper <- function(x) {
  # Simply get and return a percentile
  output <- quantile(x, 97.5 / 100, na.rm = T)
  return(output)
}

cfb <- function(v) {

  # Calculate the Coffey-Feingold-Bromberg metric for
  # a numeric vector `v` - useful as a normed measure
  # of variability for a set of proportions.

  # This implementation assumes equal weights

  # Reference:
  # Coffey, M. P., Feingold, M., & Bromberg, J. (1988).
  # A normed measures of variability among proportions.
  # Computational Statistics & Data Analysis, 7(2), 127-141.
  # https://doi.org/10.1016/0167-9473(88)90088-6

  # calculate mean (u) and sample size (n)
  u <- mean(v, na.rm = T)
  if (is.na(u)) return(NA)
  n = length(v)

  # define numerator (h)
  h = var(v)

  # calculate denominator (max h)
  # (assuming equal weights)
  r <- n*u - floor(n*u)
  h_max <- u*(1-u) - r*(1-r)/n

  # return H statistic (sqrt(h/h_max))
  return(sqrt(h / h_max))

}

#################################################################################
### Takes in raked or raw draw-level estimates and makes stat summary rasters
## Inputs:
# draw_level_cell_pred: Cells by Draws matrix which is output from predict_mbg() or from rake_predictions()
# mask: Should be the simple_raster
# return_as_raster: If TRUE returns as raster, else as table
# summary_stat: ie mean, cirange, quantile, sd
## Outputs: Summary table or raster of the cell_pred table put in
#################################################################################
make_cell_pred_summary    <- function(draw_level_cell_pred,
                                      mask                 = simple_raster,
                                      return_as_raster     = TRUE,
                                      summary_stat         = 'mean',
                                      ...){

  # make summary
  summ <- apply(draw_level_cell_pred, 1, summary_stat, ...)

  # put it in a raster
  if(return_as_raster){
    yrs = dim(draw_level_cell_pred)[1]/length(cellIdx(mask))
    #message(sprintf('Making a RasterBrick with %i layers',yrs))
    summ <- insertRaster(mask,  matrix(summ,  ncol = yrs))
  }

  return(summ)

}


#' @title Make admin prediction summary
#' @description Take an admin pred object and an sp_hierarchy list and generate a sorted,
#' cleaned table for a given set of summary statistics
#'
#' @param admin_pred the admin draws object (e.g. admin_1, etc)
#' @param sp_hierarchy_list spatial hierarchy object from the admin draws RData files
#' @param summary_stat summary statistic to apply
#' @param ... any other arguments to pass to the `summary_stats` functions
#'            (note that currently will be passed to all functions)
#' @return data table of summary stats and admin info
#' @examples
#' # load("indicator_raked_admin_draws_eb_bin0.RData") # Replace with your admin draws
#' make_admin_pred_summary(admin_2,
#'                         sp_hierarchy_list,
#'                         c("mean", "cirange", "upper", "lower"))

make_admin_pred_summary    <- function(admin_pred,
                                       sp_hierarchy_list,
                                       summary_stats = 'mean',
                                       ...){

  ### Get set up
  str_match <- stringr::str_match

  # Split up your input data
  ad_code <- subset(admin_pred, select = grep("ADM[0-9]_CODE", names(admin_pred)))
  year <- subset(admin_pred, select = "year")
  pop <- subset(admin_pred, select = "pop")
  draws <- subset(admin_pred, select = grep("V[0-9]*", names(admin_pred)))

  # Get the admin level
  ad_code_name <- names(ad_code)
  ad_level <- as.numeric(str_match(ad_code_name,"ADM([0-9])_CODE")[,2])

  # Get all admin levels
  all_admin_levels <- as.numeric(str_match(names(sp_hierarchy_list), "ADM([0-9])_CODE")[,2])
  all_admin_levels <- unique(all_admin_levels)[!is.na(unique(all_admin_levels))]

  ### Make summary
  summ_list <- lapply(summary_stats, function(ss) {
    apply(draws, 1, ss, ...)})

  if (length(summ_list) > 1) {
    output_df <- as.data.table(cbind(year, ad_code, do.call(cbind, summ_list)))
    names(output_df)[grep("V[0-9]*", names(output_df))] <- summary_stats
  } else if (length(summ_list) == 1) {
    output_df <- as.data.table(cbind(year, ad_code, summ_list[[1]]))
    names(output_df)[grep("V[0-9]*", names(output_df))] <- summary_stats
  }

  ### Add on identifiers

  # Drop smaller admin levels
  drop_levels <- all_admin_levels[all_admin_levels > ad_level]

  if (length(drop_levels) > 0) {
    for (d in drop_levels) {
      drop_cols <- names(sp_hierarchy_list)[grepl(paste0("ADM",d), names(sp_hierarchy_list))]
      sp_hierarchy_list <- subset(sp_hierarchy_list, select = !(names(sp_hierarchy_list) %in% drop_cols))
    }
  }
  sp_hierarchy_list <- unique(sp_hierarchy_list)

  output_df <- merge(sp_hierarchy_list, output_df, all.y = T, all.x = F)

  # Clean up & sort
  order_vars <- c(paste0("ADM", 0:ad_level, "_NAME"), "year")
  setorderv(output_df, order_vars)

  # Get col order
  ad_col_order <- sort(names(output_df)[grep("AD*_*", names(output_df))])
  cols_to_order <- c(ad_col_order, "region", "year", summary_stats)
  cols_to_order <- cols_to_order[cols_to_order %in% names(output_df)]
  if (length(cols_to_order) < length(names(output_df))) {
    other_cols <- names(output_df)[!(names(output_df) %in% cols_to_order)]
    cols_to_order <- c(cols_to_order, other_cols)
  }

  setcolorder(output_df, cols_to_order)

  return(output_df)

}

#' @title summarize_admins
#' @description Function to summarize admin_pred objects
#'
#' This is a wrapper for `make_admin_pred_summary()`
#'
#' @param ind indicator
#' @param ig indicator_group
#' @param summstats Summary statistics (functions) to compute.
#'                  Order will be the order they are written in csv
#'                  This is passed to `make_admin_pred_summary()`
#' @param raked Raked (T), unraked (F), or both (`c(T,F)`)?
#' @param ad_levels Admin levels to summarize (0, 1, and 2 by default)
#' @return Writes csv files to `<FILEPATH>`
#' @examples
#' summarize_admins(summstats = c("mean", "lower", "upper", "cirange"),
#'                  ad_levels = c(0,1,2),
#'                  raked     = c(T,F))
#' @rdname summarize_admins
summarize_admins <- function(ind = indicator,
                             ig = indicator_group,
                             summstats = c("mean", "lower", "upper", "cirange"),
                             raked = c(T,F),
                             ad_levels = c(0,1,2),
                             file_addin = NULL,
                             ...) {

  sharedir       <- sprintf('<FILEPATH>/%s/%s',ig,ind)
  input_dir <- paste0(sharedir, "/output/", run_date, "/")
  output_dir <- paste0(input_dir, "/pred_derivatives/admin_summaries/")
  dir.create(output_dir, recursive = T, showWarnings = F)

  # Convert raked to character
  rr <- character()
  if (T %in% raked) rr <- c(rr, "raked")
  if (F %in% raked) rr <- c(rr, "unraked")

  # If file_addin present, use it
  if (!is.null(file_addin)) file_addin <- paste0("_", file_addin)
  if (is.null(file_addin)) file_addin <- ""

  # Summarize and save admin preds
  for (rake in rr) {
    load(paste0(input_dir, ind, "_", rake, "_admin_draws_eb_bin0_0.RData"))
    sp_hierarchy_list <- mutate_if(sp_hierarchy_list, is.factor, as.character)
    sp_hierarchy_list <- mutate_at(sp_hierarchy_list, grep('_CODE', names(sp_hierarchy_list), value = T), as.numeric)

    for (ad in ad_levels) {
      message(paste0("Summarizing ", ind, ": admin ", ad, " (", rake, ")"))
      ad_summary_table <- make_admin_pred_summary(admin_pred = get(paste0("admin_", ad)),
                                                  sp_hierarchy_list,
                                                  summary_stats = summstats,
                                                  ...)
      fwrite(ad_summary_table,
             file = paste0(output_dir, ind, "_admin_", ad, "_", rake, file_addin, "_summary.csv"))
    }
  }
}

#################################################################################
### Saves post-estimation output files in the proper directory
## Inputs:
## Outputs:
#################################################################################
save_post_est   <- function(x,
                            filetype,
                            filename,
                            indic = indicator,
                            custom_output_folder = NULL){

  output_dir <- paste0('<FILEPATH>/', indicator_group, '/', indic, '/output/', run_date)
  
  if(!is.null(custom_output_folder)) output_dir <- custom_output_folder

  dir.create(output_dir, showWarnings = FALSE)

  filetype = tolower(filetype)
  if(!filetype %in% c('rdata','raster','csv'))
    stop('filetype argument has to be either rdata or raster or csv')


  if(filetype=='raster')
    writeRaster(
      x,
      file = paste0(output_dir, '/', indic,'_',filename),
      format='GTiff',
      overwrite = TRUE
    )

  if(filetype=='rdata')
    save(
      x,
      file = paste0(output_dir, '/', indic,'_',filename,'.RData'),
      compress = TRUE
    )

  if(filetype=='csv')
    write.csv(
      x,
      file = paste0(output_dir, '/', indic,'_',filename,'.csv')
    )


}


load_cell_preds <- function(indicator_group,
                            indicator,
                            rd = run_date,
                            region,
                            agebin,
                            u5m=FALSE,
                            other='',
                            ageasindic=TRUE){


  if(u5m){
    if(ageasindic==FALSE){
      load(paste0('<FILEPATH>/',indicator_group,'/',indicator,'/output/',rd,'/',
                  indicator,'_cell_draws_eb_bin',agebin,'_',region,'_0',other,'NA.RData')) # the 0 are no holdout
    } else {
      load(paste0('<FILEPATH>/',indicator_group,'/',indicator,'_age',agebin,'/output/',rd,'/',
                  indicator,'_age',agebin,'_cell_draws_eb_bin',agebin,'_',region,'_0',other,'.RData'))
    }
    cell_pred <- cptmp
  } else {
    load(paste0('<FILEPATH>/',indicator_group,'/',indicator,'/output/',rd,'/',
                indicator,'_cell_draws_eb_bin',agebin,'_',region,'_0.RData')) # the 0 are no holdouts
  }
  return(cell_pred)
}


logit_rake <- function(country_year, rake_dt, gbd) {
  rake_to_country <- strsplit(country_year, "_")[[1]][[1]]
  rake_to_year <- strsplit(country_year, "_")[[1]][[2]]
  p_N <- gbd[name == rake_to_country & year == rake_to_year, mean]
  raking_factors <- FindK(p_i = rake_dt[, p_i],
                          p_N = p_N,
                          N_i = rake_dt[, sim_N_i],
                          a = 2)
  rf <- as.data.table(rake_to_country)
  setnames(rf, 'rake_to_country', 'name')
  rf <- rf[, year := rake_to_year]
  rf <- rf[, raking_factor := raking_factors]
  return(rf)
}

G <- function(x){
  val <- logit(x)
  return(val)
}

G_Inv <- function(x){
  val <- ilogit(x)
  return(val)
}

Agg <- function(p_i, N_i){
  return(sum(p_i * N_i)/sum(N_i))
}

EvalK <- function(K, p_i, p_N, N_i, Mult=FALSE){
  if (Mult == TRUE){
    p_tilde <- G_Inv(G(p_i) * K)
  } else {
    p_tilde <- G_Inv(G(p_i) + K)
  }
  LHS <- Agg(p_tilde, N_i)
  RHS <- p_N
  SE <- (LHS-RHS)^2
  return(SE)
}

## p_i = our grid of probabilities
## p_N = national target
## sim_N_i = our grid of sample sizes (population raster)
## a = bounds for raking factor, perhaps return warning if optimal is on edge of limits
FindK <- function(p_i, p_N, N_i, a, Mult=FALSE) {
  if (Mult == TRUE){
    Limits <- c(0,a)
  } else {
    Limits <- c(-a,a)
  }
  iter <- 1
  Boundary = TRUE
  while (Boundary & iter < 10){
    Limits <- Limits * iter
    val <- optimize(EvalK, Limits, p_i = p_i, p_N = p_N, N_i = N_i, tol = 1e-20)$min
    Boundary <- (round(abs(val)) == Limits[2])
    iter <- iter + 1
  }
  return(val)
}


## ############################################################
## ############################################################
## Find approximate raw covariate weights in final INLA model
##
## ############################################################
## ############################################################
get.cov.wts <- function(rd, ## run_date
                        ind, ## indicator
                        ind_gp, ## indicator_group
                        reg,
                        age = 0,
                        holdout = 0) {
  ## ##########################################
  ## load the workspaces and objects we need ##
  ## ##########################################

  pathaddin <- paste0('_bin', age, '_', reg, '_', holdout)

  this_config <- fread(paste0('<FILEPATH>/', ind_gp, '/', ind, '/output/', rd, '/config.csv'))
  stacker_list <- this_config[V1 == 'stacked_fixed_effects', V2]
  stackers_used <- strsplit(stacker_list," ")
  stackers_used <- stackers_used[[1]][stackers_used[[1]] != "+"]

  ## load the data used to fit the stackers and reconstruct the design matrix
  load(paste0('<FILEPATH>/', ind_gp, '/', ind, '/model_image_history/', rd, pathaddin, '.RData'))
  fit.data <- as.data.frame(df)
  fit.data <- fit.data[, -(1:max(match(stackers_used, colnames(fit.data))))] # this drops all the ID and post-stacking variables, leaving just the pre-stacking ones
  if("gaul_code" %in% colnames(fit.data)) fit.data$gaul_code <- NULL
  if("period" %in% colnames(fit.data)) fit.data$period <- NULL
  X <- fit.data
  rc <- colnames(X)

  ## load the stacker fits
  load(paste0('<FILEPATH>/', ind_gp, '/', ind, '/output/', rd,
              sprintf('/child_model_list_%s_%i.RData', reg, holdout)))
  smo <- child_models
  rm(child_models)

  ## load the INLA/TMB fit
  load(paste0('<FILEPATH>/', ind_gp, '/', ind, '/output/', rd,
              sprintf('/%s_model_eb_bin%i_%s_%i.RData', ind, age, reg, holdout)))
  fit <- res_fit
  rm(res_fit)

  ## make a matrix to hold the p-values/importances
  imp.mat <- as.data.frame(matrix(ncol = length(rc), nrow = length(smo), dimnames = list(names(smo), rc)))

  ## #####################################################################################
  ## now, for each of the models, add the pvalues to the corresponding rows and columns ##
  ## #####################################################################################

  ## ~~~~~
  ## gam ~
  ## ~~~~~
  gam <- smo[['gam']]
  if(!is.null(gam) & class(gam)[1] != 'try-error'){
    smoothed   <- summary(gam)$s.table[, 4] ## smoothed table
    names(smoothed) <- substr(names(smoothed), 3, nchar(names(smoothed))-1)
    unsmoothed <- summary(gam)$p.table[, 4] ## parametric table
    all.p <- c(smoothed, unsmoothed)

    ## now match names and stick into our pval.matrix
    imp.mat["gam",] <- all.p[rc]

    ##  convert from p-val to `importance`
    imp.mat["gam", ] <- -log(imp.mat["gam", ])
    imp.mat["gam", ][is.na(imp.mat["gam", ])] <- 0
    if (sum(is.infinite(unlist(imp.mat["gam",]))) > 0) { # if a covariate has 'infinite' importance, assign it 1 and everything else 0
      imp.mat["gam",][!is.infinite(unlist(imp.mat["gam",]))] <- 0
      imp.mat["gam",][is.infinite(unlist(imp.mat["gam",]))] <- 1
    }
    imp.mat["gam", ] <- imp.mat["gam", ] / sum(imp.mat["gam", ])

  }

  ## ~~~~~
  ## gbm ~
  ## ~~~~~
  gbm <- smo[['gbm']]
  if(!is.null(gbm) & class(gbm)[1] != 'try-error'){
    rel.inf <- gbm$contributions
    imp.mat["gbm",] <- rel.inf[rc, "rel.inf"]

    ## convert to scaled importance
    imp.mat["gbm", ] <- imp.mat["gbm", ] / sum(imp.mat["gbm", ])

  }

  ## all the penalized regression DO NOT give SD or p-vals...
  ## I 'standardize' the coefs as per a suggestion in this thread:
  ## https://stats.stackexchange.com/questions/14853/variable-importance-from-glmnet

  ## ~~~~~~~
  ## lasso ~
  ## ~~~~~~~
  lasso <- smo[['lasso']]
  if(!is.null(lasso) & class(lasso)[1] != 'try-error'){
    l <- lasso$cv_1se_lambda ## this is the CV lambda from the child fit
    l.idx <- which(lasso$lambda == l)

    sds <- apply(X, 2, sd, na.rm = T)
    unscaled.coefs <- lasso$beta[, l.idx]
    unscaled.coefs <- unscaled.coefs[rc]

    ## scaled
    scaled.coefs   <- abs(unscaled.coefs) * sds

    ## put them in the matrix
    imp.mat["lasso", ] <- scaled.coefs[rc]

    ## convert to scaled importance
    imp.mat["lasso", ] <- imp.mat["lasso", ] / sum(imp.mat["lasso", ])

  }

  ## ~~~~~~~
  ## ridge ~
  ## ~~~~~~~
  ridge <- smo[['ridge']]
  if(!is.null(ridge) & class(ridge)[1] != 'try-error'){
    l <- ridge$cv_1se_lambda ## this is the CV lambda from the child fit
    l.idx <- which(ridge$lambda == l)

    sds <- apply(X, 2, sd, na.rm = T)
    unscaled.coefs <- ridge$beta[, l.idx]
    unscaled.coefs <- unscaled.coefs[rc]

    ## scaled
    scaled.coefs   <- abs(unscaled.coefs) * sds

    ## put them in the matrix
    imp.mat["ridge", ] <- scaled.coefs[rc]

    ## convert to scaled importance
    imp.mat["ridge", ] <- imp.mat["ridge", ] / sum(imp.mat["ridge", ])

  }

  ## ~~~~~~
  ## enet ~
  ## ~~~~~~
  enet <- smo[['enet']]
  if(!is.null(enet) & class(enet)[1] != 'try-error'){
    l <- enet$cv_1se_lambda ## this is the CV lambda from the child fit
    l.idx <- which(enet$lambda == l)

    sds <- apply(X, 2, sd, na.rm = T)
    unscaled.coefs <- enet$beta[, l.idx]
    unscaled.coefs <- unscaled.coefs[rc]

    ## scaled
    scaled.coefs   <- abs(unscaled.coefs) * sds

    ## put them in the matrix
    imp.mat["enet", ] <- scaled.coefs[rc]

    ## convert to scaled importance
    imp.mat["enet", ] <- imp.mat["enet", ] / sum(imp.mat["enet", ])

  }

  ## ~~~~~~~~~
  ## xgboost ~
  ## ~~~~~~~~~
  xgboost <- smo[['xgboost']]
  if(!is.null(xgboost) & class(xgboost)[1] != 'try-error'){
    load_R_packages("caret")
    # Extract row names to then assign to column names
    xg_names <- rownames(varImp(xgboost, scale = FALSE)$importance)

    
    # Extract coefficients, already correctly scaled
    scaled.coefs <- varImp(xgboost, scale = FALSE)$importance[[1]]
    

    # Assign column names
    names(scaled.coefs) <- xg_names

    # put them in the matrix
    imp.mat["xgboost", ] <- scaled.coefs[rc]
  }

  ## ##########################################
  ## Now we propagate through the INLA coefs ##
  ## ##########################################

  ## to account for different scaling in INLA we also create SDs of the covariates that go into INLA
  inla.X <- df[, paste0(names(smo), "_cv_pred"), with=F]
  if (this_config[V1 == "indicator_family", V2] == "binomial" & this_config[V1 == "stackers_in_transform_space", V2] == T) {
    inla.X <- logit(inla.X)
  }
  inla.sds   <- apply(inla.X, 2, sd, na.rm = T)

  ## get the coefficients
  if ("sdrep" %in% names(fit)) { # TMB, w/ stackers as fixed effects
    inla.coefs <- fit$sdrep$par.fixed[names(fit$sdrep$par.fixed) == "alpha_j"]
    names(inla.coefs) <- fit$fenames
    inla.coefs <- inla.coefs[stackers_used]
  } else if ("covar" %in% names(fit$summary.random)) { # INLA, w/ stackers as random effects (to sum to one)
    inla.coefs <- fit$summary.random$covar$"mean"
    names(inla.coefs) <- stackers_used
  } else { # INLA, w/ stackers as fixed effects
    inla.coefs <- fit$summary.fixed[stackers_used, "mean"]
    names(inla.coefs) <- stackers_used
  }

  ## get the scaled coefficients
  inla.coefs <- inla.coefs[names(smo)]
  scaled.inla.coefs <- inla.coefs * inla.sds

  ## make and scale the INLA weighted relative importance and add as a row
  inla.rel.imp <- apply(imp.mat, 2, function(x) sum(x*scaled.inla.coefs))
  inla.rel.imp <- abs(inla.rel.imp)
  inla.rel.imp <- inla.rel.imp/sum(inla.rel.imp)
  imp.mat <- rbind(imp.mat, "INLA COMBINED" = inla.rel.imp)

  ## add the scaled coefs as a column
  imp.mat <- cbind(imp.mat, "SCALED.INLA.COEFS" = c(scaled.inla.coefs, NA))

  ## save to general output directory
  save(imp.mat,
       file = paste0('<FILEPATH>/', ind_gp, '/', ind, '/output/', rd,
                     sprintf('/cov_wts_%s_holdout_%i.RData', reg, holdout)))

  return(imp.mat)
}


## #########################
## plot the output from cov wts in a heatmap in ggplot
## input is same as get.cov.wts() function)
## output is a ggplot object
## #########################

plot.cov.wts <- function(rd, ## run_date
                         ind, ## indicator
                         ind_gp, ## indicator_group
                         reg,
                         plot.inla.col = TRUE,
                         age = 0,
                         holdout = 0) {

  ## load a cov.wts object
  load(paste0('<FILEPATH>/', ind_gp, '/', ind, '/output/', rd,
              sprintf('/cov_wts_%s_holdout_%i.RData', reg, holdout)))
  cov.wts <- imp.mat

  ## remove the final column which has inla weight
  if(!plot.inla.col){
    cov.wts <- cov.wts[, -ncol(cov.wts)]
  }else{
    ## rescale to be -1:1
    cov.wts[, ncol(cov.wts)] <- cov.wts[, ncol(cov.wts)] / sum(cov.wts[, ncol(cov.wts)], na.rm = TRUE)
  }

  ## melt
  cw.m <- as.data.table(reshape2::melt(as.matrix(cov.wts)))
  colnames(cw.m)[1:3] <- c('Model', "Covar", "Imp")
  ## cw.m <- na.omit(cw.m)

  ## reorder factors
  cw.m$Model <- factor(cw.m$Model, levels = c("INLA COMBINED", sort(setdiff(unique(cw.m$Model), "INLA COMBINED"))))

  ## setup the plot

  ## pdf("~/test2.pdf", width = 10, height = 3)
  base_size <- 9
  p <- ggplot(cw.m, aes(Covar, Model)) +
    geom_tile(aes(fill = Imp), colour = "white") +
    scale_fill_gradient2(low = "red", mid = 'white', high = "steelblue", midpoint = 0) +
    theme_grey(base_size = base_size) +
    labs(x = "Covariate", y = "Model") +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme(axis.text.x = element_text(size = base_size * 0.8,
                                     angle = 270,
                                     hjust = 0,
                                     colour = "grey50"))
  ## print(p)
  ## dev.off()
  return(p)

}



prep_postest <- function(indicator,
                         indicator_group,
                         run_date,
                         save_objs) {

  # Save a list of objects in a standard location for parallel scripts to pull from
  main_dir <- paste0('<FILEPATH>/', indicator_group, '/', indicator, '/output/', run_date, '/')
  temp_dir <- paste0(main_dir, "temp_post_est/")
  temp_file <- paste0(temp_dir, "post_est_temp_objs.RData")
  dir.create(temp_dir, showWarnings = F)
  save(list = save_objs, file = temp_file)
}

post_load_combine_save <- function(regions = strata,
                                   summstats = c("mean", "cirange", "upper", "lower"),
                                   raked = c("raked", "unraked"),
                                   rf_table = TRUE,
                                   run_summ = TRUE,
                                   indic = indicator,
                                   ig = indicator_group,
                                   sdir = sharedir,
                                   proj = FALSE,
                                   proj_folder = NULL) {
  
  message(paste0("indic: ", indic))
  message(paste0("ig: ", ig))
  
  rake_addin <- character()
  if ("unraked" %in% raked) {
    lookup_dir <- paste0(sprintf("<FILEPATH>/%s/%s/output/%s/", ig, indic, run_date))
    ur <- length(grep(paste0(indic, ".*unraked.*raster.tif"), list.files(lookup_dir)))
    if(proj) ur <- length(grep(paste0(indic, ".*unraked_PROJ.*raster_PROJ.tif"), list.files(lookup_dir)))
    if (ur > 0) rake_addin <- c(rake_addin, unraked = "_unraked")
    if (ur == 0) rake_addin <- c(rake_addin, unraked = "")
  }
  
  if ("raked" %in% raked) {
    rake_addin <- c(rake_addin, raked = "_raked")
  }
  
  # loop through and combine all rasters
  message("\nCombining rasters...")
  for (rake in rake_addin) {
    message(names(rake_addin)[which(rake_addin == rake)])
    rr <- rake
    for (ss in summstats) {
      message(paste0("  ", ss))
      rlist <- list()
      for (reg in regions) {
        message(paste0("    ", reg))
        rlist[[reg]] <-
          brick(ifelse(proj,
                       sprintf("<FILEPATH>/%s/%s/output/%s/%s_%s%s_%s_raster_PROJ.tif", ig, indic, run_date, indic, reg, rake, ss),
                       sprintf("<FILEPATH>/%s/%s/output/%s/%s_%s%s_%s_raster.tif", ig, indic, run_date, indic, reg, rake, ss)
          ))
      }
      if (length(rlist) > 1) rlist <- do.call(raster::merge, unname(rlist)) else rlist <- rlist[[1]]
      if (ss == "cirange") ssname <- "range" else ssname <- ss # naming convention
      save_post_est(
        rlist, "raster",
        ifelse(!proj, 
               paste0(ssname, rr, "_raster"),
               paste0(ssname, rr, "_raster_PROJ")),
        indic
      )
    }
  }
  
  # do rf also
  if (rf_table) {
    message("RF table")
    rflist <- list()
    for (reg in regions) {
      rflist[[reg]] <-
        if(proj) {
          read.csv(sprintf("<FILEPATH>/%s/%s/output/%s/%s_%s_rf_PROJ.csv", ig, indic, run_date, indic, reg))  
        } else {
          read.csv(sprintf("<FILEPATH>/%s/%s/output/%s/%s_%s_rf.csv", ig, indic, run_date, indic, reg))  
        }
    }
    if(!proj) {
      save_post_est(do.call(rbind.fill, rflist), "csv", "rf", indic)
    } else {
      save_post_est(do.call(rbind.fill, rflist), "csv", "rf_PROJ", indic)
    }
    
  }
  
  # make a run summary graph
  if (run_summ) {
    graph_run_summary(
      run_date = run_date,
      indicator_group = ig,
      indicator = indic
    )
  }
}


clean_after_postest <- function(indicator,
                                indicator_group,
                                run_date,
                                strata,
                                delete_region_rasters = F) {

  # Delete intermediate files that we no longer need
  main_dir <- paste0('<FILEPATH>/', indicator_group, '/', indicator, '/output/', run_date, '/')
  temp_dir <- paste0(main_dir,'temp_post_est/')

  # Deleting - be careful!
  unlink(temp_dir, recursive = T)

  # If desired, insert code here to delete other temporary objects (eg. region-specific rasters)
  grep_string <- paste0(indicator, "_(", paste(strata, collapse = "|"), ").*_raster.tif")
  region_rasters <- grep(grep_string,list.files(main_dir), value=T) %>%
    paste0(main_dir, .)
  if(delete_region_rasters == T) unlink(region_rasters)
}



