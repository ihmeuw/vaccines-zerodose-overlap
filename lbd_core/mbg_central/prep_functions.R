make_time_stamp <- function(time_stamp) {
  
  run_date <- gsub("-","_",Sys.time())
  run_date <- gsub(":","_",run_date)
  run_date <- gsub(" ","_",run_date)
  
  if(time_stamp==FALSE) run_date <- 'scratch'
  
  return(run_date)
  
}

## Load parameters from config file into memory
#   Arguments:
#     repo            = Location where you've cloned "mbg" repository.
#     indicator_group = Category of indicator, i.e. "education"
#     indicator       = Specific outcome to be modeled within indicator category, i.e. "edu_0"
################ THIS FUNCTION HAS BEEN DEPRECATED IN FAVOR OF SET_UP_CONFIG ###################
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#'
#' @note THIS FUNCTION HAS BEEN DEPRECATED IN FAVOR OF SET_UP_CONFIG
#'
#' @param repo Location where you've cloned "mbg" repository.
#' @param indicator_group Category of indicator, i.e. "education"
#' @param indicator Specific outcome to be modeled within indicator category, i.e. "edu_0"
#' @param config_name PARAM_DESCRIPTION, Default: NULL
#' @param covs_name PARAM_DESCRIPTION, Default: NULL
#' @param post_est_only PARAM_DESCRIPTION, Default: FALSE
#' @param run_date PARAM_DESCRIPTION, Default: ''
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @rdname load_config
#' @export
load_config <- function(repo, indicator_group, indicator, config_name=NULL, covs_name = NULL, post_est_only=FALSE, run_date = '') {
  
  # Pull from config .csv file
  if (is.null(config_name)) {
    ## If new model run, pull config from <FILEPATH> repo
    if(post_est_only==FALSE) config <- fread(paste0(repo, '/', indicator_group, '/config_', indicator, '.csv'), header=FALSE)
    ## If running analysis on existing model, use config from that model's outputs folder
    if(post_est_only==TRUE) config <- fread(paste0('<FILEPATH>/', indicator_group, '/', indicator, '/output/', run_date, '/config.csv'))
  } else {
    config <- fread(paste0(repo, '/', indicator_group, '/', config_name, '.csv'), header=FALSE)
  }
  
  # If a covariate .csv file exists, use that instead
  if (!is.null(covs_name)) {
    
    # Grab fixed effects & measures (and gbd fixed effects & measures) from CSV if present
    covs <- fread(paste0(repo, '/', indicator_group, '/covariates/',covs_name, '.csv'), header = TRUE)
    # After update to data.table 1.11.4, "T" and "F" are not read in as
    # logical, but as characters, which we need to remedy here. We are
    # assuming that the "covs.csv" has an "include" and "gbd" column here
    covs[, gbd := as.logical(gbd)]
    covs[, include := as.logical(include)]
    covs <- subset(covs, include == T) # Use only those where include flag set to true
    fe <- subset(covs, gbd == F)
    gbd <- subset(covs, gbd == T)
    gbd[measure != "output", measure := "covariate"] # FIXME: This is a hack for backwards compatability -- basically it assumes you meant 'covariate' if you specified anything other than 'outcome' (eg, mean or NA)
    fixed_effects <- paste0(fe$covariate, collapse = " + ")
    fixed_effects_measures <- paste0(fe$measure, collapse = " + ")
    gbd_fixed_effects <- paste0(gbd$covariate, collapse = " + ")
    gbd_fixed_effects_measures <- paste0(gbd$measure, collapse = " + ")
    if(!("constraint" %in% names(covs))){
      fixed_effects_constraints <- paste0("c(", paste(rep(0, nrow(fe)), collapse=", "), ")")
      gbd_fixed_effects_constraints <- paste0("c(", paste(rep(0, nrow(gbd)), collapse=", "), ")")
    }
    else{
      fixed_effects_constraints <- paste0("c(", paste(unname(fe$constraint), collapse=", "), ")")
      gbd_fixed_effects_constraints <- paste0("c(", paste(unname(gbd$constraint), collapse=", "), ")")
    }
    
    # Remove any other versions from original config
    config <- subset(config, !(V1 %in% c("fixed_effects", "fixed_effects_measures", "gbd_fixed_effects", "gbd_fixed_effects_measures")))
    config <- config %>%
      rbind(., list("fixed_effects", fixed_effects)) %>%
      rbind(., list("fixed_effects_measures", fixed_effects_measures)) %>%
      rbind(., list("fixed_effects_constraints", fixed_effects_constraints)) %>%
      rbind(., list("gbd_fixed_effects", gbd_fixed_effects)) %>%
      rbind(., list("gbd_fixed_effects_measures", gbd_fixed_effects_measures)) %>%
      rbind(., list("gbd_fixed_effects_constraints", gbd_fixed_effects_constraints))
  }
  
  # Assign all the covariates to the environment
  for (param in config[, V1]) {
    message(paste("Assigning config value", param))
    assign(param, config[V1==param, V2], envir=globalenv())
  }
  
  return(config)
}

##### Overloading load_config to point to set_up_config in misc_functions
load_config <- function(...) {
  warning("load_config() and check_config() will be deprecated in favor of set_up_config()")
  set_up_config(...)
  
}



## Create directory structure
#   Arguments:
#     indicator_group = Category of indicator, i.e. "education"
#     indicator       = Specific outcome to be modeled within indicator category, i.e. "edu_0"
create_dirs <- function(indicator_group, indicator) {
  
  dir.create(paste0('<FILEPATH>/', indicator_group))
  dir.create(paste0('<FILEPATH>/', indicator_group, '/', indicator))
  
  indicator_dir <- paste0('<FILEPATH>/', indicator_group, '/', indicator)
  
  for(dir in c('output','model_image_history')) {
    dir.create(paste0(indicator_dir,'/',dir), showWarnings = FALSE)
  }
  
}

## Make template rasters (pull from central analysis folders based off the area to model specified in config)
#   Arguments:
#     simple = Single polygon that defines boundaries of the entire area you want to model over.
#   Returns: Empty raster over modeling area. To be used for cropping covariates quickly and projecting model.
get_template_raster <- function(simple) {
  
  message('Creating rasters of admin units')
  root <- '<FILEPATH>'
  
  # Centrally controlled folder of analysis shapefiles
  analysis_raster_dir <- paste0(root,'<FILEPATH>')
  
  # Load empty raster for largest analysis area, mask/crop to selected analysis area
  template_raster <- raster(paste0(analysis_raster_dir, 'stage2_analysis.tif'))
  template_raster <- mask(crop(template_raster,simple),simple)
  
  return(template_raster)
  
}


## Load input data from required location
#   Arguments:
#     indicator = Specific outcome to be modeled within indicator category, i.e. "edu_0"
#     simple    = Single polygon that defines boundaries of the entire area you want to model over.
#   Returns: Input data subset to modeling area.
load_input_data <- function(indicator, agebin = 0, removeyemen = FALSE, pathaddin = "",
                            withdate=FALSE, date='', years='five_year',range=5, update_run_date = FALSE,
                            withtag=FALSE, datatag='', use_share=FALSE, yl = year_list, region=NULL,
                            poly_ag = use_global_if_missing("poly_ag"),
                            zcol_ag = use_global_if_missing("zcol_ag"), simple=NULL) {
  
  # Ensure str_match loaded
  str_match <- stringr::str_match
  
  if(withdate){
    if(date=='')
      rd=run_date
    if(date!='')
      rd=date
  } else {
    rd = run_date
  }
  
  # Load input data by indicator
  root <- '<FILEPATH>'
  if(use_share==FALSE) load_dir <- paste0(root,'/WORK/11_geospatial/10_mbg/input_data/')
  if(use_share==TRUE) load_dir  <- '<FILEPATH>/input_data/'
  
  if(!withdate & !withtag) filename <- paste0(load_dir, indicator)
  if(withtag)              filename <- paste0(load_dir, indicator, datatag)
  if(withdate)             filename <- paste0(root,'/WORK/11_geospatial/10_mbg/input_data/dated/',rd,'/', indicator)
  
  # try to see if an RDS exists, if so use that, if not use a csv
  if(file.exists(paste0(filename,'.RDS'))){
    message('READING INPUT DATA FROM RDS FILE')
    d <- readRDS(paste0(filename,'.RDS'))
  } else {
    message('READING INPUT DATA FROM CSV FILE')
    d <- read.csv(paste0(filename,'.csv'))
  }
  
  d$latitude  <- as.numeric(as.character(d$latitude))
  d$longitude <- as.numeric(as.character(d$longitude))
  message(nrow(d))
  
  # Remove odd long/lats for point data (when not using polygon resampling)
  if(!"point" %in% names(d)) {
    # assume only point data
    d$point = 1
  }
  
  d=d[d$latitude<=90 | (d$point == 0 & as.logical(poly_ag)),]
  d=d[d$latitude>=-90 | (d$point == 0 & as.logical(poly_ag)),]
  d=d[d$longitude<=180 | (d$point == 0 & as.logical(poly_ag)),]
  d=d[d$longitude>=-180 | (d$point == 0 & as.logical(poly_ag)),]
  d <- subset(d, !is.na(latitude) | (d$point == 0 & as.logical(poly_ag)))
  d <- subset(d, latitude!=0 | (d$point == 0 & as.logical(poly_ag)))
  message(nrow(d))
  
  # Check for necessary columns
  if(!(indicator %in% names(d))) stop(paste0("Your input data does not contain a column for your indicator: ", indicator))
  
  d <- as.data.table(d)
  
  # Change all "country" assignments to national level (in case subnational in the input data)
  if (nrow(d[grepl("[A-Z]*_[.]*", country),]) > 0) {
    subnat_countries <- unique(d[grepl("[A-Z]*_[.]*", country), country])
    warning(paste0("Changing subnational to national country codes for the following: ",
                   paste0(subnat_countries, collapse = ",")))
    d[grepl("[A-Z]*_[.]*", country), country := str_match(country,"([A-Z]*)_[.]*")[,2]]
  }
  
  d$keep <- F
  # Subset to within modeling area that includes a buffer, if simple polygon is provided
  if(!is.null(simple)){
    message('subset based on simple polygon')
    d$rowid <- 1:nrow(d)
    d$keep <- F
    dpoint <- d[, c("rowid", "longitude", "latitude")]
    coordinates(dpoint) <- c("longitude", "latitude")
    proj4string(dpoint) <- proj4string(simple)
    dpoint$keep <- !is.na(over(dpoint, as(simple, "SpatialPolygons")))
    d$keep[dpoint$rowid] <- dpoint$keep
  }
  
  #If lacking simple polygon for data, then subset to within regions
  if(is.null(simple) & !is.null(region)) {
    adm0_list <- get_adm0_codes(region, shapefile_version = modeling_shapefile_version)
    # get GAUL to iso mapping
    loc_codes <- get_location_code_mapping(shapefile_version=modeling_shapefile_version)
    loc_codes <- loc_codes[ADM_CODE %in% adm0_list, ]
    regs      <- loc_codes$ihme_lc_id
    d[country %in% regs, "keep"] <- T
  }  
  if(is.null(simple) & is.null(region)) {
    warning("Missing region information. Will keep all data")
    d$keep <- T
  }
  
  
  message(paste0(round(mean(d$keep), 2)*100, '% of input data in specified template'))
  d <- d[d$keep, ]
  
  if(agebin!=0)   d = d[age%in%agebin,]
  if(removeyemen) d = d[country!='Yemen' & country!='YEM',]
  
  # remap any years as needed
  if(years=='five_year') {
    d <- d[year >= 1998 & year <= 2002, year := 2000]
    d <- d[year >= 2003 & year <= 2007, year := 2005]
    d <- d[year >= 2008 & year <= 2012, year := 2010]
    d <- d[year >= 2013 & year <= 2017, year := 2015]
  }
  
  if (nrow(subset(d, year < min(yl))) > 0) {
    warning(paste0("Dropping all data before min(year_list) = ", min(yl), "..."))
    d <- subset(d, year >= min(yl))
  }
  if (nrow(subset(d, year > max(yl))) > 0) {
    warning(paste0("Dropping all data after max(year_list) = ", max(yl), "..."))
    d <- subset(d, year <= max(yl))
  }
  
  # add in a weight column if it does not exist
  if(!"weight" %in% names(d)) {
    warning("A 'weight' column does not exists so one is added with all 1s.")
    d[,weight := 1]
  }
  
  # creaste a weighted SS to base QTs on
  if(sum(c('N','weight') %in% colnames(d)) == 2) d[,weighted_n := N*weight]
  
  # check that there is some non aggregate data
  if(!is.null(zcol_ag)) {
    if(all( (d$point==0 & as.logical(poly_ag)) | !is.na(d[[zcol_ag]]))) {
      #stop("There is only aggregate data. Currently the pipeline does not support modeling without some disaggregated data.")
    }
  } else {
    if(all(d$point == 0 & as.logical(poly_ag))) {
      #stop("There is only aggregate data. Currently the pipeline does not support modeling without some disaggregated data.")
    }
  }
  
  # Save a copy
  if(update_run_date == TRUE) {
    if(dir.exists(paste0('<FILEPATH>/', indicator_group, '/', indicator, "/output/", run_date)) == TRUE) {
      existing_dir <- paste0('<FILEPATH>/', indicator_group, '/', indicator, "/output/", run_date)
      new_try <- existing_dir
      index <- 0
      while(dir.exists(new_try)) {
        index <- index + 1
        new_try <- paste0(existing_dir, '_', index)
      }
      run_date <- paste0(run_date, '_', index)
      dir.create(new_try, showWarnings = FALSE)
      run_date_dir <- new_try
    }
    if(dir.exists(paste0('<FILEPATH>/', indicator_group, '/', indicator, "/output/", run_date)) == FALSE) {
      run_date_dir <- paste0('<FILEPATH>/', indicator_group, '/', indicator, "/output/", run_date)
      dir.create(run_date_dir, showWarnings = FALSE)
    }
    write.csv(d, file=paste0(run_date_dir, "/input_data", pathaddin, ".csv"))
    return(list(d, run_date))
  }
  
  if(update_run_date == FALSE) {
    if(agebin==0){
      run_date_dir <- paste0('<FILEPATH>/', indicator_group, '/', indicator, "/output/", run_date)
      dir.create(run_date_dir, showWarnings = FALSE)
      write.csv(d, file=paste0(run_date_dir, "/input_data", pathaddin, ".csv"))
    } else {
      run_date_dir <- paste0('<FILEPATH>/', indicator_group, '/', indicator, "_age",agebin,"/output/", run_date)
      dir.create(run_date_dir, showWarnings = FALSE)
      write.csv(d, file=paste0(run_date_dir, "/input_data", pathaddin, ".csv"))
    }
    return(d)
  }
  
}

## Read in shapefile and dissolve to single polygon for creating one big mesh (no admin1 boundaries)
#   gaul_list = any
#   buffer = distance to buffer around object
#   tolerance = tol in the gSimplify function. Higher # = coarser polygon
#   use_premade = shouold premade versions of regions / Africa be used?
#   raking = pull raking shapefile
#   shapefile_version = string specifying which shapefile version to pull

load_simple_polygon <- function(gaul_list, buffer, tolerance = 0.2,
                                subset_only = F, makeplots = F, use_premade = F,
                                custom_shapefile_path = NULL, custom_shapefile = NULL,
                                raking = F, shapefile_version = 'current', use_IND_states = F, sub_reg = reg,
                                model_just_IND = F) {
  
  # Logic check
  if (!is.null(custom_shapefile_path) & !is.null(custom_shapefile)) stop("You cannot specify both a custom shapefile and a custom shapefile path")
  
  # If using custom shapefiles, don't use premade polys
  if (!is.null(custom_shapefile_path) | !is.null(custom_shapefile)) use_premade <- F
  
  
  
  
  # Otherwise, make a new simple_poly
  # count the vertices of a SpatialPolygons object, with one feature
  vertices <- function(x) sum(sapply(x@polygons[[1]]@Polygons, function(y) nrow(y@coords)))
  
  # ~~~~~~~~~~~~~~~~~
  # load data
  
  if (is.null(custom_shapefile_path) & is.null(custom_shapefile)) {
    message("Opening master shapefile...")
    master_shape <- readOGR(get_admin_shapefile(admin_level = 0, raking = raking,
                                                version = shapefile_version))
    master_shape@data$ADM0_CODE <- as.numeric(as.character(master_shape@data$ADM0_CODE))
    subset_shape <- master_shape[master_shape@data$ADM0_CODE %in% gaul_list, ]
    
    
  } else if (!is.null(custom_shapefile_path) & is.null(custom_shapefile)) {
    message("Opening custom shapefile...")
    master_shape <- readOGR(custom_shapefile_path)
    subset_shape <- master_shape
  } else if (is.null(custom_shapefile_path) & !is.null(custom_shapefile)) {
    master_shape <- custom_shapefile
    subset_shape <- master_shape
  }
  
  if(subset_only==TRUE) {
    return(list(subset_shape=subset_shape,spoly_spdf=NULL))
  }
  
  if(subset_only==FALSE) {
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    message('Making a super-low-complexity map outline for INLA mesh creation')
    
    message(paste0('Full polygon vertices: ', vertices(subset_shape)))
    
    # Merge everything together (where intersecting)
    af <- gUnaryUnion(subset_shape)
    
    # Initial simplification
    af_simple <- gSimplify(af, tol = tolerance, topologyPreserve = TRUE)
    
    # Remove tiny features
    # Get all sub polygons and their areas
    polys <- af_simple@polygons[[1]]@Polygons
    areas <- sapply(polys, function(x) x@area)
    
    # If there is more than one sub polygon, remove the ditzels (many single-country subsets are a single polygon,
    #   like Uganda, which would break these few lines)
    if(length(areas)>1) {
      # find top 5% by area
      big_idx <- which(areas > quantile(areas, 0.95))
      
      # convert back into a spatialPolygons object
      spoly <- SpatialPolygons(list(Polygons(polys[big_idx], ID = 1)))
    }
    if(length(areas)==1) {
      spoly <- af_simple
      big_idx <- 1
    }
    
    # Buffer slightly
    spoly <- gBuffer(spoly, width = buffer)
    
    # simplify again to reduce vertex count
    spoly <- gSimplify(spoly, tol = tolerance, topologyPreserve = TRUE)
    
    # Get list of original polygons
    polys2 <- af@polygons[[1]]@Polygons
    
    # Check if all are within the simple polygon
    check_if_in_spoly <- function(the_poly, compare_to, the_proj = projection(master_shape)) {
      the_poly <- SpatialPolygons(list(Polygons(list(the_poly), ID = 1)))
      projection(the_poly) <- the_proj
      projection(compare_to) <- the_proj
      
      if(suppressWarnings(gIsValid(the_poly)) == F) return(TRUE) #Ignore invalid polygons
      
      poly_intersect <- rgeos::gIntersection(the_poly, compare_to)
      
      if(is.null(poly_intersect)) {
        return(FALSE)
      } else {
        return(ifelse((raster::area(poly_intersect) == raster::area(the_poly)), TRUE, FALSE))
      }
    }
    
    over_list <- sapply(polys2, function(x) check_if_in_spoly(x, compare_to = spoly))
    
    if (all(over_list) == FALSE) {
      # Add back in polygons if missed by above procedure (e.g. islands dropped)
      big_idx <- unique(c(big_idx, which(over_list == F)))
      spoly <- SpatialPolygons(list(Polygons(polys[big_idx], ID = 1)))
      spoly <- gBuffer(spoly, width = buffer)
      spoly <- gSimplify(spoly, tol = tolerance, topologyPreserve = TRUE)
    }
    
    # Now check again with new spoly
    over_list <- sapply(polys2, function(x) check_if_in_spoly(x, compare_to = spoly))
    
    # If still not all enclosed, tolerance probably too high. Return warning
    if (all(over_list) == FALSE) {
      number_false = length(over_list[over_list == F])
      number_total = length(over_list)
      warning(paste0(number_false, " of ", number_total, " polygons are NOT enclosed within your simple polygon. \n",
                     "Adjust your buffer and tolerance values."))
    }
    
    # Return results
    message(paste0('Simplified vertices: ', vertices(spoly)))
    
    # plot to check it encloses all of the important bits
    if(makeplots) plot(spoly)
    if(makeplots) plot(af, add = TRUE, border = grey(0.5))
    
    # turn into an SPDF
    spoly_spdf <- SpatialPolygonsDataFrame(spoly,
                                           data = data.frame(ID = 1),
                                           match.ID = FALSE)
    
    # add projection information
    projection(spoly_spdf) <- projection(master_shape)
    
    return(list(subset_shape=subset_shape,spoly_spdf=spoly_spdf))
  }
}


#' @title Buils spatial FEM mesh
#'
#' @description Build a finite-elements triangulation mesh to use in
#'   SPDE GP approximation
#'
#' @param d data.frame or data.table of observational data with
#'   'longitude' and 'latitude' columns
#'
#' @param simple simple_polygon defining modeling domain
#'
#' @param max_edge String of R vector notation containing two positive
#'   reals. e.g. "c(0.25, 5)" First number represents the max triangle
#'   edge length of the mesh within the modeling domain. Second number
#'   represents max triangle edge length in the exterior buffer
#'   portion of the mesh. Increase the first number for testing
#'   purposes (for comp speedup but it comes with increased
#'   approximation error and mesh artifacts in final
#'   predictions). Units are in lat-long degrees. Used if s2mesh=FALSE.
#'
#'   NOTE: the first number is also passed to the minimum allowed
#'   triangle edge length (the `cutoff`) within the interior of the
#'   domain. Since this first number is now specifying the min and the
#'   max, the meshing algorithm cannot usually obey both
#'   specifications completely. In reality, this seems to mean (from
#'   some empirical tests) that the meshing algorithm obeys the
#'   minimum allowed triangle edge length but some edges may be longer
#'   than the 'max' interior edge length.
#'
#' @param mesh_offset string of 2 element numeric vector in R
#'   notation. e.g. "c(1, 5)". Describes the automatic extension
#'   distance from the 'simple' boundary. The entries are the inner
#'   and outer extension distances, respectively. Units are in
#'   lat-long distance. Used if s2mesh=FALSE.
#'
#' @param plot_mesh Logical. Should a plot of the mesh be generated?
#'   Not currently implemented for s2mesh
#'
#' @param s2mesh Logical. Should the mesh be created on the surface of
#'   a sphere? If TRUE, s2params is used to specify mesh parameters
#'   instead of max_edge and mesh_offset
#'
#' @param s2params string of 3 element numeric vector in R
#'   notation. e.g. "c(25, 500, 1000)". The entries describe the
#'   minimum triangle edge length allowed, hos far to extend the mesh
#'   beyond the 'simple' boundary, and the maximum allowed triangle
#'   edge length, respectively. Units are in kilometers. Used only if
#'   s2mesh=TRUE.
#'
#' @return an 'inla.mesh' object
#'

build_space_mesh <- function(d, simple, max_edge, mesh_offset,
                             plot_mesh = F, s2mesh = FALSE,
                             s2params = NULL){
  
  if(!as.logical(s2mesh)){ ## build mesh on R2 plane
    message(paste0('Creating spatial mesh, max edge parameter: ', max_edge))
    max.edge <- eval(parse(text=max_edge))
    mesh_offset <- eval(parse(text=mesh_offset))
    mesh_s <- inla.mesh.2d(
      boundary = inla.sp2segment(simple),
      loc = cbind(d$longitude,d$latitude),
      max.edge = max.edge,
      offset = mesh_offset,
      cutoff = max.edge[1]
    )
    
    if (plot_mesh) {
      plot(mesh_s, asp=1)
      points(d$longitude, d$latitude, col=d$year)
    }
    
  } else{ ## build mesh on sphere surface
    
    if(is.null(s2params)){
      stop("You've chosen to build an s2 mesh but haven't provided parameters to do so (i.e. s2params=NULL)!")
    }
    
    message(paste0('Creating spatial SPHERICAL mesh, min edge, max edge, extension kms are: ', s2params))
    s2params <- eval(parse(text = s2params))
    
    ## convert data locs to 3d coords on the unit-sphere
    true.radius.of.earth = 6371
    s3 <- lonlat3D(d$longitude, d$latitude)
    
    ## convert boundary of simple_polygon to 3d coords on unit-sphere
    boundary <- inla.sp2segment(simple)
    boundary.loc <- lonlat3D(boundary$loc[, 1], boundary$loc[, 2])
    
    ## make a s2 domain mesh using data locs and boundary locs
    all.loc <- rbind(s3, boundary.loc)
    mesh_s <- inla.mesh.create(loc=all.loc,
                               cutoff = s2params[1] / true.radius.of.earth, ## minimum triangle edge allowed
                               extend = list(offset = s2params[2] / true.radius.of.earth), ## how far to extend mesh
                               refine=list(max.edge = s2params[3] / true.radius.of.earth)) ## max triangle edge allowed
    
    ## TODO add 3d mesh plotting
    ## an example of how to do this can be found in the code examples here:
    ## http://www.r-inla.org/examples/case-studies/simpson2011
    if(plot_mesh){
      ## plot 3d mesh
    }
  }
  
  return(mesh_s)
  
}

## Create temporal mesh (defaulting to the four period U5M approach for now, come back and make more flexible later)
build_time_mesh <- function(periods=1:4) {
  mesh_t <- inla.mesh.1d(
    loc = c(periods),
    degree = 1,
    boundary = rep("free", 2)
  )
  return(mesh_t)
}



#' @title Rasterize with border checks
#'
#' @description Rasterizing using a shapefile and a template raster, such that
#' we account for any pixels that are on the border of \code{field} units, which
#' could have been lost due to raster::rasterize only evaluating on centroids
#'
#' @param shapes SpatialPolygonDataFrame.. Input shapefile
#'
#' @param template_raster SpatialPolygonDataFrame.. The reference raster (usually WorldPop)
#'
#' @param field String The field with appropriate administrative unit information (usually ADM0_CODE)
#'
#' @param link_table String or data.table. If data.table it is used as-is. If String: either an absolute
#'   file path to an RDS file OR a short name for the administrative shape file e.g., "2019_02_27" or "current".
#'
#' @param id_raster String or raster object. If link table is a data.table,
#'   id_raster should be a raster RDS file. If link table is an absolute path,
#'   id_raster should also be an absolute path. Otherwise if link_table is a
#'   relative path, id_raster will be inferred.
#'
#' @return A raster with border and insides properly filled
#'
#' @details rasterize_check_coverage has three distinct use cases based off of the value of link_table
#'
#' 1. \code{link_table} is NULL. In this case rasterize_check_coverage will behave identically to raster::rasterize
#'
#' 2. \code{link_table} is a String referring to relase of admin shapefiles ("current" or e.g., "2019_02_27"). In this case
#'    \code{field} should be "ADM0_CODE", "ADM1_CODE" or "ADM2_CODE". This will load the lbd_standard_link.rds file,
#'    from the related admin shapefile directory, aggregate area_fraction as necessary to match the level of \code{field},
#'    and then apply those values to pixels in the space defined by \code{shapes}.
#'
#' 3. \link{link_table} is a data.table OR a String absolute path to a RDS file containing a data.table. This will use the
#'    provided \code{link_table} to assign values to the result raster similarly to use case #2.
#'
#' Note that for both use cases 2 and 3 all pixel_id coordinates must be in the same raster space. This is currently the
#' area defined by cropping the world raster to the pixels occupied by stage 1 and stage 2 countries.
#'
#' @export
#'
rasterize_check_coverage <- function(shapes, template_raster, field, ...,  link_table = modeling_shapefile_version, id_raster = NULL, use_IND_states = F) {
  # backwards-compatible behavior - just call rasterize()
  if (is.null(link_table)) return(raster::rasterize(shapes, template_raster, field = field, ...))
  
  # Validate arguments
  is_admin_link_table <- FALSE
  if (is.data.table(link_table)) {  # link table provided
    is_admin_link_table <- TRUE
    
    if (is.null(id_raster)) {
      stop("Link table given without ID raster.")
    }
  } else if (R.utils::isAbsolutePath(link_table)) {  # absolute path to link table provided
    link_table <- readRDS(link_table)
    
    if (is.null(id_raster)) {
      stop("Link table path given without ID raster path.")
    } else {
      id_raster <- readRDS(id_raster)
    }
  } else if (is_admin_shapefile_string(link_table)) {  # compute path to link table
    is_admin_link_table <- TRUE
    # load link table with pre-computed ownership percentages for each pixel cell
    link_table_file <- file.path(get_admin_shape_dir(link_table), "lbd_standard_link.rds")
    id_raster_file <- file.path(get_admin_shape_dir(link_table), "lbd_standard_id_raster.rds")
    link_table <- readRDS(link_table_file)
    id_raster <- readRDS(id_raster_file)
  } else {
    stop("link_table argument was neither a data.table, an admin shapefile string, or an absolute path to a RDS file.")
  }
  
  if (! field %in% names(link_table)) {
    msg <- paste("WARNING: rasterize_check_coverage called with field", field,
                 "which is not present in link_table. Defaulting to raster::rasterize()")
    message(msg)
    return(raster::rasterize(shapes, template_raster, field = field, ...))
  }
  
  # aggregate link table generically for admin 0/1/2
  # Note: we need `with=FALSE` because `field` is passed as a parameter (not a hard-coded string)
  table <- link_table[,c("pixel_id", field, "area_fraction"), with = FALSE]
  if (is_admin_link_table && field != "ADM2_CODE") {
    # sum rows; area_fraction now represents the total area coverage by ADM0/1_CODE instead of ADM2_CODE
    table <- table[, .(area_fraction = sum(area_fraction)), by = c("pixel_id", field)]
  }
  # subset table so that we have 1 entry per pixel_id - the value of `field` with the maximum
  # area_fraction value for that pixel_id
  # https://stackoverflow.com/a/24558696
  pixel_owner <- table[table[, .I[which.max(area_fraction)], by = pixel_id]$V1]

  
  pixel_owner <- pixel_owner[order(pixel_id)]

  
  # create reference pixel owner from id raster
  reference_pixel_owner <- id_raster
  raster::values(reference_pixel_owner) <- NA
  
  # subset to only those pixels owned by a shape we're interested in
  owned_pixels <- pixel_owner[pixel_owner[[field]] %in% shapes[[field]]]
  reference_pixel_owner[owned_pixels$pixel_id] <- owned_pixels[[field]]
  
  result <- raster::crop(reference_pixel_owner, template_raster, snap = "near")
  if (raster::ncell(result) != raster::ncell(template_raster)) {
    message <- paste("Error in creating result raster. Should have created a raster of shape",
                     paste(dim(result), collapse=","),
                     "but instead created a raster of shape",
                     paste(dim(template_raster), collapse=","))
    stop(message)
  }
  return(result)
}


#' @title Build a simple raster and associated pop raster
#' @description Builds a rasterized version of subset_shape to define
#'   the raster version of the modeling domain. Also returns a
#'   population raster for the modeling domain
#' @param subset_shape A \code{SpatialPolygonsDataFrame} to be
#'   rasterized. Usually output from
#'   \code{\link[load_simple_polygon()]}.
#' @param field Name of data entry in subset_shape to use as values in
#'   simple_raster, Default: NULL
#' @param raking Logical. Use raking shapefiles? Default: F
#' @param link_table String or data.table. If data.table it is used as-is, if string it is the shapefile version to specify correct link table
#' @param id_raster Raster. If building simple raster on link table, provide id_raster used for building the link table
#' @param pop_measure String. name of population measure to return in pop_raster. Default: 'total'
#' @param pop_release String. worldpop version date. If NULL, take most recent Default: NULL
#' @param pop_start_year Integer. first year of worldpop to load. Default: 2000
#' @param pop_end_year Integer. last year of worldpop to load. Default: 2018
#' @param shapefile_version String. which shapefile version is used.
#' @return A named list cotaining a 'simple_raster' and a 'pop_raster'
#' @details Rasterizes subset_shape, giving pixels values of
#'   \code{field}. Also loads and crops and returns a worldpop raster
#'   for the same domain with annual layers.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   ## Load simple polygon template to model over
#' adm_list           <- get_adm0_codes('wssa', shapefile_version = 'current')
#' simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list, buffer = 1, tolerance = 0.4,
#'   shapefile_version = 'current')
#' subset_shape        <- simple_polygon_list[[1]]
#' simple_polygon      <- simple_polygon_list[[2]]
#'
#' ## Load list of raster inputs (pop and simple)
#' raster_list        <- build_simple_raster_pop(subset_shape)
#' simple_raster      <- raster_list[['simple_raster']]
#' pop_raster         <- raster_list[['pop_raster']]
#' }
#' }
#' @seealso
#'  \code{\link[raster]{mask}}
#' @rdname build_simple_raster_pop
#'
#' @export
#'
#' @importFrom raster mask unique merge
#'
build_simple_raster_pop <- function(subset_shape,
                                    field=NULL,
                                    raking=F,
                                    link_table = modeling_shapefile_version,
                                    id_raster = NULL,
                                    pop_measure = 'total',
                                    pop_release = NULL,
                                    pop_start_year = 2000,
                                    pop_end_year = 2018,
                                    use_IND_states=F) {
  
  
  if (is.null(field)) {
    if ("GAUL_CODE" %in% names(subset_shape@data)) field <- "GAUL_CODE"
    if ("ADM0_CODE" %in% names(subset_shape@data)) field <- "ADM0_CODE"
  }
  
  if(raking) {
    field <- 'loc_id'
    # no 'loc_id' field in the link table, so we can't use it
    link_table <- NULL
  }
  
  ## if unspecified, get the most recent worldpop release
  if (is.null(pop_release)) {
    helper <- CovariatePathHelper$new()
    pop_rast_path  <- helper$covariate_paths(covariates = 'worldpop',
                                             measures = pop_measure,
                                             releases = pop_release)
    pop_release <- helper$newest_covariate_release(pop_rast_path)
  }
  
  # To ensure correct "snap" method is used, first convert to a template raster that is masked to population
  pop_rast <- brick(paste0("<FILEPATH>/", pop_release, "/1y/worldpop_total_1y_2010_00_00.tif"))
  template_raster <- raster::crop(pop_rast, raster::extent(subset_shape), snap = "out")
  
  # load in the population raster given the measure and the release
  cropped_pop <- load_worldpop_covariate(template_raster,
                                         covariate = 'worldpop',
                                         pop_measure = pop_measure,
                                         pop_release = pop_release,
                                         start_year = pop_start_year,
                                         end_year = pop_end_year,
                                         interval = 12)[['worldpop']]
  
  ## Fix rasterize
  if(use_IND_states==F){
    initial_raster <- rasterize_check_coverage(subset_shape, cropped_pop, field = field, link_table = link_table, id_raster = id_raster)
  }else{
    initial_raster1 <- rasterize_check_coverage(subset_shape[subset_shape$ADM0_NAME != 'India',], cropped_pop, field = field, link_table = link_table, id_raster = id_raster, use_IND_states = F)
    initial_raster2 <- rasterize_check_coverage(subset_shape[subset_shape$ADM0_NAME == 'India',], cropped_pop, field = 'ADM1_CODE', link_table = link_table, id_raster = id_raster, use_IND_states = use_IND_states)
    initial_raster <- merge(initial_raster1, initial_raster2)
  }
  
  if (length(subset(subset_shape, !(get(field) %in% unique(initial_raster)))) != 0 & use_IND_states == F) {
    rasterized_shape <-
      raster::merge(
        rasterize_check_coverage(subset(subset_shape, !(get(field) %in% unique(initial_raster))),
                                 cropped_pop,
                                 field = field,
                                 link_table = link_table,
                                 id_raster = id_raster),
        initial_raster)
  }
  if (length(subset(subset_shape, !(get(field) %in% unique(initial_raster)))) == 0 | use_IND_states == T) {
    rasterized_shape <- initial_raster
  }
  masked_pop <- raster::mask(x = cropped_pop, mask = rasterized_shape)
  
  raster_list <- list()
  raster_list[["simple_raster"]] <- rasterized_shape
  raster_list[["pop_raster"]] <- masked_pop
  
  return(raster_list)
}


## #############################################################################
## GET ADMIN CODES AND RELATED FUNCTIONS
## #############################################################################


#' @title Load the GAUL lookup table
#'
#' @description Loads the most recent version of the lookup table that links
#'   ADM0 codes with other identifiers such as GBD location IDs, MBG modeling
#'   regions, and ISO codes, among others.
#'
#' @return Returns a data.table of the ADM0 lookup table.
#'
load_adm0_lookup_table <- function() {
  lookup_table <- load_stage_list()
  # Set ISO codes as lowercase for easy lookup
  lookup_table$iso3 <- tolower(lookup_table$iso3)
  return(lookup_table)
}



#' @title Pull custom modeling regions
#'
#' @description Define modeling regions that are not simple combinations of the
#'   default MBG regions (in other words, regions that are not combinations of
#'   four-letter MBG regions such as "wssa" or "seas+ocea" or ISO-3 codes such
#'   as 'ZAF' or 'CHN').
#'
#' @param custom_region character vector of custom modeling regions
#'
#' @return Returns a named list of custom modeling regions with associated
#'    "standard" (non-custom) modeling regions that can be directly interpreted
#'    by get_adm0_codes().
#'
pull_custom_modeling_regions <- function(custom_regions){
  custom_regions <- tolower(custom_regions)
  
  # FULL LIST OF ALL REFERENCE REGIONS
  # If you need to add a new custom region, add it to this list
  ref_reg_list <- list(
    'africa'          = 'noaf+essa+wssa+cssa+sssa-yem',
    'middle_east'     = 'mide+stan-pak',
    'eastern_europe'  = "blr+est+lva+ltu+mda+ukr",
    'latin_america'   = 'caca+trsa+ansa',
    'south_asia'      = 'soas+chn_d2+pak',
    'central_america' = 'caca',
    'south_america'   = 'ansa+trsa',
    # se_asia was historically inclusive of East Asia + SE Asia
    'se_asia'         = 'eaas+seas+ocea+png',
    
    #data coverage regions
    'africa_dcp' = 'noaf+essa+wssa+cssa+sssa+yem',
    'middle_east_dcp' = 'mide+stan-yem-pak',
    'latin_america_dcp' = 'latin_america+cub',
    'south_asia_dcp' = 'south_asia-mdv-syc',
    'se_asia_dcp' = 'eaas+seas+png+idn+phl+tls+mys+twn',
    
    'stage1' = 'noaf+essa+wssa+cssa+sssa-yem',
    # ONLY stage 2 countries (not inclusive of Stage 1)
    'stage2' = 'ansa+caca+stan+eaas+mide+ocea+soas+seas+trsa+yem',
    'stage3' = 'all-stage1-stage2',
    
 'vax_seas' = 'seas+ocea-asm-fji-kir-wsm-ton',
   'vax_cssa' = 'cssa',
    'vax_essa' = 'essa+syc',
    'vax_wssa' = 'wssa',
    # #########################################################################################
   #########################################################################################
   
  )
  # Warn if there are any custom regions not in the reference list
  missing_regions <- custom_regions[ !(custom_regions %in% names(ref_reg_list)) ]
  if( length(missing_regions) > 0 ){
    message(paste0('WARNING: The following custom regions are not defined: ',
                   paste(missing_regions, collapse=','))
    )
  }
  # Return a named list of all custom regions
  custom_regions_list <- ref_reg_list[ names(ref_reg_list) %in% custom_regions ]
  return(custom_regions_list)
}


#' @title Get Admin0 (Country) Codes
#'
#' @description Pull Admin0 Codes for the specified countries or regions,
#'   optionally excluding certain Admin0 codes.
#'
#' @param adm0s Vector of ISO-3 codes, four-letter modeling region names, region
#'   group names as defined in get_region_groups, or 'all' (returns all Admin0
#'   codes)
#' @param strict (default `FALSE`) Causes the function to fail if an ISO code or
#'   region name is not included in the full list.
#' @param lookup_table (default `NULL`) Sets a data.frame or data.table to be
#'   used as the lookup_table. If `NULL`, the lookup table will be loaded using
#'   the `load_adm0_lookup_table()` function.
#' @param core_repo (default NULL) THIS ARGUMENT IS DEPRECATED AND WILL BE
#'   REMOVED IN THE FUTURE. Please remove it from your function calls.
#' @param adm0_type (default 'detect') Which class of admin0 codes
#'   should be pulled? Must be one of 'gaul', 'gadm', or 'detect'. If
#'   'gaul' or 'gadm', it will simply use what was specified. If
#'   'detect', the function will detect the adm0_type directly from an
#'   admin_shapefile. It will check the adm0_type specified by the
#'   shapefile_version argument if it is specified, otherwise if will
#'   detect the adm0_type using 'current'.
#' @param shapefile_version string specifying shapefile_version to be
#'   used in adm0_type determination if adm0_type is set to 'detect'
#' @param subnational_raking Logical. set to true if you want to use raking version of shapefile
#' @return Vector of numeric ADM0 codes. Warns if not all regions align with a
#'   ADM0 code.
#'
get_adm0_codes <- function(adm0s,
                           strict = FALSE,
                           lookup_table = NULL,
                           core_repo = NULL,
                           adm0_type = 'detect',
                           shapefile_version = NULL,
                           subnational_raking = FALSE
){
  
  # Check adm0_type input
  if(!(adm0_type %in% c('gaul', 'gadm', 'detect'))){
    stop("You must select either 'gaul', 'gadm', or 'detect' for adm0_type!")
  }
  # Check that data.table has been loaded
  if(!('data.table' %in% .packages())){
    stop("Please load the data.table package before running get_adm0_codes!")
  }
  # core_repo deprecation message
  if( !is.null(core_repo) ){
    message(paste0(
      "WARNING: The 'core_repo' argument has been deprecated in get_adm0_codes()",
      " -- please remove this argument."
    ))
  }
  
  # Determine adm0_type
  if(adm0_type == 'detect'){
    if(is.null(shapefile_version)){
      adm0_type <- detect_adm_shapefile_date_type(
        get_admin_shapefile(version = 'current', raking = subnational_raking)
      )$shpfile_type
    } else{
      adm0_type <- detect_adm_shapefile_date_type(
        get_admin_shapefile(version = shapefile_version, raking = subnational_raking)
      )$shpfile_type
    }
  }
  message(sprintf('Using ADM codes from %s', adm0_type))
  
  # Read in the ADM0 code data.table once at the beginning of the function
  if(is.null(lookup_table)) lookup_table <- load_adm0_lookup_table()
  
  # Check that a valid adm0_type has been defined
  # If so, define the field that will be pulled from the lookup table
  adm_type <- tolower(adm0_type)
  adm0_field_reference = list(
    gadm = 'gadm_geoid',
    gaul = 'GAUL_CODE'
  )
  if(!(adm0_type %in% names(adm0_field_reference))){
    stop(paste0(
      "Invalid adm0_type value. Valid values are: ",
      paste(names(adm0_field_reference,collapse=', '))
    ))
  }
  pull_field <- adm0_field_reference[[adm0_type]]
  
  ## PROCESS INCOMING ADM0 TEXT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # In case multiple character strings were passed in, concatenate with '+'
  adm0s <- paste(adm0s, collapse='+')
  # Set all as lowercase
  adm0s <- tolower(adm0s)
  
  # Split on - and + characters
  adm0s <- gsub(';','',adm0s)
  adm0s <- gsub('\\+',';;\\+',adm0s)
  adm0s <- gsub('-',';;-',adm0s)
  adm0s_vec <- unlist(base::strsplit(adm0s, ';;'))
  
  # Making the function past-compatible, for modelers who are still using the
  #  original 'name' modeling region for North Africa (now 'noaf').
  #  'name_historic' is now a custom region, equivalent to 'noaf-esh', in the
  #  pull_custom_modeling_regions().
  adm0s_vec <- gsub('^([+-])?name$', '\\1name_historic', adm0s_vec)
  
  # If the string was empty, end now
  if (length(adm0s_vec)==0){
    message("No GAUL codes were returned; exiting")
    return(numeric(0))
  }
  
  # The first region and all regions beginning in '+' are included
  # Strip the beginning + signs
  include <- gsub('\\+','', c(adm0s_vec[1], adm0s_vec[ grepl('^\\+',adm0s_vec) ]) )
  # All regions beginning with '-' are treated as exclusions
  # Strip any minus signs from them
  exclude <- gsub('-','', adm0s_vec[ grepl('^-',adm0s_vec) ] )
  
  ## CHECK FOR CUSTOM REGIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Define regex for all non-custom regions
  standard_code_regex <- "^([a-z]{3,4}|[a-z]{3}_d[0-9]{1,5})$"
  # Split both included and excluded data into custom and non-custom regions
  standard_inc <- include[ grepl(standard_code_regex, include) ]
  custom_inc   <- include[ !grepl(standard_code_regex, include) ]
  # Check in excluded data
  standard_exc <- exclude[ grepl(standard_code_regex, exclude) ]
  custom_exc   <- exclude[ !grepl(standard_code_regex, exclude) ]
  
  # Get definitions and ADM0 codes for included and excluded GAULs, if needed
  # `custom_inc_codes` and `custom_exc_codes` will be added/subtracted at the end
  # Check out that functional recursion
  recursive_pull <- function(x) get_adm0_codes(
    x, lookup_table=lookup_table, adm0_type=adm0_type
  )
  
  if ( length(custom_inc) > 0 ){
    custom_inc_defs  <- pull_custom_modeling_regions(custom_inc)
    custom_inc_codes <- unlist(lapply(custom_inc_defs, recursive_pull))
  } else {
    custom_inc_codes <- integer(0)
  }
  if ( length(custom_exc) > 0 ){
    custom_exc_defs  <- pull_custom_modeling_regions(custom_exc)
    custom_exc_codes <- unlist(lapply(custom_exc_defs, recursive_pull))
  } else {
    custom_exc_codes <- integer(0)
  }
  
  ## PROCESS NON-CUSTOM REGIONS (MBG MODELING REGIONS + ISOS) ~~~~~~~~~~~~~~~~~~
  # Helper function to pull all GAUL codes for a set of modeling regions or
  #  ISO codes
  utility_pull_adm0_codes <- function(region_codes, standard_regex){
    # Return none if the length of the inputs is 0
    if (length(region_codes)==0) return(integer(0))
    # Input data assertions
    if( !all(grepl(standard_regex, region_codes)) ){
      stop('All region codes must match the MBG or ISO code formats.')
    }
    isos     <- region_codes[nchar(region_codes)!=4]
    mbg_regs <- region_codes[nchar(region_codes)==4]
    if ('all' %in% isos){
      # 'all' is a special case where all GAUL codes are pulled
      pulled_adm0s <- unique(lookup_table[get(pull_field)>=0, get(pull_field)])
    } else (
      # Pull all GAUL codes matching the ISOs or regions
      pulled_adm0s <- unique(
        lookup_table[((iso3 %in% isos)|(mbg_reg %in% mbg_regs)) & (get(pull_field)>=0),
                     get(pull_field)]
      )
    )
    # Warn if any iso codes or MBG regions aren't in the lookup table
    missing_isos <- isos[ !(isos %in% lookup_table[,iso3]) & (isos!='all') ]
    if (length(missing_isos) > 0) message(paste0("WARNING: Missing these ISOs: ",
                                                 paste(missing_isos,collapse=',')))
    missing_regs <- mbg_regs[ !(mbg_regs %in% lookup_table[,mbg_reg]) ]
    if (length(missing_regs) > 0) message(paste0("WARNING: Missing these MBG regions: ",
                                                 paste(missing_regs,collapse=',')))
    return(as.integer(pulled_adm0s))
  }
  
  # Pull GAUL codes for standard include and exclude regions
  standard_inc_codes <- utility_pull_adm0_codes(standard_inc, standard_code_regex)
  standard_exc_codes <- utility_pull_adm0_codes(standard_exc, standard_code_regex)
  
  ## COMBINE TO CREATE FINAL ADM0 CODE SET ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Get all (standard + custom) codes to include
  all_include_codes <- unique( c(standard_inc_codes, custom_inc_codes) )
  # Get all (standard + custom) codes to exclude
  all_exclude_codes <- unique( c(standard_exc_codes, custom_exc_codes) )
  # Final codes = codes to include - codes to exclude
  return_codes <- all_include_codes[ !(all_include_codes %in% all_exclude_codes) ]
  
  return(return_codes)
}


#' @title Get GAUL Codes (TO BE DEPRECATED)
#'
#' @description Pull Admin0 codes for the specified countries or regions,
#'   optionally excluding certain Admin0 codes. The only difference from the
#'   standard get_adm0_codes() function is that the default adm0_type is 'gaul'
#'   rather than 'gadm'.
#'
#' @note THIS FUNCTION IS BEING DEPRECATED - ALWAYS USE get_adm0_codes() INSTEAD.
#'
get_gaul_codes <- function(
  adm0s,
  strict = FALSE,
  lookup_table = NULL,
  core_repo = NULL,
  adm0_type = 'gaul',
  shapefile_version = NULL
){
  get_adm0_codes(
    adm0s        = adm0s,
    strict       = strict,
    lookup_table = lookup_table,
    core_repo    = core_repo,
    adm0_type    = adm0_type,
    shapefile_version = shapefile_version
  )
}



add_gauls_regions <- function(df, simple_raster, shapefile_version = 'current') {
  
  # Extract GAUL_CODE from simple_raster using lat/longs
  df$ADM0_CODE <- raster::extract(simple_raster, df[ , c('longitude', 'latitude'), with=F])
  
  # Add names of regions by GAUL_CODE
  for(r in Regions) {
    df <- df[ADM0_CODE %in% get_adm0_codes(r, shapefile_version = shapefile_version), region := r]
  }
  
  ##Check that merge of GAUL_CODE worked properly
  df <- df[!is.na(region), good_records := 1]
  message(paste0(length(df$good_records) - length(df[good_records==1, N]), ' out of ', length(df$good_records), ' could not have GAUL/region extracted properly. Probably coastal points, need to fix.'))
  df <- df[good_records==1, ]
  
}


gaul_convert <- function(countries, from = "iso3", verbose = F, shapefile_version = 'current') {
  
  # Purpose: Convert a vector of countries (ihme_loc_id format) to vector of GAUL codes
  # Inputs:
  #         countries: vector of countries in ihme_loc_id format
  #         from: format of input
  #               options: "iso3" = "ihme_loc_id" = "ihme_lc_id"
  #                        "name" (the loc_name or loc_nm_short)
  #               ("iso3" is treated as "ihme_loc_id" for backwards compatability)
  #
  # Outputs: a vector of gaul codes
  
  # load reference table

    j_root <- "<FILEPATH>"

  
  str_match <- stringr::str_match
  
  # Catch if already passed gaul codes
  if(class(countries) =="numeric") return (countries)
  if(all(grepl("^[[:digit:]]+$", countries))) return(countries)
  
  gaul_table <- get_location_code_mapping(shapefile_version = shapefile_version)
  
  # convert input & output to lower case for easier matching
  #lowercase_cols <- c("short_name", "official_name", "iso3", "iso2", "uni", "undp")
  #gaul_table[, (lowercase_cols) := lapply(.SD, tolower), .SDcols = lowercase_cols,]
  
  ## ## convert input & output to lower case for easier matching
  if (verbose == T) message("\nlowercasing columns. the columns that get lowered are:")
  for(i in 1:ncol(gaul_table)){
    if(any(class(gaul_table[[i]]) %in% c("factor", "character"))) {
      if (verbose == T) message(sprintf("On column: %s", colnames(gaul_table)[i]))
      gaul_table[[i]] <- tolower(gaul_table[[i]])
    }
  }
  
  # Lowercase & ensure character for input
  countries <- tolower(countries)
  countries <- as.character(countries)
  
  ## Catch if a subnational in XXX_##### IHME syntax
  ## This returns national-level gaul codes only
  if(any(grepl("_", countries))) {
    countries[grepl("_", countries)] <- str_match(countries[grepl("_", countries)], "(.*)_.*")[,2]
  }
  
  if (from == "iso3" | from == "ihme_loc_id" | from == "ihme_lc_id") {
    
    gaul_code <- sapply(countries, function(x) gaul_table[ihme_lc_id == x, GAUL_CODE]) %>% as.numeric
    
  } else if(from == "name") {
    
    # Matching only national-level for now
    # Drop undefined & subnational rows
    gaul_table_nat <- subset(gaul_table, GAUL_CODE != -1)
    gaul_table_nat <- subset(gaul_table_nat, level == 3)
    
    gaul_code <- sapply(countries, function(x) gaul_table_nat[loc_nm_sh == x, GAUL_CODE]) %>% as.numeric
    
    #check to see if this matched all of the provided items in the vector; use partial / fuzzy matching if not
    if(length(gaul_code[is.na(gaul_code)]) > 0) {
      
      # Create a table to fill in
      table_matching <- cbind(countries, gaul_code) %>% as.data.table
      names(table_matching) <- c("country", "gaul_code")
      table_matching$gaul_code <- as.numeric(table_matching$gaul_code)
      
      approx_matched <- table_matching[is.na(gaul_code), country]
      
      # Indicate that approximate matching took place
      
      message("\nNot all country names provided were found in the lookup table.")
      message("Attempting to match names provided with those in lookup table.")
      message(paste0("Approximate matching attempted for: ", paste(approx_matched, collapse = ', '), "\n"))
      
      approx_match <- function(country) {
        # First, try matching to long form of name
        gaul_code <- gaul_table_nat[grep(country, gaul_table_nat$loc_name),]$GAUL_CODE
        
        # If that doesn't work, grep within the short name
        if (length(gaul_code) == 0) gaul_code <- gaul_table_nat[grep(country, gaul_table_nat$loc_nm_sh),]$GAUL_CODE
        
        # If that doesn't work, grep within the long name
        if (length(gaul_code) == 0) gaul_code <- gaul_table_nat[grep(country, gaul_table_nat$loc_name),]$GAUL_CODE
        
        # Could fill in other matching here if desired
        
        # Warn if nonspecific
        if (length(gaul_code) > 1) warning(paste0("\"", country, "\" matches multiple country names in the lookup table. Please be more specific."))
        
        # Finally, if no matches, return NA
        if (length(gaul_code) != 1) gaul_code <- NA
        
        return(as.numeric(gaul_code))
      }
      
      # Try approximate matching
      table_matching[is.na(gaul_code)]$gaul_code <- sapply(table_matching[is.na(gaul_code)]$country,approx_match)
      
      not_matched <- table_matching[is.na(gaul_code)]$country
      
      # Error checking
      if(length(not_matched) > 0) {
        warning(paste0("Some countries could not be matched:\n", paste(not_matched, collapse=', ')))
      }
      gaul_code <- table_matching$gaul_code %>% as.numeric
    }
    
  } else {
    # Error catching for non-supported country type
    stop("\nPlease enter a valid country code type")
  }
  
  if(length(gaul_code[is.na(gaul_code)]) > 0){
    # Error catching for failure to match all country codes
    message(paste0("CAUTION! Returning NA values.\nMatches not found for all country codes in input list.\n",
                   "Please check your input values"))
  }
  return(gaul_code)
}

## Make map of period indices to run any time periods in your data (like annual instead of 5-year)
make_period_map <- function(modeling_periods) {
  data_period <- sort(modeling_periods)
  period_ids <- seq(data_period)
  period_map <- as.data.table(data_period)
  period_map <- period_map[, period_id := period_ids]
  return(period_map)
}

interpolate_gbd <- function(gbd) {
  new_gbd <- list()
  for(this_year in c(2000:2015)) {
    if(this_year %in% 2000:2004) copied_data <- gbd[year == 2000, ]
    if(this_year %in% 2005:2009) copied_data <- gbd[year == 2005, ]
    if(this_year %in% 2010:2014) copied_data <- gbd[year == 2010, ]
    if(this_year %in% 2015:2015) copied_data <- gbd[year == 2015, ]
    copied_data <- copied_data[, year := this_year]
    new_gbd[[as.character(this_year)]] <- copied_data
  }
  new_gbd <- rbindlist(new_gbd)
  new_gbd <- new_gbd[order(name, year)]
  new_gbd <- new_gbd[!(year %in% c(2000,2005,2010,2015)), mean := NA]
  library(zoo)
  for(country in unique(new_gbd[, name])) {
    new_gbd <- new_gbd[name == country, mean := na.approx(new_gbd[name == country, mean])]
  }
  return(new_gbd)
}

make_regions_map <- function(regions, subset_shape, shapefile_version = 'current') {
  col_list <- c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf','#999999')
  i <- 1
  plot(subset_shape)
  for(reg in regions) {
    shapes <- subset_shape[subset_shape$GAUL_CODE %in% get_adm0_codes(reg, shapefile_version = shapefile_version), ]
    plot(shapes, add=TRUE, col=col_list[i])
    i <- i + 1
  }
}


merge_with_ihme_loc <-function(d, re=Regions, shapefile_version = 'current'){
  message(nrow(d))
  gaul_to_loc_id <- get_location_code_mapping(shapefile_version = shapefile_version)
  d <- d[, ihme_lc_id := as.character(country)]
  
  d <- merge(d, gaul_to_loc_id, by='ihme_lc_id',all.x=T)
  
  message(nrow(d))
  for(r in re) {
    d <- d[GAUL_CODE %in% get_adm0_codes(r, shapefile_version = shapefile_version), region := r]
  }
  return(d)
}





#' @title Set root path for J drive
#'
#' @description A function to set the "root" based on the OS used.
#' @author Rebecca Stubbs
#'
set_root<-function(){
  root<<-ifelse(Sys.info()[1]=="<FILEPATH>", "<FILEPATH>", "<FILEPATH>") # Setting "root" as global variable
}




#' @title Crop set and mask
#'
#' @description  Crop the raster to the extent of a simple_raster object,
#' set the extent, and then mask it based on the simple_raster object.
#'
#' @author Rebecca Stubbs
#'
#' @param raster_object A raster object (raster, Brick, etc) of a specific region
#' @param simple_raster A simple raster object that serves as the template for that region
#' @return A raster/Brick/Stack, with the extents and masks of the simple raster.
crop_set_mask<-function(raster_object,template_raster){
  raster_object<-raster::crop(raster_object, extent(template_raster))
  raster_object<-raster::setExtent(raster_object, template_raster)
  raster_object<-raster::mask(raster_object,template_raster)
  return(raster_object)
}


#' @title Vectorize covariate constraints
#'
#' @description Creates a vector of constraints based on the MBG config
#'
#' @author Neal Marquez
#'
#' @param config data.frame, config data frame
#'
#' @return vector of constraints with names for the covariates
#'
#'
covariate_constraint_vectorize <- function(config){
  fixed_effects <-
    unlist(strsplit(unlist(config[unlist(config[,1]) == "fixed_effects",2]), " \\+ "))
  gbd_fixed_effects <-
    unlist(strsplit(unlist(config[unlist(config[,1]) == "gbd_fixed_effects",2]), " \\+ "))
  fixed_effects_constraints <-
    eval(parse(text=unlist(config[unlist(config[,1]) == "fixed_effects_constraints",2])))
  gbd_fixed_effects_constraints <-
    eval(parse(text=unlist(config[unlist(config[,1]) == "gbd_fixed_effects_constraints",2])))
  
  if(is.null(gbd_fixed_effects_constraints) | length(gbd_fixed_effects_constraints) == 0){
    gbd_fixed_effects_constraints <- rep(0, length(gbd_fixed_effects))
  } else if(length(gbd_fixed_effects_constraints) == 1){
    gbd_fixed_effects_constraints <- rep(gbd_fixed_effects_constraints, length(gbd_fixed_effects))
  }
  
  if(is.null(fixed_effects_constraints) | length(fixed_effects_constraints) == 0){
    fixed_effects_constraints <- rep(0, length(fixed_effects))
  } else if(length(fixed_effects_constraints) == 1){
    fixed_effects_constraints <- rep(fixed_effects_constraints, length(fixed_effects))
  }
  
  constraints <- c(fixed_effects_constraints, gbd_fixed_effects_constraints)
  names(constraints) <- c(fixed_effects, gbd_fixed_effects)
  
  return(constraints)
}


