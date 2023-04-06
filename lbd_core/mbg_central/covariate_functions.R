
#' @title Load GBD Covariates
#' @description A faster version of load_gbd_covariates. I've left the old one for backwards capability (dccasey 8/23/2017)
#' @param covs A character vector listing GBD covariates/outputs to extract. For covariates, this
#'     should be indicated by covariate_name_short, while for outputs, this should be indicated by
#'     acause. Usually fulfilled by gbd_fixed_effects.
#' @param measures A character vector coresponding to 'covs' listing the type of GBD quantity for
#'     each item. Options are 'covariate' and 'outcome'. Usually fulfilled by gbd_fixed_effects_measures.
#' @param year_ids A numeric vector of year_ids. Usually fulfilled by year_list.
#' @param age_ids A string of age_ids. Usually fulfilled by gbd_fixed_effects_age.
#' @param template A raster layer of the buffered modelling area. usually it will be cov_layers[[1]][[1]].
#'     If NULL, a default template is loaded using load_and_crop_covariates_annual()
#' @param use_subnationals Logical. If true, the function will replace admin0 with a subnational units
#'     where possible. Use with caution because it's not been tested outside of Africa. It might not
#'     work for countries with multiple levels of subnational units (e.g. UK or India).
#' @param simple_polygon simple_polygon object used for the modeling region. made in load_simple_polygon
#' @param interval_mo number of months in a time unit. usually 12 to correspond 'annual' year_loadxs
#' @param year_list vector of years. If NULL, defaults to global value "year_list".
#' @return A list of covariates
load_gbd_covariates = function(covs, measures, year_ids, age_ids,
                               template, use_subnationals = F,
                               simple_polygon, interval_mo,
                               year_list = use_global_if_missing("year_list")){

  # check to see if the template is class raster, if not it is most
  # likely NULL since we pass in cov_layers[[1]][[1]] by default and
  # that is only created if we loaded in geospatial covs. otherwise,
  # we load in a geospatial cov to use as template
  if(class(template) != 'RasterLayer'){
    message('Loading in raster template for GBD covs since template argument was not a RasterLayer')
    template <-  load_and_crop_covariates_annual(covs            = 'evi',
                                                 measures        = 'median',
                                                 simple_polygon  = simple_polygon,
                                                 start_year      = min(year_ids),
                                                 end_year        = max(year_ids),
                                                 interval_mo     = as.numeric(interval_mo))[[1]][[1]]
  }
  
  # Load the analysis shapefile
  ## since this is fixed and is not related to shapefile_version, we
  ## also fix the shapefile_version passed to load_gbd_data inside
  ## fetch_gbd_covs to be one that worked for this shapeset
  world_shape <- readRDS('<FILEPATH>/GBD2016_analysis_final.rds')

  # If we are not using subnationals, keep only national units; otherwise remove the parent units
  if (!use_subnationals) {
    world_shape <- world_shape[world_shape$level==3,]
  } else {
    world_shape <- world_shape[!world_shape$loc_id %in% unique(world_shape$parent_id),]
  }

  world_shape <- crop(world_shape, template)
  # we must skip using the link_table as it is not relevant to this shapefile
  afras <- rasterize_check_coverage(world_shape, template, 'GAUL_CODE', fun = 'last', link_table = NULL)
  afras <- crop_set_mask(afras, template)
  
  # Check to make sure all of the relevant gauls are in the file
  if (!all(world_shape$GAUL_CODE %in% unique(afras))) {
    afras <- rasterize_check_coverage(world_shape[!world_shape$GAUL_CODE %in% unique(afras),], afras, 'GAUL_CODE', fun = 'first', update = T)
    afras <- crop_set_mask(afras, template)
  }


  # Loop over requested covariates
  fetch_gbd_cov = function(name, measure, afras) {

    # Load country-level results
    message("  Loading: ", name)
    gbd <- load_gbd_data(gbd_type = measure,
                         gbd_name          = name,
                         gaul_list         = unique(afras),
                         measure_id        = 5,
                         age_group_id      = age_ids,
                         metric_id         = 3,
                         year_ids          = year_ids,
                         return_by_age_sex = 'no',
                         gbd_round_id      = 6,
                         decomp_step = 'step1',
                         shapefile_version = shapefile_version,
                         collapse_age_sex  = TRUE)


    if (nrow(gbd) != nrow(unique(gbd[, list(name, year)]))) stop(paste0(name, "is not unique by location-year"))

    # For each gaul code and year, update the values
    blank = brick(lapply(year_list, function(x) afras * NA))

    for (yyy in 1:length(year_ids)) {
      for (ggg in unique(afras)) {
        blank[[yyy]][which(raster::getValues(afras) == ggg)] = gbd[name == ggg & year == year_list[yyy], mean]
      }
    }

    names(blank) <- rep(name, times = dim(blank)[3])

    return(blank)
  }

  all_gbd_layers <- lapply(1:length(covs), function(x) fetch_gbd_cov(covs[x], measures[x], afras))
  names(all_gbd_layers) <- covs

  return(all_gbd_layers)
}



#' @title Checks for covariate issues
#'
#' @description This function looks for issues with covariates that will cause issues down the line and tries to fix them.
#' The two issues looked for here are uniform covariates and pixel coverage. Uniform Covariates: If one of the extracted variables
#' does not vary over the data area. Pixel Coverage: If one of the covariates is missing in large parts of your model region,
#' these areas will be NA in your results later on. This often happens if you are using a modeled covariate in a new or
#' partially new region which was not previously modelled. This function takes in all parallel model objects that would need to
#' change if a covariate was removed for the code below in the parallel mode to work.
#'
#' @param cc object to change: cs_covs
#' @param afe object to change: all_fixed_effects
#' @param afeb object to change: all_fixed_effects_brt
#' @param acl object to change: all_cov_layers
#' @param tc object to change: the_covs
#' @param check_uniform boolean, do a check for uniform covariates, defaults to TRUE
#' @param check_pixelcount boolean, do a check for covariates with too few pixel (ie too low of geographic coverage in the region), defaults to TRUE
#' @param check_pixelcount_thresh for a pixelcount check, what proportion of the maximum observed pixel coverage is needed to keep the covariate in? must be between 0 and 1. defaults to 0.95
#'
#' @return the necessary objects as a named list that will later get assigned
#' to the environment within the parallel_model.R script
#' @export
check_for_cov_issues   <- function(cc                      = cs_covs,
                                   afe                     = all_fixed_effects,
                                   afeb                    = all_fixed_effects_brt,
                                   fe                      = fixed_effects,
                                   tc                      = the_covs,
                                   acl                     = all_cov_layers,
                                   check_uniform           = TRUE,
                                   check_pixelcount        = TRUE,
                                   check_pixelcount_thresh = 0.95){
  
  # check that threshold is between 0 and 1
  if(check_pixelcount_thresh < 0 | check_pixelcount_thresh > 1){
    stop('check_pixelcount_thresh must be between 0 and 1')
  }
  
  # make a few useful objects to track names and dropped covs
  covs <- as.character(cc$cs_df$name)
  covs <- covs[!grepl('gaul_code', covs)]
  dropcovs <- c()
  
  # run through extracted data and find any non-varying covariates
  if(check_uniform == TRUE){
    message('Checking for uniform covariates ... ')
    for(covname in covs){
      if (all( abs(na.omit(cc$covs[[covname]]) - mean(na.omit(cc$covs[[covname]]))) == 0 )){
        message(sprintf('WARNING: %s did not vary in your data and is being removed as a covariate', covname))
        dropcovs <- c(dropcovs,covname)
      }
    }
    if(length(dropcovs) == 0){
      message('  All clear')
    }
  }
  
  # check to see if any of the covariates have significant area missing, if so drop them 
  # this typically happens if you use a modelled surface as a covariate in a region that was missing one or more countries
  # count the number of pixels in each covariate, make sure they are within 95% non-NA of the max covered covariate
  if(check_pixelcount == TRUE){
    message('Checking for pixel counts ... ') 
    px_cnt <- c()
    for(covname in covs){
      dimtoget <- max(dim(acl[[covname]])[3]) # grab just the most recent (or synoptic) layer, assuming they are all the same
      px_cnt   <- c(px_cnt, length(cellIdx(acl[[covname]][[dimtoget]]))) 
    }
    threshdrops <- covs[px_cnt/max(px_cnt) < check_pixelcount_thresh]
    if(length(threshdrops) == 0){
      message('  All clear')
    } else {
      message(sprintf('WARNING: the following covariates had < %i per cent pixel coverage in this region and are being dropped:\n %s', 
                      round(check_pixelcount_thresh,2)*100, paste(threshdrops,collapse=', ')))
      dropcovs <- c(dropcovs, threshdrops)
    }
  }
  
  # make sure dropcovs is unique
  dropcovs <- unique(dropcovs)
  
  # if any non-varying covariates were detected, then remove them from any objects
  if(length(dropcovs) > 0){
    for(dc in dropcovs){
      
      # drop from the cs_covs object, which itself is a list
      cc$covs[, (dc) := NULL]
      cc$cs_df <- cc$cs_df[!as.character(cc$cs_df$name) %in% dc,]
      
      # drop from the covariate layers list of rasters
      acl <- replace(acl, dc, NULL)
    }
  } else {
    message('No non-varying covariates were detected in the data. Yay!')
  }

  # redo the all fixed effects pseudo-formula strings
  fe   <- paste0(format_covariates(fe)  [!format_covariates(fe)   %in% dropcovs], collapse = ' + ')
  afe  <- paste0(format_covariates(afe) [!format_covariates(afe)  %in% dropcovs], collapse = ' + ')
  afeb <- paste0(format_covariates(afeb)[!format_covariates(afeb) %in% dropcovs], collapse = ' + ')
  tc   <- tc[!tc %in% dropcovs]
  
  # return the objects as a named list that will later get assigned to the environment within the parallel_model.R script
  return(list( cs_covs               = cc,
               all_fixed_effects     = afe,
               all_fixed_effects_brt = afeb,
               fixed_effects         = fe,
               the_covs              = tc,
               all_cov_layers        = acl))
}


#' @title Covariate loader for standard Model-based Geostatistics.
#'
#' @description Loads covariate data in bulk and returns suitable raster/brick objects.
#'
#' @details
#' Covariates are stored in a 00_MBG_STANDARD directory (see
#' \code{MbgStandardCovariateLoader$public_fields$cov_dir} or an instances $cov_dir attribute)
#' as .tif files and loaded as raster objects (either as a layer or a brick of layers).
#'
#' @param start_year A numeric indicating the earliest year of data to attempt to retrieve.
#' @param end_year A numeric indicating the latest year of data to attempt to retrieve.
#' @param interval A numeric number of months that the data is provided in. Annnual data
#'  would have an interval of 12
#' @param covariate_config A data.table with covariate loading information. This must contain
#'  columns "covariate", "measure", and "release". All three should be string values with
#'  the "release" being a timestamp in the form of YYYY_MM_DD e.g., "2019_06_13".
#'
#' @examples
#' config.table <- data.table(
#'   covariate = c("access", "evi"),
#'   measure = c("mean", "median"),
#'   release = c("2019_06_10", "2019_06_10")
#' )
#'
#' start_year <- 1998
#' end_year <- 2017
#' interval_mo <- 60
#' template_raster <- suppressWarnings(empty_world_raster())
#'
#' loader <- MbgStandardCovariateLoader$new(
#'   start_year = start_year,
#'   end_year = end_year,
#'   interval = interval_mo,
#'   covariate_config = config.table
#' )
#' # all loaded covariates will be cropped and masked to the provided template raster.
#' lcc <- loader$get_covariates(template_raster)
#'
#' @rdname MbgStandardCovariateLoader
#' @export
MbgStandardCovariateLoader <- R6::R6Class("MBGStandardCovariateLoader",
  public = list(
    # init function
    initialize = function(start_year, end_year, interval, covariate_config, cov_dir = NULL) {
      private$start_year <- start_year
      private$end_year <- end_year
      private$interval <- interval
      private$covariate_config <- covariate_config
      private$path_helper <- CovariatePathHelper$new()
      if (!is.null(cov_dir)) {
        private$path_helper$cov_dir <- cov_dir
      }
      self$validate() # performs additional assignments
    },
    validate = function() {
      if (!dir.exists(private$path_helper$cov_dir)) {
        stop(sprintf("Covariate directory %s does not exist or is not accessible!", private$path_helper$cov_dir))
      }
      # validate interval
      if (!private$interval %in% private$valid_intervals) {
        stop(sprintf(
          "Only intervals %s supported. If you need monthly contact the LBD Core team",
          paste(private$valid_intervals, collapse = "/")
        ))
      }

      interval.years <- private$interval / 12
      private$all_periods <- seq(private$start_year, private$end_year, interval.years)

      private$duration <- paste0(interval.years, "y")

      measure.dirs <- private$path_helper$covariate_paths(
        covariates = private$covariate_config$covariate,
        measures = private$covariate_config$measure,
        releases = private$covariate_config$release
      )
      if (!all(dir.exists(measure.dirs))) {
        private$error_for_missing_covariates()
      }
      private$measure_dirs <- measure.dirs
      # TODO: can we check for duration/"synoptic" dirs and error faster IFF not present?
    },
    # public functions
    get_covariates = function(template_raster) {
      "returns a list of covariates, just like load_and_crop_covariates_annual"
      cov.list <- list()
      for (i in 1:nrow(private$covariate_config)) {
        covariate <- private$covariate_config[i, covariate]
        duration <- private$best_duration_dir(private$measure_dirs[i], private$duration)
        if (duration$is.synoptic) {
          message(sprintf("Loading %s which is synoptic", covariate))
          rast <- private$load_synoptic_covariate(duration$dir, i)
        } else {
          message(sprintf("Loading %s which is not synoptic", covariate))
          rast <- private$load_covariate(duration$dir, i)
        }
        cropped <- raster::crop(rast, raster::extent(template_raster))
        cov.list[[covariate]] <- raster::mask(cropped, template_raster)
      }
      return(cov.list)
    }
  ), # end public
  # private functions
  private = list(
    path_helper = NULL,
    start_year = NULL,
    end_year = NULL,
    interval = NULL, # in months
    valid_intervals = c(12, 24, 60),
    duration = NULL, # e.g., '5y'/ '2y' / '1y'
    covariate_config = NULL, # data.table
    measure_dirs = NULL,
    all_periods = NULL,

    best_duration_dir = function(measure.dir, duration) {
      # ideal: data exists for requested covariate/measure/duration
      dir <- file.path(measure.dir, duration)
      if (dir.exists(dir)) {
        return(list(dir = dir, is.synoptic = FALSE))
      }

      meta <- private$path_helper$covariate_metadata_from_path(measure.dir)

      # acceptable: data exists for covariate/measure as synoptic data
      dir <- file.path(measure.dir, "synoptic")
      if (dir.exists(dir)) {
        message(sprintf("%s measure (%s / %s) is synoptic only", meta$measure, meta$covariate, meta$release))
        return(list(dir = dir, is.synoptic = TRUE))
      }

      # error: no data available
      measure <- basename(measure.dir)
      err.msg <- paste(
        duration, "duration for measure", measure, "for covariate", covariate,
        "does not exist and is not synoptic"
      )
      stop(err.msg)
    },
    load_covariate = function(dir, i) {
      rasters <- list()
      covariate <- private$covariate_config[i, covariate]
      measure <- private$covariate_config[i, measure]

      periods <- private$get_periods(dir, private$all_periods)
      n_periods <- length(private$all_periods)
      for (period.index in 1:n_periods) {
        period <- private$all_periods[period.index]
        if (period %in% periods$missing) {
          best.period <- private$get_closest_period(period, periods$present)
          msg <- sprintf(
            "WARNING! We are substituting in %s data from period: %i to use as if it were for period: %i",
            covariate, best.period, period
          )
          message(msg)
        } else {
          best.period <- period
        }
        best.file <- sprintf("%s_%s_%s_%i_00_00.tif", covariate, measure, private$duration, best.period)
        best.path <- file.path(dir, best.file)
        rasters[[period.index]] <- raster::raster(best.path) # BREAKING THINGS
      }
      result <- raster::stack(rasters[1:period.index])
      names(result) <- rep(paste0(covariate, ".", 1:n_periods))
      return(result)
    },
    load_synoptic_covariate = function(dir, i) {
      covariate <- private$covariate_config[i, covariate]
      measure <- private$covariate_config[i, measure]
      path <- file.path(
        dir,
        paste(covariate, measure, "synoptic.tif", sep = "_")
      )

      if (!file.exists(path)) {
        err.msg <- paste("Searched for the following file and it does not exist:", path)
        stop(err.msg)
      }
      result <- raster::raster(path)
      names(result) <- covariate
      return(result)
    },
    get_periods = function(dir, expected_periods) {
      # background: files have very explicit filenames e.g., "cruststmn_median_5y_2000_00_00"
      # this is COVARIATE_MEASURE_DURATION_YEAR_MONTH_DAY

      # get files with duration (e.g., "5y") in name. Others should be ignored
      files <- list.files(dir)
      files <- files[grep(private$duration, files)]
      # strip extension, get unique values
      base_names <- unique(unlist(lapply(files, private$filename_without_extension)))
      # extract YEAR (third to last value) and convert to numeric
      periods <- as.numeric(sort(unlist(lapply(
        strsplit(base_names, split = "_", fixed = TRUE),
        function(pieces) {
          pieces[length(pieces) - 2]
        }
      ))))

      missing.periods <- setdiff(expected_periods, periods)
      if (length(missing.periods) > 0) {
        message("WARNING! You are trying to load a raster covariate but the following years are missing:")
        message(paste(missing.periods, collapse = " ", sep = ""))
        message("WARNING! We will map adjacent nearby years to these missing periods to fill in your dataset")
      }
      return(list(present = periods, missing = missing.periods))
    },
    get_closest_period = function(desired, available) {
      distance <- abs(desired - available)
      available[which.min(distance)]
    },
    filename_without_extension = function(f) {
      # returns filename without extension
      # unlike tools::file_path_sans_exit this will remove ALL extensions, not the first
      # e.g., foo.bar.baz becomes foo, not foo.bar
      strsplit(f, ".", fixed = TRUE)[[1]][1]
    },
    error_for_missing_covariates = function() {
      # TODO: should we build one big error message instead of a staged one?
      #       right now a user might have to run this 3 times to find all their errors:
      #       missing covariates, missing releases, missing measures

      # test for requested COVARIATES which do not exist
      cov.dirs <- private$path_helper$covariate_paths(covariates = private$covariate_config$covariate)
      if (!all(dir.exists(cov.dirs))) {
        missing.index <- which(!dir.exists(cov.dirs))
        covariates <- private$covariate_config[missing.index, covariate]
        msg <- paste(
          "You have selected some covariates in fixed_effects which do not exist:",
          paste0(covariates, collapse = ", ")
        )
        stop(msg)
      }

      # test for requested MEASURES which do not exist
      measure.dirs <- private$path_helper$covariate_paths(
        covariates = private$covariate_config$covariate,
        measures = private$covariate_config$measure
      )
      if (!all(dir.exists(measure.dirs))) {
        missing.index <- which(!dir.exists(measure.dirs))
        covariates <- private$covariate_config[missing.index, covariate]
        measures <- private$covariate_config[missing.index, measure]
        msg <- paste(
          "The following measures for covariates do not exist:",
          paste(measures, " (", covariates, ")", sep = "", collapse = "; ")
        )
        stop(msg)
      }

      # test for requested RELEASES which do not exist
      release.dirs <- private$path_helper$covariate_paths(
        covariates = private$covariate_config$covariate,
        measures = private$covariate_config$measure,
        releases = private$covariate_config$release
      )
      if (!all(dir.exists(release.dirs))) {
        missing.index <- which(!dir.exists(release.dirs))
        covariates <- private$covariate_config[missing.index, covariate]
        measures <- private$covariate_config[missing.index, measure]
        releases <- private$covariate_config[missing.index, release]
        msg <- paste(
          "The following releases for covariate / measure do not exist:",
          paste(releases, " (", covariates, " / ", measures, ")", sep = "", collapse = "; ")
        )
        stop(msg)
      }
    }
  ) # end private
)


#' @title Helper object for dealing with covariate paths
#' @description CovariatePathHelper turns covariate/measure/release data into paths via \code{covariate_paths()}
#'  and converts it back via \code{covariate_metadata_from_path}. This is used internally by
#'  \code{MbgStandardCovariateLoader}.
#' @examples
#' \dontrun{
#'
#' covariates <- c("access", "evi")
#' measures <- c("mean", "median")
#' releases <- c("2019_06_10", "2019_06_10")
#' helper <- CovariatePathHelper$new()
#' paths <- helper$covariate_paths(covariates = covariates,
#'                                 measures = measures,
#'                                 releases = releases)
#' metadata <- helper$covariate_metadata_from_path(paths[1])
#' }
#' @rdname CovariatePathHelper
#' @export
CovariatePathHelper <- R6::R6Class("CovariatePathHelper",
  public = list(
    cov_dir = "<FILEPATH>",

    covariate_paths = function(covariates, measures = NULL, releases = NULL) {
      if (is.null(measures)) {
        file.path(self$cov_dir, covariates)
      } else if (is.null(releases)) {
        file.path(self$cov_dir, covariates, measures)
      } else {
        file.path(self$cov_dir, covariates, measures, releases)
      }
    },
    covariate_metadata_from_path = function(path) {
      # remove cov_dir from string (will now begin with "/") and split on "/" (accounting for OS platform)
      pieces <- strsplit(sub(self$cov_dir, "", path), .Platform$file.sep, fixed = TRUE)[[1]]
      # first value is "", subsequent values are interesting
      result <- list(
        covariate = pieces[2],
        measure = pieces[3],
        release = pieces[4]
      )
      return(result)
    },
    newest_covariate_release = function(covariate_paths) {
      # USE.NAMES = FALSE causes the returned vector to only support numeric indexing
      return(sapply(covariate_paths, private$newest_covariate_release_single, USE.NAMES = FALSE))
    }
  ), # end public
  private = list(
    newest_covariate_release_single = function(path) {
      all.dirs <- list.dirs(path, full.names = FALSE, recursive = FALSE)
      release.dirs <- all.dirs[grep("^\\d{4}_\\d{2}_\\d{2}$", all.dirs)]
      return(max(release.dirs))
    }
  ) # end private
)

#' @title Load worldpop covariate raster and return
#' @description Loads a covariate raster (worldpop by default, can be overridden) and returns.
#' @param template_raster the raster template which all returned data will match.
#' @param covariate the covariate to load. Defaults to "worldpop"
#' @param pop_measure the covariate measure to load.
#' @param pop_release the covariate measure release to use.
#' @param start_year the first year to locate data for. Defaults the minimum value in the configuration value \code{year_list}
#' @param end_year the last year to locate data for. defaults the maximum value in the configuration value \code{year_list}
#' @param interval the number of months between data readings. Defaults to the global interval_mo
#' @examples
#' \dontrun{
#' worldpop <- load_worldpop_covariate(simple_polygon, measure = 'a0004t', release = '2019_06_10')
#'
#' raked_worldpop <- load_worldpop_covariate(simple_polygon,
#'                                           covariate = 'worldpop_raked',
#'                                           measure = 'a0004t',
#'                                           release = '2019_06_10',
#'                                           start_year = 2000,
#'                                           end_year = 2017,
#'                                           interval = 12)
#' }
#' @export
#' @rdname load_worldpop_covariate
#' @return list with named value containing your covariate data as a raster.
load_worldpop_covariate <- function(template_raster,
                                    covariate = "worldpop",
                                    pop_measure,
                                    pop_release,
                                    start_year = min(year_list),
                                    end_year = max(year_list),
                                    interval = interval_mo) {
  worldpop_config <- data.table(covariate = c(covariate),
                                measure = c(pop_measure),
                                release = c(pop_release))

  loader <- MbgStandardCovariateLoader$new(start_year = start_year,
                                           end_year = end_year,
                                           interval = interval,
                                           covariate_config = worldpop_config)
  return(loader$get_covariates(template_raster))
}

#' @title Reads covariate configuration file.
#' @description Read covariate configuration file providing sensible default values \code{header = TRUE, fill = TRUE}
#' @param path_or_text the file path or literal CSV text to read in.
#' @param ... any additional arguments to pass to \code{data.table::fread}
#' @return configuration as a data.table
#' @rdname read_covariate_config
#' @export
read_covariate_config <- function(path_or_text, ...) {
  data.table::fread(path_or_text,
    fill = TRUE, # fill blank fields in rows with uneven length
    header = TRUE, # first line is a header
    ...
  )
}


#' @title Update fixed effect covariate configuration with optional/missing data
#' @description Updates the fixed effect configuration with \emph{release} values for each
#' covariate/measure pair if it is not provided by using the most recent available release.
#' @param fixed_effect_config the fixed effect configuration, a data.table
#' @return NULL (the fixed_effect_config is modified in-place)
#' @rdname update_fixed_effect_config_with_missing_release
#' @export
update_fixed_effect_config_with_missing_release <- function(fixed_effect_config) {
  # add "release" column if not present.
  # fill it with the default fill value of data.table::fread
  if (!"release" %in% names(fixed_effect_config)) {
    # use explicit NA type so subsequent updating with character works
    # https://stackoverflow.com/a/15554528
    fixed_effect_config[, release := NA_character_]
    indices.to.update <- 1:nrow(fixed_effect_config)
  } else {
    indices.to.update <- which(fixed_effect_config$release == "")
  }
  if (length(indices.to.update) == 0) {
    return(NULL)
  }
  helper <- CovariatePathHelper$new()
  covariate.paths <- helper$covariate_paths(
    covariates = fixed_effect_config[indices.to.update, covariate],
    measures = fixed_effect_config[indices.to.update, measure])
  releases <- helper$newest_covariate_release(covariate.paths)
  fixed_effect_config[indices.to.update, release := releases]
  return(NULL)
}
