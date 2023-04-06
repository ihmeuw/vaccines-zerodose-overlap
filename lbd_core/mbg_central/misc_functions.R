


waitformodelstofinish <- function(sleeptime=100,
                                  path =  paste0('<FILEPATH>', indicator_group, '/', indicator, '/output/', run_date),
                                  rd   = run_date,
                                  lv   = loopvars,
                                  showfiles = TRUE,
                                  showcluster = FALSE){
  
  n_finished <- length(grep("fin_", list.files(path)))
  
  lv <- data.table(lv)
  names(lv) <- c("reg", "holdout")
  lv$reg <- as.character(lv$reg)
  
  lv[, file := paste0(path, "/", "fin__bin0_", reg, "_", holdout)]
  
  while(n_finished != nrow(lv)){
    
    n_finished <- length(grep("fin_", list.files(path)))
    
    message('\n====================================================================================')
    message(sprintf('=====================      Run Date: %s      ======================',rd))
    message(paste0('\nAt ',Sys.time(),' .... ',n_finished,' Models have written output.'))
    if(showfiles){
      message('\nCurrently missing models:')
      for(i in 1:nrow(lv))
        if(file.exists(lv[i, file]) == F)
          message(paste('Region =',lv[i,1],'| Holdout =',lv[i,2]))
    }
    n_cluster <- system("qstat | grep job_ | wc | awk '{print $1}'", intern = T)
    message(paste0('\nFuthermore, there are still ', n_cluster, ' jobs running on the cluster.'))
    if(showcluster){
      system("qstat -r | grep jobname | grep -oP \"(?<=job_)(.*)\"")
    }
    message('\n====================================================================================')
    message('====================================================================================')
    message("\n")
    Sys.sleep(sleeptime)
  }
  
  unlink(lv$file) # Clean up by deleting extra files once done with the loop
}



#' @title Combine aggregation
#' @description Combine aggregation objects across region
#' @param run_date, indicator, indicator_group for this run
#' @param ages: single value or vector of ages
#' @param regions: vector of regions used
#' @param holdouts: vector of holdouts used, e.g. just 0 or c(1,2,3,4,5,0)
#' @param raked: vector of raked values, e.g. just T, just F, or c(T,F)
#' @param dir_to_search: which directory to search in (defaults to share directory)
#' @param delete_region_files: logical. Should we delete the region-specific intermediate files?
#' @param merge_hierarchy_list: logical. Do you want to merge the sp_hierarchy_list onto your admin tables?
#' @param check_for_dupes PARAM_DESCRIPTION, Default: F
#' @return rdata files for each combo of age/holdout/raked
#'   each with admin_0, admin_1, admin_2 data table objects & the sp_hierarchy_list object
#'   that maps them to names of admin units
combine_aggregation <- function(rd = run_date,
                                indic = indicator,
                                ig = indicator_group,
                                ages,
                                regions,
                                holdouts,
                                raked,
                                dir_to_search = NULL,
                                delete_region_files = T,
                                merge_hierarchy_list = F,
                                check_for_dupes = F) {
  
  # Combine aggregation objects across region

  
  # Args:
  #   run_date, indicator, indicator_group for this run
  #   ages: single value or vector of ages
  #   regions: vector of regions used
  #   holdouts: vector of holdouts used, e.g. just 0 or c(1,2,3,4,5,0)
  #   raked: vector of raked values, e.g. just T, just F, or c(T,F)
  #   dir_to_search: which directory to search in (defaults to share directory)
  #   delete_region_files: logical. Should we delete the region-specific intermediate files?
  #   merge_hierarchy_list: logical. Do you want to merge the sp_hierarchy_list onto your admin tables?
  
  # Outputs:
  #   rdata files for each combo of age/holdout/raked
  #   each with admin_0, admin_1, admin_2 data table objects & the sp_hierarchy_list object
  #   that maps them to names of admin units
  
  if (is.null(dir_to_search)) {
    dir_to_search <- paste0("<FILEPATH>",ig,"/",indic,"/output/",run_date,"/")
  }
  
  message("Combining aggregation results...")
  
  for(rake in raked) {
    for (holdout in holdouts) {
      for (age in ages) {
        message(paste0("\nWorking on age: ", age, " | holdout: ", holdout, " | raked: ", rake))
        
        # Set up lists
        ad0 <- list()
        ad1 <- list()
        ad2 <- list()
        sp_h <- list()
        
        
        for (reg in regions) {
          message(paste0("  Region: ", reg))
          
          load(paste0(dir_to_search, indic, "_", ifelse(rake, "raked", "unraked"),
                      "_admin_draws_eb_bin", age, "_", reg, "_", holdout, ".RData"))
          
          if(merge_hierarchy_list == T) {
            # Prepare hierarchy list for adm0
            ad0_list <- subset(sp_hierarchy_list, select = c("ADM0_CODE", "ADM0_NAME", "region")) %>% unique
            
            # Prepare hierarchy list for adm1
            ad1_list <- subset(sp_hierarchy_list,
                               select = c("ADM0_CODE", "ADM1_CODE", "ADM0_NAME", "ADM1_NAME", "region")) %>%
              unique
            
            # Merge
            admin_0 <- merge(ad0_list, admin_0, by = "ADM0_CODE", all.y = T)
            admin_1 <- merge(ad1_list, admin_1, by = "ADM1_CODE", all.y = T)
            admin_2 <- merge(sp_hierarchy_list, admin_2, by = "ADM2_CODE", all.y = T)
            rm(ad0_list, ad1_list)
          }
          if(check_for_dupes){
            adms <- get_adm0_codes(reg)
            sp_hier <- get_sp_hierarchy()
            include_ad0 <- sp_hier$ADM0[ADM0_CODE %in% adms, ADM0_CODE]
            include_ad1 <- sp_hier$ADM1[ADM0_CODE %in% adms, ADM1_CODE]
            include_ad2 <- sp_hier$ADM2[ADM0_CODE %in% adms, ADM2_CODE]
            
            ad0[[reg]] <- admin_0[ADM0_CODE %in% include_ad0]
            ad1[[reg]] <- admin_1[ADM1_CODE %in% include_ad1]
            ad2[[reg]] <- admin_2[ADM2_CODE %in% include_ad2]
            sp_h[[reg]] <- sp_hierarchy_list
          } else{
            ad0[[reg]] <- admin_0
            ad1[[reg]] <- admin_1
            ad2[[reg]] <- admin_2
            sp_h[[reg]] <- sp_hierarchy_list
          }
          
          rm(admin_0, admin_1, admin_2, sp_hierarchy_list)
        }
        
        # Get to long format & save
        message("  Combining...")
        admin_0 <- rbindlist(ad0)
        admin_1 <- rbindlist(ad1)
        admin_2 <- rbindlist(ad2)
        sp_hierarchy_list <- rbindlist(sp_h)
        
        message("  Saving combined file...")
        save(admin_0, admin_1, admin_2, sp_hierarchy_list,
             file = paste0(dir_to_search, indic, "_",
                           ifelse(rake, "raked", "unraked"),
                           "_admin_draws_eb_bin", age, "_",
                           holdout, ".RData"))
      }
    }
  }
  
  if (delete_region_files == T) {
    # Make sure all full files are written
    combos <- expand.grid(ifelse(raked, "raked", "unraked"), ages, holdouts)
    files_to_check <- sapply(1:nrow(combos), function(i) {
      paste0(dir_to_search, indic, "_", combos[i, 1],
             "_admin_draws_eb_bin", combos[i, 2], "_", combos[i,3], ".RData")
    })
    
    if (all(file.exists(files_to_check))) {
      message("All anticipated combined files were created successfully.  Deleting intermediate files...")
      combos <- expand.grid(ifelse(raked, "raked", "unraked"), ages, regions, holdouts)
      files_to_delete <- sapply(1:nrow(combos), function(i) {
        paste0(dir_to_search, indic, "_", combos[i, 1],
               "_admin_draws_eb_bin", combos[i, 2], "_", combos[i,3], "_", combos[i, 4], ".RData")
      })
      unlink(files_to_delete)
    } else {
      warning("Did not delete intermediate files - not all output files created successfully!")
    }
  }
  
  # Finally, delete the "fin" files
  fin_files_to_delete <- list.files(dir_to_search, pattern = "fin_agg_", full.names=T)
  unlink(fin_files_to_delete)
}


#' @title Get singularity
#'
#' @description \code{get_singularity} determines which Singularity image to use. The
#' default is the 'default' keyword. Image names without the full path to the
#' file defined are assumed to exist as the default location:
#'    \code{<FILEPATH>}
#' In either case it will test to make sure the file exists and exit if it does
#' not. The default image is hardcoded into the shell script used to launch
#' Singularity containers:
#'   \code{lbd_core/mbg_central/share_scripts/shell_sing.sh}
#'
#' @param image A string that defines which Singularity image to launch
#'   [default = 'default']. If the 'default' keyword is passed in or left blank,
#'   the default keyword will be returned. Either the full path to the image
#'   may be provided or only the Singularity image name. In the latter case,
#'   the image is assumed to live in the default image location:
#'   \code{<FILEPATH>}
#'
#' @return When image = 'default', 'default' is returned. When this keyword is
#'   passed to the shell_sing.sh script through `qsub` it will use the default
#'   Singularity image hardcoded into it. Otherwise, if the function is
#'   successful at verifying that the Singularity image file specified exists,
#'   it will return the full path to that image.
#'
#' @seealso This function is used by:
#'   \code{\link{parallelize}}
#'   \code{\link{make_qsub}}
#'   \code{\link{make_qsub_share}}
#'   \code{\link{make_qsub_postest}}
#'   \code{\link{submit_aggregation_script}}
#'
get_singularity <- function(image = 'default') {
  if(image == 'default') {        # use default image
    sing_image <- 'default'
  } else if(grepl('/', image)) {  # user supplied path to image
    sing_image <- image
  } else {                        # image at default location
    sing_image <- paste0('<FILEPATH>/', image)
  }
  # If something other than the default image is being used, let's make sure
  # the image file actually exists:
  if(!sing_image == 'default' & !file.exists(sing_image)) {
    stop(paste0("Could not locate Singularity image: ", sing_image))
  }
  return(sing_image)
}



insertRaster <- function (raster, new_vals, idx = NULL) {
  
  
  # calculate cell index if not provided
  if (is.null(idx)) idx <- cellIdx(raster)
  
  # check the index makes superficial sense
  stopifnot(length(idx) == nrow(new_vals))
  stopifnot(max(idx) <= ncell(raster))
  
  # create results raster
  n <- ncol(new_vals)
  raster_new <- raster::brick(replicate(n,
                                        raster[[1]],
                                        simplify = FALSE))
  names(raster_new) <- colnames(new_vals)
  
  # update the values
  for(i in 1:n) {
    raster_new[[i]][idx] <- new_vals[, i]
  }
  
  return (raster_new)
  
}
condSim <- function (vals, weights = NULL, group = NULL, fun = NULL, ...) {
  # given a matrix of pixel-level prevalence samples `prev`
  # where each rows are pixels and columns are draws, a vector
  # of corresponding pixel populations `pop`, and an optional pixel
  # grouping factor `group`, return draws for the total deaths in each
  # group, or overall if groups are not specified
  
  # get dimensions of vals
  ncell <- nrow(vals)
  ndraw <- ncol(vals)
  
  # capture function as a string
  fun_string <- deparse(substitute(fun))
  
  # check fun accepts a
  
  # check dimensions of weights and group, set to 1 if not specified
  if (is.null(weights)) {
    weights <- rep(1, ncell)
  } else {
    if (length(weights) != ncell) {
      stop (sprintf('number of elements in weights (%i) not equal to number of cells in vals (%i)',
                    length(weights),
                    ncell))
    }
  }
  
  if (is.null(group)) {
    group <- rep(1, length(weights))
  } else {
    if (length(group) != ncell) {
      stop (sprintf('number of elements in group (%i) not equal to number of cells in vals (%i)',
                    length(group),
                    ncell))
    }
  }
  
  # otherwise, get the levels in group and create a matrix of results
  levels <- unique(na.omit(group))
  nlevel <- length(levels)
  
  ans <- matrix(NA,
                ncol = ndraw,
                nrow = nlevel)
  rownames(ans) <- levels
  
  # loop through levels in group, getting the results
  for (lvl in 1:nlevel) {
    
    # get an index o pixels in the level
    idx <- which(group == levels[lvl])
    
    # by default, calculate a weighted sum
    if (is.null(fun)) {
      
      # get draws and add to results
      # exception for if area has 1 cell, transpose matrix so it conforms (RB)
      if(all(dim(t(vals[idx, ]))==c(1,ndraw))){
        ans[lvl, ] <- weights[idx] %*% t(vals[idx, ])
      } else {
        ans[lvl, ] <- weights[idx] %*% vals[idx, ]
      }
      
    } else {
      
      # otherwise, apply function to each column
      ans[lvl, ] <- apply(vals[idx, ], 2, fun, weights = weights[idx], ...)
      
    }
    
  }
  
  # if only one level, make this a vector
  if (nlevel == 1) ans <- as.vector(ans)
  
  # return result
  return (ans)
  
}



# given the admin level, a GAUL code and the admin1 and 2 shapefiles,
# return an SPDF with the relevant area
getPoly <- function (level, code, admin1, admin2) {
  
  # get admin level
  if (level == 1) {
    admin <- admin1
  } else {
    admin <- admin2
  }
  
  # find the reight area
  idx <- match(code, admin$GAUL_CODE)
  
  # if it's valid
  if (length(idx) == 1) {
    ans <- admin[idx, ]
  } else {
    # handle errors
    warning (paste0("something's wrong on row ", i))
    ans <- NULL
  }
  
  # return result
  return (ans)
}





getCols <- function (df, period = 1) {
  # subset results matrix
  df[, grep(sprintf('^p%s_', period),
            colnames(df))]
}

rnormMatrix <- function (n, mean, sd) {
  # sample random normals with matrix parameters
  # returna an array as a result
  
  # coerce to matrix
  mean <- as.matrix(mean)
  sd <- as.matrix(sd)
  
  # get & check dimensions
  ncol <- ncol(mean)
  nrow <- nrow(mean)
  stopifnot(ncol(sd) == ncol)
  stopifnot(nrow(sd) == nrow)
  
  # convert to vector
  mean <- as.vector(mean)
  sd <- as.vector(sd)
  
  # sample
  draws <- rnorm(n * length(mean), mean, sd)
  
  # reshape
  ans <- array(draws, dim = c(nrow, ncol, n))
  
  return (ans)
}


# functions to centre and scale matrices, columnwise

getCentreScale <- function (x, exclude = NULL, na.rm = TRUE) {
  # get dataframe of centreing and scaling values to convert x
  # to the standard normal. exclude is an optional character vector
  # giving column names to exclude from scaling
  
  # get means and SDs for all columns
  df <- data.frame(name = colnames(x),
                   mean = colMeans(x, na.rm = na.rm),
                   sd = apply(x, 2, sd, na.rm = na.rm))
  rownames(df) <- NULL
  
  # replace any zero standard deviations with 1
  # to avoid divide-by-zero errors
  df$sd[df$sd == 0] <- 1
  
  # if any named covariates are to be excluded, set mean to 0 and sd to 1
  if (!is.null(exclude)) {
    idx <- match(exclude, df$name)
    df$mean[idx] <- 0
    df$sd[idx] <- 1
  }
  
  return (df)
}

centreScale <- function (x, df, inverse = FALSE) {
  # apply pre-calculated centreing/scaling to matrix x,
  # with fixed dataframe of means/sds df
  # or uncentre/unscale if inverse = TRUE
  
  # get the centreing/scaling dataframe if not available
  if (is.null(df))
    df <- getCentreScale(x)
  
  # get index to match up values with column names
  names <- colnames(x)
  idx <- match(names, df$name)
  
  if (any(is.na(idx))) {
    stop ('could not match up column names with the values in df')
  }
  
  df <- df[idx, ]
  
  # apply transformations
  if (!inverse) {
    # move to standard normal
    
    # centre
    x <- sweep(x, MARGIN = 2, STATS = df$mean, FUN = '-')
    # scale
    x <- sweep(x, MARGIN = 2, STATS = df$sd, FUN = '/')
    
  } else {
    # inverse case, move from standard normal to original
    
    # unscale
    x <- sweep(x, MARGIN = 2, STATS = df$sd, FUN = '*')
    # uncentre
    x <- sweep(x, MARGIN = 2, STATS = df$mean, FUN = '+')
    
  }
  
  return (x)
  
}



elogit <- function (y, n) log ( (y + 0.5) / (n - y + 0.5) )



# These 3 functions are used as a summstats argument via a configuration option. DO NOT DELETE
cirange = function(x){
  z=quantile(x,probs=c(.025,.975),na.rm=T)
  return(z[2]-z[1])
}
lower = function(x) quantile(x,probs=.025,na.rm=T)
upper = function(x) quantile(x,probs=.975,na.rm=T)






# Timer Functions ---------------------------------------------------------

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


graph_run_summary <- function(run_date,
                              indicator,
                              indicator_group,
                              return_graph = F) {
  
  # Function to graph run_summary files
  # Requires creation of a run_summary .csv file in your run_date directory
  
  if (Sys.info()["sysname"] == "Linux") {
    j_root <- "<FILEPATH>"
    #package_lib <- paste0(j_root,'/temp/geospatial/packages') # Library for all MBG versioned packages.
    #.libPaths(package_lib)
  }
  require(data.table)
  require(magrittr)
  require(ggplot2)
  require(RColorBrewer)
  dir <- paste0("<FILEPATH>", indicator_group, "/", indicator, "/output/", run_date, "/")
  
  file <- list.files(dir, pattern = "run_summary.*.csv")
  
  # Catch if file does not exist
  if (length(file) == 0) {
    message("No run summary file found to graph... exiting function.")
    return(NULL)
  }
  
  # Otherwise continue and create file name
  file <- paste0(dir, file)
  
  df <- read.csv(file, stringsAsFactors = F) %>% as.data.table
  
  grab_time_hours <- function(x) {
    
    v_time <- unlist(strsplit(x, " "))
    hours <- v_time[1]
    hours <- substr(hours, 0, nchar(hours) - 1) %>% as.numeric
    
    minutes <- v_time[2]
    minutes <- substr(minutes, 0, nchar(minutes) - 1) %>% as.numeric
    
    seconds <- v_time[1]
    seconds <- substr(seconds, 0, nchar(seconds) - 1) %>% as.numeric
    
    hours <- hours + minutes/60 + seconds/(60*60)
    return(hours)
    
  }
  
  df$time <- sapply(df$time, grab_time_hours)
  df$step <- factor(df$step, levels = c("Stacking - GAM", "Stacking - GBM", "Stacking - lasso", "Stacking - ridge", "Stacking - enet",
                                        "MBG - fit model", "MBG - predict model", "Cross-validation",
                                        "Stacking - all", "MBG - all", "Entire script"))
  
  summary_steps <- c("Stacking - all", "MBG - all", "Entire script")
  
  df_graph <- subset(df, !(step %in% summary_steps))
  df_graph[, holdout := as.character(holdout != 0)]
  
  g_plot <- ggplot(data = df_graph, aes(x = step, y = time, color = region)) +
    stat_summary(fun.y = mean, geom = "line", aes(group = region, color = region)) +
    geom_point(aes(color=region, shape=holdout)) +
    scale_color_brewer(palette = "Set1") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Step",
         y = "Time (hours)",
         title = paste0("Run Date: ", run_date),
         color = "Region",
         shape = "Holdout?")
  
  png(filename = paste0(dir, "run_summary_", indicator, "_", run_date, ".png"),
      type = "cairo",
      units = "in",
      width = 8,
      height = 4.5,
      pointsize = 12,
      res = 300)
  
  print(g_plot)
  
  dev.off()
  
  if (return_graph == T) return(g_plot)
  
}





check_config <- function(cr = core_repo) {
  
  # TODO: update package to use data()
  # data(must_haves)
  must_haves <- read.csv(paste0(cr, '/mbg_central/share_scripts/common_inputs/config_must_haves.csv'), header = F, stringsAsFactors = F)$V1
  
  message("\nRequired covariates: ")
  for(confs in must_haves){
    if (exists(confs)) {
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == 'use_gp') {
      message("You are missing a 'use_gp' argument in your config. Defaulting it to TRUE")
      use_gp <<- TRUE
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == 'use_stacking_covs') {
      message("You are missing a 'use_stacking_covs' argument in your config. Defaulting it to TRUE")
      use_stacking_covs <<- TRUE
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == 'use_raw_covs') {
      message("You are missing a 'use_raw_covs' argument in your config. Defaulting it to FALSE")
      use_raw_covs <<- FALSE
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == 'fit_with_tmb') {
      message("You are missing a 'fit_with_tmb' argument in your config. Defaulting it to FALSE")
      fit_with_tmb <<- FALSE
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == 'gbd_fixed_effects_measures') {
      message("You are missing a 'gbd_fixed_effects_measures' argument in your config. Defaulting it to 'covariate' for all elements of gbd_fixed_effects")
      gbd_fixed_effects_measures <<- paste(rep("covariate", length(strsplit(gbd_fixed_effects, split = " \\+ ")[[1]])), collapse = " + ")
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == 'gbd_fixed_effects_age') {
      message("You are missing a 'gbd_fixed_effects_age' argument in your config. Defaulting to '2 3 4 5'")
      gbd_fixed_effects_age <<- '2 3 4 5'
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == 'z_list') {
      message("You are missing a 'z_list' argument in your config. Defaulting it to 0")
      z_list <<- 0
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == 'zcol') {
      message("You are missing a 'zcol' argument in your config. Defaulting it to z_column_default_blank")
      zcol <<- 'z_column_default_blank'
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == 'summstats') {
      message("You are missing a 'summstats' argument in your config. Defaulting to c('mean','lower','upper','cirange')")
      summstats <<- c('mean','lower','upper','cirange')
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == 'scale_gaussian_variance_N') {
      message("You are missing a 'scale_gaussian_variance_N' argument in your config. Defaulting to TRUE")
      scale_gaussian_variance_N <<- TRUE
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == "nugget_prior") {
      message("You are missing a 'nugget_prior' argument in your config. Defaulting to 'list(prior = 'loggamma', param = c(2, 1))'")
      nugget_prior <<- "list(prior = 'loggamma', param = c(2, 1))"
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == "ctry_re_prior") {
      message("You are missing a 'ctry_re_prior' argument in your config. Defaulting to 'list(prior = 'loggamma', param = c(2, 1))'")
      ctry_re_prior <<- "list(prior = 'loggamma', param = c(2, 1))"
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == "nid_re_prior") {
      message("You are missing a 'nid_re_prior' argument in your config. Defaulting to 'list(prior = 'loggamma', param = c(2, 1))'")
      nid_re_prior <<- "list(prior = 'loggamma', param = c(2, 1))"
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == "use_nid_res") {
      message("You are missing a 'use_nid_res' argument in your config. Defaulting to FALSE")
      use_nid_res <<- FALSE
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == "rho_prior") {
      message("You are missing a 'rho_prior' argument in your config. Defaulting to 'list(prior = 'normal', param = c(0, 0.1502314))'")
      rho_prior <<- "list(prior = 'normal', param = c(0, 1/(2.58^2)))"
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == "gp_sum_to_zero") {
      message("You are missing a 'gp_sum_to_zero' argument in your config. Defaulting to FALSE")
      gp_sum_to_zero <<- FLASE
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == "use_s2_mesh") {
      message("You are missing a 'use_s2_mesh' argument in your config. Defaulting to FALSE")
      use_s2_mesh <<- FALSE
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == "s2_mesh_params") {
      message("You are missing a 's2_mesh_params' argument in your config. Defaulting to c(50, 500, 1000)")
      s2_mesh_params <<- "c(25, 500, 1000)"
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == "sparse_ordering") {
      message("You are missing a 'sparse_ordering' argument in your config. Defaulting to TRUE")
      sparse_ordering <<- TRUE
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == "modeling_shapefile_version") {
      message("You are missing a 'modeling_shapefile_version' argument in your config. Defaulting to 'current'")
      modeling_shapefile_version <<- "current"
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == "raking_shapefile_version") {
      message("You are missing a 'raking_shapefile_version' argument in your config. Defaulting to 'current'")
      raking_shapefile_version <<- "current"
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == "subnational_raking") {
      message("You are missing a 'subnational_raking' argument in your config. Defaulting to TRUE")
      subnational_raking <<- TRUE
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == "check_cov_pixelcount") {
      message("You are missing a 'check_cov_pixelcount' argument in your config. Defaulting to FALSE")
      check_cov_pixelcount <<- FALSE
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == "gbd_fixed_effects_constraints") {
      message("You are missing a 'gbd_fixed_effects_constraints' argument in your config. Defaulting to FALSE")
      gbd_fixed_effects_constraints <<- "c(0)"
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == "fixed_effects_constraints") {
      message("You are missing a 'fixed_effects_constraints' argument in your config. Defaulting to FALSE")
      fixed_effects_constraints <<- "c(0)"
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == "memory") {
      message("You are missing a 'memory' argument in your config. Defaulting to 10G")
      memory <<- 10
      
    } else if (confs == "singularity_version") {
      message("You are missing a 'singularity_version' argument in your config. Defaulting to 'default'")
      singularity_version <<- "default"
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == "queue") {
      message("You are missing a 'queue' argument in your config. Defaulting to 'long.q', unless you have use_geos_nodes to TRUE, which will override this to geospatial.q")
      queue <<- "long.q"
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == "run_time") {
      message("You are missing a 'run_time' argument in your config. Defaulting to 16 days ('16:00:00:00')")
      run_time <<- "16:00:00:00"
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == "countries_not_to_rake") {
      message("You are missing a 'countries_not_to_rake' argument in your config. Defaulting to ESH+GUF")
      countries_not_to_rake <<- "ESH+GUF"
      message(paste0("  ", confs, ": ", get(confs)))
      
    } else if (confs == "countries_not_to_subnat_rake") {
      message("You are missing a 'countries_not_to_subnat_rake' argument in your config. Defaulting to PHL+NGA+PAK+ETH+KEN")
      countries_not_to_subnat_rake <<- "PHL+NGA+PAK+ETH+KEN"
      message(paste0("  ", confs, ": ", get(confs)))
      
    } else if (confs == "rake_countries") {
      message("You are missing a 'rake_countries' argument in your config. Defaulting to TRUE")
      rake_countries <<- TRUE
      message(paste0("  ", confs, ": ", get(confs)))
      
    } else if (confs == "use_space_only_gp") {
      message("You are missing a 'use_space_only_gp' argument in your config. Defaulting to FALSE")
      use_space_only_gp <<- FALSE
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == "st_gp_int_zero") {
      message("You are missing a 'st_gp_int_zero' argument in your config. Defaulting to FALSE")
      st_gp_int_zero <<- FALSE
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == "s_gp_int_zero") {
      message("You are missing a 's_gp_int_zero' argument in your config. Defaulting to FALSE")
      s_gp_int_zero <<- FALSE
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == "use_time_only_gmrf") {
      message("You are missing a 'use_time_only_gmrf' argument in your config. Defaulting to FALSE")
      use_time_only_gmrf <<- FALSE
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == "time_only_gmrf_type") {
      message("You are missing a 'time_only_gmrf_type' argument in your config. Defaulting to FALSE")
      time_only_gmrf_type <<- "rw2"
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else if (confs == "spde_prior") {
      message("You are missing a 'spde_prior' argument in your config. Defaulting to 'list(type='pc')'")
      spde_prior <<- "list(type='pc')"
      message(paste0('  ', confs, ': ', get(confs)))
      
    } else {
      stop(paste0(confs, " is missing, add it to your config"))
    }
  }
  
  
  ## Test for subnational random effect
  if(exists("use_subnat_res", envir = .GlobalEnv)) {
    stopifnot(exists("subnat_country_to_get", envir = .GlobalEnv))
    # stopifnot(length(eval(parse(text = subnat_country_to_get))) == 1)
  } else {
    use_subnat_res <<- FALSE
    subnat_country_to_get <<- FALSE
  }
  
  
  message("\nAdditional config arguments: ")
  extras <- config$V1[!(config$V1 %in% must_haves)]
  for (extra in extras) message(paste0('  ', extra, ': ', get(extra)))
  
  ## print out shapefile info
  m.sf.info <- detect_adm_shapefile_date_type(shpfile_path = get_admin_shapefile(version = modeling_shapefile_version))
  r.sf.info <- detect_adm_shapefile_date_type(shpfile_path = get_admin_shapefile(version = raking_shapefile_version))
  message("\n\n\nSHAPEFILE VERSION INFORMATION: ")
  message(sprintf("\n--MODELING SHAPEFILE VERSION: %s -- which contains %s codes", m.sf.info$shpfile_date, toupper(m.sf.info$shpfile_type)))
  message(sprintf("\n--RAKING SHAPEFILE VERSION:   %s -- which contains %s codes\n", r.sf.info$shpfile_date, toupper(r.sf.info$shpfile_type)))
}


#' @title Easy eval-parse
#' @description Allows for easily eval-parsing through a config dataset
#' @param data The data.table 
#' @param column The column with the string call
#' @return Evaluated call
#' @export
#' @rdname ez_evparse
ez_evparse <- function(data, column) {
  return(eval(parse(text = data[, column, with = FALSE])))
}






#' @title Set up config
#' @description Setting up configuration variables for an MBG run
#' @param repo Location where you've cloned the MBG repository for your indicator.
#' @param core_repo Location where you've cloned the lbd_core repository. Not necessary in the package version.
#' @param indicator_group Category of indicator, i.e. "education"
#' @param indicator Specific outcome to be modeled within indicator category, i.e. "edu_0"
#' @param config_name Name of configuration file in the indicator folder, Default: NULL
#' @param config_file Full path to configuration file that overrides \code{config_name}, Default: NULL
#' @param covs_name Name of covariates configuration file, Default: NULL
#' @param covs_file Full path to covariates configuration file that overrides \code{covs_name}, Default: NULL
#' @param post_est_only Set up only for post estimation? Default: FALSE
#' @param run_date Run date, Default: ''
#' @param push_to_global_env Should the config parameters be pushed to the global environment? Default: TRUE
#' @param run_tests Run the assertion tests? This will run the tests and error out if there's an 
#' inconsistent config parameter. Default: TRUE
#' @param return_list Return a list result or just the config? Default FALSE
#' @return Depends on return_list. If FALSE (default) returns just the MBG config (a list). If True, returns a
#' named list of configs, where "config" is the usual MBG config, and "fixed_effects_config" is the config info
#' of the fixed effects
#' @details DETAILS
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   config <- load_config(repo = core_repo,
#'     indicator_group = indicator_group,
#'     indicator = indicator,
#'     config_name = 'config_training',
#'     covs_name = 'covs_training')
#' }
#' }
#' @rdname set_up_config
#' @importFrom assertthat is.flag is.string is.number
#' @export
set_up_config <- function(repo, 
                          core_repo = repo,
                          indicator_group, 
                          indicator, 
                          config_name = NULL, 
                          config_file = NULL, 
                          covs_name = NULL, 
                          covs_file = NULL,
                          post_est_only = FALSE, 
                          run_date = "", 
                          push_to_global_env = TRUE,
                          run_tests = TRUE,
                          return_list = FALSE) {
  
  ###### Block 1: Equivalent to load_config ######
  
  print("[1/6] Load the configs")
  
  ####### Logic checking for model config ####### 
  ## Make sure only one of config_name or config_file are not null
  if (!is.null(config_name) & !is.null(config_file)) {
    stop("You must specify just one of config_name or config_file, not both", call. = FALSE)
  }
  
  ## Pull config from indicator repo
  if (is.null(config_name) & is.null(config_file)) {
    message("Pulling config from default folder, since config_name and config_file are NULL")
    ## If new model run, pull config from /share repo
    if (post_est_only == FALSE) 
      config <- data.table::fread(paste0(repo, "/",  "config_", indicator, ".csv"), header = FALSE)
    ## If running analysis on existing model, use config from that model's outputs folder
    if (post_est_only == TRUE) 
      config <- data.table::fread(paste0("<FILEPATH>", indicator_group, "/", indicator, "/output/", run_date, "/config.csv"))
  }
  
  ## Pull by specific config name
  if (!is.null(config_name) & is.null(config_file)) {
    message("Pulling config from specified name")
    config <- data.table::fread(paste0(repo, "/", config_name, ".csv"), header = FALSE)
  }
  ## Pull specified config file
  if (is.null(config_name) & !is.null(config_file)) {
    message("Pulling config from specified filepath")
    config <- data.table::fread(config_file, header = FALSE)
  }
  
  ####### Logic checking for covariates config ####### 
  ## Make sure only one of covs_name or covs_file are not null
  if (!is.null(covs_name) & !is.null(covs_file)) {
    stop("You must specify just one of covs_name or covs_file, not both", call. = FALSE)
  }
  
  ## Covs not pulled
  if (is.null(covs_name) & is.null(covs_file)) {
    message("Not pulling covs since covs_name and covs_file are NULL")
    covs <- NULL
  }
  
  ## Pull by specific covs name
  if (!is.null(covs_name) & is.null(covs_file)) {
    message("Pulling covs from specified name")
    covs <- read_covariate_config(paste0(repo, "/covariates/", covs_name, ".csv"))
  }
  
  ## Pull specified covs file
  if (is.null(covs_name) & !is.null(covs_file)) {
    message("Pulling covs from specified filepath")
    covs <- read_covariate_config(covs_file)
  }
  
  ## For parsimony, let's make sure that the config column names are V1 and V2
  config <- data.table(config)
  if(colnames(config)[1] != "V1" & colnames(config)[2] != "V2") {
    warning("Renaming config column names to V1 and V2. Please verify that 'config' is properly built")
    colnames(config) <- c("V1", "V2")
  }
  
  
  # If a covariate .csv file exists, use that instead
  if (!is.null(covs)) {
    
    # Grab fixed effects & measures (and gbd fixed effects & measures) from CSV if present
    
    # After update to data.table 1.11.4, 'T' and 'F' are not read in as logical, 
    ## but as characters, which we need to remedy here. 
    ## We are assuming that the 'covs.csv' has an 'include' and 'gbd' column here
    covs[, `:=`(gbd, as.logical(gbd))]
    covs[, `:=`(include, as.logical(include))]
    covs <- subset(covs, include == T)  # Use only those where include flag set to true
    fe <- subset(covs, gbd == F)
    update_fixed_effect_config_with_missing_release(fe)
    gbd <- subset(covs, gbd == T)
    gbd[measure != "output", `:=`(measure, "covariate")]  # FIXME: This is a hack for backwards compatability -- basically it assumes you meant 'covariate' if you specified anything other than 'outcome' (eg, mean or NA)
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
    
    # Remove any other versions from original config and 
    # override with covariates config outputs
    all_varz <- c(
      "fixed_effects", "fixed_effects_measures", "fixed_effects_constraints",
      "gbd_fixed_effects", "gbd_fixed_effects_measures", "gbd_fixed_effects_constraints"
    )
    for(varz in all_varz) {
      if(!(varz %in% colnames(config))) {
        config <- rbindlist(list(config, data.table(V1 = varz, V2 = get(varz))))
      } else {
        config[V1 == varz, V2:= get(varz)]
      }
    }
    
  }
  
  
  ###### Block 2: Add fields in config that are not in the default set ######
  
  print("[2/6] Add fields that are in the default config set but not in user's config")
  
  ## Load in the default config dataset
  if (.in.package()) {
    data("default_config_values", package = packageName())
  } else {
    default_config_values <- data.table::fread(file.path('<FILEPATH>/lbd_core/mbg_central/share_scripts/common_inputs/config_values.csv'), header = TRUE, stringsAsFactors = FALSE)
  }
  
  ## Now, go through each of the values in `config_values` and 
  ## add on all the fields that are not in the user-specified config
  config <- set_default_config_values(config, default_config_values)
  
  
  ###### Block 3: Extra parameters in config ######
  
  print("[3/6] Add fields that are in user's config but not in the default config set")
  message("\nAdditional covariates: ")
  extras <- config$V1[!(config$V1 %in% names(default_config_values))]
  for (extra in extras) {
    message(paste0("  ", extra, ": ", config[V1 == extra, V2] ))
  }
  
  ###### Block 4: Print out shapefile info from config. Resolve 'current' to fixed version date ######
  
  print("[4/6] Print out shapefile info from config")
  
  ## get the shapefile info
  m.sf.info <- detect_adm_shapefile_date_type(shpfile_path = get_admin_shapefile(version = config[V1 == 'modeling_shapefile_version', V2]))
  r.sf.info <- detect_adm_shapefile_date_type(shpfile_path = get_admin_shapefile(version = config[V1 == 'raking_shapefile_version', V2]))
  
  ## replace shapefile version in config (and env variable) with the actual version date
  ## if a specific date was already set, nothing changes.
  ## if 'current' had been selected, then it will be replaced by the version date that is currently symlinked to 'current'
  config[V1 == 'modeling_shapefile_version', V2 := m.sf.info$shpfile_date]
  config[V1 == 'raking_shapefile_version', V2 := r.sf.info$shpfile_date]
  
  ## print out the shapefile info
  message("\n\n\nSHAPEFILE VERSION INFORMATION: ")
  message(sprintf("\n--MODELING SHAPEFILE VERSION: %s -- which contains %s codes", m.sf.info$shpfile_date, toupper(m.sf.info$shpfile_type)))
  message(sprintf("\n--RAKING SHAPEFILE VERSION:   %s -- which contains %s codes\n", r.sf.info$shpfile_date, toupper(r.sf.info$shpfile_type)))
  
  
  ###### Block 5: Run tests on all the configuration variables loaded ######
  if(run_tests) {
    print("[5/6] Running simple type-assertion tests on config parameters")
    if (.in.package()) {
      data("config_tests", package = packageName())
    } else {
      config_tests <- data.table::fread(paste0('<FILEPATH>/mbg_central/share_scripts/common_inputs/config_tests.csv'), header = TRUE, stringsAsFactors = FALSE)
    }
    
    ## Test for params only in the config_tests list of params
    for (param in sort(config[, V1])) {
      cat(paste0("Testing config parameter: ", param, " "))
      if(param %in% config_tests$variable) {
        test_call_1 <- config_tests[variable == param, test_call]
        test_call_2 <- config_tests[variable == param, extra_test1]
        test_call_3 <- config_tests[variable == param, extra_test2]
        
        if(test_call_1 != "") {
          ## For a string in the config file, the eval-parse combo will
          ## fail to evaluate it, and so we build in this exception for that
          tryCatch(
            get(test_call_1)(ez_evparse(config[V1 == param, ], "V2")),
            error = function(e) {
              if(attributes(e)$class[[1]] == 'simpleError') {
                if (test_call_1 == "is.string") {
                  message(paste0("Assertion on ", param, " errored out because it's tested as a string. Please check for the real type manually"))
                } else {
                  stop(sprintf("%s errored with message: %s", test_call_1, geterrmessage()))
                }
              }
            }
          )
        }
        if(test_call_2 != ""  ) {
          tryCatch(
            assertthat::assert_that(eval(parse(text = test_call_2))),
            error = function(e) {
              stop(paste0("The following test failed: ", test_call_2) )
            }
          )
        }
        if(test_call_3 != ""  ) {
          tryCatch(
            assertthat::assert_that(eval(parse(text = test_call_3))),
            error = function(e) {
              stop(paste0("The following test failed: ", test_call_3) )
            }
          )
        }
        cat(" OK. \n")
      }
    }
    
    ## Stop if using z or poly aggregation strategies without TMB
    if(as.logical(config[V1 == "poly_ag", "V2"]) | config[V1 == "zcol_ag", "V2"] != "NULL") {
      if(!as.logical(config[V1 == "fit_with_tmb", "V2"])) {
        stop("Must use TMB when incorporating polygon or aggregated z data")
      }
      if(as.logical(config[V1 == "makeholdouts", "V2"])) {
        stop("There is aggregated data and functionality for making holdouts is 
             not yet implemented. Set makeholdouts to FALSE.")
      }
      if(as.logical(config[V1 == "test", "V2"])) {
        stop("Testing with aggregated data not yet implemented. Set test to FALSE.")
      }
    }
    
  } else {
    warning("[5/6] Skipping over type-assertion")
  }
  
  
  
  
  ###### Final Block :  Assign all the covariates to the environment if desired ######
  
  if(push_to_global_env) {
    print("[6/6] Pushing config parameters into global environment")
    for (param in config[, V1]) { 
      assign(param, config[V1 == param, V2], envir = globalenv())
    }
    if (!is.null(covs)) {
      assign("fixed_effects_config", fe, envir = globalenv())
      assign("gbd_fixed_effects_config", gbd, envir = globalenv())
    }
    # Processing of z config arguments
    if (zcol_ag != "NULL") {
      if (exists("z_ag_mat_file")) {
        assign("z_ag_mat", read.csv(z_ag_mat_file, header=F), envir = globalenv())
      } else {
        assign("z_ag_mat", NULL, envir = globalenv())
        assign("zcol_ag_id", NULL, envir = globalenv())
      }
      if (exists("z_map_file")) {
        assign("z_map", read_z_mapping(z_map_file), envir = globalenv())
      } else {
        assign("z_map", NULL, envir = globalenv())
      }
    } else {
      assign("zcol_ag", NULL, envir = globalenv())
      assign("z_map", NULL, envir = globalenv())
      assign("z_ag_mat", NULL, envir = globalenv())
      assign("zcol_ag_id", NULL, envir = globalenv())
    }
  } else {
    print("[6/6] Config parameters not passed into global environment")
  }
  
  ## Return the config data.table
  message("Saving out config...")
  if (return_list) {
    return(list("config" = config, "fixed_effects_config" = fe))
  } else {
    return(config)
  }
}






## load_from_parallelize() ################################################
#' @title Load from parallelize
#' @description This function takes the two things passed in a qsub created by
#' \code{parallelize()} - the temp file name and which row in loopvars the current iteration
#' of the child script should load from - and loads the appropriate \code{save_objs} and
#' \code{expand_vars} from \code{parallelize()} into the environment of the child script.
#'
#' Note that this is meant to be run from the \code{child script}; by default
#' both \code{fname} and \code{rownumber} should be loaded appropriately from
#' \code{commandArgs()}
#'
#' @param fname filename of the temp file created by \code{parallelize()}
#' @param rownumber which row of the `lv` object should this particular
#'                  instance of the child script load from
#' @return nothing; assigns objects to child script global environment.
#' @examples
#' \dontrun{
#' # Note: this is within the CHILD SCRIPT, not the master script
#'
#' # Ensure that you've first loaded your functions
#' # (need `load_from_parallelize()` loaded before you can use this)
#' # A good place to put this is right after you've sourced all the
#' # mbg_central function scripts.  Then simply run the
#' # function to set up your environment:
#'
#' load_from_parallelize()
#' }
load_from_parallelize <- function(fname = as.character(commandArgs()[6]),
                                  rownumber = as.numeric(commandArgs()[7])) {
  
  message(paste0("fname: ", fname))
  message(paste0("rownumber: ", rownumber))
  tmp_dir <- "<FILEPATH>"
  
  # Load in the temporary object
  load(paste0(tmp_dir, fname, ".RData"), envir = .GlobalEnv, verbose = T)
  
  lv <- data.frame(lv)
  
  # Load loopvars
  this_lv <- lv[rownumber, -which(names(lv) == "jobname"), drop = FALSE]
  
  # Assign loopvars
  for (n in names(this_lv)) {
    assign(n, this_lv[1, which(names(this_lv) == n)], envir = .GlobalEnv)
    message(paste0(n, ": ", get(n)))
  }
}

#' @title  Clean up path
#' @description Sometimes extra '/''s are added to file paths here and there when they are
#' constructed which makes it difficult to compare a file path (string) to
#' another. This function will clean up an existing file path removing
#' additional '/''s to make it possbile to compare them reliably
#'
#' @param messy_path a file path (character vector). This assumes that the path
#'    is constructed with one or more '/' as a separator and not '\'
#'
#' @return a 'clean' path with a single '/' as a separator and no trailing '/'
#'
clean_path <- function(messy_path) {
  # convert all '/' to white space separating the diretory names between
  dir_names  <- strsplit(messy_path, '/')[[1]]
  # paste together the diretory names with '/' separators ignoring empty space
  clean_path <- paste(dir_names[which(dir_names != "")], collapse = '/')
  # make sure the tack on any prepended '/' if the original path began with a '/'
  if(substr(messy_path, start = 1, stop = 1) == '/') clean_path <- paste0('/', clean_path)
  return(clean_path)
}

#' @title Get Git Status
#'
#' @description Given a path to a directory, determine if it is a git repository, and if so,
#' return a string of git hash, branch, date of last commit, and status plus the
#' diff if requested
#'
#' @param repo The git repository directory
#' @param repo_name The name of the repo to print in the header
#' @param show_diff Print the diff
#'
#' @return character vector of a report intended to be printed with cat
#'
get_git_status <- function(repo, repo_name, show_diff = FALSE) {
  # start the report with the name of the repo in a header
  report <- paste("\n\n**********", repo_name, "**********\n", repo,
                  collapse = "\n ")
  # first check and see if this is a git repository and give a warning if it isn't
  if(!dir.exists(paste0(repo, '/.git'))) {
    report <- paste0(c(report,
                       "\n WARNING: 'repo' does not appear to be a git repository",
                       "\n          Cannot print git hash\n"))
  } else {
    # Collect some git commands for the report
    repo_git_cmd <- paste0("cd ", repo, "; git")
    branch       <- system(paste(repo_git_cmd, "branch | grep \'*\' | cut -d ' ' -f2-"),
                           intern = TRUE)
    commit       <- system(paste(repo_git_cmd, "log -1 --pretty=oneline"),
                           intern = TRUE)
    commit_date  <- system(paste(repo_git_cmd, "log -1 --format=%cd"),
                           intern = TRUE)
    status       <- system(paste(repo_git_cmd, "status --long"), intern = TRUE)
    if(show_diff) diff <- system(paste(repo_git_cmd, "diff HEAD"), intern = TRUE)
    # finish constructing the report
    report <- paste(c(report,
                      "\n******* Branch *******", branch,
                      "\n******* Commit *******", commit,
                      "\n** Commit Date/Time **", commit_date,
                      "\n******* Status *******", status, "\n"),
                    collapse = "\n ")
    if(show_diff) report <- paste(c(report,
                                    "********Diff********", diff, "\n"),
                                  collapse = "\n ")
  }
  return(report)
}


#' @title Record Git Status
#' @description Retrieve information about current git status (eg, branch, commit,
#' uncommitted changes) for core repository and optionally a separate
#' indicator repository. A check can also be made to see if your core
#' repository is in sync with the LBD core repo
#' (\code{<FILEPATH>/lbd_core}) if a fork is being used.
#'
#' Can be used at a minimum to print the git hash of the code being used
#' for posterity.
#'
#' @param core_repo file path to the lbd_core repo
#' @param indic_repo file path to an indicator-specific repo (optional)
#' @param show_diff logical. If there are uncomitted changes, should the
#'     output from git diff be shown?
#' @param check_core_repo logical. Will check whatever has been set as
#'     'core_repo' is the default LBD core code master repo and will give
#'     messages and warnings if not.
#' @param file file path to a text file where output should be saved. This
#'     is optional. If no file path is provided, the output will instead be
#'     printed to the screen.
#' @return Git status
record_git_status <- function(core_repo,
                              indic_repo = NULL,
                              show_diff = FALSE,
                              check_core_repo = TRUE,
                              file = NULL) {
  
  # the core code repo directory
  lbd_core_repo <- '<FILEPATH>/lbd_core'
  
  # if a file is specified, start a sink to record output
  if (!is.null(file)) sink(file)
  
  # print out the core_repo status
  cat(get_git_status(repo = core_repo, repo_name = 'Core repo', show_diff = FALSE))
  
  # if the user wants to make sure their repo is up-to-date with LBD master
  if(check_core_repo) {
    check_repo_report <- "\n********** REPO CHECK **********\n"
    # The two repo paths are the same
    if(clean_path(core_repo) == lbd_core_repo) {
      check_repo_report <- paste0(c(check_repo_report,
                                    "'core_repo' set to default LBD core code master repo: '",
                                    core_repo, "'\n"))
    } else {
      # Get the git hashes of the standard LBD core code repo and user repo and compare
      lbd_core_repo_hash <- system(paste0('cd ', lbd_core_repo, '; git rev-parse HEAD'),
                                   intern = TRUE)
      core_repo_hash     <- system(paste0('cd ', core_repo, '; git rev-parse HEAD'),
                                   intern = TRUE)
      # The two repos paths are not the same, but the hash matches (separate, up-to-date clones)
      if(lbd_core_repo_hash == core_repo_hash) {
        check_repo_report <- paste0(c(check_repo_report,
                                      "Current 'core_repo' clone is up-to-date with LBD core code master repo:\n'",
                                      lbd_core_repo, "' == '", core_repo, "'\n"))
        # The two repo paths are not the same and the hash doesn't match, repos are out of sync
      } else {
        # print out the LBD core code master repo information for reference against current
        # repo being used as a warning.
        check_repo_report <- paste0(c(check_repo_report,
                                      "WARNING: Current 'core_repo' clone is out of sync with the LBD core code master repo:\n'",
                                      lbd_core_repo, "' != '", core_repo, "'\n",
                                      "\n\n** LBD Core Code MASTER Repo Git Info **\n"))
        check_repo_report <- paste0(c(check_repo_report,
                                      get_git_status(repo = lbd_core_repo,
                                                     repo_name = 'LBD Core Repo',
                                                     show_diff = FALSE)))
      }
    }
    cat(check_repo_report)
  }
  # run git log and git status on the indicator repo
  if (!is.null(indic_repo)) {
    cat(get_git_status(repo = indic_repo, repo_name = 'Indicator repo', show_diff = FALSE))
  }
  # if a file is specified, end the sink recording output
  if (!is.null(file)) sink()
}




#' @title Plots, in 3d, a s2-manifold or r2-manifold mesh
#'
#' @description This function allows 3d drawing and imaging of a mesh.inla
#' object. it can plot either a mesh constructed on the R2 or S2
#' manifold. NOTE that it requires the 'rgl' R package and it spawns
#' interactive graphics windows. it is untested on the cluster and is
#' meant for use on local machines
#'
#' @param mesh an inla.mesh object
#'
#' @param draw.edges Logical. Draw the edges between the vertices?
#'
#' @param draw.segments Logical. Draw the segments that bound the mesh
#'   object?
#'
#' @param draw.plane Logical. Draw a planar shape to aid in displaying
#'   curvature of mesh?
#'
#' @param node.cols Numeric vector with length equal to the number of
#'   mesh vertices (mesh$n). e.g. pass in the posterior mean of the
#'   spatial random effects heights to visualize the fitted GP.
#'
#' @param col.type String taking value of either 'bw' or 'col' and
#'   determining whether the color of the surface should be drawn in
#'   black and white or in color. Only used if a node.cols vector is
#'   passed in.
#'
#' @param window.dims 2d numeric vector describing the width and
#'   height of the plotting in pixels
#'
#' @param plot.mirror Logical. Should a mirror image of the mesh be
#'   added to the plot (could help visualize mesh if printing to a
#'   static image)
#'
#' @param returns nothing but spawns an interactive plot window
#'
#' @examples
#' \dontrun{
#' # plot a mesh. add color to the background that results from
#' # the linear interpolation of randomly generated nodal (basis
#' # height) values
#' draw.s2.mesh(mesh_s, draw.edges = T, draw.segments = T, col.type = 'col',
#'              node.cols = rnorm(n = mesh_s$n), draw.plane = F)
#'
#' ## take a snapshot to save to file
#' fig.path <- '/path/to/outputdir/'
#' rgl.snapshot(file.path(fig.path, "mesh.png"), top=TRUE)
#'
#' ## shut down the graphics window
#' rgl.close()
#' }
draw.mesh <- function(mesh, draw.edges=TRUE, draw.segments=TRUE,
                      draw.plane = F, node.cols = NULL,
                      col.type = 'bw', window.dims = c(840, 400),
                      plot.mirror = FALSE){
  
  require('rgl') ## this is an interactive R graphics. won't work on the cluster
  
  window.dims = c(50, 50, 50 + window.dims[1], 50 + window.dims[2])
  
  if(is.null(node.cols)){
    node.cols <- rep(1, mesh$n)
  }
  
  if(col.type == 'bw')  cp <- function (n,...) { return (grey.colors(n,0.95,0.05,...))}
  if(col.type == 'col') cp <- grDevices::colorRampPalette(c("darkblue", "blue", "cyan",
                                                            "yellow", "red", "darkred"))
  
  mesh0 = inla.mesh.create(loc=cbind(0,0), extend=list(offset=1.1,n=4))
  
  mesh01 = mesh0
  mesh02 = mesh0
  mesh1 = mesh
  mesh2 = mesh
  mesh02$loc[,1] = mesh02$loc[,1]*(-1)
  mesh02$loc[,3] = mesh02$loc[,3]*(-1)
  mesh2$loc[,1] = mesh2$loc[,1]*(-1)
  mesh2$loc[,3] = mesh2$loc[,3]*(-1)
  
  mesh01$loc[,1] = mesh01$loc[,1]-1.1
  mesh02$loc[,1] = mesh02$loc[,1]+1.1
  mesh1$loc[,1] = mesh1$loc[,1]-1.1
  mesh2$loc[,1] = mesh2$loc[,1]+1.1
  
  rgl::open3d(windowRect=window.dims)
  if(draw.plane){
    plot(mesh01, rgl=TRUE, col="white", color.palette=cp,
         draw.vertices=FALSE, draw.edges=FALSE, add=TRUE)
    if(plot.mirror){
      plot(mesh02, rgl=TRUE, col="white", color.palette=cp,
           draw.vertices=FALSE, draw.edges=FALSE, add=TRUE)
    }
  }
  plot(mesh1, rgl=TRUE, col=node.cols, color.palette=cp,
       draw.vertices=FALSE, draw.edges=draw.edges, add=TRUE,
       draw.segments=draw.segments)
  if(plot.mirror){
    plot(mesh2, rgl=TRUE, col=node.cols, color.palette=cp,
         draw.vertices=FALSE, draw.edges=draw.edges, add=TRUE,
         draw.segments=draw.segments)
  }
  
  rgl::view3d(0,0,fov=0,zoom=0.4)
  rgl::rgl.bringtotop()
}


#' @title lonlat3D
#' @description  This function takes in a vector of longitude and a vector of
#' latitude and returns coordinates on the S2 sphere (globe living in
#' 3D) in (x, y, z) coords on a sphere with radius 1
#'
#' @param lon numeric vector of longitude coords
#' @param lat numeric vector of latitude coords
#'
#' @return 3 column numeric matrix where each row is a (x,y,z) of the
#'   transformed (long, lat) coords

lonlat3D <- function(lon,lat){
  cbind(cos((lon/180)*pi)*cos((lat/180)*pi),
        sin((lon/180)*pi)*cos((lat/180)*pi),
        sin((lat/180)*pi))
}


#' @title Get output regions
#' @description This function takes in a directory where mbg modeling has stored
# outputs (*cell_pred* objects) and infers the regions specified in
# the model
#'
#' @description Determines modeling regions from written output dir objects
#'
#' @param in_dir directory path containing completed mbg cell_pred objects
#'
#' @return A vector string of region names

get_output_regions <- function(in_dir) {
  return(unique(stringr::str_match(list.files(in_dir, pattern = paste0('_cell_draws_eb_')),
                                   '_cell_draws_[^_]+_[^_]+_(.*)_[0-9].RData')[,2]))
}



#' Return value if provided in function, else return global of same name.
#'
#' Returns the value named \code{name} from the calling function's environment. If that
#' results in an error (e.g., beause the value was not provided) then return an identically
#' named value from the global environment. If no such value exists in the global
#' environment then error.
#'
#' @param name character name of the value to return.
#'
#' @return the value.
#'
#' @export
use_global_if_missing <- function(name) {
  tryCatch(
    {
      return(get(name, pos = parent.frame(1)))
    },
    error = function(e) {
      if (name %in% names(.GlobalEnv)) {
        return(.GlobalEnv[[name]])
      } else {
        stop(sprintf("Variable %s not provided to function and not available in global environment", name))
      }
    }
  )
}

#' @title in a package?
#' @description Predicate: is this code running in a package?
#'
#' @return TRUE if code is running in a package, FALSE otherwise.
#' @export
.in.package <- function() {
  !is.null(utils::packageName())
}


#' @title Returns the stage master list
#'
#' @description Loads from J drive in lbd_core code and from a pre-saved data file in the
#' lbd.mbg package.
#'
#' @export
load_stage_list <- function() {
  if (.in.package()) {
    data("stage_master_list")
    stage_master_list
  } else {
    data.table::fread("<FILEPATH>/stage_master_list.csv")
  }
}
