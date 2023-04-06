#####################################################################
## Generic parallel script for running MBG models                  ##
## <AUTHOR>   ##
#####################################################################

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~ SETUP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## grab arguments

library(argparse)
parser <- ArgumentParser()

parser$add_argument("--reg",   help="region",  type="character")
parser$add_argument("--age",      help="age", default=0, type="integer")
parser$add_argument("--run_date",  help="run_date",  type="character")
parser$add_argument("--test",    help="testing", default=0, type="integer")
parser$add_argument("--holdout",     help="which holdout", default=0, type="integer")
parser$add_argument("--indicator",      help="which vaccine", type="character")
parser$add_argument("--indicator_group",      help="vaccine", default="vaccine", type="character")

 seed                     <- 121317
args      <- parser$parse_args()
list2env(args, environment())

print(reg)
print(age)
print(run_date)
print(test)
print(holdout)
print(indicator)
print(indicator_group)

## make a pathaddin that get used widely
pathaddin <- paste0('_bin',age,'_',reg,'_',holdout)

## load an image of the main environment
load(paste0('<FILEPATH>', indicator_group, '/', indicator, '/model_image_history/pre_run_tempimage_', run_date, pathaddin,'.RData'))

args      <- parser$parse_args()
list2env(args, environment())

print(reg)
print(age)
print(run_date)
print(test)
print(holdout)
print(indicator)
print(indicator_group)



pathaddin <- paste0('_bin',age,'_',reg,'_',holdout)
outputdir <- file.path('<FILEPATH>',indicator_group,indicator,'output',run_date,'/')
dir.create(outputdir, showWarnings = FALSE)

## print run options
message("options for this run:\n")
for(arg in c('reg','age','run_date','test','holdout',
             'indicator','indicator_group','pathaddin','outputdir'))
  message(paste0(arg,':\t',get(arg),'\t // type: ',class(get(arg))))

# print out session info so we have it on record
sessionInfo()

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
lapply(package_list, require, character.only = TRUE)
library(TMB)
load_mbg_functions(core_repo)

source(paste0('/ihme/code/geospatial/', Sys.info()["user"], '/lbd_hiv/mbg/hiv_prev_disagg/3_functions/read_inla_prior.r'))

## Throw a check for things that are going to be needed later
message('Looking for things in the config that will be needed for this script to run properly')
check_config()

# We need to be in the singularity image, and specifically the LBD one if using TMB
if(!is_singularity()) {
  stop('YOU MUST USE THE SINGULARITY IMAGE TO FIT YOUR MODELS.')
}

if(as.logical(fit_with_tmb) & !is_lbd_singularity()) {
  stop('YOU MUST USE THE LBD SINGULARITY IMAGE IF YOU WANT TO FIT YOUR MODEL USING TMB.')
}

## Print the core_repo hash and check it
message("Printing git hash for 'core_repo'")# and checking against LBD Core Code master repo")

record_git_status(core_repo = core_repo, check_core_repo = FALSE)

## Make sure this inla patch is implemented
#if(grepl("geos", Sys.info()[4])) 
#  INLA:::inla.dynload.workaround()
INLA:::inla.binary.install()

inla.setOption("enable.inla.argument.weights", TRUE)

#stop("stop here and check")
## cores to use
cores_to_use <- Sys.getenv("SGE_HGR_fthread")
message(paste("Model set to use", cores_to_use, "cores"))


#Load pardiso license
inla.setOption("pardiso.license", "<FILEPATH>")

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~ Prep MBG inputs/Load Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PID <- Sys.getpid()
tic("Entire script") # Start master timer


  
  message('You have chosen to not skip directly to inla.')
  
  ## some set up
  if (class(year_list) == "character") year_list <- eval(parse(text=year_list))
  if (class(z_list)    == "character") z_list    <- eval(parse(text=z_list))
  
  if(!exists('interacting_gp_1_effects') & use_gp == TRUE) interacting_gp_1_effects <- "c('space', 'year')"
  if (class(interacting_gp_1_effects)    == "character")   interacting_gp_1_effects <- eval(parse(text=interacting_gp_1_effects))
  
  if(length(trim(strsplit(gbd_fixed_effects, "\\+")[[1]]))==0) gbd_fixed_effects <- ''
  
  if(!exists('use_IND_states')) use_IND_states <- FALSE
  if(!(reg %like% 'vax_soas_sub')) use_IND_states <- FALSE
  
  if(!exists('reweight')) reweight <- FALSE
  message(paste0('reweight: ', reweight))
  if(!exists('set_as_points')) set_as_points <- FALSE
  
  
  ## Load simple polygon template to model over
  gaul_list           <- get_adm0_codes(reg, shapefile_version = modeling_shapefile_version)
  simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list, buffer = 1, tolerance = 0.4,
                                             shapefile_version = modeling_shapefile_version, use_IND_states = use_IND_states)
  subset_shape        <- simple_polygon_list[[1]]
  simple_polygon      <- simple_polygon_list[[2]]
  
  ## Load list of raster inputs (pop and simple)
  raster_list        <- build_simple_raster_pop(subset_shape, pop_release=pop_release, use_IND_states = use_IND_states)
  simple_raster      <- raster_list[['simple_raster']]
  pop_raster         <- raster_list[['pop_raster']]
  
  ## Load input data based on stratification and holdout, OR pull in data as normal and run with the whole dataset if holdout == 0.
  ## For holdouts, we have depreciated val level, so each val level must be recorded in a different run date
  if(holdout==0) {
    message('Holdout == 0 so loading in full dataset using load_input_data()')
    if(!exists("loadinputdata_from_rundate")) {
      df <- load_input_data(indicator   = gsub(paste0('_age',age),'',indicator),
                            agebin      = age,
                            removeyemen = FALSE,
                            pathaddin   = pathaddin,
                            years       = yearload,
                            withtag     = as.logical(withtag),
                            datatag     = datatag,
                            use_share   = as.logical(use_share),
                            yl          = year_list,
                            region      = reg,
                            simple = simple_polygon) #NULL) #Make this null if you don't want buffer data in the model
      
      if(indicator=='dpt12_cond'){
        df[,dpt12_cond:=round(dpt12_cond, 8)]
        df[,N:=round(N, 8)]
      }
      if(indicator=='mcv1_cov'){
        df[,mcv1_cov:=round(mcv1_cov, 8)]
        df[,N:=round(N, 8)]
      }
      if(indicator=='bcg1_cov'){
        df[,bcg1_cov:=round(bcg1_cov, 8)]
        df[,N:=round(N, 8)]
      }
      if(indicator=='dpt3_cov'){
        df[,dpt3_cov:=round(dpt3_cov, 8)]
        df[,N:=round(N, 8)]
      }
      if(indicator=='polio3_cov'){
        df[,polio3_cov:=round(polio3_cov, 8)]
        df[,N:=round(N, 8)]
      }
      if(indicator=='dpt3_timeliness_ratio'){
        df[,dpt3_timeliness_ratio:=round(dpt3_timeliness_ratio, 8)]
        df[,N:=round(N, 8)]
      }
      #!#!#!#!#!#!#!#
      #Drop all ages below 1y and above
      df <- df[age_bin > 3 | (is.na(age_bin) & age_bin_agg!='1:3')]
      
      fwrite(df, paste0(outputdir,'/input_data_',reg,'_subsetted.csv'))
      
      #Merge the z and agebin columns
      df[is.na(age_bin) & !is.na(z), age_bin:=z]
      df[,z:=NULL]
      
      
      save(df, file=paste0(outputdir,"processed_input_data", pathaddin,".rda"))
    } else {
      load(paste0('<FILEPATH>/',indicator_group,"/",indicator,'/output/',loadinputdata_from_rundate,"/","processed_input_data", pathaddin,".rda"))
      save(df, file=paste0(outputdir,"processed_input_data", pathaddin,".rda")) # save to current rundate dir
    }
  } else {
    message(paste0('Holdout != 0 so loading holdout data only from holdout ',holdout))
    message('Please be sure you have a list object called stratum_ho in your environment.')
    
    stratum_ho <- readRDS(paste0('<FILEPATH>/vaccine/', indicator, '/output/', run_date, '/stratum.rds'))
    
    ## if stratifies by age then make sure loads correctly
    if(age!=0) df <- as.data.table(stratum_ho[[paste('region',reg,'_age',age,sep='__')]])
    if(age==0) df <- as.data.table(stratum_ho[[paste('region',reg,sep='__')]])
    df <- df[fold != holdout, ]
    df$first_entry <- 1
    df$agg_weight <- 1
  }
  
  if(as.logical(use_subnat_res)) {
    
    #get adm0 codes for countries to get subnat REs
    if("all" %in% subnat_country_to_get){
      countries_to_get_subnat_res <- get_adm0_codes(reg, shapefile_version = modeling_shapefile_version)
    } else {
      countries_to_get_subnat_res <- get_adm0_codes(subnat_country_to_get, shapefile_version = modeling_shapefile_version)[get_adm0_codes(subnat_country_to_get, shapefile_version = modeling_shapefile_version) %in% get_adm0_codes(reg, shapefile_version = modeling_shapefile_version)]
    }
    
    #load and subset standard admin1 shape to countries with subnational REs 
    subnat_full_shp    <- readRDS(get_admin_shapefile( admin_level = 1, raking = F, suffix = '.rds', version = modeling_shapefile_version ))
    subnat_shapefile <- raster::subset(subnat_full_shp, 
                                       ADM0_CODE %in% countries_to_get_subnat_res)
    
    simple_polygon_list2 <- load_simple_polygon(gaul_list = NULL, 
                                                buffer = 1, tolerance = 0.4, 
                                                custom_shapefile = subnat_shapefile)
    subset_shape2        <- simple_polygon_list2[[1]]
    
    ## Load list of raster inputs (pop and simple)
    raster_list2        <- build_simple_raster_pop(subset_shape2, field = "ADM1_CODE")
    simple_raster2      <- raster_list2[['simple_raster']]
    
    #simple_raster2 has to be the same size as the simple raster for predict to work correctly
    simple_raster2 <- raster::extend(simple_raster2, extent(simple_raster))
    simple_raster2 <- raster::crop(simple_raster2, extent(simple_raster))
    
    ## Merge ADM0/1 codes to df
    adm1_subset_lox <- over(SpatialPoints(df[,.(long = longitude, lat = latitude)], 
                                          CRS(proj4string(subnat_shapefile))), subnat_shapefile)
    df[, subnat_re_ADM1_CODE := as.numeric(as.character(adm1_subset_lox$ADM1_CODE))]
    df[, subnat_re_ADM0_CODE := as.numeric(as.character(adm1_subset_lox$ADM0_CODE))]
    
    #create new ADM1 columns for each country in subnat_country_to_get so data can be fit separately
    for(i in 1:length(unique(na.omit(df$subnat_re_ADM0_CODE)))) {
      df[subnat_re_ADM0_CODE == unique(na.omit(df$subnat_re_ADM0_CODE))[i], (paste0("SUBNAT", i)) := subnat_re_ADM1_CODE]
    }
  } else {
    simple_raster2 <- NULL
  }
  
  ## if testing, we only keep 1000 or so observations
  if(as.logical(test)){
    test_pct <- as.numeric(test_pct)
    
    message(paste0('Test option was set on and the test_pct argument was found at ',test_pct,'% \n
                 ... keeping only ', round(nrow(df)*(test_pct/100),0),' random rows of data.'))
    
    set.seed(seed)
    #seed <- increment_seed(seed)
    df <- df[sample(nrow(df),round(nrow(df)*(test_pct/100),0)),]
    
    message('Also, making it so we only take 100 draws')
    samples <- 100
  }
  
  ## if there is another weight column, multiply it with weight now
  if(exists('other_weight')) if(other_weight!='') {
    message(paste0('Multiplying weight and ',other_weight))
    df[['weight']] <- df[['weight']]*df[[other_weight]]
  }
  
  ## Some built in data checks that cause known problems later on
  if(indicator_family=='binomial' & any(df[,get(indicator)]/df$N > 1))
    stop('You have binomial data where k > N. Check your data before proceeding')
  if(any(df[['weight']] %in% c(Inf,-Inf) | any(is.na(df[['weight']] ))))
    stop('You have illegal weights (NA,Inf,-Inf). Check your data before proceeding')
  
  ## Save distribution of data for this region
  png(paste0(outputdir, reg, '.png'))
  if(indicator_family=='binomial') hist(df[df$first_entry==1, get(indicator)]/df$N[df$first_entry==1]) else hist(df[df$first_entry==1, get(indicator)])
  dev.off()
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~ Pull Covariates ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  ## Define modeling space. In years only for now.
  if(yearload=='annual') period_map <- make_period_map(modeling_periods = c(min(year_list):max(year_list)))
  if(yearload=='five-year') period_map <- make_period_map(modeling_periods = seq(min(year_list), max(year_list), by = 5))
  
  ## Make placeholders for covariates
  cov_layers <- gbd_cov_layers <- NULL
  
  ## Pull all covariate bricks/layers
  # if (nrow(fixed_effects_config) > 0) {
  message('Grabbing raster covariate layers')
  
  if(indicator == 'mcv1_cov'|indicator == 'mcv1_disagg'){
    try(reg_opt <- fread(paste0('<FILEPATH>/vaccine/covariates/mcv1_covs_', reg, '.csv')))
  }
  if(indicator == 'dpt3_timeliness_ratio'){
    reg_opt <- fread(paste0('<FILEPATH>/vaccine/covariates/timeliness_covs.csv'))
  }
  if(indicator == 'polio3_cov'){
    reg_opt <- fread(paste0('<FILEPATH>/vaccine/covariates/polio_covs_', reg, '.csv'), na.strings="")
  }
  if(indicator %in% c('dpt12_cond', 'dpt3_cov')){
    reg_opt <- fread(paste0('<FILEPATH>/vaccine/covariates/dpt_covs_', reg, '.csv'), na.strings="")
  }
  if(indicator == 'bcg1_cov'){
    reg_opt <- fread(paste0('<FILEPATH>/vaccine/covariates/bcg_covs_', reg, '.csv'), na.strings="")
  }  
  
  if(!(exists('reg_opt'))) reg_opt <- fread(paste0('<FILEPATH>/vaccine/covariates/mcv1_covs_vax_ansa.csv'))
  
  reg_opt <- reg_opt[which(reg_opt$gbd=='F'),]
  
  fixed_effects_config <- reg_opt
  fixed_effects_config$release <- ifelse(is.na(fixed_effects_config$release), "", fixed_effects_config$release)
  
  effects <- reg_opt$covariate
  measures <- reg_opt$measure
  loader <- MbgStandardCovariateLoader$new(start_year = min(year_list),
                                           end_year = max(year_list),
                                           interval = as.numeric(interval_mo),
                                           covariate_config = fixed_effects_config)
  cov_layers <- loader$get_covariates(simple_polygon)
  #}
  
  ## Pull country level gbd covariates
  if (nchar(gbd_fixed_effects) > 0) {
    message('Grabbing GBD covariates')

    effects <- trim(strsplit(gbd_fixed_effects, "\\+")[[1]])
    measures <- trim(strsplit(gbd_fixed_effects_measures, "\\+")[[1]])
    gbd_cov_layers <- load_gbd_covariates(covs     = effects,
                                          measures = measures,
                                          year_ids = year_list,
                                          age_ids  = gbd_fixed_effects_age,
                                          template = cov_layers[[1]][[1]],
                                          simple_polygon = simple_polygon,
                                          interval_mo = interval_mo)
  }
  
  ## Combine all covariates
  all_cov_layers <- c(cov_layers, gbd_cov_layers)
  
  ## regenerate all fixed effects equation from the cov layers
  all_fixed_effects <- paste(names(all_cov_layers), collapse = " + ")
  
  ## Make stacker-specific formulas where applicable
  all_fixed_effects_brt <- all_fixed_effects
  
  ## Set Up Country Fixed Effects
  if(use_child_country_fes == TRUE | use_inla_country_fes == TRUE) {
    message('Setting up country fixed effects')
    fe_gaul_list <- unique(c(gaul_convert(unique(df[, country]),
                                          shapefile_version =
                                            modeling_shapefile_version),
                             gaul_list))
    fe_template  <- cov_layers[[1]][[1]]
    simple_polygon_list <- load_simple_polygon(gaul_list   = fe_gaul_list,
                                               buffer      = 0.4,
                                               subset_only = TRUE,
                                               shapefile_version = modeling_shapefile_version)
    fe_subset_shape     <- simple_polygon_list[[1]]
    id_raster <- NULL
    gaul_code <- rasterize_check_coverage(fe_subset_shape, fe_template, field = 'ADM0_CODE') #!#!#!#!#!#!
    gaul_code <- setNames(gaul_code,'gaul_code')
    gaul_code <- create_categorical_raster(gaul_code)

    
    ## update covlayers and add country fixed effects to the
    all_cov_layers = update_cov_layers(all_cov_layers, gaul_code)
    all_fixed_effects_cfes = paste(all_fixed_effects,
                                   paste(names(gaul_code)[1:length(names(gaul_code))],
                                         collapse = " + "), sep=" + ")
    
    ## update specific stacker formulas (for now we just want country effects in BRT)
    all_fixed_effects_brt <- all_fixed_effects_cfes
  }
  
  ## Add these to the fixed effects if we want them in stacking
  if(use_child_country_fes == TRUE) {
    gaul_fes <- paste(names(gaul_code)[2:length(names(gaul_code))], collapse = " + ")
    all_fixed_effects = paste(all_fixed_effects, gaul_fes, sep = " + ")
  }
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~ Stacking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  tic("Stacking - all") ## Start stacking master timer
  
  ## Figure out which models we're going to use
  child_model_names <- stacked_fixed_effects        %>%
    gsub(" ", "", .)          %>%
    strsplit(., "+", fixed=T) %>%
    unlist
  message(paste0('Child stackers included are: ',paste(child_model_names,collapse=' // ')))
  
  the_covs <- format_covariates(all_fixed_effects)
  
  df[, period := year - 1999]
  df[, period_id := year - 1999]
  
  ## copy the dataset to avoid unintended namespace conflicts
  the_data <- copy(df)
  
  ## only use data where we know what age group or point
  ag_data <- the_data[the_data$agg_weight!=1, ]
  the_data <- the_data[the_data$agg_weight==1, ]
  
  ## shuffle the data into six folds
  set.seed(seed)
  #seed <- increment_seed(seed)
  the_data <- the_data[sample(nrow(the_data)),]
  the_data[,fold_id := cut(seq(1,nrow(the_data)),breaks=as.numeric(n_stack_folds),labels=FALSE)]
  
  ## add a row id column
  the_data[, a_rowid := seq(1:nrow(the_data))]
  
  ## extract covariates to the points and subset data where its missing covariate values
  cs_covs <- extract_covariates(the_data,
                                all_cov_layers,
                                id_col              = "a_rowid",
                                return_only_results = TRUE,
                                centre_scale        = TRUE,
                                period_var          = 'year',
                                period_map          = period_map)
  
  # A check to see if any of the variables do not vary across the data. This could break model later so we check and update some objects
  covchecklist <- check_for_cov_issues(check_pixelcount = check_cov_pixelcount,
                                       check_pixelcount_thresh = ifelse(exists("pixelcount_thresh"), as.numeric(pixelcount_thresh), 0.95))
  for(n in names(covchecklist)){
    assign(n, covchecklist[[n]])
  }
  
  # plot covariates as a simple diagnostic here
  pdf(sprintf('%s/raw_covariates_%s.pdf',outputdir,pathaddin), height=12, width=12)
  for(covname in names(all_cov_layers)){
    plot(all_cov_layers[[covname]],main=covname,maxpixel=1e6)
  }
  dev.off()
  
  ## Check for data where covariate extraction failed
  rows_missing_covs <- nrow(the_data) - nrow(cs_covs[[1]])
  if (rows_missing_covs > 0) {
    pct_missing_covs <- round((rows_missing_covs/nrow(the_data))*100, 2)
    warning(paste0(rows_missing_covs, " out of ", nrow(the_data), " rows of data ",
                   "(", pct_missing_covs, "%) do not have corresponding ",
                   "covariate values and will be dropped from child models..."))
    if (rows_missing_covs/nrow(the_data) > 0.1) {
      stop(paste0("Something has gone quite wrong: more than 10% of your data does not have ",
                  "corresponding covariates.  You should investigate this before proceeding."))
    }
  }
  
  the_data <- merge(the_data, cs_covs[[1]], by = "a_rowid", all.x = F, all.y = F)
  
  ## store the centre scaling mapping
  covs_cs_df  <-  cs_covs[[2]]
  
  if(as.logical(use_raw_covs) | as.logical(use_stacking_covs)) {
    ## this will drop rows with NA covariate values
    the_data    <- na.omit(the_data, c(indicator, 'N', the_covs))
  }
  
  ## stop if this na omit demolished the whole dataset
  if(nrow(the_data) == 0) stop('You have an empty df, make sure one of your covariates was not NA everywhere.')
  
  ## seperated out into a different script
  if(as.logical(use_stacking_covs)){
    message('Fitting Stackers')
    
    # Run the child stacker models 
    child_model_run <- run_child_stackers(models = child_model_names, input_data = the_data)
    
    # Bind the list of predictions into a data frame
    child_mods_df <- do.call(cbind, lapply(child_model_run, function(x) x[[1]]))
    
    ## combine the children models with the_data
    the_data  <- cbind(the_data, child_mods_df)
    
    ## Rename the child model objects into a named list
    child_model_objs <- setNames(lapply(child_model_run, function(x) x[[2]]), child_model_names)
    
    #make sure the child stacking models are saved & accessible
    save(child_model_objs, file=paste0(outputdir,reg,'_parent_region_stacking_models.RData'))
    
    
    
    ## return the stacked rasters
    stacked_rasters <- make_stack_rasters(covariate_layers = all_cov_layers, #raster layers and bricks
                                          period           = min(period_map[, period_id]):max(period_map[, period_id]),
                                          child_models     = child_model_objs,
                                          indicator_family = indicator_family,
                                          centre_scale_df  = covs_cs_df)
    
    ## plot stackers
    pdf(paste0(outputdir, 'stacker_rasters', pathaddin, '.pdf'))
    for(i in 1:length(stacked_rasters))
      plot(stacked_rasters[[i]],main=names(stacked_rasters[[i]]),maxpixel=ncell(stacked_rasters[[i]]))
    dev.off()
    
    message('Stacking is complete')
  } ## if(use_stacking_covs)
  
  ## add aggregate data back in, with stacking predictions from the full model
  if (nrow(ag_data) > 0) {
    ag_data[, a_rowid := 1:.N + max(the_data$a_rowid)]
    if(as.logical(use_stacking_covs)) {
      ag_stackers <- extract_covariates(ag_data,
                                        stacked_rasters,
                                        id_col              = "a_rowid",
                                        return_only_results = TRUE,
                                        centre_scale        = FALSE,
                                        period_var          = "year",
                                        period_map          = period_map)
      ag_stackers <- ag_stackers[, c("a_rowid", child_model_names, child_model_names), with = F]
      
      stacker_names <- c(paste0(child_model_names, "_full_pred"), paste0(child_model_names, "_cv_pred"))
      setnames(ag_stackers, c("a_rowid", stacker_names))
      
      ag_data <- merge(ag_data, ag_stackers, by = "a_rowid")
      
      if(any(is.na(ag_data[, ..stacker_names]))) {
        stop("There are NAs in predictions from stackers for aggregated data. 
               Please contact the core code team if you encounter this problem.")
      }
    } else {
      ag_covs <- extract_covariates(ag_data,
                                    all_cov_layers,
                                    id_col              = "a_rowid",
                                    return_only_results = TRUE,
                                    centre_scale        = FALSE,
                                    period_var          = 'year',
                                    period_map          = period_map)
      cov_names <- names(all_cov_layers)
      cs_ag_covs <- centreScale(ag_covs[, ..cov_names], df = covs_cs_df)
      ag_covs <- cbind(ag_covs[,- ..cov_names], cs_ag_covs)
      
      ag_data <- merge(ag_data, ag_covs, by = "a_rowid")
      
      if(as.logical(use_raw_covs)) {
        if(any(is.na(ag_data[, ..cov_names]))) {
          stop("There are NAs in covariates for aggregated data. 
               Please contact the core code team if you encounter this problem.")
        }
      }
    }
    
    the_data <- rbind(the_data, ag_data, fill = T)
  }
  
  
  
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~ Final Pre-MBG Processing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  ## set the fixed effects to use in INLA based on config args
  if(!as.logical(use_stacking_covs) & !as.logical(use_raw_covs) & !as.logical(use_inla_country_fes)){
    all_fixed_effects <- ''
  }
  if(as.logical(use_stacking_covs) & !as.logical(use_raw_covs) & !as.logical(use_inla_country_fes)){
    all_fixed_effects <- stacked_fixed_effects
  }
  if(!as.logical(use_stacking_covs) & as.logical(use_raw_covs) & !as.logical(use_inla_country_fes)){
    all_fixed_effects <- paste(fixed_effects, gbd_fixed_effects, sep = " + ") ## from config
  }
  if(!as.logical(use_stacking_covs) & !as.logical(use_raw_covs) & as.logical(use_inla_country_fes)){
    all_fixed_effects <- paste(names(gaul_code)[2:length(names(gaul_code))], collapse = " + ")
  }
  if(as.logical(use_stacking_covs) & as.logical(use_raw_covs) & !as.logical(use_inla_country_fes)){
    all_fixed_effects <- paste(stacked_fixed_effects, fixed_effects, sep = " + ")
  }
  if(as.logical(use_stacking_covs) & !as.logical(use_raw_covs) & as.logical(use_inla_country_fes)){
    all_fixed_effects <- paste(stacked_fixed_effects, paste(names(gaul_code)[2:length(names(gaul_code))], collapse = " + "), sep = " + ")
  }
  if(!as.logical(use_stacking_covs) & as.logical(use_raw_covs) & as.logical(use_inla_country_fes)){
    all_fixed_effects <- paste(fixed_effects, paste(names(gaul_code)[2:length(names(gaul_code))], collapse = " + "), sep = " + ")
  }
  if(as.logical(use_stacking_covs) & as.logical(use_raw_covs) & as.logical(use_inla_country_fes)){
    all_fixed_effects <- paste(stacked_fixed_effects, fixed_effects, paste(names(gaul_code)[2:length(names(gaul_code))], collapse = " + "), sep = " + ")
  }
  
  ## copy things back over to df
  df <- copy(the_data)
  
  ## remove the covariate columns so that there are no name conflicts when they get added back in
  df <- df[,paste0(the_covs) := rep(NULL, length(the_covs))]
  
  ## Double-check that gaul codes get dropped before extracting in save_mbg_input()
  df <- df[, grep('gaul_code_*', names(df), value = T) := rep(NULL, length(grep('gaul_code_*', names(df), value = T)))]
  
  ## create a full raster list to carry though to the shiny/next steps
  if(as.logical(use_stacking_covs)){
    cov_list      <- c(unlist(stacked_rasters),unlist(all_cov_layers))
    child_mod_ras <- cov_list[child_model_names]
  }else{
    cov_list <- unlist(all_cov_layers)
    child_model_names <- ''
  }
  
  toc(log = T) ## End stacking master timer
  
  ## make sure this inla patch is implemented if running on geos
  if(grepl("geos", Sys.info()[4])) INLA:::inla.dynload.workaround()
  
  set.seed(seed)
  #seed <- increment_seed(seed)
  
  #In case config hasn't been updated to account for multiple mesh options
  if(!exists('s2_mesh_params_int')) s2_mesh_params_int <- s2_mesh_params
  
  ## Build spatial mesh over modeling area--for interacting space terms
  mesh_int <- build_space_mesh(d           = df,
                               simple      = simple_polygon,
                               max_edge    = mesh_s_max_edge,
                               mesh_offset = mesh_s_offset,
                               s2mesh = as.logical(use_s2_mesh),
                               s2params = s2_mesh_params_int)
  
  ##build spatial mesh for space-only GP
  mesh_s   <- build_space_mesh(d           = df,
                               simple      = simple_polygon,
                               max_edge    = mesh_s_max_edge,
                               mesh_offset = mesh_s_offset,
                               s2mesh = as.logical(use_s2_mesh),
                               s2params = s2_mesh_params)
  
  ## Build temporal mesh (standard for now)
  if (length(unique(year_list)) == 1) {
    mesh_t <- NULL
  } else { 
    mesh_t <- build_time_mesh(periods = eval(parse(text = mesh_t_knots)))
  }
  
  
  ## ## For raw covs, don't want to center-scale (as that will happen in `build_mbg_data_stack()`)

  if (as.logical(use_raw_covs) == TRUE) {
    centre_scale_covs <- FALSE
  } else {
    centre_scale_covs <- TRUE
  }
  
  df[, N := round(N)]
  df[, (indicator) := get(indicator) %>% round]
  
  
  ## Save all inputs for MBG model into correct location on <FILEPATH>
  save_mbg_input(indicator         = indicator,
                 indicator_group   = indicator_group,
                 df                = df,
                 simple_raster     = simple_raster,
                 mesh_s            = mesh_s,
                 mesh_int          = mesh_int,
                 mesh_t            = mesh_t,
                 cov_list          = cov_list,
                 pathaddin         = pathaddin,
                 run_date          = run_date,
                 child_model_names = child_model_names,
                 all_fixed_effects = all_fixed_effects,
                 period_map        = period_map,
                 centre_scale      = centre_scale_covs)
  


## reload data an prepare for MBG
load(paste0('<FILEPATH>/', indicator_group, '/', indicator, '/model_image_history/', run_date, pathaddin, '.RData'))

# Bound GBM to 0-1 if desired
if (exists("gbm_bounded_0_1") & use_stacking_covs==T) {
  if (as.logical(gbm_bounded_0_1) == T) {
    message("Truncating Stacker values > 0.999 to 0.999")
    values(cov_list[["gbm"]])[values(cov_list[["gbm"]]) >0.999 & !is.na(values(cov_list[["gbm"]]))] <- 0.999
    values(cov_list[["gam"]])[values(cov_list[["gam"]]) >0.999 & !is.na(values(cov_list[["gam"]]))] <- 0.999
    values(cov_list[["lasso"]])[values(cov_list[["lasso"]]) >0.999 & !is.na(values(cov_list[["lasso"]]))] <- 0.999
    stacker_cols <- grep(paste0("(", paste(child_model_names, collapse="|"), ")(.*_pred)"), names(df), value=T)
    replace_one <- function(x) {
      x[x>0.999 & !is.na(x)] <- 0.999
      return(x)
    }
    df[, (stacker_cols) := lapply(.SD, replace_one), .SDcols = stacker_cols]
  }
}

## convert stackers to transform space, if desired
## NOTE: we do this here to ensure that the stacker rasters are saved in prevalence/untransformed space
## this is useful for diagnostics and other code that was built expecting the untransformed rasters
if (as.logical(stackers_in_transform_space) & indicator_family == 'binomial' & as.logical(use_stacking_covs)){
  message('Converting stackers to logit space')
  
  ## transform the rasters
  for (ii in child_model_names) {
    
    ## Preserve variable names in the raster first
    tmp_rastvar <- names(cov_list[[ii]])
    
    ## Logit
    cov_list[[ii]] <- logit(cov_list[[ii]])
    
    ## Reassign names
    names(cov_list[[ii]]) <- tmp_rastvar
    rm(tmp_rastvar)
  }
  
  ## transform the stacker values that are in df
  stacker_cols <- grep(paste0("(", paste(child_model_names, collapse="|"), ")(.*_pred)"), names(df), value=T)
  df[, (stacker_cols) := lapply(.SD, logit), .SDcols = stacker_cols]
  
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~ Run MBG ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


tic("MBG - all") ## Start MBG master timer

## for stacking, overwrite the columns matching the model_names so that we can trick inla into being our stacker
if(as.logical(use_stacking_covs)){
  df[, paste0(child_model_names) := lapply(child_model_names, function(x) get(paste0(x,'_cv_pred')))]
}

if(!as.logical(fit_with_tmb)){
  ## Generate MBG formula for INLA call (will run but not used by TMB)
  mbg_formula <- build_mbg_formula_with_priors(fixed_effects               = all_fixed_effects,
                                               add_nugget                  = use_inla_nugget,
                                               nugget_prior                = nugget_prior,
                                               add_ctry_res                = use_country_res,
                                               ctry_re_prior               = ctry_re_prior,
                                               temporal_model_theta1_prior = rho_prior,
                                               no_gp                       = !as.logical(use_gp),
                                               use_space_only_gp           = as.logical(use_space_only_gp),
                                               use_time_only_gmrf          = as.logical(use_time_only_gmrf),
                                               # use_sz_gp                   = as.logical(use_sz_gp),
                                               # use_tz_gp                   = as.logical(use_tz_gp),
                                               time_only_gmrf_type         = time_only_gmrf_type,
                                               coefs.sum1                  = coefs_sum1,
                                               subnat_RE                   = use_subnat_res,
                                               subnat_country_to_get       = subnat_country_to_get,
                                               subnat_re_prior             = subnat_re_prior,
                                               timebycountry_RE            = use_timebyctry_res,
                                               adm0_list                    = gaul_list)
}

## For INLA we need to add data for missing time points to ensure we get predictions
##  for all relevant time points. The 0 observations do not contribute to the 
##  model fitting but they prevent INLA from auto-removing 
##  random effects that (conditionally) have no data impacting their fit
if (!as.logical(fit_with_tmb)) {
  if(use_timebyctry_res) {
    ## If we are using a time only effect by country then we need to make sure 
    ##  all year effects are estimated for each country.
    df$adm0code <- gaul_convert(df$country)
    for(adm0_code in gaul_list) {
      dfsub <- df[df$adm0code == adm0_code, ]
      missing_years <- setdiff(year_list, dfsub$year)
      if (length(missing_years) > 0) {
        fake_data <- dfsub[1:length(missing_years), ]
        fake_data[, year := missing_years]
        fake_data[, c(indicator, 'N', 'weight') := 0]
        fake_data[, period := NULL]
        fake_data <- merge(fake_data, period_map)
        df <- rbind(df, fake_data)
      }
    }
  } else {
    ## If not, we only need to make sure we have an observation for each missing
    ##  year (country irrelevant)
    missing_years <- setdiff(year_list, df$year)
    if (length(missing_years) > 0) {
      fake_data <- df[1:length(missing_years), ]
      fake_data[, year := missing_years]
      fake_data[, c(indicator, 'N', 'weight') := 0]
      fake_data[, period := NULL]
      fake_data <- merge(fake_data, period_map)
      df <- rbind(df, fake_data)
    }
  }
}


if(skiptoinla==TRUE){
  df <- process_input_data(df,
                           pop_release                = pop_release,
                           interval_mo                = interval_mo,
                           modeling_shapefile_version = modeling_shapefile_version,
                           poly_ag                    = as.logical(poly_ag),
                           zcol                       = 'z_column_default_blank',
                           zcol_ag                    = zcol_ag,
                           zcol_ag_id                 = zcol_ag_id,
                           z_map                      = z_map,
                           z_ag_mat                   = z_ag_mat)
}

####Try a different weighting system?
if(reweight==T){
  message('reweighting before inla')
  df[,N:=N*weight]
  df[,(indicator):=get(indicator)*weight]
  df[,weight:=1]
}

# get covariate constraints
cov_constraints <- covariate_constraint_vectorize(config)

## Create SPDE INLA stack
df <- subset(df, N > 0 )

## If needed, add fake data to make sure INLA estimates all years
missing_years <- setdiff(year_list, df$year)
if (length(missing_years) > 0 & !as.logical(fit_with_tmb)) {
  fake_data <- df[1:length(missing_years), ]
  fake_data[, year := missing_years]
  fake_data[, c(indicator, 'N', 'weight') := 0]
  fake_data[, period := NULL]
  fake_data <- merge(fake_data, period_map)
  df <- rbind(df, fake_data)
}

#Make sure the zcol is accounted for even if not yet in config
if(!exists('zcol')) {
  zcol = 'agebin'
  df[[zcol]] <- 0
}

if(!('period' %in% names(df))) df[,period:=period_id]

message(paste0('this is the use_error_iid_re: '), as.logical(use_error_iid_re))

input_data <- build_mbg_data_stack(df            = df, # note that merge (if using TMB) will return data in a different (but internally consistent) order, just different than df
                                   fixed_effects = all_fixed_effects,
                                   mesh_s        = mesh_s,
                                   mesh_int = mesh_int,
                                   mesh_t        = mesh_t, # not currently implemented with tmb
                                   spde_prior    = spde_prior,
                                   use_ctry_res  = use_country_res,
                                   use_subnat_res  = use_subnat_res, 
                                   #use_nid_res   = use_nid_res,
                                   use_gp = as.logical(use_gp),
                                   st_gp_int_zero = as.logical(st_gp_int_zero), 
                                   use_space_only_gp = as.logical(use_space_only_gp),
                                   s_gp_int_zero = as.logical(s_gp_int_zero),
                                   use_time_only_gmrf = as.logical(use_time_only_gmrf),
                                   use_timebyctry_res = use_timebyctry_res,
                                   adm0_list = gaul_list,
                                   use_age_only_gmrf = as.logical(use_age_only_gmrf),
                                   use_nugget    = use_inla_nugget, # implemented with tmb
                                   exclude_cs    = child_model_names, # raw covs will get center scaled here though (see notes above)
                                   coefs.sum1    = coefs_sum1, # not currenlty implemented tmb
                                   tmb           = fit_with_tmb,
                                   scale_gaussian_variance_N = scale_gaussian_variance_N,
                                   shapefile_version = modeling_shapefile_version, 
                                   zl            = z_list, # if this is not zero and tmb==TRUE, it will trigger 3rd kronecker and fixed effects
                                   zcol          = zcol, # must not be null if z_list is present
                                   cov_constraints = cov_constraints,
                                   use_error_iid_re = as.logical(use_error_iid_re))   

## combine all the inputs, other than cs_df these are not used if you are using TMB
stacked_input  <- input_data[[1]]
spde_int <- spde      <- input_data[[2]] ## used for space-time gp
cs_df          <- input_data[[3]]
spde.sp        <- input_data[[4]] ## used for space only (time stationary) gp

## Generate other inputs necessary
outcome <- df[[indicator]] # N+_i - event obs in cluster
N       <- df$N                  # N_i - total obs in cluster
weights <- df$weight

## catch in case there is no weight column
if(is.null(weights)){
  weights = rep(1,nrow(df))
}


# save RDS file of input data for replication
saveRDS(object = input_data, ## save this here in case predict dies
        file = sprintf('/ihme/geospatial/mbg/%s/%s/output/%s/%s_data_input_list_%s_holdout_%i_agebin_%i.RDS',
                       indicator_group, indicator, run_date, ifelse(fit_with_tmb,'tmb','inla'), reg, holdout, age))

input_data<-readRDS( ## save this here in case predict dies
  file = sprintf('/ihme/geospatial/mbg/%s/%s/output/%s/%s_data_input_list_%s_holdout_%i_agebin_%i.RDS',
                 indicator_group, indicator, run_date, ifelse(fit_with_tmb,'tmb','inla'), reg, holdout, age))

################################################################################################################
tic("MBG - fit model") ## Start MBG - model fit timer

## Set the number of cores to be equal to input;
## If missing, then revert to cores_to_use value
if(Sys.getenv("OMP_NUM_THREADS") != "") {
  setompthreads(Sys.getenv("OMP_NUM_THREADS"))
} else {
  print("Threading information not found; setting cores_to_use as the input OpenMP threads.")
  setompthreads(cores_to_use)
}

if(Sys.getenv("MKL_NUM_THREADS") != "") {
  setmklthreads(Sys.getenv("MKL_NUM_THREADS"))
} else {
  print("Threading information not found; setting cores_to_use as the input MKL threads.")
  setmklthreads(cores_to_use)
}


## Fit MBG model
if(!as.logical(skipinla)) {
  if(fit_with_tmb == FALSE) {
    message('Fitting model with R-INLA')
    
    model_fit <- fit_mbg(indicator_family = indicator_family,
                         stack.obs        = stacked_input,
                         spde             = spde_int,
                         cov              = outcome,
                         N                = N,
                         int_prior_mn     = intercept_prior,
                         f_mbg            = mbg_formula,
                         run_date         = run_date,
                         keep_inla_files  = keep_inla_files,
                         cores            = 1, #cores_to_use,
                         wgts             = weights,
                         intstrat         = intstrat,
                         fe_sd_prior      = 1 / 9, ## this actually sets precision!. prec=1/9 -> sd=3
                         sparse_ordering  = as.logical(sparse_ordering))    
    
  } else {
    message('Fitting model with TMB')
    message(sprintf('%i Data points and %i mesh nodes',nrow(df),length(input_data$Parameters$Epsilon_int_1)))
    
    # run the model
    system.time(
      model_fit <- fit_mbg_tmb( lbdcorerepo     = core_repo,
                                cpp_template    = 'mbg_tmb_model',
                                tmb_input_stack = input_data,
                                control_list    = list(trace=1, eval.max=500, iter.max=300, abs.tol=1e-20),
                                optimizer       = 'nlminb', # TODO add optimx
                                ADmap_list      = NULL,
                                sparse_ordering = as.logical(sparse_ordering),
                                int_gp_1_effs   = interacting_gp_1_effects) 
    )
  }
  
  saveRDS(object = model_fit, ## save this here in case predict dies
          file = sprintf('<FILEPATH>/%s/%s/output/%s/%s_model_fit_pre_preds_%s_holdout_%i_agebin_%i.RDS',
                         indicator_group, indicator, run_date, ifelse(fit_with_tmb,'tmb','inla'), reg, holdout, age))
  
  
}else{
  ## skipped fitting INLA so just load model and move to predict
  model_fit <- readRDS( file = sprintf('<FILEPATH>/%s/%s/output/%s/%s_model_fit_pre_preds_%s_holdout_%i_agebin_%i.RDS',
                                       indicator_group, indicator, run_date, ifelse(fit_with_tmb,'tmb','inla'), reg, holdout, age))
  
}

toc(log = T) ## End MBG - model fit timer

tic("MBG - predict model") ## Start MBG - model predict timer

## Run predict_mbg on chunks of 50 samples (to avoid memory issues)
message('Making predictions in 50 draw chunks.')

max_chunk <- 50
samples   <- as.numeric(samples)

## Create vector of chunk sizes
chunks <- rep(max_chunk, samples %/% max_chunk)
if (samples %% max_chunk > 0) chunks <- c(chunks, samples %% max_chunk)

if(fit_with_tmb == FALSE) {
  pm <- lapply(chunks, function(samp) {
    predict_mbg(res_fit       = model_fit,
                cs_df         = cs_df,
                mesh_s        = mesh_int,
                mesh_t        = mesh_t,
                cov_list      = cov_list,
                samples       = samp,
                simple_raster = simple_raster,
                transform     = transform,
                coefs.sum1    = coefs_sum1,
                pred_gp       = as.logical(use_gp),
                yl            = year_list,
                use_space_only_gp = as.logical(use_space_only_gp),
                use_time_only_gmrf = as.logical(use_time_only_gmrf),
                use_timebyctry_res = as.logical(use_timebyctry_res),
                shapefile_version = modeling_shapefile_version,
                simple_raster_subnats = simple_raster2,
                subnat_country_to_get = subnat_country_to_get)[[3]]
  })
  
  ## Make cell preds and a mean raster
  cell_pred <- do.call(cbind, pm)
  mean_ras  <- insertRaster(simple_raster,matrix(rowMeans(cell_pred),ncol = max(period_map$period)))
  toc(log = T) # Stop MBG - model predict timer
  
  
  
  message('Saving preds')
  save_mbg_preds(config     = config,
                 time_stamp = time_stamp,
                 run_date   = run_date,
                 mean_ras   = mean_ras,
                 sd_ras     = NULL,
                 res_fit    = model_fit,
                 cell_pred  = cell_pred,
                 df         = df,
                 pathaddin  = pathaddin)
  
  
  # plot the mean raster
  pdf(paste0(outputdir,'mean_rasterXX', pathaddin, '.pdf'))
  plot(mean_ras,maxpixel=1e6)
  dev.off()
  
} else { 
  

  
  
  outputdir_chunks <- file.path(outputdir, 'prediction_chunks')
  dir.create(outputdir_chunks, showWarnings = FALSE)
  
  for (chunk in 1:length(chunks)) {
    pm<-  predict_mbg_tmb(samples          = chunks[chunk],
                          seed                 = NULL,
                          tmb_input_stack      = input_data,
                          model_fit_object     = model_fit,
                          fes                  = all_fixed_effects, # TODO use input_data or model_fit object for this (in case its changed due to checks)
                          sr                   = simple_raster,
                          yl                   = year_list,
                          zl                   = z_list,
                          covs_list            = cov_list,  #instead of the temporary covs_list
                          clamp_covs           = clamp_covs,
                          cov_constraints = covariate_constraint_vectorize(config),
                          use_space_only_gp = as.logical(use_space_only_gp),
                          use_time_only_gmrf = as.logical(use_time_only_gmrf),
                          use_age_only_gmrf = as.logical(use_age_only_gmrf),
                          int_gp_1_effect =            as.logical(use_gp),
                          #use_covs         = FALSE,
                          int_gp_1_effs   = interacting_gp_1_effects,
                          use_sz_gp = as.logical(use_sz_gp),
                          use_tz_gp = as.logical(use_tz_gp),
                          mesh_int = mesh_int,
                          mesh_s = mesh_s)  
    
    for(z in 1:length(z_list)) {
      if (length(z_list)>1){
        pm_zx <- pm[[z]] 
      }  else {
        pm_zx<- pm
      }
      
      pathaddin <- paste0('_bin', z_list[z], '_chunk',chunk, '_', reg,'_',holdout)
      message(paste0('Saving draws for chunk ', chunk, ' for agebin ', z_list[z]))
      save(
        pm_zx,
        file     = paste0(outputdir_chunks, '/', indicator,'_cell_draws_eb',pathaddin,'.RData'),
        compress = TRUE
      )
      
      rm(pm_zx)
      gc()
    }
    rm(pm)
    gc()
  }
  
  for (z in z_list) {
    cell_pred<-NULL
    message('Reload and combine temporary chunk files agebin ', z)
    for (chunk in 1:length(chunks)) {
      pathaddin <- paste0('_bin',z, '_chunk',chunk, '_', reg,'_',holdout)
      load(paste0(outputdir_chunks, '/', indicator,'_cell_draws_eb',pathaddin,'.RData'))
      
      cell_pred <- cbind(cell_pred, pm_zx)
      
      rm(pm_zx)
      gc()
    }
    
    pathaddin <- paste0('_bin',z, '_',reg,'_',holdout) # new pathaddin
    
    # make a mean raster
    library(matrixStats)
    mean_ras  <- insertRaster(simple_raster,matrix(rowMeans(cell_pred),ncol = max(period_map$period)))
    sd_ras    <- insertRaster(simple_raster,matrix(  rowSds(cell_pred),ncol = max(period_map$period)))
    
    # save z specific objects
    writeRaster(
      mean_ras,
      file      = paste0(outputdir, '/', indicator,'_prediction_eb',pathaddin),
      overwrite = TRUE
    )
    
    save(
      cell_pred,
      file = (paste0(outputdir, "/", indicator, "_cell_draws_eb", pathaddin, ".RData")),
      compress = TRUE
    )
    
    #Create mean & sd pdfs
    pdf(paste0(outputdir,'mean_raster', pathaddin, '.pdf'))
    plot(mean_ras,main='mean',maxpixel=1e6)
    plot(sd_ras,main='sd',maxpixel=1e6)
    dev.off()
    
    #Clean up
    rm(cell_pred)
    rm(mean_ras)
    rm(sd_ras)
    gc()
    
  }
  
  #Remove temporary chunk files
  message('Remove temporary chunk files')
  for (z in z_list) {
    for (chunk in 1:length(chunks)) {
      pathaddin <- paste0('_bin',z,'_chunk', chunk, '_', reg,'_',holdout)
      file.remove(paste0(outputdir_chunks, '/', indicator,'_cell_draws_eb',pathaddin,'.RData'))
    }
  }
}     



##################################################################
message("wrapping up")
## timer stuff
toc(log = T) # End master timer

## Format timer
ticlog   <- tic.log(format = F)
df_timer <- generate_time_log(ticlog)
df_timer[, region := reg]
df_timer[, holdout := holdout]
setcolorder(df_timer, c("region", "holdout", "step", "time"))

## Pull run time for this run
run_time_all <- df_timer[step == "Entire script", time]

## Write to a run summary csv file in the output directory
output_file <- paste0(outputdir, "run_summary_", indicator, "_", run_date, ".csv")

## Pull in file contents from other reg/holdouts (if exists)
if (file.exists(output_file)) {
  file_contents <- read.csv(output_file, stringsAsFactors = F) %>% as.data.table
  df_timer      <- rbind(file_contents, df_timer)
}


# Write a an empty file to indicate done with this parallel script
write(NULL, file = paste0(outputdir, "/fin_", pathaddin))

## Write CSV
write.csv(df_timer,file = output_file, row.names = FALSE)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FIN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~