### rake and aggregate covariates and dtp1


## clear environment
rm(list=ls())

## Set repo location and indicator group
user               <- Sys.info()['user']
core_repo          <- '<FILEPATH>/lbd_core/'
vax_repo         <- '<FILEPATH>/vaccine/'

## sort some directory stuff
commondir      <- sprintf('<FILEPATH>/lbd_core/mbg_central/share_scripts/common_inputs')
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
load_mbg_functions(core_repo)

library(rgdal)
library(raster)
library(sp)
library(data.table)
library(rgeos)
library(sf)
library(dplyr)
library(doParallel)

# Custom load indicator-specific functions
source(paste0(vax_repo,'functions/misc_vaccine_functions.R'))
interval_mo <- 12

library(argparse)
parser <- ArgumentParser()

parser$add_argument("--region",   help="modeling region",  type="character")
parser$add_argument("--shapefile_version",   help="shapefile version",  type="character")
parser$add_argument("--cov_end_year",   help="last year in covariate estimation year list",  type="integer")
parser$add_argument("--dpt_end_year",   help="last year in dpt estimation year list",  type="integer")
parser$add_argument("--pop_measure_cov",   help="population group used for covariate",  type="character")
parser$add_argument("--run_date",   help="overlap run date",  type="character")
parser$add_argument("--dpt3_run_date",   help="run date for dpt3 results",  type="character")
parser$add_argument("--dpt12_run_date",   help="run date for dpt12_cond results",  type="character")


args      <- parser$parse_args()
list2env(args, environment())

cov_year_list <- 2000:cov_end_year
year_list<-2000:dpt_end_year
pop_release <- '2020_03_20'
#####################################################################
# load the cell id to admin units link
modeling_shapefile_version <- shapefile_version

ad0_codes <- get_adm0_codes(region, shapefile_version = shapefile_version)

simple_polygon_list <- load_simple_polygon(gaul_list = ad0_codes, buffer = 1, tolerance = 0.4, shapefile_version = shapefile_version)
subset_shape        <- simple_polygon_list[[1]]
# Load simple raster by region
raster_list        <- build_simple_raster_pop(subset_shape, pop_measure = 'a0004t') 
simple_raster      <- raster_list[['simple_raster']]
remove(raster_list)
gc()

#load itning
itn<- lapply(cov_year_list, function(x){
  f<-raster(paste0('<FILEPATH>/2020_GBD2021_Africa_ITN_Coverage_2000/2020_GBD2021_Africa_ITN_Coverage_',x,'.tif'))
})

if(length(year_list)!=length(cov_year_list)){
  m<-length(year_list) -length(cov_year_list)
  for(x in 1:m){
    itn[[length(itn)+1]] <- itn[[1]]
  }
}
itn<-brick(itn)

#Subset itning down to the region
itn<-crop(itn,simple_raster)
itn<-extend(itn, extent(simple_raster))
itn<-raster::mask(itn,simple_raster)

pixel_id <- seegSDM:::notMissingIdx(simple_raster)

#get link table
link_table <- get_link_table(simple_raster, shapefile_version = shapefile_version)

#####################################################################
# collect and load the under 5 population data
simple_polygon <- simple_raster

covdt_dtp3 <- load_populations_cov(region, pop_measure='a0004t', measure = 'count', simple_polygon, simple_raster, year_list, interval_mo, pixel_id = pixel_id)
covdt_dtp0 <- load_populations_cov(region, pop_measure='a0000t', measure = 'count', simple_polygon, simple_raster, year_list, interval_mo, pixel_id = pixel_id)
covdt_cov <- load_populations_cov(reg, pop_measure_cov, measure = 'count', simple_polygon, simple_raster, year_list, interval_mo, pixel_id = pixel_id)


scalars <- fread(file =paste0("<FILEPATH>/vaccine/dpt3_cov/output/",dpt3_run_date,"/dpt3_cov_", region, "_pop_rf.csv"))
scalars[,V1:=NULL]
## load the raking factors
fractional_rf <- fread(paste0("<FILEPATH>/vaccine/dpt3_cov/output/",dpt3_run_date,"/dpt3_cov_", region, "_rf.csv"))
try(fractional_rf[,location_id:=loc],silent=T)
fractional_rf[,V1:=NULL]

#####################################################################
# Prepping the cell_pred and link table to be linked and then merging them
link <- prep_link_table(
  link_table = link_table,
  simple_raster = simple_raster,
  pixel_id = pixel_id
)

cell_ids <- link_table[[2]]

# getting the connector for sub-national raking
connector <- get_gbd_locs(
  rake_subnational = T,
  reg = region,
  shapefile_version = shapefile_version
)

# getting the connector for sub-national raking
nat_connector <- get_gbd_locs(
  rake_subnational = F,
  reg = region,
  shapefile_version = shapefile_version
)

# merge the connector on to the link table
link <- sub_nat_link_merge(
  rake_subnational=T,
  link,
  connector,
  nat_connector
)

remove(link_table, simple_polygon_list, subset_shape)
gc()



load(paste0("<FILEPATH>/dpt3_cov/output/",dpt3_run_date,"/dpt3_cov_cell_draws_eb_bin0_", region, "_0.RData"))

try(cell_pred<-raked_cell_pred)
try(remove(raked_cell_pred))
gc()

# set cell pred as a data table, and rename things
cell_pred <- prep_cell_pred(cell_pred = cell_pred,
                            cell_ids  = cell_ids,
                            pixel_id  = pixel_id,
                            covdt     = covdt_dtp3)


cell_pred = cbind(cell_pred, covdt_dtp0[,list(dtp0_pop=pop)])
cell_pred = cbind(cell_pred, covdt_cov[,list(cov_pop=pop)])


remove(covdt_cov,covdt_dtp0)
gc()

#Merge in itn
# #merge itn percent and counts first to deal w some NAs
itnv<-as.vector(itn)
sv<-as.vector(simple_raster)
itnv<-itnv[!is.na(sv)]
# itnv[is.na(itnv)] <- 0

remove(itn)
gc()

cell_pred_itn <- data.table(rep(1,(length(year_list) * length(unique(link$ID)))))
cell_pred_itn[, cell_pred_id := .I] #cell_pred_itn object ID
cell_pred_itn[,cell_id := rep(cell_ids, times = nrow(cell_pred_itn) / length(cell_ids))]  #cell id references the africa map
cell_pred_itn[,pixel_id := rep(pixel_id, times = nrow(cell_pred_itn) / length(pixel_id))] #pixel id references the regional map
cell_pred_itn[,itn := itnv] #add itn percent data


cell_pred = cbind(cell_pred, cell_pred_itn[,list(itn)])

overs <- paste0("V", 1:500)

#Align missing pixels from the cell pred and new covariate
cell_pred[is.na(itn) & !is.na(V2),(overs):=NA]
cell_pred[!is.na(itn) & is.na(V2),itn:=NA]
cell_pred[is.na(itn),c('pop','cov_pop','dtp0_pop'):=NA]

# merge cell_pred on the link
cell_pred = merge(link, cell_pred, by.x = 'ID', by.y = 'cell_id',allow.cartesian = TRUE)

# merge on the link

remove(cell_pred_itn, itnv, sv)
gc()

#################################################################################################

## merge scalars onto the cell_pred
cell_pred <- merge(cell_pred, scalars, all.x=T, by = c("location_id", "year"))
cell_pred <- merge(cell_pred, fractional_rf, all.x=T, by = c("location_id", "year"))
try(cell_pred$rf <- as.numeric(as.character(cell_pred$raking_factor)),silent=T)
cell_pred$rf <- as.numeric(as.character(cell_pred$rf))

############################################################
# adding the raking factors and scaling the populations

message("adding raking factors")
# convert to fractional populations
cell_pred <- cell_pred[, cov_pop := cov_pop * area_fraction]
cell_pred <- cell_pred[, dtp0_pop := dtp0_pop * area_fraction]
cell_pred <- cell_pred[, dtp3_pop := pop * area_fraction * pop_scalar]


# multiply the cell_pred by the area fraction for the dedupe function (so that each cell will add to 1 and the constituent rates are weighted by area)
overs <- paste0("V", 1:500)

cell_pred <- cell_pred[, (overs) := lapply(overs, function(x) invlogit(logit(get(x)) + rf))]
cell_pred <- cell_pred[, (overs) := lapply(overs, function(x) get(x) * dtp3_pop) ]
cell_pred <- cell_pred[, itn := lapply('itn', function(x) get(x) * cov_pop) ]
############################################################
# creating a raked counts aggregations
message("creating a raked counts aggregations")

overs <- c(overs,'dtp3_pop','itn','cov_pop','dtp0_pop')

# do calculations!
admin_2_dtp3 <- cell_pred[, lapply(.SD, sum, na.rm=T), .SDcols=c(overs), by = c("year", "ADM2_CODE", "ADM0_CODE")]
admin_1_dtp3 <- cell_pred[, lapply(.SD, sum, na.rm=T), .SDcols=c(overs), by = c("year", "ADM1_CODE", "ADM0_CODE")]
admin_0_dtp3 <- cell_pred[, lapply(.SD, sum, na.rm=T), .SDcols=c(overs), by = c("year", "ADM0_CODE")]

overs <- paste0("V", 1:500)

admin_0_dtp3 <- admin_0_dtp3[, (overs) := lapply(.SD, function(x) x / dtp3_pop), .SDcols=(overs)]
admin_1_dtp3 <- admin_1_dtp3[, (overs) := lapply(.SD, function(x) x / dtp3_pop), .SDcols=(overs)]
admin_2_dtp3 <- admin_2_dtp3[, (overs) := lapply(.SD, function(x) x / dtp3_pop), .SDcols=(overs)]


admin_0_dtp3[,mean:=apply(.SD,1,  mean, na.rm=T), .SDcols=overs]
admin_1_dtp3[,mean:=apply(.SD,1,  mean, na.rm=T), .SDcols=overs]
admin_2_dtp3[,mean:=apply(.SD,1,  mean, na.rm=T), .SDcols=overs]


admin_0_dtp3 <- admin_0_dtp3[, itn_pct := lapply(.SD, function(x) x / cov_pop), .SDcols='itn']
admin_1_dtp3 <- admin_1_dtp3[, itn_pct := lapply(.SD, function(x) x / cov_pop), .SDcols='itn']
admin_2_dtp3 <- admin_2_dtp3[, itn_pct := lapply(.SD, function(x) x / cov_pop), .SDcols='itn']

admin_0_dtp3<-admin_0_dtp3[,list(dtp3_pop,cov_pop,dtp0_pop,dtp3=mean,itn_pct,itn,
                                 year,ADM0_CODE)]
admin_1_dtp3<-admin_1_dtp3[,list(dtp3_pop,cov_pop,dtp0_pop,dtp3=mean,itn_pct,itn,
                                 year,ADM0_CODE,ADM1_CODE)]
admin_2_dtp3<-admin_2_dtp3[,list(dtp3_pop,cov_pop,dtp0_pop,dtp3=mean,itn_pct,itn,
                                 year,ADM0_CODE,ADM2_CODE)]
###############################################################################################################
#Repeat aggregation for dpt12_cond
###############################################################################################################
remove(cell_pred,scalars,fractional_rf)
gc()

scalars <- fread(file =paste0("<FILEPATH>/vaccine/dpt12_cond/output/",dpt12_run_date,"/dpt12_cond_", region, "_pop_rf.csv"))
scalars[,V1:=NULL]
## load the raking factors
fractional_rf <- fread(paste0("<FILEPATH>/vaccine/dpt12_cond/output/",dpt12_run_date,"/dpt12_cond_", region, "_rf.csv"))
try(fractional_rf[,location_id:=loc],silent=T)
fractional_rf[,V1:=NULL]




load(paste0("<FILEPATH>/vaccine/dpt12_cond/output/",dpt12_run_date,"/dpt12_cond_cell_draws_eb_bin0_", region, "_0.RData"))

try(cell_pred<-raked_cell_pred)
try(remove(raked_cell_pred))
gc()

# set cell pred as a data table, and rename things
cell_pred <- prep_cell_pred(cell_pred = cell_pred,
                            cell_ids  = cell_ids,
                            pixel_id  = pixel_id,
                            covdt     = covdt_dtp3)

# merge cell_pred on the link
cell_pred = merge(link, cell_pred, by.x = 'ID', by.y = 'cell_id',allow.cartesian = TRUE)

## merge scalars onto the cell_pred
cell_pred <- merge(cell_pred, scalars, all.x=T, by = c("location_id", "year"))
cell_pred <- merge(cell_pred, fractional_rf, all.x=T, by = c("location_id", "year"))
try(cell_pred$rf <- as.numeric(as.character(cell_pred$raking_factor)),silent=T)
cell_pred$rf <- as.numeric(as.character(cell_pred$rf))

############################################################
# adding the raking factors and scaling the populations

message("adding raking factors")
# convert to fractional populations
cell_pred <- cell_pred[, dtp3_pop := pop * area_fraction * pop_scalar]


# multiply the cell_pred by the area fraction for the dedupe function (so that each cell will add to 1 and the constituent rates are weighted by area)
overs <- paste0("V", 1:500)

cell_pred <- cell_pred[, (overs) := lapply(overs, function(x) invlogit(logit(get(x)) + rf))]
cell_pred <- cell_pred[, (overs) := lapply(overs, function(x) get(x) * dtp3_pop) ]

#make sure pop NAs line up with pred NAs
cell_pred[is.na(V2),c('pop','dtp3_pop'):=NA]

# creating a raked counts aggregations
message("creating a raked counts aggregations")

overs <- c(overs,'dtp3_pop')

# do calculations!
admin_2_dtp12 <- cell_pred[, lapply(.SD, sum, na.rm=T), .SDcols=c(overs), by = c("year", "ADM2_CODE", "ADM0_CODE")]
admin_1_dtp12 <- cell_pred[, lapply(.SD, sum, na.rm=T), .SDcols=c(overs), by = c("year", "ADM1_CODE", "ADM0_CODE")]
admin_0_dtp12 <- cell_pred[, lapply(.SD, sum, na.rm=T), .SDcols=c(overs), by = c("year", "ADM0_CODE")]

overs <- paste0("V", 1:500)

admin_0_dtp12 <- admin_0_dtp12[, (overs) := lapply(.SD, function(x) x / dtp3_pop), .SDcols=(overs)]
admin_1_dtp12 <- admin_1_dtp12[, (overs) := lapply(.SD, function(x) x / dtp3_pop), .SDcols=(overs)]
admin_2_dtp12 <- admin_2_dtp12[, (overs) := lapply(.SD, function(x) x / dtp3_pop), .SDcols=(overs)]


admin_0_dtp12[,mean:=apply(.SD,1,  mean, na.rm=T), .SDcols=overs]
admin_1_dtp12[,mean:=apply(.SD,1,  mean, na.rm=T), .SDcols=overs]
admin_2_dtp12[,mean:=apply(.SD,1,  mean, na.rm=T), .SDcols=overs]

admin_0_dtp12<-admin_0_dtp12[,list(dtp12=mean,year,ADM0_CODE)]
admin_1_dtp12<-admin_1_dtp12[,list(dtp12=mean,year,ADM0_CODE,ADM1_CODE)]
admin_2_dtp12<-admin_2_dtp12[,list(dtp12=mean,year,ADM0_CODE,ADM2_CODE)]

admin_0_dtp12
admin_0_dtp3

###########################################################################################################################
#Combine dtp12_cond and dtp3 to get dtp1, and then dtp0
###########################################################################################################################
admin_0<-merge(admin_0_dtp3,admin_0_dtp12, by=c('ADM0_CODE','year'))
admin_1<-merge(admin_1_dtp3,admin_1_dtp12, by=c('ADM0_CODE','ADM1_CODE','year'))
admin_2<-merge(admin_2_dtp3,admin_2_dtp12, by=c('ADM0_CODE','ADM2_CODE','year'))

#calculate dtp1, subtract to get dtp0
admin_0[,dtp0_pct:= 1 - (dtp12*(1-dtp3)+dtp3)]
admin_1[,dtp0_pct:= 1 - (dtp12*(1-dtp3)+dtp3)]
admin_2[,dtp0_pct:= 1 - (dtp12*(1-dtp3)+dtp3)]

admin_0[,dtp0:=dtp0_pct*dtp0_pop]
admin_1[,dtp0:=dtp0_pct*dtp0_pop]
admin_2[,dtp0:=dtp0_pct*dtp0_pop]

#Remove fake years
admin_0 <- admin_0[year %in% cov_year_list]
admin_1 <- admin_1[year %in% cov_year_list]
admin_2 <- admin_2[year %in% cov_year_list]

## save aggregations
save_dir <- paste0('<FILEPATH>/zerodose_overlap_analyses/',run_date,'/itn/')
dir.create(save_dir, recursive = T)

write.csv(admin_0, paste0(save_dir, region, '_admin_0.csv'), row.names = FALSE)
write.csv(admin_1, paste0(save_dir, region, '_admin_1.csv'), row.names = FALSE)
write.csv(admin_2, paste0(save_dir, region, '_admin_2.csv'), row.names = FALSE)


message('done!')
