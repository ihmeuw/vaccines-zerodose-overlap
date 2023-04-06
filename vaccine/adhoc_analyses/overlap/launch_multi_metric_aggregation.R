library(data.table)
library(sf)

###############################################
#
run_align    <- FALSE
run_combines <- TRUE

run_date          <- '2023_03_03'
shapefile_version <- '2022_05_09'
dpt3_run_date     <- '2023_02_15'
dpt12_run_date    <- '2023_02_15'

dpt_end_year      <- 2021

stunting_end_year <- 2019
u5m_end_year      <- 2017
ort_end_year      <- 2017
itn_end_year      <- 2019
lf_end_year       <- 2018

pop_measure_stunting <- 'a0004t'
pop_measure_u5m      <- 'a0004t'
pop_measure_ort      <- 'a0004t'
pop_measure_itn      <- 'total'
pop_measure_lf       <- 'total'

reg_mem <-data.table(region=c('vax_cssa','vax_essa','vax_seas','vax_wssa'),
                     memory=c(115,100,200,130))

stunting_regions<-u5m_regions<-ort_regions<-lf_regions<-c('vax_cssa','vax_wssa','vax_essa','vax_seas')
itn_regions       <-c('vax_cssa','vax_wssa','vax_essa')


if(run_align){
  #####################################################################################
#Launching for which regions

dir.create(paste0('<FILEPATH>/zerodose_overlap_analyses/',run_date,'/'))
dir.create(paste0('<FILEPATH>/zerodose_overlap_analyses/',run_date,'/stunting/'))
dir.create(paste0('<FILEPATH>/zerodose_overlap_analyses/',run_date,'/stunting/errors/'))
dir.create(paste0('<FILEPATH>/zerodose_overlap_analyses/',run_date,'/stunting/output/'))

#Launch aggregation script by region
for(r in stunting_regions) {

  
  mem <- reg_mem[region==r]$memory

sbatch<- paste0('<ADDRESS> ',
                '-s /<FILEPATH>/vaccine/adhoc_analyses/overlap/aggregate_stunting_rasters.R '
                ,
                paste('--region',r,
                      '--shapefile_version', shapefile_version,
                      '--cov_end_year', stunting_end_year,
                      '--dpt_end_year',dpt_end_year,
                      '--pop_measure_cov',pop_measure_stunting,
                      '--run_date',run_date,
                      '--dpt3_run_date',dpt3_run_date,
                      '--dpt12_run_date',dpt12_run_date)
)

system(sbatch)

}

###########################################################################################
#U5M
dir.create(paste0('<FILEPATH>/zerodose_overlap_analyses/',run_date,'/'))
dir.create(paste0('<FILEPATH>/zerodose_overlap_analyses/',run_date,'/u5m/'))
dir.create(paste0('<FILEPATH>/zerodose_overlap_analyses/',run_date,'/u5m/errors/'))
dir.create(paste0('<FILEPATH>/zerodose_overlap_analyses/',run_date,'/u5m/output/'))


for(r in u5m_regions) {

  
  mem <- reg_mem[region==r]$memory
  sbatch<- paste0('<ADDRESS> ',
                  '-s /<FILEPATH>/vaccine/adhoc_analyses/overlap/aggregate_u5m_rasters.R '
                  ,
                  paste('--region',r,
                        '--shapefile_version', shapefile_version,
                        '--cov_end_year', u5m_end_year,
                        '--dpt_end_year',dpt_end_year,
                        '--pop_measure_cov',pop_measure_u5m,
                        '--run_date',run_date,
                        '--dpt3_run_date',dpt3_run_date,
                        '--dpt12_run_date',dpt12_run_date)
  )
  
  
  system(sbatch)
  
}
###############################################################################
#ORT
dir.create(paste0('<FILEPATH>/zerodose_overlap_analyses/',run_date,'/'))
dir.create(paste0('<FILEPATH>/zerodose_overlap_analyses/',run_date,'/ort/'))
dir.create(paste0('<FILEPATH>/zerodose_overlap_analyses/',run_date,'/ort/errors/'))
dir.create(paste0('<FILEPATH>/zerodose_overlap_analyses/',run_date,'/ort/output/'))

for(r in ort_regions) {

  mem <- reg_mem[region==r]$memory
  sbatch<- paste0('<ADDRESS> ',
                  '-s /<FILEPATH>/vaccine/adhoc_analyses/overlap/aggregate_ort_rasters.R '
                  ,
                  paste('--region',r,
                        '--shapefile_version', shapefile_version,
                        '--cov_end_year', ort_end_year,
                        '--dpt_end_year',dpt_end_year,
                        '--pop_measure_cov',pop_measure_ort,
                        '--run_date',run_date,
                        '--dpt3_run_date',dpt3_run_date,
                        '--dpt12_run_date',dpt12_run_date)
  )
  
  
  system(sbatch)
  
}

###############################################################################
#itn
dir.create(paste0('<FILEPATH>/zerodose_overlap_analyses/',run_date,'/'))
dir.create(paste0('<FILEPATH>/zerodose_overlap_analyses/',run_date,'/itn/'))
dir.create(paste0('<FILEPATH>/zerodose_overlap_analyses/',run_date,'/itn/errors/'))
dir.create(paste0('<FILEPATH>/zerodose_overlap_analyses/',run_date,'/itn/output/'))


for(r in itn_regions) {

  mem <- reg_mem[region==r]$memory
  sbatch<- paste0('<ADDRESS> ',
                  '-s /<FILEPATH>/vaccine/adhoc_analyses/overlap/aggregate_itn_rasters.R '
                  ,
                  paste('--region',r,
                        '--shapefile_version', shapefile_version,
                        '--cov_end_year', itn_end_year,
                        '--dpt_end_year',dpt_end_year,
                        '--pop_measure_cov',pop_measure_itn,
                        '--run_date',run_date,
                        '--dpt3_run_date',dpt3_run_date,
                        '--dpt12_run_date',dpt12_run_date)
  )
  
  
  system(sbatch)
  
}

#LF

dir.create(paste0('<FILEPATH>/zerodose_overlap_analyses/',run_date,'/'))
dir.create(paste0('<FILEPATH>/zerodose_overlap_analyses/',run_date,'/lf/'))
dir.create(paste0('<FILEPATH>/zerodose_overlap_analyses/',run_date,'/lf/errors/'))
dir.create(paste0('<FILEPATH>/zerodose_overlap_analyses/',run_date,'/lf/output/'))

for(r in lf_regions) {
 
  mem <- reg_mem[region==r]$memory
  sbatch<- paste0('<ADDRESS> ',
                  '-s /<FILEPATH>/vaccine/adhoc_analyses/overlap/aggregate_lf_rasters.R '
                  ,
                  paste('--region',r,
                        '--shapefile_version', shapefile_version,
                        '--cov_end_year', lf_end_year,
                        '--dpt_end_year',dpt_end_year,
                        '--pop_measure_cov',pop_measure_lf,
                        '--run_date',run_date,
                        '--dpt3_run_date',dpt3_run_date,
                        '--dpt12_run_date',dpt12_run_date)
  )
  
  
  system(sbatch)
  
}

}
######################################################################
#Combine regional data for each metric
######################################################################
#Add location name information
if(run_combines){
shp_adm2 <- as.data.table(read_sf(dsn=paste0('<FILEPATH>/',shapefile_version,'/lbd_standard_admin_2.shp')))
shp_adm1 <- as.data.table(read_sf(dsn=paste0('<FILEPATH>/',shapefile_version,'/lbd_standard_admin_1.shp')))
shp_adm0 <- as.data.table(read_sf(dsn=paste0('<FILEPATH>/',shapefile_version,'/lbd_standard_admin_0.shp')))

shp_adm2[,geometry:=NULL]
shp_adm1[,geometry:=NULL]
shp_adm0[,geometry:=NULL]


save_dir <- paste0('<FILEPATH>/zerodose_overlap_analyses/',run_date,'/')
dir.create(save_dir, recursive = T)
#########################################################################################
#stunting
#Combine regions, add country & district names
print('combine stunting')
a0<-rbindlist(lapply(stunting_regions,function(r){
  a0<-fread(paste0(save_dir,'/stunting/',r,'_admin_0.csv'))
}))
a1<-rbindlist(lapply(stunting_regions,function(r){
  a1<-fread(paste0(save_dir,'/stunting/',r,'_admin_1.csv'))
}))
a2<-rbindlist(lapply(stunting_regions,function(r){
  a2<-fread(paste0(save_dir,'/stunting/',r,'_admin_2.csv'))
}))

a0<-merge(a0,shp_adm0,by='ADM0_CODE')
a1<-merge(a1,shp_adm1,by=c('ADM0_CODE','ADM1_CODE'))
a2<-merge(a2,shp_adm2,by=c('ADM0_CODE','ADM2_CODE'))

a2[,c('geo_id', 'ad2_id', 'ad0_parent', 'ad1_parent'):=NULL]

write.csv(a0, paste0(save_dir, '/stunting/all_admin_0.csv'), row.names = FALSE)
write.csv(a1, paste0(save_dir, '/stunting/all_admin_1.csv'), row.names = FALSE)
write.csv(a2, paste0(save_dir, '/stunting/all_admin_2.csv'), row.names = FALSE)

####################################################################################
#U5M

#Combine regions, add country & district names

print('combine u5m')
a0<-rbindlist(lapply(u5m_regions,function(r){
  a0<-fread(paste0(save_dir,'/u5m/',r,'_admin_0.csv'))
}))
a1<-rbindlist(lapply(u5m_regions,function(r){
  a1<-fread(paste0(save_dir,'/u5m/',r,'_admin_1.csv'))
}))
a2<-rbindlist(lapply(u5m_regions,function(r){
  a2<-fread(paste0(save_dir,'/u5m/',r,'_admin_2.csv'))
}))

a0<-merge(a0,shp_adm0,by='ADM0_CODE')
a1<-merge(a1,shp_adm1,by=c('ADM0_CODE','ADM1_CODE'))
a2<-merge(a2,shp_adm2,by=c('ADM0_CODE','ADM2_CODE'))

a2[,c('geo_id', 'ad2_id', 'ad0_parent', 'ad1_parent'):=NULL]

write.csv(a0, paste0(save_dir, '/u5m/all_admin_0.csv'), row.names = FALSE)
write.csv(a1, paste0(save_dir, '/u5m/all_admin_1.csv'), row.names = FALSE)
write.csv(a2, paste0(save_dir, '/u5m/all_admin_2.csv'), row.names = FALSE)


####################################################################################
#ort
#Combine regions, add country & district names

print('combine ort')
a0<-rbindlist(lapply(ort_regions,function(r){
  a0<-fread(paste0(save_dir,'/ort/',r,'_admin_0.csv'))
}))
a1<-rbindlist(lapply(ort_regions,function(r){
  a1<-fread(paste0(save_dir,'/ort/',r,'_admin_1.csv'))
}))
a2<-rbindlist(lapply(ort_regions,function(r){
  a2<-fread(paste0(save_dir,'/ort/',r,'_admin_2.csv'))
}))


a0<-merge(a0,shp_adm0,by='ADM0_CODE')
a1<-merge(a1,shp_adm1,by=c('ADM0_CODE','ADM1_CODE'))
a2<-merge(a2,shp_adm2,by=c('ADM0_CODE','ADM2_CODE'))

a2[,c('geo_id', 'ad2_id', 'ad0_parent', 'ad1_parent'):=NULL]

write.csv(a0, paste0(save_dir, '/ort/all_admin_0.csv'), row.names = FALSE)
write.csv(a1, paste0(save_dir, '/ort/all_admin_1.csv'), row.names = FALSE)
write.csv(a2, paste0(save_dir, '/ort/all_admin_2.csv'), row.names = FALSE)


####################################################################################
#itns
#Combine regions, add country & district names

print('combine itn')
a0<-rbindlist(lapply(itn_regions,function(r){
  a0<-fread(paste0(save_dir,'/itn/',r,'_admin_0.csv'))
}))
a1<-rbindlist(lapply(itn_regions,function(r){
  a1<-fread(paste0(save_dir,'/itn/',r,'_admin_1.csv'))
}))
a2<-rbindlist(lapply(itn_regions,function(r){
  a2<-fread(paste0(save_dir,'/itn/',r,'_admin_2.csv'))
}))


a0<-merge(a0,shp_adm0,by='ADM0_CODE')
a1<-merge(a1,shp_adm1,by=c('ADM0_CODE','ADM1_CODE'))
a2<-merge(a2,shp_adm2,by=c('ADM0_CODE','ADM2_CODE'))

a2[,c('geo_id', 'ad2_id', 'ad0_parent', 'ad1_parent'):=NULL]

write.csv(a0, paste0(save_dir, '/itn/all_admin_0.csv'), row.names = FALSE)
write.csv(a1, paste0(save_dir, '/itn/all_admin_1.csv'), row.names = FALSE)
write.csv(a2, paste0(save_dir, '/itn/all_admin_2.csv'), row.names = FALSE)



####################################################################################
#lf
#Combine regions, add country & district names

print('combine lf')
a0<-rbindlist(lapply(lf_regions,function(r){
  a0<-fread(paste0(save_dir,'/lf/',r,'_admin_0.csv'))
}))
a1<-rbindlist(lapply(lf_regions,function(r){
  a1<-fread(paste0(save_dir,'/lf/',r,'_admin_1.csv'))
}))
a2<-rbindlist(lapply(lf_regions,function(r){
  a2<-fread(paste0(save_dir,'/lf/',r,'_admin_2.csv'))
}))

a0<-merge(a0,shp_adm0,by='ADM0_CODE')
a1<-merge(a1,shp_adm1,by=c('ADM0_CODE','ADM1_CODE'))
a2<-merge(a2,shp_adm2,by=c('ADM0_CODE','ADM2_CODE'))

a2[,c('geo_id', 'ad2_id', 'ad0_parent', 'ad1_parent'):=NULL]

write.csv(a0, paste0(save_dir, '/lf/all_admin_0.csv'), row.names = FALSE)
write.csv(a1, paste0(save_dir, '/lf/all_admin_1.csv'), row.names = FALSE)
write.csv(a2, paste0(save_dir, '/lf/all_admin_2.csv'), row.names = FALSE)

}
