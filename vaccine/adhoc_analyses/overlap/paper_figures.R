library(data.table)
library(ggplot2)
library(sf)
library(rgdal)
library(Hmisc)
library(gridExtra)
#Required for using st_intersection
sf::sf_use_s2(FALSE)
library(scales)
library(plyr)
require(DescTools)
library(viridis)
library(tictoc)
library(cowplot)


rm(list=ls())

run_date<- '2023_03_03'
shapefile_version <- '2022_05_09'

#Which plots to run?
make_global_high_priority_plots   <- TRUE
make_national_high_priority_plots <- TRUE
make_global_all_priority_plots    <- TRUE
make_national_all_priority_plots  <- TRUE
make_auc_plots                    <- TRUE #includes step plots
make_multi_metric_overlap_plots   <- TRUE
make_combo_plots                  <- TRUE #can only work if national high priority plots also == TRUE

dir.create(paste0('<FILEPATH>/zerodose_overlap_analyses/',run_date,'/results/paper/global_priority_plots/'),recursive = T)
dir.create(paste0('<FILEPATH>/zerodose_overlap_analyses/',run_date,'/results/paper/national_priority_plots/'),recursive = T)

## Set repo location and indicator group
user               <- Sys.info()['user']
core_repo          <- '<FILEPATH>/lbd_core/'
vax_repo           <- '<FILEPATH>/vaccine/'

## sort some directory stuff
commondir      <- sprintf('<FILEPATH>/lbd_core/mbg_central/share_scripts/common_inputs')
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
load_mbg_functions(repo=core_repo)

# Custom load indicator-specific functions
source(paste0(vax_repo,'functions/misc_vaccine_functions.R'))

#Get the red/purple/blue color scale for overlap figures
colors <-    c('#E3DAE3', '#DDA3A6', '#D07C80', '#DA434A',
               '#B6CBE1','#BDAFC4',  '#bc7c8f',  '#ae3a4e',
               '#9BBDDF', '#89a1c8', '#806a8a','#77324c',
               '#5AA1E9', '#4885c1', '#435786','#3f2949' )
#names for calling folders
covariate_names_abbrev <- c('stunting','u5m','ort','itn',
                            'lf')
#names for labels
covariate_names_full <- c('stunting','mortality','ORT','ITNs',
                          'LF')
#prefixes to add to labels
covariate_additions <-c('Children with','Under-5','Missed','Missed',
                        '')
#What year are we looking at for which variable
cov_years <- c(2019,2017,2017,2019,
               2018)

#read in shapefiles for later
shp_adm2_0 <- read_sf(dsn=paste0('<FILEPATH>/',shapefile_version,'/lbd_standard_admin_2.shp'))
shp_adm0   <- read_sf(dsn=paste0('<FILEPATH>/',shapefile_version,'/lbd_standard_admin_0_stage_1_2.shp'))
shp_adm0   <- shp_adm0[,1:4]

#tables to write to for summarizing across metrics
s<-data.table()
sc<-data.table()

o<-data.table()
oc<-data.table()

cc<-data.table()
nn<-data.table()

mc0<-data.table()

shp_adm0[shp_adm0$ADM0_NAME=='Democratic Republic of the Congo',]$ADM0_NAME <- 'D.R.C.'
shp_adm2_0[shp_adm2_0$ADM0_NAME=='Democratic Republic of the Congo',]$ADM0_NAME <- 'D.R.C.'

countries <- c('Nigeria','Angola','Ethiopia','Indonesia','D.R.C.')
shp_adm0 <- shp_adm0[shp_adm0$ADM0_NAME %in% countries,]
shp_adm2_0 <- shp_adm2_0[shp_adm2_0$ADM0_NAME %in% countries,]
###########################################################################################
#Loop over metrics
for(cov_num in 1:length(covariate_names_full)){
  
  remove(g,g1,legend,legend1,df,df_c,a)
  cov_name<-covariate_names_abbrev[cov_num]
  print(cov_name)
  #read in data for metric
  df<-fread(paste0('<FILEPATH>/zerodose_overlap_analyses/',run_date,'/',covariate_names_abbrev[cov_num],'/all_admin_2.csv'))
  df[ADM0_NAME=='Democratic Republic of the Congo',]$ADM0_NAME <- 'D.R.C.'
  
  df<-df[(ADM0_NAME %in% countries)]
  
  #subset shapefile to countries used here
  shp_adm2<-shp_adm2_0[shp_adm2_0$ADM0_NAME %in% unique(df$ADM0_NAME),]
  
  df<-df[dtp0_pop!=0]
  
  #Rename covariates for uniformity in coding
  df$cov <- df[[covariate_names_abbrev[cov_num]]]
  df[[covariate_names_abbrev[cov_num]]] <- NULL
  df$cov_pct <- df[[paste0(covariate_names_abbrev[cov_num],'_pct')]]
  df[[paste0(covariate_names_abbrev[cov_num],'_pct')]] <- NULL
  
  #If covariate is ORS or Bednets, need to calculate inverted values to get unmet need
  if(cov_name=='ort')df[,cov_pct:=cov_pct/100]
  if(cov_name %in% c('ort','itn','ri_mcv1')){
    df[,cov:=(cov/(cov_pct))-cov]
    df[,cov_pct:=1-cov_pct]
  }
  if(cov_name=='itn')df[is.na(cov) & cov_pct==1,cov:=cov_pop]
  
  df_all_years<-copy(df)
  
  df<-df[year==cov_years[cov_num]]
  
  #Get district order for cov counts (sco)
  setorderv(df,'cov',order=-1)
  df[,sco:=1:.N]
  
  #Get district order for cov pct (spo)
  setorderv(df,'cov_pct',order=-1)
  df[,spo:=1:.N]
  
  #Get district order for zerodose counts (zco)
  setorderv(df,'dtp0',order=-1)
  df[,zco:=1:.N]
  
  #Get district order for zerodose pct (zpo)
  setorderv(df,'dtp0_pct',order=-1)
  df[,zpo:=1:.N]
  ######################################################################################################################
  #Make a temp legend for high priority plots
  #Make fake dataset for legend
  fake<-data.table(ADM0_NAME=c('Ethiopia','Nigeria','Angola'),
                   comb=c('Both','No-DTP only',paste0(paste(covariate_additions[cov_num], covariate_names_full[cov_num]),' only')),
                   color=as.factor(c('#3f2949','#DA434A','#5AA1E9')),order=c(3,1,2))
  fake<-merge(shp_adm0,fake,by='ADM0_NAME')
  f1<-ggplot()+
    geom_sf(data=fake,aes(fill=reorder(comb,order)))+
    scale_fill_manual(values=c('#DA434A','#5AA1E9','#3f2949'))+
    labs(fill='')
  
  f1 <- cowplot::get_legend(f1)
  ##############################################################################################
  #Global high priority plots
  if(make_global_high_priority_plots){
    print('making global high priority plots')
    tic()
    #Plot quantile overlap: percentage space
    df[, zd_pct_quart := cut(dtp0_pct, breaks = wtd.quantile(dtp0_pct, dtp0_pop, c(0, 0.25, 0.5, 0.75, 1), na.rm = T), labels = F, include.lowest = T)]
    df[, cov_pct_quart := cut(cov_pct, breaks = wtd.quantile(cov_pct, cov_pop, c(0, 0.25, 0.5, 0.75, 1), na.rm = T), labels = F, include.lowest = T)]
    
    # #Fix some NAs
    df[is.na(zd_pct_quart) & zpo==1,zd_pct_quart:=4]
    df[is.na(cov_pct_quart) & spo==1,cov_pct_quart:=4]
    df[is.na(zd_pct_quart) & zpo==max(df$zpo),zd_pct_quart:=1]
    df[is.na(cov_pct_quart) & spo==max(df$spo),cov_pct_quart:=1]
    df[dtp0_pct==max(df$dtp0_pct),zd_pct_quart:=4]
    df[cov_pct==max(df$cov_pct),cov_pct_quart:=4]
    
    
    
    levels <- CJ(cov_quart = unique(df$cov_pct_quart),
                 zd_quart = unique(df$zd_pct_quart))
    levels[, comb := factor(paste(cov_quart, zd_quart))]
    
    
    
    
    
    df[, comb := factor(paste(cov_pct_quart, zd_pct_quart), levels = levels(levels$comb))]
    df[!(comb %in% c('1 4','2 4','3 4','4 4',
                     '4 1','4 2', '4 3')),comb:=NA]
    df[comb == '4 4', comb:='Both']
    df[comb %in% c('1 4','2 4','3 4'),comb:='No-DTP only']
    df[comb %in% c('4 1','4 2','4 3'),comb:=paste0(paste(covariate_additions[cov_num], covariate_names_full[cov_num]),' only')]
    df[comb=='Both',order:=3]
    df[comb=='No-DTP only',order:=1]
    df[comb==paste0(paste(covariate_additions[cov_num], covariate_names_full[cov_num]),' only'),order:=2]
    
    #merge with shapefile data
    df1<-merge(shp_adm2,df, by=c('ADM0_CODE','ADM2_CODE','ADM1_CODE','ADM0_NAME','ADM2_NAME','ADM1_NAME'))
    
    
    #Get overlap %
    global_prev_overlap<-nrow(df[comb =='Both'])/nrow(df[comb %in% c ('Both','No-DTP only')])
    print('global prev overlap: ')
    print(global_prev_overlap)
    
    #Plot quantile overlap: count space
    df[, zd_quart := cut(dtp0, breaks = wtd.quantile(dtp0, dtp0_pop, c(0, 0.25, 0.5, 0.75, 1), na.rm = T), labels = F, include.lowest = T)]
    df[, cov_quart := cut(cov, breaks = wtd.quantile(cov,cov_pop,  c(0, 0.25, 0.5, 0.75, 1), na.rm = T), labels = F, include.lowest = T)]
    
    df[is.na(zd_quart),zd_quart:=1]
    df[is.na(cov_quart),cov_quart:=1]
    df[dtp0==max(df$dtp0),zd_quart:=4]
    df[cov==max(df$cov),cov_quart:=4]
    
    levels <- CJ(cov_quart = unique(df$cov_quart),
                 zd_quart = unique(df$zd_quart))
    levels[, comb := factor(paste(cov_quart, zd_quart))]
    
    
    
    
    df[, comb := factor(paste(cov_quart, zd_quart), levels = levels(levels$comb))]
    df[,cr:=NULL]
    df[!(comb %in% c('1 4','2 4','3 4','4 4',
                     '4 1','4 2', '4 3')),comb:=NA]
    df[comb == '4 4', comb:='Both']
    df[comb %in% c('1 4','2 4','3 4'),comb:='No-DTP only']
    df[comb %in% c('4 1','4 2','4 3'),comb:=paste0(paste(covariate_additions[cov_num], covariate_names_full[cov_num]),' only')]
    df[comb=='Both',order:=3]
    df[comb=='No-DTP only',order:=1]
    df[comb==paste0(paste(covariate_additions[cov_num], covariate_names_full[cov_num]),' only'),order:=2]
    
    
    df$comb<-droplevels(df$comb)
    
    
    #Get overlap %
    global_counts_overlap<-nrow(df[comb =='Both'])/nrow(df[comb %in% c ('Both','No-DTP only')])
    print('global counts overlap: ')
    print(global_counts_overlap)
    
    
    #merge with shapefile data
    df2<-merge(shp_adm2,df, by=c('ADM0_CODE','ADM2_CODE','ADM1_CODE','ADM0_NAME','ADM2_NAME','ADM1_NAME'))
    
    g1<-list()
    #plot
    for(z in 1:length(countries)){
      c=countries[z]
      g1[[z]]<-ggplot()+
        geom_sf(data=df1[df1$ADM0_NAME==c & !is.na(df1$comb),],aes(fill=reorder(comb,order)),color=NA)+
        geom_sf(data=df1[df1$ADM0_NAME==c & is.na(df1$comb),],fill='dark gray',color=NA)+
        scale_fill_manual(values=c('#DA434A','#5AA1E9','#3f2949')[sort(unique(df1[df1$ADM0_NAME==c & !is.na(df1$comb),]$order))])+
        geom_sf(data=shp_adm0[shp_adm0$ADM0_NAME==c,],fill=NA,color='black',size=0.13)+
        theme_void()+
        theme(text=element_text(size=8))+
        labs(fill='')
      
    }
    
    
    #save legend for country-specific plots
    legend<-cowplot::get_legend(g1[[1]]) #Take from nigeria
    
    for(z in c(1:length(countries))){
      g1[[z]] <-g1[[z]]+theme(legend.position = 'none')
    }
    
    a<- plot_grid(g1[[1]],g1[[2]],g1[[3]],g1[[4]],g1[[5]],f1, labels = c(countries,''), label_size = 7, nrow = 2)
    # a<- plot_grid(plotlist = g1, labels = c(countries), label_size = 7, nrow = 2)
    # a<- plot_grid(a,legend, labels = c('',''), label_size = 7, ncol = 2,rel_widths = c(4,1))
    #Combine and save both maps
    png(paste0('<FILEPATH>/zerodose_overlap_analyses/',run_date,'/results/paper/global_priority_plots/global_overlap_high_priority_prev_',covariate_names_abbrev[cov_num],'.png'),
        width = 8.5, height = 5, res=800, units = 'in')
    print(a)
    dev.off()
    
    g1<-list()
    #plot
    for(z in 1:length(countries)){
      c=countries[z]
      g1[[z]]<-ggplot()+
        geom_sf(data=df2[df2$ADM0_NAME==c & !is.na(df2$comb),],aes(fill=reorder(comb,order)),color=NA)+
        geom_sf(data=df2[df2$ADM0_NAME==c & is.na(df2$comb),],fill='dark gray',color=NA)+
        scale_fill_manual(values=c('#DA434A','#5AA1E9','#3f2949'))+
        geom_sf(data=shp_adm0[shp_adm0$ADM0_NAME==c,],fill=NA,color='black',size=0.13)+
        theme_void()+
        theme(text=element_text(size=8))+
        labs(fill='')
      
    }
    
    
    
    for(z in c(1:length(countries))){
      g1[[z]] <-g1[[z]]+theme(legend.position = 'none')
    }
    
    a<- plot_grid(g1[[1]],g1[[2]],g1[[3]],g1[[4]],g1[[5]],legend, labels = c(countries,''), label_size = 7, nrow = 2)
    
    # a<- plot_grid(plotlist = g1, labels = c(countries), label_size = 7, nrow = 2)
    # a<- plot_grid(a,legend, labels = c('',''), label_size = 7, ncol = 2,rel_widths = c(4,1))
    #Combine and save both maps
    png(paste0('<FILEPATH>/zerodose_overlap_analyses/',run_date,'/results/paper/global_priority_plots/global_overlap_high_priority_counts_',covariate_names_abbrev[cov_num],'.png'),
        width = 8.5, height = 5, res=800, units = 'in')
    print(a)
    dev.off()
    
    #Fill global data tables for combining metrics data later
    s<-rbind(s, data.table(overlap=global_prev_overlap,type='Prevalence',cov=covariate_names_abbrev[cov_num]))
    s<-rbind(s, data.table(overlap=global_counts_overlap,type='Counts',cov=covariate_names_abbrev[cov_num]))
    toc()
  } #make high priority national plot
  ###################################################################################################################
  #Get country-specific priorities
  if(make_national_high_priority_plots){
    
    print('making national high priority plots')
    tic()
    df[,zd_pct_quart:=NULL]
    df[,cov_pct_quart:=NULL]
    df[,zd_quart:=NULL]
    df[,cov_quart:=NULL]
    df[,comb:=NULL]
    df[,cr:=NULL]
    df[,cr.x:=NULL]
    df[,cr.y:=NULL]

    for(c in sort(unique(df$ADM0_NAME))){
      print(c)
      df_c<-df[ADM0_NAME==c]
      
      try(df_c[, zd_pct_quart := cut(dtp0_pct, breaks = wtd.quantile(dtp0_pct, dtp0_pop, c(0, 0.25, 0.5, 0.75, 1), na.rm = T), labels = F, include.lowest = T)],silent=T)
      # #if pop-weighted means dont work, use unweighted
      if(is.null(df_c$zd_pct_quart)) 
        df_c[, zd_pct_quart := cut(dtp0_pct, breaks = quantile(dtp0_pct, c(0, 0.25, 0.5, 0.75, 1), na.rm = T), labels = F, include.lowest = T)]
      try(df_c[, cov_pct_quart := cut(cov_pct, breaks = wtd.quantile(cov_pct, cov_pop, c(0, 0.25, 0.5, 0.75, 1), na.rm = T), labels = F, include.lowest = T)],silent=T)
      if(is.null(df_c$cov_pct_quart)) 
        try(df_c[, cov_pct_quart := cut(cov_pct, breaks = quantile(cov_pct, c(0, 0.25, 0.5, 0.75, 1), na.rm = T), labels = F, include.lowest = T)],silent=T)
      
      df_c[is.na(zd_pct_quart) & dtp0_pct==max(df_c$dtp0_pct),zd_pct_quart:=4]
      df_c[is.na(cov_pct_quart) & cov_pct==max(df_c$cov_pct),cov_pct_quart:=4]

      # #Correct some NAs
      df_c[is.na(zd_pct_quart) & zpo==1,zd_pct_quart:=4]
      df_c[is.na(cov_pct_quart) & spo==1,cov_pct_quart:=4]
      df_c[is.na(zd_pct_quart) & zpo==max(df_c$zpo),zd_pct_quart:=1]
      df_c[is.na(cov_pct_quart) & spo==max(df_c$spo),cov_pct_quart:=1]
      df_c[dtp0_pct==max(df_c$dtp0_pct),zd_pct_quart:=4]
      df_c[cov_pct==max(df_c$cov_pct),cov_pct_quart:=4]
      
      
      
      levels <- CJ(cov_quart = unique(df_c$cov_pct_quart),
                   zd_quart = unique(df_c$zd_pct_quart))
      levels[, comb := factor(paste(cov_quart, zd_quart))]
      
      
      
      
      df_c[, comb := factor(paste(cov_pct_quart, zd_pct_quart), levels = levels(levels$comb))]
      df_c[!(comb %in% c('1 4','2 4','3 4','4 4',
                         '4 1','4 2', '4 3')),comb:=NA]
      df_c[comb == '4 4', comb:='Both']
      df_c[comb %in% c('1 4','2 4','3 4'),comb:='No-DTP only']
      df_c[comb %in% c('4 1','4 2','4 3'),comb:=paste0(paste(covariate_additions[cov_num], covariate_names_full[cov_num]),' only')]
      df_c[comb=='Both',order:=3]
      df_c[comb=='No-DTP only',order:=1]
      df_c[comb==paste0(paste(covariate_additions[cov_num], covariate_names_full[cov_num]),' only'),order:=2]
      
      
      natl_prev_overlap = nrow(df_c[comb=='Both'])/nrow(df_c[comb %in% c('Both','No-DTP only')])
      
      df_c1 <- copy(df_c)
      df_c1$comb<-droplevels(df_c1$comb)
      df_c1$color<-as.factor(df_c1$comb)
      df_c1[comb=='Both',color:='#3f2949']
      df_c1[comb=='No-DTP only',color:='#DA434A']
      df_c1[!is.na(comb) & !(comb %in% c('Both','No-DTP only')),color:='#5AA1E9']
      df_c1$color<-droplevels(df_c1$color)
      df_c1<-merge(shp_adm2[shp_adm2$ADM0_NAME==c,],df_c1,by=c('ADM0_CODE','ADM2_CODE','ADM1_CODE','ADM0_NAME','ADM2_NAME','ADM1_NAME'))
      
      df_c[,comb:=NULL]
      
      #Calculate quantile overlap: count space
      try(df_c[, zd_quart := cut(dtp0, breaks = wtd.quantile(dtp0, dtp0_pop, c(0, 0.25, 0.5, 0.75, 1), na.rm = T), labels = F, include.lowest = T)],silent=T)
      # #if pop-weighted means dont work, use unweighted
      if(is.null(df_c$zd_quart)) 
        df_c[, zd_quart := cut(dtp0, breaks = quantile(dtp0, c(0, 0.25, 0.5, 0.75, 1), na.rm = T), labels = F, include.lowest = T)]
      try(df_c[, cov_quart := cut(cov, breaks = wtd.quantile(cov, cov_pop, c(0, 0.25, 0.5, 0.75, 1), na.rm = T), labels = F, include.lowest = T)],silent=T)
      if(is.null(df_c$cov_quart)) 
        try(df_c[, cov_quart := cut(cov, breaks = quantile(cov, c(0, 0.25, 0.5, 0.75, 1), na.rm = T), labels = F, include.lowest = T)],silent=T)

      df_c[is.na(zd_quart) & dtp0_pop < 1, zd_quart:=1]
      df_c[is.na(cov_quart) & cov_pop < 1,cov_quart:=1]
      df_c[is.na(zd_quart) & dtp0==max(df_c$dtp0),zd_quart:=4]
      df_c[is.na(cov_quart) & cov==max(df_c$cov),cov_quart:=4]
      df_c[dtp0==max(df_c$dtp0),zd_quart:=4]
      df_c[cov==max(df_c$cov),cov_quart:=4]
      
      df_c[is.na(zd_quart),zd_quart:=1]
      df_c[is.na(cov_quart),cov_quart:=1]
      
      levels <- CJ(cov_quart = unique(df_c$cov_quart),
                   zd_quart = unique(df_c$zd_quart))
      levels[, comb := factor(paste(cov_quart, zd_quart))]
      
      
      
      
      df_c[, comb := factor(paste(cov_quart, zd_quart), levels = levels(levels$comb))]
      df_c[!(comb %in% c('1 4','2 4','3 4','4 4',
                         '4 1','4 2', '4 3')),comb:=NA]
      df_c[comb == '4 4', comb:='Both']
      df_c[comb %in% c('1 4','2 4','3 4'),comb:='No-DTP only']
      df_c[comb %in% c('4 1','4 2','4 3'),comb:=paste0(paste(covariate_additions[cov_num], covariate_names_full[cov_num]),' only')]
      df_c[comb=='Both',order:=3]
      df_c[comb=='No-DTP only',order:=1]
      df_c[comb==paste0(paste(covariate_additions[cov_num], covariate_names_full[cov_num]),' only'),order:=2]
      
      
      natl_counts_overlap = nrow(df_c[comb=='Both'])/nrow(df_c[comb %in% c('Both','No-DTP only')])
      
      sc<-rbind(sc, data.table(overlap=natl_prev_overlap,type='Prevalence',cov=covariate_names_abbrev[cov_num],country=c))
      sc<-rbind(sc, data.table(overlap=natl_counts_overlap,type='Counts',cov=covariate_names_abbrev[cov_num],country=c))
      
      
      df_c2 <- copy(df_c)
      df_c2$comb<-droplevels(df_c2$comb)
      df_c2$color<-as.factor(df_c2$comb)
      df_c2[comb=='Both',color:='#3f2949']
      df_c2[comb=='No-DTP only',color:='#DA434A']
      df_c2[!is.na(comb) & !(comb %in% c('Both','No-DTP only')),color:='#5AA1E9']
      df_c2$color<-droplevels(df_c2$color)
      df_c2<-merge(shp_adm2[shp_adm2$ADM0_NAME==c,],df_c2,by=c('ADM0_CODE','ADM2_CODE','ADM1_CODE','ADM0_NAME','ADM2_NAME','ADM1_NAME'))
      
      
      g1<-ggplot()+
        geom_sf(data=shp_adm0[shp_adm0$ADM0_NAME==c,],fill='light gray',color=NA)+
        geom_sf(data=df_c1[is.na(df_c1$comb),],fill='dark gray',color='white',size=0.1)
      for(i in 1:length(levels(df_c1$color))){
        g1<-g1+
          geom_sf(data=df_c1[df_c1$comb==levels(df_c1$comb)[i] & !is.na(df_c1$comb),],
                  fill=levels(df_c1$color)[i],color='white',size=0.1)
      }
      g1<-g1+
        geom_sf(data=shp_adm0[shp_adm0$ADM0_NAME==c,],fill=NA,color='black',size=0.13)+
        theme_void()+
        theme(text=element_text(size=8))+
        labs(title='(A) By prevalence (%)',fill='')
      
      
      g2<-ggplot()+
        geom_sf(data=shp_adm0[shp_adm0$ADM0_NAME==c,],fill='light gray',color=NA)+
        geom_sf(data=df_c2[is.na(df_c2$comb),],fill='dark gray',color='white',size=0.1)
      for(i in 1:length(levels(df_c2$color))){
        g2<-g2+
          geom_sf(data=df_c2[df_c2$comb==levels(df_c2$comb)[i] & !is.na(df_c2$comb),],
                  fill=levels(df_c2$color)[i],color='white',size=0.1)
      }
      g2<-g2+
        geom_sf(data=shp_adm0[shp_adm0$ADM0_NAME==c,],fill=NA,color='black',size=0.13)+
        theme_void()+
        theme(text=element_text(size=8))+
        labs(title='(B) By counts',fill='')
      
      
      png(paste0('<FILEPATH>/zerodose_overlap_analyses/',run_date,'/results/paper/national_priority_plots/high_priority_overlap_by_country_',
                 covariate_names_abbrev[cov_num],'_',c,'.png'),
          width = 9, height = 5.4,units='in',res=600)
      
      print(grid.arrange(g1,g2, f1,ncol=3, widths=c(0.4,0.4,0.2)))
      
      dev.off()
    }
    
    
  } #make high priority national plots
  #############################################################################################################
  #Maps for full priority: global
  if(make_global_all_priority_plots){
    
    print('making global all priority plots')
    tic()
    #Plot quantile overlap: percentage space
    df[, zd_pct_quart := cut(dtp0_pct, breaks = wtd.quantile(dtp0_pct, dtp0_pop, c(0, 0.25, 0.5, 0.75, 1), na.rm = T), labels = F, include.lowest = T)]
    df[, cov_pct_quart := cut(cov_pct, breaks = wtd.quantile(cov_pct, cov_pop, c(0, 0.25, 0.5, 0.75, 1), na.rm = T), labels = F, include.lowest = T)]
    
    #Fix some NAs
    df[is.na(zd_pct_quart) & zpo==1,zd_pct_quart:=4]
    df[is.na(cov_pct_quart) & spo==1,cov_pct_quart:=4]
    df[is.na(zd_pct_quart) & zpo==max(df$zpo),zd_pct_quart:=1]
    df[is.na(cov_pct_quart) & spo==max(df$spo),cov_pct_quart:=1]
    
    zd_pct_breaks<-round(wtd.quantile(df$dtp0_pct, df$dtp0_pop, c(0, 0.25, 0.5, 0.75, 1), na.rm = T)*100)
    cov_pct_breaks<-round(wtd.quantile(df$cov_pct, df$cov_pop, c(0, 0.25, 0.5, 0.75, 1), na.rm = T)*100)
    
    levels <- CJ(cov_quart = unique(df$cov_pct_quart),
                 zd_quart = unique(df$zd_pct_quart))
    levels[, comb := factor(paste(cov_quart, zd_quart))]
    levels[,cr:=colors]
    
    
    legend <- ggplot(levels) +
      geom_raster(aes(x = factor(zd_quart), y = factor(cov_quart), fill = comb), show.legend = F) +
      scale_fill_manual(values = levels$cr) +
      scale_x_discrete(labels = c(paste0(zd_pct_breaks[1],'-',zd_pct_breaks[2]),
                                  paste0(zd_pct_breaks[2],'-',zd_pct_breaks[3]),
                                  paste0(zd_pct_breaks[3],'-',zd_pct_breaks[4]),
                                  paste0(zd_pct_breaks[4],'-',zd_pct_breaks[5])), expand = c(0, 0)) +
      scale_y_discrete(labels = c(paste0(cov_pct_breaks[1],'-',cov_pct_breaks[2]),
                                  paste0(cov_pct_breaks[2],'-',cov_pct_breaks[3]),
                                  paste0(cov_pct_breaks[3],'-',cov_pct_breaks[4]),
                                  paste0(cov_pct_breaks[4],'-',cov_pct_breaks[5])), expand = c(0, 0)) +
      coord_equal() +
      labs(x = ' \nNo-DTP (%)', y = paste(covariate_additions[cov_num],
                                          covariate_names_full[cov_num],' (%)\n '), title = NULL) +
      theme_minimal() +
      theme(axis.line = element_blank(), axis.ticks = element_blank(),
            panel.grid = element_blank(), axis.text.y = element_text(angle = 90, hjust = 0.5),
            plot.margin = unit(c(0, 0, 0, 0), "in"), panel.border = element_rect(fill = NA, color = 'black'))+
      theme(text=element_text(size=6))+
      theme(axis.text.y=element_text(angle=45, hjust=1))+
      theme(axis.text.x=element_text(angle=45, hjust=1))
    
    
    df[, comb := factor(paste(cov_pct_quart, zd_pct_quart), levels = levels(levels$comb))]
    df<-merge(df,levels[,list(comb,cr)],by='comb')
    
    #merge with shapefile data
    df1<-merge(shp_adm2,df, by=c('ADM0_CODE','ADM2_CODE','ADM1_CODE','ADM0_NAME','ADM2_NAME','ADM1_NAME'))
    g1<-list()
    for(z in 1:length(countries)){
      c=countries[z]
      g<-ggplot()+
        geom_sf(data=shp_adm0[shp_adm0$ADM0_NAME==c,],fill='gray',color=NA)+
        theme_void()+
        theme(text=element_text(size=8))
      for(i in 1:16){
        g<-g+
          geom_sf(data=df1[df1$ADM0_NAME==c & df1$comb==levels$comb[i],],fill=levels$cr[i],color=NA, show.legend = F)
      }
      g<-g+
        geom_sf(data=shp_adm0[shp_adm0$ADM0_NAME==c,],fill=NA,color='black',size=0.13)
      g1[[z]]<-g
    }
    
    a<- plot_grid(g1[[1]],g1[[2]],g1[[3]],g1[[4]],g1[[5]],legend, labels = c(countries,''), label_size = 7, nrow = 2)
    
    # a<- plot_grid(plotlist = g1, labels = c(countries), label_size = 7, nrow = 2)
    # a<- plot_grid(a,legend, labels = c('',''), label_size = 7, ncol = 2,rel_widths = c(4,1))
    png(paste0('<FILEPATH>/zerodose_overlap_analyses/',run_date,'/results/paper/global_priority_plots/global_overlap_all_priority_prev_',covariate_names_abbrev[cov_num],'.png'),
        width = 8.5, height = 5, res=800, units = 'in')
    print(a)
    dev.off()
    #Plot quantile overlap: count space
    df[, zd_quart := cut(dtp0, breaks = wtd.quantile(dtp0, dtp0_pop, c(0, 0.25, 0.5, 0.75, 1), na.rm = T), labels = F, include.lowest = T)]
    df[, cov_quart := cut(cov, breaks = wtd.quantile(cov,cov_pop,  c(0, 0.25, 0.5, 0.75, 1), na.rm = T), labels = F, include.lowest = T)]
    
    df[is.na(zd_quart),zd_quart:=1]
    df[is.na(cov_quart),cov_quart:=1]
    
    levels <- CJ(cov_quart = unique(df$cov_quart),
                 zd_quart = unique(df$zd_quart))
    levels[, comb := factor(paste(cov_quart, zd_quart))]
    levels[,cr:=colors]
    
    zd_breaks<-wtd.quantile(df$dtp0, df$dtp0_pop, c(0, 0.25, 0.5, 0.75, 1), na.rm = T)
    cov_breaks<-wtd.quantile(df$cov, df$cov_pop, c(0, 0.25, 0.5, 0.75, 1), na.rm = T)
    
    try(zd_breaks[zd_breaks<10] <- round_any(zd_breaks[zd_breaks<10],1),silent=T)
    try(zd_breaks[zd_breaks<100 & zd_breaks>=10] <- round_any(zd_breaks[zd_breaks<100 & zd_breaks>=10],5),silent=T)
    try(zd_breaks[zd_breaks<1000 & zd_breaks>=100] <- round_any(zd_breaks[zd_breaks<1000 & zd_breaks>=100],10),silent=T)
    #try(zd_breaks[zd_breaks>=1000] <- round_any(zd_breaks[zd_breaks>=1000],100),silent=T)
    try(zd_breaks[zd_breaks>=1000] <- round_any(zd_breaks[zd_breaks>=1000],1000),silent=T)
    try(zd_breaks[zd_breaks>=1000] <- paste0(zd_breaks[zd_breaks>=1000]/1000, 'k'))
    
    
    try(cov_breaks[cov_breaks<10] <- round_any(cov_breaks[cov_breaks<10],1),silent=T)
    try(cov_breaks[cov_breaks<100 & cov_breaks>=10] <- round_any(cov_breaks[cov_breaks<100 & cov_breaks>=10],5),silent=T)
    try(cov_breaks[cov_breaks<1000 & cov_breaks>=100] <- round_any(cov_breaks[cov_breaks<1000 & cov_breaks>=100],10),silent=T)
    #try(cov_breaks[cov_breaks>=1000] <- round_any(cov_breaks[cov_breaks>=1000],100),silent=T)
    try(cov_breaks[cov_breaks>=1000] <- round_any(cov_breaks[cov_breaks>=1000],1000),silent=T)
    
    cov_breaks_less_mill<-cov_breaks[cov_breaks<1000000]
    cov_breaks_greater_mill<-cov_breaks[cov_breaks>=1000000]
    
    try(cov_breaks_less_mill[cov_breaks_less_mill>=1000 & cov_breaks_less_mill<1000000] <- paste0(cov_breaks_less_mill[cov_breaks_less_mill>=1000 & cov_breaks_less_mill<1000000]/1000,'k'))
    try(cov_breaks_greater_mill[cov_breaks_greater_mill>=1000000] <- paste0(round(cov_breaks_greater_mill[cov_breaks_greater_mill>=1000000]/1000000,1),'m'))
    cov_breaks<-c(cov_breaks_less_mill,cov_breaks_greater_mill)
    
    legend1 <- ggplot(levels) +
      geom_raster(aes(x = factor(zd_quart), y = factor(cov_quart), fill = comb), show.legend = F) +
      scale_fill_manual(values = levels$cr) +
      scale_x_discrete(labels = c(paste0(zd_breaks[1],'-',zd_breaks[2]),
                                  paste0(zd_breaks[2],'-',zd_breaks[3]),
                                  paste0(zd_breaks[3],'-',zd_breaks[4]),
                                  paste0(zd_breaks[4],'-',zd_breaks[5])), expand = c(0, 0)) +
      scale_y_discrete(labels = c(paste0(cov_breaks[1],'-',cov_breaks[2]),
                                  paste0(cov_breaks[2],'-',cov_breaks[3]),
                                  paste0(cov_breaks[3],'-',cov_breaks[4]),
                                  paste0(cov_breaks[4],'-',cov_breaks[5])), expand = c(0, 0)) +
      coord_equal() +
      labs(x = ' \nNo-DTP (counts)', y = paste(covariate_additions[cov_num],
                                               covariate_names_full[cov_num],'(counts)\n '), title = NULL) +
      theme_minimal() +
      theme(axis.line = element_blank(), axis.ticks = element_blank(),
            panel.grid = element_blank(), axis.text.y = element_text(angle = 90, hjust = 0.5),
            plot.margin = unit(c(0, 0, 0, 0), "in"), panel.border = element_rect(fill = NA, color = 'black'))+
      theme(text=element_text(size=6))+
      theme(axis.text.y=element_text(angle=45, hjust=1))+
      theme(axis.text.x=element_text(angle=45, hjust=1))
    
    
    df[, comb := factor(paste(cov_quart, zd_quart), levels = levels(levels$comb))]
    df[,cr:=NULL]
    df<-merge(df,levels[,list(comb,cr)],by='comb')
    
    #merge with shapefile data
    df2<-merge(shp_adm2,df, by=c('ADM0_CODE','ADM2_CODE','ADM1_CODE','ADM0_NAME','ADM2_NAME','ADM1_NAME'))
    g1<-list()
    for(z in 1:length(countries)){
      c=countries[z]
      g<-ggplot()+
        geom_sf(data=shp_adm0[shp_adm0$ADM0_NAME==c,],fill='gray',color=NA)+
        theme_void()+
        theme(text=element_text(size=8))
      for(i in 1:16){
        g<-g+
          geom_sf(data=df2[df2$ADM0_NAME==c & df2$comb==levels$comb[i],],fill=levels$cr[i],color=NA, show.legend = F)
      }
      g<-g+
        geom_sf(data=shp_adm0[shp_adm0$ADM0_NAME==c,],fill=NA,color='black',size=0.13)
      g1[[z]]<-g
    }
    
    
    
    
    a<- plot_grid(g1[[1]],g1[[2]],g1[[3]],g1[[4]],g1[[5]],legend1, labels = c(countries,''), label_size = 7, nrow = 2)
    
   png(paste0('<FILEPATH>/zerodose_overlap_analyses/',run_date,'/results/paper/global_priority_plots/global_overlap_all_priority_counts_',covariate_names_abbrev[cov_num],'.png'),
        width = 8.5, height = 5, res=800, units = 'in')
    print(a)
    dev.off()
  } #make all priority global plots
  ########################################################################################################################################################
  #Country-specific all-quartile overlap plots
  if(make_national_all_priority_plots){
    
    print('making national all priority plots')
    tic()
    #Calculate and plot country-specific quantiles
    df[,zd_pct_quart:=NULL]
    df[,cov_pct_quart:=NULL]
    df[,zd_quart:=NULL]
    df[,cov_quart:=NULL]
    df[,comb:=NULL]
    df[,cr:=NULL]
    df[,cr.x:=NULL]
    df[,cr.y:=NULL]

    for(c in sort(unique(df$ADM0_NAME))){
      
      
      df_c<-df[ADM0_NAME==c]
      
      
      try(df_c[, zd_pct_quart := cut(dtp0_pct, breaks = wtd.quantile(dtp0_pct, dtp0_pop, c(0, 0.25, 0.5, 0.75, 1), na.rm = T), labels = F, include.lowest = T)],silent=T)
      # #if pop-weighted means dont work, use unweighted
      if(is.null(df_c$zd_pct_quart))
        df_c[, zd_pct_quart := cut(dtp0_pct, breaks = quantile(dtp0_pct, c(0, 0.25, 0.5, 0.75, 1), na.rm = T), labels = F, include.lowest = T)]
      try(df_c[, cov_pct_quart := cut(cov_pct, breaks = wtd.quantile(cov_pct, cov_pop, c(0, 0.25, 0.5, 0.75, 1), na.rm = T), labels = F, include.lowest = T)],silent=T)
      if(is.null(df_c$cov_pct_quart))
        try(df_c[, cov_pct_quart := cut(cov_pct, breaks = quantile(cov_pct, c(0, 0.25, 0.5, 0.75, 1), na.rm = T), labels = F, include.lowest = T)],silent=T)
      
       df_c[is.na(zd_pct_quart) & dtp0_pct==max(df_c$dtp0_pct),zd_pct_quart:=4]
      df_c[is.na(cov_pct_quart) & cov_pct==max(df_c$cov_pct),cov_pct_quart:=4]
      
      
      if(cov_num==1) df_c[ADM2_CODE==1015038 & is.na(cov_pct_quart), cov_pct_quart:=2]
      
      zd_pct_breaks<-wtd.quantile(df_c$dtp0_pct, df_c$dtp0_pop, c(0, 0.25, 0.5, 0.75, 1), na.rm = T)*100
      cov_pct_breaks<-wtd.quantile(df_c$cov_pct, df_c$cov_pop, c(0, 0.25, 0.5, 0.75, 1), na.rm = T)*100
      
      try(zd_pct_breaks[zd_pct_breaks<0.01] <- round(zd_pct_breaks[zd_pct_breaks<0.01],3),silent=T)
      try(zd_pct_breaks[zd_pct_breaks<0.1 & zd_pct_breaks>=0.01] <- round(zd_pct_breaks[zd_pct_breaks<0.1& zd_pct_breaks>=0.01],2),silent=T)
      try(zd_pct_breaks[zd_pct_breaks>=0.1 & zd_pct_breaks<1] <- round(zd_pct_breaks[zd_pct_breaks>=0.1 & zd_pct_breaks<1],1),silent=T)
      try(zd_pct_breaks[zd_pct_breaks>=1.0] <- round(zd_pct_breaks[zd_pct_breaks>=1],1),silent=T)
      
      try(cov_pct_breaks[cov_pct_breaks<0.01] <- round(cov_pct_breaks[cov_pct_breaks<0.01],3),silent=T)
      try(cov_pct_breaks[cov_pct_breaks<0.1 & cov_pct_breaks>=0.01] <- round(cov_pct_breaks[cov_pct_breaks<0.1 & cov_pct_breaks>=0.01],2),silent=T)
      try(cov_pct_breaks[cov_pct_breaks>=0.1 & cov_pct_breaks<1] <- round(cov_pct_breaks[cov_pct_breaks>=0.1 & cov_pct_breaks<1],1),silent=T)
      try(cov_pct_breaks[cov_pct_breaks>=1.0] <- round(cov_pct_breaks[cov_pct_breaks>=1],1),silent=T)
      
      
      levels <- CJ(cov_quart = 1:4,
                   zd_quart = 1:4)
      levels[, comb := factor(paste(cov_quart, zd_quart))]
      levels[,cr:=colors]
      
      
      legend <- ggplot(levels) +
        geom_raster(aes(x = factor(zd_quart), y = factor(cov_quart), fill = comb), show.legend = F) +
        scale_fill_manual(values = levels$cr) +
        scale_x_discrete(labels = c(paste0(0,'-',zd_pct_breaks[2]),
                                    paste0(zd_pct_breaks[2],'-',zd_pct_breaks[3]),
                                    paste0(zd_pct_breaks[3],'-',zd_pct_breaks[4]),
                                    paste0(zd_pct_breaks[4],'-',zd_pct_breaks[5])), expand = c(0, 0)) +
        scale_y_discrete(labels = c(paste0(0,'-',cov_pct_breaks[2]),
                                    paste0(cov_pct_breaks[2],'-',cov_pct_breaks[3]),
                                    paste0(cov_pct_breaks[3],'-',cov_pct_breaks[4]),
                                    paste0(cov_pct_breaks[4],'-',cov_pct_breaks[5])), expand = c(0, 0)) +
        coord_equal() +
        labs(x = ' \nNo-DTP (%)', y = paste(covariate_additions[cov_num], covariate_names_full[cov_num],
                                            ' (%)\n '), title = NULL) +
        theme_minimal() +
        theme(axis.line = element_blank(), axis.ticks = element_blank(),
              panel.grid = element_blank(), axis.text.y = element_text(angle = 90, hjust = 0.5),
              plot.margin = unit(c(0.2,0.2,0.2,0.2), "in"),
              panel.border = element_rect(fill = NA, color = 'black'))+
        theme(text=element_text(size=8))+
        theme(axis.text.y=element_text(angle=45, hjust=1))+
        theme(axis.text.x=element_text(angle=45, hjust=1))
      
      
      df_c[, comb := factor(paste(cov_pct_quart, zd_pct_quart), levels = levels(levels$comb))]
      df_c<-merge(df_c,levels[,list(comb,cr)],by='comb')
      
      
      g0<-ggplot()+
        theme_classic()+
        theme(legend.position = 'none',
              plot.margin = unit(c(0.1,0.1,0.1,0.5), "in"))+
        labs(x='No-DTP (%)',y=paste(covariate_additions[cov_num], covariate_names_full[cov_num],' (%)'))
      for(i in 1:16){
        g0<-g0+
          geom_point(data=df_c[df_c$comb==levels$comb[i],],aes(x=dtp0_pct*100,y=cov_pct*100,size=dtp0_pop),color=levels$cr[i])
      }
      g0<-g0+lims(x=c(0,max(df_c$dtp0_pct)*100),y=c(0,max(df_c$cov_pct)*100)) 
      
      #merge with shapefile data
      df_c1<-merge(shp_adm2,df_c, by=c('ADM0_CODE','ADM2_CODE','ADM1_CODE','ADM0_NAME','ADM2_NAME','ADM1_NAME'))
      
      g<-ggplot()+
        geom_sf(data=shp_adm0[shp_adm0$ADM0_NAME==c,],fill='gray',color=NA)+
        theme_void()+
        labs(title=paste0('(A) No-DTP and ',paste(tolower(covariate_additions[cov_num]), covariate_names_full[cov_num]),
                          ' quartiles (%): ', c))+
        theme(text=element_text(size=8))
      for(i in 1:16){
        g<-g+
          geom_sf(data=df_c1[df_c1$comb==levels$comb[i],],fill=levels$cr[i],color='light gray',size=0.02, show.legend = F)
      }
      g<-g+
        geom_sf(data=shp_adm0[shp_adm0$ADM0_NAME==c,],fill=NA,color='black',size=0.13)
      
      
      try(df_c[, zd_quart := cut(dtp0, breaks = wtd.quantile(dtp0, dtp0_pop, c(0, 0.25, 0.5, 0.75, 1), na.rm = T), labels = F, include.lowest = T)],silent=T)
      # #if pop-weighted means dont work, use unweighted
      if(is.null(df_c$zd_quart))
        df_c[, zd_quart := cut(dtp0, breaks = quantile(dtp0, c(0, 0.25, 0.5, 0.75, 1), na.rm = T), labels = F, include.lowest = T)]
      try(df_c[, cov_quart := cut(cov, breaks = wtd.quantile(cov, cov_pop, c(0, 0.25, 0.5, 0.75, 1), na.rm = T), labels = F, include.lowest = T)],silent=T)
      if(is.null(df_c$cov_quart))
        try(df_c[, cov_quart := cut(cov, breaks = quantile(cov, c(0, 0.25, 0.5, 0.75, 1), na.rm = T), labels = F, include.lowest = T)],silent=T)
      #Hotfix for Guyana ORT
      if(is.null(df_c$cov_quart) & c=='Guyana' & cov_num==3)
        try(df_c[, cov_quart := cut(cov, breaks = wtd.quantile(cov,cov_pop, c(0, 0.25, 0.5, 1), na.rm = T), labels = F, include.lowest = T)+1],silent=T)
      if(is.null(df_c$cov_quart) & c=='Mongolia' & cov_num==3)
        try(df_c[, cov_quart := cut(cov, breaks = wtd.quantile(cov,cov_pop, c(0, 0.5, 0.75, 1), na.rm = T), labels = F, include.lowest = T)+1],silent=T)
      
      df_c[is.na(zd_quart) & dtp0_pop < 1, zd_quart:=1]
      df_c[is.na(cov_quart) & cov_pop < 1,cov_quart:=1]
      df_c[is.na(zd_quart) & dtp0==max(df_c$dtp0),zd_quart:=4]
      df_c[is.na(cov_quart) & cov==max(df_c$cov),cov_quart:=4]
      
      
      levels <- CJ(cov_quart = 1:4,
                   zd_quart = 1:4)
      levels[, comb := factor(paste(cov_quart, zd_quart))]
      levels[,cr:=colors]
      
      zd_breaks<-wtd.quantile(df_c$dtp0, df_c$dtp0_pop, c(0, 0.25, 0.5, 0.75, 1), na.rm = T)
      cov_breaks<-wtd.quantile(df_c$cov, df_c$cov_pop, c(0, 0.25, 0.5, 0.75, 1), na.rm = T)
      
      
      
      try(zd_breaks[zd_breaks<10] <- round_any(zd_breaks[zd_breaks<10],1),silent=T)
      try(zd_breaks[zd_breaks<100 & zd_breaks>=10] <- round_any(zd_breaks[zd_breaks<100 & zd_breaks>=10],5),silent=T)
      try(zd_breaks[zd_breaks<1000 & zd_breaks>=100] <- round_any(zd_breaks[zd_breaks<1000 & zd_breaks>=100],10),silent=T)
      #try(zd_breaks[zd_breaks>=1000] <- round_any(zd_breaks[zd_breaks>=1000],100),silent=T)
      try(zd_breaks[zd_breaks>=1000] <- round_any(zd_breaks[zd_breaks>=1000],1000),silent=T)
      try(zd_breaks[zd_breaks>=1000] <- paste0(zd_breaks[zd_breaks>=1000]/1000,'k'))
      
      
      try(cov_breaks[cov_breaks<10] <- round_any(cov_breaks[cov_breaks<10],1),silent=T)
      try(cov_breaks[cov_breaks<100 & cov_breaks>=10] <- round_any(cov_breaks[cov_breaks<100 & cov_breaks>=10],5),silent=T)
      try(cov_breaks[cov_breaks<1000 & cov_breaks>=100] <- round_any(cov_breaks[cov_breaks<1000 & cov_breaks>=100],10),silent=T)
      #try(cov_breaks[cov_breaks>=1000] <- round_any(cov_breaks[cov_breaks>=1000],100),silent=T)
      try(cov_breaks[cov_breaks>=1000] <- round_any(cov_breaks[cov_breaks>=1000],1000),silent=T)
      cov_breaks_less_mill<-cov_breaks[cov_breaks<1000000]
      cov_breaks_greater_mill<-cov_breaks[cov_breaks>=1000000]
      
      try(cov_breaks_less_mill[cov_breaks_less_mill>=1000 & cov_breaks_less_mill<1000000] <- paste0(cov_breaks_less_mill[cov_breaks_less_mill>=1000 & cov_breaks_less_mill<1000000]/1000,'k'))
      try(cov_breaks_greater_mill[cov_breaks_greater_mill>=1000000] <- paste0(round(cov_breaks_greater_mill[cov_breaks_greater_mill>=1000000]/1000000,1),'m'))
      cov_breaks<-c(cov_breaks_less_mill,cov_breaks_greater_mill)
      
      
      legend1 <- ggplot(levels) +
        geom_raster(aes(x = factor(zd_quart), y = factor(cov_quart), fill = comb), show.legend = F) +
        scale_fill_manual(values = levels$cr) +
        scale_x_discrete(labels = c(paste0(0,'-',zd_breaks[2]),
                                    paste0(zd_breaks[2],'-',zd_breaks[3]),
                                    paste0(zd_breaks[3],'-',zd_breaks[4]),
                                    paste0(zd_breaks[4],'-',zd_breaks[5])), expand = c(0, 0)) +
        scale_y_discrete(labels = c(paste0(0,'-',cov_breaks[2]),
                                    paste0(cov_breaks[2],'-',cov_breaks[3]),
                                    paste0(cov_breaks[3],'-',cov_breaks[4]),
                                    paste0(cov_breaks[4],'-',cov_breaks[5])), expand = c(0, 0)) +
        coord_equal() +
        labs(x = ' \nNo-DTP (counts)', y = paste(covariate_additions[cov_num], covariate_names_full[cov_num],'(counts)\n '), title = NULL) +
        theme_minimal() +
        theme(axis.line = element_blank(), axis.ticks = element_blank(),
              panel.grid = element_blank(), axis.text.y = element_text(angle = 90, hjust = 0.5),
              plot.margin = unit(c(0.2,0.2,0.2,0.2), "in"),
              panel.border = element_rect(fill = NA, color = 'black'))+
        theme(text=element_text(size=8))+
        theme(axis.text.y=element_text(angle=45, hjust=1))+
        theme(axis.text.x=element_text(angle=45, hjust=1))
      
      
      df_c[, comb := factor(paste(cov_quart, zd_quart), levels = levels(levels$comb))]
      df_c[,cr:=NULL]
      df_c<-merge(df_c,levels[,list(comb,cr)],by='comb')
      
      #merge with shapefile data
      df_c2<-merge(shp_adm2,df_c, by=c('ADM0_CODE','ADM2_CODE','ADM1_CODE','ADM0_NAME','ADM2_NAME','ADM1_NAME'))
      
      #plot
      g1<-ggplot()+
        geom_sf(data=shp_adm0[shp_adm0$ADM0_NAME==c,],fill='gray',color=NA)+
        theme_void()+
        labs(title=paste0('(B) No-DTP and ',paste(tolower(covariate_additions[cov_num]), covariate_names_full[cov_num]),' quartiles (counts): ', c))+
        theme(text=element_text(size=8))
      for(i in 1:16){
        g1<-g1+
          geom_sf(data=df_c2[df_c2$comb==levels$comb[i],],fill=levels$cr[i],color='light gray',size=0.02, show.legend = F)
      }
      g1<-g1+
        geom_sf(data=shp_adm0[shp_adm0$ADM0_NAME==c,],fill=NA,color='black',size=0.13)
      
      
      g2<-ggplot()+
        theme_classic()+
        theme(legend.position = 'none',
              plot.margin = unit(c(0.1,0.1,0.1,0.2), "in"))+
        labs(x='No-DTP (counts)',y=paste(covariate_additions[cov_num], covariate_names_full[cov_num],'(counts)'))
      
      for(i in 1:16){
        g2<-g2+
          geom_point(data=df_c[df_c$comb==levels$comb[i],],aes(x=dtp0,y=cov,size=dtp0_pop),color=levels$cr[i])
      }
      g2<-g2+
        scale_y_continuous(labels = function(l) ifelse(l <= 9999, l, comma(l)) ,limits=c(0, max(df_c$cov)))+
        scale_x_continuous(labels = function(l) ifelse(l <= 9999, l, comma(l)),limits=c(0, max(df_c$dtp0)))+
        theme(axis.text.y=element_text(angle=45, hjust=1))+
        theme(axis.text.x=element_text(angle=45, hjust=1))
      
      png(paste0('<FILEPATH>/zerodose_overlap_analyses/',run_date,'/results/paper/national_priority_plots/all_overlap_by_country_',covariate_names_abbrev[cov_num],
                 '_',c,'.png'),
          width = 10, height = 7,units='in',res=600)
      
      print(grid.arrange(g,legend,g0,g1,legend1,g2, widths=c(0.40,0.25,0.35),nrow=2))
      dev.off()
      
    }
    
    
    toc()
  } #make all priority national plots
  ####################################################################################################
  #AUC curves
  if(make_auc_plots){
    
    print('making auc plots')
    tic()
    nd_all<-data.table()
    for(c in unique(df_all_years$ADM0_NAME)){
      for(y in c(2000, cov_years[cov_num])){
        newdata <- df_all_years[ADM0_NAME==c & year==y]
        
        #Find fewest districts with 25% of dtp0 kids
        newdata <- newdata[order(-dtp0),]
        
        newdata$cumsum <- cumsum(newdata$dtp0)
        
        total0 <- sum(newdata$dtp0, na.rm=TRUE)
        newdata$dtp0_percent <- newdata$cumsum / total0
        
        total1 <- sum(newdata$cov, na.rm=TRUE)
        
        
        
        nd<-newdata[,list(year=y,cov_count_value=cov,dtp0_count_value=dtp0,dtp0,
                          ADM0_NAME,ADM2_NAME,dtp0_total=total0, cov_total=total1, dtp0_pop)]
        
        nd_all<-rbind(nd_all,nd)
        
      }
    }
    nd_all[,dtp0_percent:=dtp0_count_value/dtp0_total*100]
    nd_all[,cov_percent:=cov_count_value/cov_total*100]
    
    for(c in unique(nd_all$ADM0_NAME)){
      for(y in c(2000,cov_years[cov_num])){
        nd_all[ADM0_NAME==c & year==y,dtp0_cumsum_percent:=cumsum(dtp0_percent)]
        nd_all[ADM0_NAME==c & year==y,cov_cumsum_percent:=cumsum(cov_percent)]
      }
    }
    
    nd_all[,dtp0_cumsum_percent:=round(dtp0_cumsum_percent,4)]
    nd_all[,cov_cumsum_percent:=round(cov_cumsum_percent,4)]
    
    for(c in unique(nd_all$ADM0_NAME)){
      for(y in c(2000,cov_years[cov_num])){
        x<- AUC(c(0,nd_all[ADM0_NAME==c & year==y]$dtp0_cumsum_percent/100), c(0,nd_all[ADM0_NAME==c & year==y]$cov_cumsum_percent/100), method='spline',from=0,to=1,
                subdivisions=2500)
        nd_all[ADM0_NAME==c & year==y,auc:=x]
      }
    }
    
    
    nd_map<-unique(nd_all[year==cov_years[cov_num],list(ADM0_NAME,auc)])
    nd_map1<-unique(nd_all[,list(ADM0_NAME,auc,year)])
    setorderv(nd_map,'auc',order=-1)
    nd_map[,rank_order:=1:.N]
    
    nd_all<-merge(nd_all[year==cov_years[cov_num]],nd_map[,list(ADM0_NAME,rank_order)],by='ADM0_NAME')
    

    for(c in sort(unique(nd_all$ADM0_NAME))){
      x<- round(AUC(c(0,nd_all[ADM0_NAME==c]$dtp0_cumsum_percent/100), c(0,nd_all[ADM0_NAME==c]$cov_cumsum_percent/100), method='spline',from=0,to=1,
                    subdivisions=2500),3)
      rank<-nd_all[ADM0_NAME==c]$rank_order[1]
      a<-  ggplot(data=nd_all[ADM0_NAME==c],aes(x=dtp0_cumsum_percent,y=cov_cumsum_percent))+geom_point(aes(size=dtp0_pop),shape=1)+
        theme_bw()+
        lims(x=c(0,100),y=c(0,100))+
        labs(title=paste0(c,': cumulative proportions reached with no-DTP targeting'),
             subtitle=paste0('AUC = ',x),
             x='Percent children with no-DTP (%)',y=paste0('Percent ', paste(tolower(covariate_additions[cov_num]), covariate_names_full[cov_num]), ' (%)'),size='< 1 population')+
        # geom_step(lty='dashed')+
        scale_size_continuous(labels=function(l) ifelse(l <= 9999, l, comma(l)))+
        geom_abline(slope=1,intercept = 0,color='red')+
        #    annotate('text',label="Equal targeting",x=30,y=50)+
        annotate('text',label=paste0("Greater targeting for ",paste(tolower(covariate_additions[cov_num]), covariate_names_full[cov_num])),x=30,y=95)+
        annotate('text',label=paste0("Worse targeting for ",paste(tolower(covariate_additions[cov_num]), covariate_names_full[cov_num])),x=75,y=5)+
        #geom_segment(x = 41, y = 50,xend = 45, yend = 50,lineend = "butt",linejoin = "bevel",size = .5,arrow = arrow(length = unit(0.1, "inches")))+
        geom_segment(x = 1, y = 93,xend = 1, yend = 97,lineend = "butt",linejoin = "bevel",size = .5,arrow = arrow(length = unit(0.1, "inches")))+
        geom_segment(x = 47, y = 7,xend = 47, yend = 3,lineend = "butt",linejoin = "bevel",size = .5,arrow = arrow(length = unit(0.1, "inches")))
      png(paste0('<FILEPATH>/zerodose_overlap_analyses/',run_date,'/results/paper/national_priority_plots/auc_step_plot_',covariate_names_abbrev[cov_num],'_',c,'.png'),
          width = 8, height = 7,units='in',res=600)
      
      print(a)
      dev.off()
    }
    
    
    nd_map1[,cov:=covariate_names_abbrev[cov_num]]
    
    nn<-rbind(nn, nd_map1)
    toc()
  } #make auc plots
  
  if(make_multi_metric_overlap_plots){
    print('making data for high priority overlap multi metric')
    tic()
    df[,zd_pct_quart:=NULL]
    df[,cov_pct_quart:=NULL]
    df[,zd_quart:=NULL]
    df[,cov_quart:=NULL]
    df[,comb:=NULL]
    df[,cr:=NULL]
    df[,cr.x:=NULL]
    df[,cr.y:=NULL]
    
    for(c in sort(unique(df$ADM0_NAME))){
      print(c)
      df_c<-df[ADM0_NAME==c]
      
      try(df_c[, zd_pct_quart := cut(dtp0_pct, breaks = wtd.quantile(dtp0_pct, dtp0_pop, c(0, 0.25, 0.5, 0.75, 1), na.rm = T), labels = F, include.lowest = T)],silent=T)
      # #if pop-weighted means dont work, use unweighted
      if(is.null(df_c$zd_pct_quart)) 
        df_c[, zd_pct_quart := cut(dtp0_pct, breaks = quantile(dtp0_pct, c(0, 0.25, 0.5, 0.75, 1), na.rm = T), labels = F, include.lowest = T)]
      try(df_c[, cov_pct_quart := cut(cov_pct, breaks = wtd.quantile(cov_pct, cov_pop, c(0, 0.25, 0.5, 0.75, 1), na.rm = T), labels = F, include.lowest = T)],silent=T)
      if(is.null(df_c$cov_pct_quart)) 
        try(df_c[, cov_pct_quart := cut(cov_pct, breaks = quantile(cov_pct, c(0, 0.25, 0.5, 0.75, 1), na.rm = T), labels = F, include.lowest = T)],silent=T)
      
      
      #hardcode a bug out
      
     df_c[is.na(zd_pct_quart) & dtp0_pct==max(df_c$dtp0_pct),zd_pct_quart:=4]
      df_c[is.na(cov_pct_quart) & cov_pct==max(df_c$cov_pct),cov_pct_quart:=4]
      if(cov_name=='stunting') df_c[ADM2_CODE==1015038 & is.na(cov_pct_quart), cov_pct_quart:=2]
      
      # #Fix some NAs
      df_c[is.na(zd_pct_quart) & zpo==1,zd_pct_quart:=4]
      df_c[is.na(cov_pct_quart) & spo==1,cov_pct_quart:=4]
      df_c[is.na(zd_pct_quart) & zpo==max(df_c$zpo),zd_pct_quart:=1]
      df_c[is.na(cov_pct_quart) & spo==max(df_c$spo),cov_pct_quart:=1]
      df_c[dtp0_pct==max(df_c$dtp0_pct),zd_pct_quart:=4]
      df_c[cov_pct==max(df_c$cov_pct),cov_pct_quart:=4]
      
      
      #Calculate quantile overlap: count space
      try(df_c[, zd_quart := cut(dtp0, breaks = wtd.quantile(dtp0, dtp0_pop, c(0, 0.25, 0.5, 0.75, 1), na.rm = T), labels = F, include.lowest = T)],silent=T)
      # #if pop-weighted means dont work, use unweighted
      if(is.null(df_c$zd_quart)) 
        df_c[, zd_quart := cut(dtp0, breaks = quantile(dtp0, c(0, 0.25, 0.5, 0.75, 1), na.rm = T), labels = F, include.lowest = T)]
      try(df_c[, cov_quart := cut(cov, breaks = wtd.quantile(cov, cov_pop, c(0, 0.25, 0.5, 0.75, 1), na.rm = T), labels = F, include.lowest = T)],silent=T)
      if(is.null(df_c$cov_quart)) 
        try(df_c[, cov_quart := cut(cov, breaks = quantile(cov, c(0, 0.25, 0.5, 0.75, 1), na.rm = T), labels = F, include.lowest = T)],silent=T)
     
      df_c[is.na(zd_quart) & dtp0_pop < 1, zd_quart:=1]
      df_c[is.na(cov_quart) & cov_pop < 1,cov_quart:=1]
      df_c[is.na(zd_quart) & dtp0==max(df_c$dtp0),zd_quart:=4]
      df_c[is.na(cov_quart) & cov==max(df_c$cov),cov_quart:=4]
      df_c[dtp0==max(df_c$dtp0),zd_quart:=4]
      df_c[cov==max(df_c$cov),cov_quart:=4]
      
      df_c[is.na(zd_quart),zd_quart:=1]
      df_c[is.na(cov_quart),cov_quart:=1]
      
      df_c<-df_c[,list(ADM0_NAME,ADM2_CODE,zd_pct_quart,cov_pct_quart,zd_quart,cov_quart)]
      df_c <- df_c[zd_pct_quart==4|cov_pct_quart==4|zd_quart==4|cov_quart==4]
      df_c[,cov:=covariate_names_abbrev[cov_num]]
      
      mc0<-rbind(mc0,df_c)
      
    }
  }
} #####loop over metrics
message('done with metric-specific plots. now run multi-metric plots')
if(make_combo_plots){
  tic()
  #############################################################################
  #Make some summarizing plots across metrics
  #############################################################################
  #Maps for how many metrics are high priority
  mc<-copy(mc0)
  
  mc[,zd_pct_quart:=as.character(zd_pct_quart)]
  mc[zd_pct_quart==4,zd_pct_quart:='high']
  mc[zd_pct_quart!='high',zd_pct_quart:='low']
  mc[,cov_pct_quart:=as.character(cov_pct_quart)]
  mc[cov_pct_quart==4,cov_pct_quart:='high']
  mc[cov_pct_quart!='high',cov_pct_quart:='low']
  
  mc[,zd_quart:=as.character(zd_quart)]
  mc[zd_quart==4,zd_quart:='high']
  mc[zd_quart!='high',zd_quart:='low']
  mc[,cov_quart:=as.character(cov_quart)]
  mc[cov_quart==4,cov_quart:='high']
  mc[cov_quart!='high',cov_quart:='low']
  
  mc1<-mc[cov=='stunting' & zd_pct_quart=='high']
  mc1[,cov:='no-DTP']
  mc1[,cov_pct_quart:=zd_pct_quart]
  mc1[,zd_pct_quart:=NULL]
  mc1[,zd_quart:=NULL]
  mc1[,cov_quart:=NULL]
  
  mc2<-mc[cov=='stunting' & zd_quart=='high']
  mc2[,cov:='no-DTP']
  mc2[,cov_quart:=zd_quart]
  mc2[,zd_quart:=NULL]
  mc2[,zd_pct_quart:=NULL]
  mc2[,cov_pct_quart:=NULL]
  
  mc3<-merge(mc1,mc2, by=c('ADM0_NAME','ADM2_CODE','cov'),all=T)
  mc3[is.na(cov_pct_quart),cov_pct_quart:='low']
  mc3[is.na(cov_quart),cov_quart:='low']
  
  mc<-mc[,list(ADM0_NAME,ADM2_CODE,cov_pct_quart,cov_quart,cov)]
  mc<-mc[cov_pct_quart=='high'|cov_quart=='high']
  
  mc<-rbind(mc,mc3)
  
  mc4<-as.data.table(table(mc[cov_pct_quart=='high',list(ADM2_CODE)]))
  mc4<-merge(shp_adm2_0,mc4,by='ADM2_CODE',all=T)
  
  
  g1<-list()
  #plot
  for(z in 1:length(countries)){
    c=countries[z]
    g1[[z]]<-ggplot()+
      geom_sf(data=mc4[mc4$ADM0_NAME==c,],aes(fill=N),color=NA)+
      geom_sf(data=shp_adm0[shp_adm0$ADM0_NAME==c,],fill=NA,color='black',size=0.13)+
      geom_sf(data=mc4[mc4$ADM0_NAME==c & mc4$ADM2_CODE %in% 
                         mc[ADM0_NAME==c & cov_pct_quart=='high' &
                              cov=='no-DTP']$ADM2_CODE,],color="#E7DADA",size=1,fill=NA)+
      #   geom_sf(data=mc4[mc4$ADM0_NAME==c & is.na(N),],fill='dark gray',color=NA)+
      scale_fill_viridis()+
      theme_void()+
      theme(text=element_text(size=8))+
      labs(fill='Number of\nindicators',title='(A) By prevalence')
    
  }
  
  mc5<-as.data.table(table(mc[cov_quart=='high',list(ADM2_CODE)]))
  mc5<-merge(shp_adm2_0,mc5,by='ADM2_CODE',all=T)
  
  
  g2<-list()
  #plot
  for(z in 1:length(countries)){
    c=countries[z]
    g2[[z]]<-ggplot()+
      geom_sf(data=mc5[mc5$ADM0_NAME==c,],aes(fill=N),color=NA)+
      geom_sf(data=shp_adm0[shp_adm0$ADM0_NAME==c,],fill=NA,color='black',size=0.13)+
      geom_sf(data=mc5[mc5$ADM0_NAME==c & mc5$ADM2_CODE %in% 
                         mc[ADM0_NAME==c & cov_quart=='high' &
                              cov=='no-DTP']$ADM2_CODE,],fill=NA,color="#E7DADA",size=1)+
      #   geom_sf(data=mc5[mc5$ADM0_NAME==c & is.na(N),],fill='dark gray',color=NA)+
      scale_fill_viridis()+
      theme_void()+
      theme(text=element_text(size=8))+
      labs(fill='Number of\nindicators',title='(B) By counts')
    
  }
  
  #Combine and save both maps
  
  for(z in 1:length(countries)){
    png(paste0('<FILEPATH>/zerodose_overlap_analyses/',run_date,'/results/paper/national_priority_plots/high_priority_metric_count_map_',countries[z],'.png'),
        width = 8, height = 5,units='in',res=600)
    grid.arrange(g1[[z]],g2[[z]],ncol=2)
    dev.off()
  }
  
  
  
  
  
  
  
  #############################################################################  

  sc[,cov:=as.factor(cov)]
  levels(sc$cov)<-c('Missed ITNs','LF','Missed ORT','Stunting','U5M')
  sc[cov=='Stunting',order:=1]
  sc[cov=='U5M',order:=2]
  sc[cov=='Missed ORT',order:=3]
  sc[cov=='Missed ITNs',order:=4]
  sc[cov=='LF',order:=5]
  sc[type=='Prevalence',order2:=1]
  sc[type=='Counts',order2:=2]
  
  
  
  nn[,cov:=as.factor(cov)]
  levels(nn$cov)<-c('Missed ITNs','LF','Missed ORT','Stunting','U5M')
  nn[cov=='Stunting',order:=1]
  nn[cov=='U5M',order:=2]
  nn[cov=='Missed ORT',order:=3]
  nn[cov=='Missed ITNs',order:=4]
  nn[cov=='LF',order:=5]
  
 
  
  library(ggpattern)  
  

  cbPalette <- c("#E69F00","#56B4E9","#009E73",
                 "#F0E442","#0072B2","#D55E00","#CC79A7")
  
  
  c<-ggplot(data=sc)+
    geom_bar_pattern(position='dodge',pattern_spacing=0.01,color='black',stat='identity',
                     aes(fill=reorder(cov,order), y=overlap*100,x=country,
                         pattern=reorder(type,order2)))+
    theme_bw()+
    scale_fill_manual(values=cbPalette)+
    scale_pattern_manual(values=c('none','stripe')) +
    labs(x='',y='Percent overlap (%)',fill='Indicator',pattern='Prioritized by:')+
    guides(fill = guide_legend(override.aes = list(pattern='none') ),
           pattern = guide_legend(override.aes = list(fill = "white")))
  
  png(paste0('<FILEPATH>/zerodose_overlap_analyses/',run_date,'/results/paper/national_priority_plots/prop_overlap_all_countries_by_metrics_bars.png'),
      width = 10, height = 5, res=800, units = 'in')
  print(c)
  dev.off()
  
 
  
  
  #################################################################################
  #auc
  l_all<-nn[ADM0_NAME %in% c(countries)]
  l_all[,country:=as.factor(ADM0_NAME)]
  l_all[,ADM0_NAME:=NULL]
  l_all[,rank_order:=NULL]
  
  
  l_all<-rbind(l_all,data.table(cov='Missed ITNs',country='Indonesia',auc=NA,order=3,
                                year=c(2000,cov_years[which(covariate_names_abbrev=='itns')])))
  
  l_all[,year:=as.character(year)]
  l_all[year!='2000',year:='Recent']
  
  l_all$country <- droplevels(l_all$country)
  
  
  c<-ggplot(data=l_all)+
    geom_bar_pattern(position='dodge',pattern_spacing=0.01,color='black',
                     stat='identity',
                     aes(fill=reorder(cov,order), y=auc,x=country,
                         pattern=year))+
    theme_bw()+
    geom_hline(yintercept=0.5,lty='dashed')+
    scale_fill_manual(values=cbPalette)+
    scale_pattern_manual(values=c('none','stripe')) +
    labs(x='',y='AUC',fill='Indicator',pattern='Year')+
    guides(fill = guide_legend(override.aes = list(pattern='none') ),
           pattern = guide_legend(override.aes = list(fill = "white")))
  
  # # #Combine and save both maps
  png(paste0('<FILEPATH>/zerodose_overlap_analyses/',run_date,'/results/paper/national_priority_plots/auc_all_countries_by_metrics_bars.png'),
      width = 10, height = 5, res=800, units = 'in')
  print(c)
  dev.off()
  
  
  
  
  
  
  
  toc()
}
message('done with all!')