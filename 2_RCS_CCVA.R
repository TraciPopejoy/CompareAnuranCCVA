## ----setup, include=FALSE-----------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(sf); library(tidyverse); library(scales) # load necessary libraries.
library(terra); library(exactextractr)
library(cowplot) #for plotting maps later


## ----data prep, eval=T--------------------------------------------------------------------------
library(Hmisc)
# Set data paths -----
# data used
PATH_iucn_ranges <-'/home/tracidubose/RCS_Anuran_input/ANURA/' #IUCN range maps downloaded 10/2020
PATH_HUC12 <- '/home/tracidubose/huc12_wgs84_sf_02010323.rds'
PATH_FocalSpecies_OccurrenceData <- "/home/tracidubose/RCS_Anuran_input/Anuran Occurrence Data/"
PATH_out<-'b_RegionalRCS/spp_hucs/' # place to store occupied hucs with climate values

# Set spatial coordinate reference system (CRS).
crs.geo <- st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84")
crs.albers <- st_crs("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs")

# Read in fixed data -----
# Focal species & Spatial Extent determined in 0_SettingTaxSpaceExtents.Rmd
anuran.taxa<-c("Acris gryllus", "Acris crepitans",
               "Anaxyrus terrestris", "Anaxyrus fowleri", "Anaxyrus quercicus", 
               "Dryophytes avivoca",  "Dryophytes andersonii", "Dryophytes cinereus", 
               "Dryophytes chrysoscelis", "Dryophytes femoralis", "Dryophytes gratiosus",
               "Dryophytes squirellus", "Gastrophryne carolinensis",
               "Lithobates capito", "Lithobates grylio", "Lithobates okaloosae", 
               "Lithobates heckscheri", "Lithobates sphenocephalus", "Lithobates areolatus", 
               "Lithobates virgatipes", "Pseudacris ornata", "Pseudacris feriarum", 
               "Lithobates sevosus",
               "Pseudacris brachyphona", "Pseudacris fouquettei", "Scaphiopus holbrookii")
#IUCN maps for native ranges
iucn_anuran<-read_sf(PATH_iucn_ranges) %>%
  filter(binomial %in% anuran.taxa) %>% 
  st_transform(crs.albers) 
# prism file location
PATH_prism<-'/home/tracidubose/prism_files/'
PATH_climout<- '/home/tracidubose/b_RegionalRCS/climate_outliers/'

calcmode <- function(a){  
  vector <- unique(a)  
  vector[which.max(tabulate(match(a, vector)))]} 

# climate outliers QA/QC
climate_outliers<-read.csv('RCS_rerun1207/climate_outlier_pts.csv')

# spatial extent
huc12s<-readRDS(PATH_HUC12) %>%
  mutate(huc2=substr(huc12, 1,2)) %>% 
  filter(huc2 %in% c("03","06","08"))
spat.ext<-readRDS('b_RegionalRCS/spatial_extent.rds') %>% st_transform(crs.albers) %>%
  st_buffer(.1)
ggplot() + geom_sf(data=spat.ext) + theme_void()


## ----load prism, eval=T-------------------------------------------------------------------------
library(prism)
clim.variables<-c('ppt','Tmax','Tmin')
spat.big<-spat.ext %>% st_buffer(100)
for(cm in clim.variables){
  options(prism.path = paste0("//home/tracidubose/prism_files/", cm, "/"))
  PRISM_data <- pd_stack(prism_archive_ls()) 
  spat.big<- st_transform(spat.big, crs(PRISM_data))
  PRISM_se<-crop(PRISM_data, st_bbox(spat.big))
  ter.prism<-rast(PRISM_se)
  names(ter.prism)<- gsub(paste0('PRISM_', str_to_lower(cm),'_stable_4kmM2_'), '',
                          gsub(paste0('PRISM_', str_to_lower(cm),'_stable_4kmM3_'), '',
                             gsub('_bil', '', names(ter.prism))))
  assign(paste0(cm,'.prism'), ter.prism)
}
rm(PRISM_data, spat.big, PRISM_se, ter.prism)


## ----extract huc means, eval=T------------------------------------------------------------------
head(Tmax.prism)[1:5]
huc.Tmax<-exact_extract(Tmax.prism, huc12s, 'mean') %>%
  bind_cols(huc12=as.character(huc12s$huc12))
head(huc.Tmax)[c(1:3,126)]
head(Tmin.prism)[1:5]
huc.Tmin<-exact_extract(Tmin.prism, huc12s, 'mean')%>%
  bind_cols(huc12=as.character(huc12s$huc12))
head(huc.Tmin)[c(1:3,126)]
head(ppt.prism)[1:5]
huc.ppt<-exact_extract(ppt.prism, huc12s, 'mean')%>%
  bind_cols(huc12=as.character(huc12s$huc12))
head(huc.ppt)[c(1:3,126)]


## ----functions, eval=T--------------------------------------------------------------------------
myfun<-function(x,index){bind_cols( rowname=as.character(index), x)}
find_area<-function(taxa, # which species to find the aoo for
                    spatial.file1, # which sp object to use for filtering/grain size
                    native_range=T, #do you want to exclude points >100km outside the native range
                    watershed=F, #are you using watersheds as grain sizes?
                    rad=1, #if buffering points, what should the buffer radius be?
                    clim.vars #list of climate variables
                    ){
  # pulls in the spatial file using get and save original crs
  if(!is.na(spatial.file1)){
    spatial.file<-get(spatial.file1)
    crs.set<-st_crs(spatial.file)}else{
      crs.set<-crs.albers
    }

  #read in sp occurrence data
  geodata<-read.csv(list.files(PATH_FocalSpecies_OccurrenceData, pattern=taxa, full.names = T)) %>%
    dplyr::select(family, genus, species, Longitude, Latitude, year, source, final.taxa) %>%
    filter(!is.na(year), year > 1894) %>% 
    st_as_sf(coords=c("Longitude","Latitude"), crs = crs.geo) %>% #convert to a sf dataframe
    rowid_to_column() # add this to align climate outliers table with this
  orig.nrow<-nrow(geodata)
  
  # restrict points to those found within IUCN native range
   if(native_range==T){
    nat_range <- iucn_anuran %>% filter(binomial %in% taxa) %>%
      st_union() %>% st_buffer(100) 
    geodata <- geodata %>% st_transform(crs.albers) %>% st_filter(nat_range)
  } 
  postNR.nrow<-nrow(geodata)
  
  #identify any climate outliers to exclude
   # create 1 km buffer
  #spp_BUF <- st_buffer(geodata, dist = 1) %>%
   # st_transform("EPSG:4269")
    # extract the values of interest from the raster for each polygon
  #clim_out<-NULL
  #for(env.r in clim.vars){
   # env.raster<-rast(paste0(PATH_prism,env.r,'_mean_usa.tif'))
    #clim_out_v<-spp_BUF %>% as_tibble() %>%
     # dplyr::select(-geometry) %>%
      #bind_cols(clim_value=exact_extract(env.raster, spp_BUF, fun='mean', progress=F)) %>%
      #mutate(mean_v=mean(clim_value),
       #    sd_v=sd(clim_value),
        #   zscore=(mean_v-clim_value)/sd_v) %>%
    #filter(abs(zscore) > 4)
  #clim_out <- bind_rows(clim_out, clim_out_v)
  #}
  #if(nrow(clim_out)!=0){
   # print(paste(length(unique(clim_out$rowid)), 'of', nrow(geodata), 'points removed as climate outliers for', env.r, 'and', taxa))
  # keep record of all excluded points
  #write.csv(clim_out, paste0(PATH_climout, taxa, '_climout.csv'))}
  
  # exclude climate outliers identified in DuBose 2022 and PRISM data
  clim_out <- climate_outliers %>% filter(species == taxa)
  dat_sp<- geodata %>% 
    filter(!(rowid %in% clim_out$rowid))
  postCL.nrow <-nrow(dat_sp)
 
  # Calculate area occupied per species. 
  if(watershed==T){
    dat_sp <- st_transform(dat_sp, st_crs(spatial.file)) # ensuring fit between spatial file and points
    dat_sp <- st_join(dat_sp, spatial.file) %>% # spatial join to identify occupied hucs
      filter(!is.na(huc12)) %>%
      # only keep unique combinations of hucs & years
      distinct(huc12, year, .keep_all=T)
    all.poly.nrow<-nrow(dat_sp) #identify how many watersheds/year combinations exist
    all.poly.n<-dat_sp %>% distinct(huc12) %>% nrow() # how many watersheds occupied
    area.type<- paste0(spatial.file1,'_WS')
    used.nrow<-nrow(dat_sp)
    min.year<-min(dat_sp$year)
    area.i<-dat_sp %>% 
      # only keep one polygon per huc
      distinct(huc12, .keep_all=T) %>% 
      pull(WSAREA) %>% sum() %>% as.numeric()
    
    # Climate Sensitivity based on the annual AOO calculated above
    CS_out<-NULL
    for(env.r in paste0('huc.',clim.vars)){
    env.raster<-get(env.r)
    CS_out_v <- dat_sp %>% 
      # get rid of a bulky column
      as_tibble() %>% select(-geometry) %>%
      # change year to one where climate averages exist
      mutate(year = case_when(year < 1925 ~ 1925,
                                year > 2018 ~ 2018,
                                T ~ year)) %>%
      # keep all occupied hucs and join all climate mean data
      left_join(env.raster, by="huc12") %>% 
      # reorder columns
      dplyr::select(huc12, year, states, starts_with('mean.')) %>%
      # translate climate data from wide to long
      pivot_longer(starts_with('mean.')) %>%
      # pull out the year represented in each row
        mutate(climate_year=gsub('mean.', '', name)) %>%
      # only keep rows that are equal to or from 29 years before
      filter(year >= climate_year, year-30 < climate_year) %>%
      group_by(huc12, year) %>%
      # calculate the moving window value (based on annual means) of each huc year combination
      dplyr::summarize(clim_sum_val=case_when(env.r == 'huc.ppt'~ mean(value, na.rm=T),
                                              env.r == 'huc.Tmax'~ max(value, na.rm=T),
                                              env.r == 'huc.Tmin'~ min(value, na.rm=T)), .groups='drop') %>%
      ungroup() 
      # across all the data, calculate the mean and sd of the hucxyear mean
      CS_out<-bind_rows(CS_out, 
                        CS_out_v %>% 
                          dplyr::summarize(clim=gsub('huc.','',env.r),
                                           Cmean=mean(clim_sum_val),
                                           Csd=sd(clim_sum_val)),
                            # calculate some stats about the years
                        data.frame(modyear=calcmode(dat_sp$year), 
                                   newyears=sum(dat_sp$year > 2018),
                                   oldyear=sum(dat_sp$year < (2020-50))))
      if(spatial.file1=='huc12s'){
      write.csv(CS_out_v, paste0(PATH_out, taxa, '_climate_info_', gsub('huc.', '', env.r),'.csv'), row.names = F)}
    }
  }else{
    # Generate 1km point buffers. You must use a projected CRS. For CRS, we use: USA Contiguous albers equal area. 
     area.type<- paste0(spatial.file1,'_',rad,'km')
     dat_sp <- st_transform(dat_sp, st_crs(spatial.file)) %>% #overwrites to the projection 
      st_filter(spatial.file) %>% # only keep points pts in the spatial extent
       st_transform(crs.albers) %>% st_buffer(rad) # buffer by the requested radius
     #print(head(dat_sp))
  # calculate total area 
     area.i<-st_union(dat_sp) %>% st_area(.) %>% as.numeric()
     all.poly.nrow<-nrow(dat_sp)
     all.poly.n<-all.poly.nrow
     min.year<-min(dat_sp$year)
     # prep the data frame for climate variable extraction  
     if(nrow(dat_sp) >= 10){
     dat_sp <- dat_sp %>%
       count(year) %>% # this creates a single polygon for each year
       rownames_to_column('rowname') %>% # this adds a rowname to join it to the climate extracted data
       bind_cols(area.i=st_area(.)) %>% # calculates the AOO at each year
       st_transform(st_crs(ppt.prism)) # transform to match climate data
#     print(head(dat_sp))
     used.nrow<-nrow(dat_sp) # keep track of years used for climate
     
    CS_out<-NULL
    for(env.r in paste0(clim.vars, '.prism')){
    env.raster<-get(env.r) # pull in the climate data
    # extract the climate data at each year polygon, returns a list
    env_vals_raw<-exact_extract(env.raster, dat_sp, coverage_area=T, progress=F, force_df=T)
    # convert the list to a dataframe and join with dat_sp by the rowname
    env_vals<-do.call('rbind',
                      mapply(myfun, env_vals_raw, seq_along(env_vals_raw), SIMPLIFY = F)) %>% 
      left_join(dat_sp %>% as_tibble() %>% select(-geometry), by='rowname') %>%
      # reorder columns to make it easier to check
      dplyr::select(rowname, year, area.i, coverage_area, n, everything()) %>%
      # convert climate data from wide formate to long format
      pivot_longer(-c(rowname, year, area.i, coverage_area, n))
    # if loop to try to keep an eye out for points that have no climate data
    if(nrow(filter(env_vals, is.na(value))) != 0){
      print(paste('rownames', paste(filter(env_vals, is.na(value)) %>% pull(rowname) %>% unique(), collapse=', '), 
                  'suspicious points for', 'ppt', 'and', taxa))}
    oldyear <- dat_sp %>% filter(year < (2020-50)) %>% pull(area.i) %>% sum() %>% as.numeric()
    most_year <-dat_sp %>% filter(max(year)==year) %>% pull(year)
    
    CS_out_v <- env_vals %>%
       # change year to one where climate averages exist
      mutate(year = case_when(year < 1925 ~ 1925,
                                year > 2018 ~ 2018,
                                T ~ year)) %>%
      mutate(climate_year=gsub('mean.', '', name)) %>%
      filter(year >= climate_year, year-30 < climate_year) %>%
      # calculate the climate average at each area
      group_by(rowname, year, coverage_area) %>% 
      dplyr::summarize(clim_sum_val=case_when(env.r == 'ppt.prism'~ mean(value, na.rm=T),
                                              env.r == 'Tmax.prism'~ max(value, na.rm=T),
                                              env.r == 'Tmin.prism'~ min(value, na.rm=T)), .groups='drop') %>%
      # cleaning Infs created by trying to calculate min / max on areas outside the raster
      mutate(across(where(is.numeric), ~na_if(., Inf)), across(where(is.numeric), ~na_if(., -Inf))) %>%
      ungroup() %>%
      dplyr::summarize(modyear=most_year,
                   oldyear=oldyear,
                   clim=gsub('.prism','',env.r),
                   Cmean=wtd.mean(clim_sum_val, weights=coverage_area, normwt = F, na.rm=T),
                   Csd=sqrt(wtd.var(clim_sum_val, weights=coverage_area, normwt = F, na.rm=T)))
    CS_out<-bind_rows(CS_out, CS_out_v)
      }
     }else{
      print(paste(nrow(dat_sp), 'in the study region'))
       area.i=0
       used.nrow<-nrow(dat_sp); min.year = NA; area.sqkm = 0; all.poly.nrow<-all.poly.n<-0
       CS_out <- data.frame(clim=clim.vars, Cmean=rep(NA,length(clim.vars)), Csd=rep(NA, length(clim.vars)))
       }
  }  
  RCS_comp <- data.frame(scientific_name = taxa, 
                     area.type = area.type,
                     area.sqkm = area.i,
                     orig.nrow=orig.nrow,
                     postNR.nrow=postNR.nrow,
                     postCL.nrow=postCL.nrow,
                     n.polygons=all.poly.nrow,
                     total.WS=all.poly.n,
                     used.nrow=used.nrow,
                     min.year=min.year,
                     stringsAsFactors = F) %>%
    bind_cols(CS_out %>%
                pivot_wider(names_from = clim,
                            values_from=c("Cmean","Csd")))
  
  return(RCS_comp)
  print(paste('done with', taxa, 'at', spatial.file1))
}



## ----subreg map, echo=F-------------------------------------------------------------------------
library(maps)
se<-st_as_sf(map('state',c('georgia', 'florida', 'south carolina', 'alabama'), 
                 fill=T, plot=F)) %>%
  st_transform(crs.albers)
al<-se %>% filter(ID == 'alabama')
fl<-se %>% filter(ID == 'florida')
sc<-se %>% filter(ID == 'south carolina')
ga<-se %>% filter(ID == 'georgia')
huc3<-huc12s %>% filter(huc2=="03") %>% st_union()
huc6<-huc12s %>% filter(huc2=="06") %>% st_union()
huc8<-huc12s %>% filter(huc2=="08") %>% st_union()
ggplot()+geom_sf(data=huc3)+geom_sf(data=huc6)+geom_sf(data=huc8)+
  geom_sf(data=spat.ext)+
  geom_sf(data=se, fill=NA)+theme_void()


## ----calculate the RCS at many scales, warning=F, eval=T----------------------------------------
RCS_comp_all<-NULL
start.time<-Sys.time()
for(spp in anuran.taxa[anuran.taxa !='Lithobates sevosus']){
  buff_out<-NULL
  for(buf.ext in c("spat.ext", "huc3","huc6", "huc8","al","fl","sc","ga")){
  buff1<-find_area(taxa=spp,
                  spatial.file1=buf.ext, rad=1,
                  clim.vars = c('ppt','Tmax','Tmin'))
  buff5<-find_area(taxa=spp,
                  spatial.file1=buf.ext, rad=5,
                  clim.vars = c('ppt','Tmax','Tmin'))
  buff_out<-bind_rows(buff_out, buff1, buff5)
  print(paste('done with',buf.ext,'for',spp))
  }
  # only have to run hucs once
  huc12all<-find_area(spp, 
          "huc12s", native_range=T, watershed=T, 
          clim.vars=c('ppt','Tmax','Tmin'))

  RCS_comp_all<-bind_rows(RCS_comp_all, huc12all, buff_out)
}
end.time<-Sys.time()
end.time-start.time


## -----------------------------------------------------------------------------------------------
# read in IUCN range map
ls_aoo <- read_sf('RCS_Anuran_input/ANURA') %>%
  filter(binomial=='Lithobates sevosus',
         grepl('Extant', legend)) %>%
  rename(year=yrcompiled) 
# HUC12 watershed grain size ----
dat_sp<-st_join(ls_aoo, huc12s) %>% # spatial join to identify occupied hucs
      filter(!is.na(huc12)) %>%
      # only keep unique combinations of hucs & years
      distinct(huc12, year, .keep_all=T)
# Climate Sensitivity based on the annual AOO calculated above
    CS_out<-NULL
    for(env.r in paste0('huc.',clim.vars)){
    env.raster<-get(env.r)
    CS_out_v <- dat_sp %>% 
      # get rid of a bulky column
      as_tibble() %>% select(-geometry) %>%
      # keep all occupied hucs and join all climate mean data
      left_join(env.raster, by="huc12") %>% 
      # reorder columns
      dplyr::select(huc12, year, states, starts_with('mean.')) %>%
      # translate climate data from wide to long
      pivot_longer(starts_with('mean.')) %>%
      # pull out the year represented in each row
        mutate(climate_year=gsub('mean.', '', name)) %>%
      # only keep rows that are equal to or from 29 years before
      filter(year >= climate_year, year-30 < climate_year) %>%
      group_by(huc12, year) %>%
      # calculate the moving window value (based on annual means) of each huc year combination
      dplyr::summarize(clim_sum_val=case_when(env.r == 'huc.ppt'~ mean(value, na.rm=T),
                                              env.r == 'huc.Tmax'~ max(value, na.rm=T),
                                              env.r == 'huc.Tmin'~ min(value, na.rm=T)), .groups='drop') %>%
      ungroup() 
      # across all the data, calculate the mean and sd of the hucxyear mean
      CS_out<-bind_rows(CS_out, 
                        CS_out_v %>% 
                          dplyr::summarize(clim=gsub('huc.','',env.r),
                                           Cmean=mean(clim_sum_val),
                                           Csd=sd(clim_sum_val)),
                            # calculate some stats about the years
                        data.frame(modyear=calcmode(dat_sp$year), 
                                   newyears=sum(dat_sp$year > 2018),
                                   oldyear=sum(dat_sp$year < (2020-50))))
    
write.csv(CS_out_v, paste0(PATH_out, 'Lithobates sevosus_climate_info_', gsub('huc.', '', env.r),'.csv'), row.names = F)}
# summarize huc12 analysis
area.i<-dat_sp %>% 
      # only keep one polygon per huc
      distinct(huc12, .keep_all=T) %>% 
      pull(WSAREA) %>% sum() %>% as.numeric()
ls.huc12.res <- data.frame(scientific_name = "Lithoates sevosus", 
                     area.type = 'huc12s_WS',
                     area.sqkm = area.i,
                     n.polygons=all.poly.nrow,
                     total.WS=dat_sp %>% distinct(huc12) %>% nrow(), # how many watersheds occupied,
                     min.year=2008,
                     stringsAsFactors = F) %>%
    bind_cols(CS_out %>%
                pivot_wider(names_from = clim,
                            values_from=c("Cmean","Csd")))

# Now do it for 1km and 5km buffers ----
# IUCN extant distribution is only within MS but extinct in AL
# only running it at entire extent because no 'points' in AL
ls.buff.res<-NULL
for(rad in c(1,5)){
 dat_sp <- st_transform(ls_aoo, st_crs(spat.ext)) %>% #overwrites to the projection 
      st_filter(spat.ext) %>% # only keep points pts in the spatial extent
       st_transform(crs.albers) %>% st_buffer(rad) # buffer by the requested radius
     #print(head(dat_sp))
  # calculate total area 
  area.i<-st_union(dat_sp) %>% st_area(.) %>% as.numeric()
  # prep the data frame for climate variable extraction  
 dat_sp <- dat_sp %>% st_transform(st_crs(ppt.prism)) # transform to match climate data
    CS_out<-NULL
    for(env.r in paste0(clim.vars, '.prism')){
    env.raster<-get(env.r) # pull in the climate data
    # extract the climate data at each year polygon, returns a list
    env_vals_raw<-exact_extract(env.raster, dat_sp, fun='mean', progress=F, force_df=T)
    CS_out_v <- env_vals %>%
      pivot_longer(-species) %>%
      mutate(climate_year=gsub('mean.', '', name)) %>%
      filter(year >= climate_year, year-30 < climate_year) %>%
      dplyr::summarize(clim=gsub('.prism','',env.r),
                       Cmean=case_when(env.r == 'ppt.prism'~ mean(value, na.rm=T),
                                              env.r == 'Tmax.prism'~ max(value, na.rm=T),
                                              env.r == 'Tmin.prism'~ min(value, na.rm=T)),
                       Csd=sd(value, na.rm=T))
    CS_out<-bind_rows(CS_out, CS_out_v)
      }
ls.buff.res<-bind_rows(ls.buff.res,
                       data.frame(scientific_name = "Lithobates sevosus", 
                     area.type = paste0('spat.ext_', rad,'km'),
                     area.sqkm = area.i,
                     min.year=2008,
                     stringsAsFactors = F) %>%
    bind_cols(CS_out %>%
                pivot_wider(names_from = clim,
                            values_from=c("Cmean","Csd")))
}

# combine it with the rest of the RCS info and write out to save it
RCS_comp_all<-bind_rows(RCS_comp_all, ls.huc12.res, ls.buff.res)
write.csv(RCS_comp_all, 'b_RegionalRCS/buffer RCS components.csv', row.names = F)


## ----code to show where climate NAs can occur, echo=F, warning=F, eval=T------------------------
ag_pts<-read.csv(list.files(PATH_FocalSpecies_OccurrenceData, pattern="Acris gryllus", full.names = T)) %>%
    dplyr::select(family, genus, species, Longitude, Latitude, year, source, final.taxa) %>%
    filter(!is.na(year), year > 1894) %>% 
    st_as_sf(coords=c("Longitude","Latitude"), crs = crs.geo) %>% #convert to a sf dataframe
    rowid_to_column()
ag_sp <- st_transform(ag_pts, st_crs(spat.ext)) %>% #overwrites to the projection 
      st_filter(spat.ext) %>% # only keep points pts in the spatial extent
       st_transform(crs.albers) %>% st_buffer(5) %>%
   count(year) %>% # this creates a single polygon for each year
       rownames_to_column('rowname') %>% # this adds a rowname to join it to the climate extracted data
       bind_cols(area=st_area(.)) %>% # calculates the AOO at each year
       st_transform(st_crs(ppt.prism))
find_area(taxa='Acris gryllus',
                  spatial.file1='spat.ext', rad=5,
                  clim.vars = c('ppt','Tmax','Tmin'))
library(tidyterra)
ggplot()+#geom_sf(data=spat.ext) +
  geom_spatraster(data=ppt.prism, aes(fill=`1954`))+
  geom_sf(data = ag_sp %>% filter(rowname==28), color='white', fill=NA)+
  coord_sf(xlim=st_bbox(dat_sp[28,])[c(1,3)],
           ylim=st_bbox(dat_sp[28,])[c(2,4)])


## ----evaluate huc3, eval=T----------------------------------------------------------------------
all.spp.huc<-NULL
# read in all the individual species x climate variable files. 
for(i in list.files(PATH_out)){
  all.spp.huc<-bind_rows(all.spp.huc,
                         read.csv(paste0('b_RegionalRCS/spp_hucs/', i), 
                                  colClasses = c('character','integer','numeric')) %>%
                           # make sure to keep species name and climate variable
                           mutate(scientific_name=gsub('_.*','',i),
                                  climate_var=gsub('_','',substr(i, nchar(i)-7, nchar(i)-4))))
}

# join with huc12 tibble 
main.huc.df<-left_join(all.spp.huc %>% 
                         # pivot from long to wide so that each row is unique species x huc
                         pivot_wider(names_from=climate_var, values_from=clim_sum_val),
          huc12s %>% as_tibble() %>% select(huc2, huc12, WSAREA, states),
          by='huc12')

huc.states<-NULL
for(stt in c("AL","FL","SC","GA")){
  st.huc<-main.huc.df %>% 
    filter(grepl(stt, states)) # only keep hucs that have state abbreviation in states column
  # calculate the area
  area.i<-st.huc %>% 
    # keep unique rows with species x huc
    distinct(scientific_name, huc12, .keep_all=T) %>%
    group_by(scientific_name) %>%
    dplyr::summarize(area.type=paste0(tolower(stt), '_WS'),
                     area.sqkm=sum(WSAREA))
  # calculate the climate niche breadths
  huc.sum<-st.huc %>% 
    # keep unique rows with species x huc x year
    distinct(scientific_name, huc12, year, .keep_all=T) %>%
    group_by(scientific_name) %>% # for each species
    summarize(area.type=paste0(tolower(stt), '_WS'),
              n.polygons=n(), # count the n huc x year
              total.WS=n_distinct(huc12), # count the unique hucs
              # calculate the sd and mean across all hucs
              across(c(ppt,Tmax,Tmin), sd, .names='Csd_{.col}'),
              across(c(ppt,Tmax,Tmin), mean, .names='Cmean_{.col}'),
              min.year=min(year),
              modyear=calcmode(year),
              oldyear=sum(year < (2020-50)))
  # join it with the huc states data frame
  huc.states<-bind_rows(huc.states,
            left_join(area.i,huc.sum, by=c('scientific_name','area.type')))
}

huc.huc<-main.huc.df %>%
  group_by(scientific_name, huc2) %>%
  dplyr::summarize(area.type=paste0('huc', unique(huc2), '_WS'),
                     area.sqkm=sum(WSAREA),
                   n.polygons=n(),
              total.WS=n_distinct(huc12),
              across(c(ppt,Tmax,Tmin), sd, .names='Csd_{.col}'),
              across(c(ppt,Tmax,Tmin), mean, .names='Cmean_{.col}'),
              min.year=min(year),
              modyear=calcmode(year),
              oldyear=sum(year < (2020-50)))

# bind these results with the rest of the RCS component results
RCS_comp_all<-bind_rows(RCS_comp_all, huc.states, huc.huc) %>%
  # only keep one record for species, area type, and area of occupancy
  distinct(scientific_name, area.type, area.sqkm, .keep_all=T) 


## ----extent pic, fig.width=4, fig.height=2.5, eval=T--------------------------------------------
ggplot()+
  geom_sf(data=huc12s %>% filter(grepl('AL',states)), aes(fill='watershed'))+
  geom_sf(data=se %>% filter(ID == 'alabama'), fill=NA, aes(color='buffer'), linewidth=2)+
  scale_color_manual('Grain Size', values=c('goldenrod','forestgreen'), 
                     aesthetics = c('color','fill'))+
  theme_void()+theme(legend.position='left')


## ----write out, eval=T--------------------------------------------------------------------------
RCS_comp_all %>% arrange(scientific_name)
write.csv(RCS_comp_all, 'b_RegionalRCS/regional_rcs_011923.csv', row.names = F)

