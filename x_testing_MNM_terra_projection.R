# Code to try to run MNMs at a future climate
# Traci P. DuBose 1/25/2023
# currently does not work

.libPaths('/home/tracidubose/R/OOD/Ubuntu-20.04-4.1.0')
print(.libPaths())
#library(devtools)
#devtools::install_github('mrke/NicheMapR')
#devtools::install_github('ilyamaclean/microclima')
library(NicheMapR)

# data download code by Mike Kerney in the Google Group
# data download described in AnuranMNMs/1_Building_InputTraits_Pts.Rmd
# https://groups.google.com/g/nichemapr/c/nOmFNhKTFrM/m/xTNL-GJ_AQAJ
PATH_ncep<-'ncep_data'
# data download described in google group message: https://groups.google.com/g/nichemapr/c/DC8dI5-0GXc
PATH_terra <- 'terra_climate'

start.time<-Sys.time()
m2<-micro_ncep(loc = c(-80.417778, 37.23),
               dstart = "01/01/2012", dfinish = "31/12/2012",
               DEP = c(0, 2.5,  5,  10,  15,  30,  50, 75, 100,  200), #cm
               minshade = 0, maxshade = 90, runshade=1,
               scenario=2,
               spatial = PATH_ncep,
               terra_source=PATH_terra,
               dem.res=1000, # requested resolution DEM from elevatr in meters
               Usrhyt = 0.01, # local height for organism
               ERR = 1.5,
               run.gads=2)
end.time<-Sys.time()
timdif<-end.time-start.time
cat(paste(timdif, 'min', '\n'))
saveRDS(m2, 'qq_m_w_scen.rds')
cat(paste0(data.frame(lat=m2$longlat[2],
                      long=m2$longlat[1],
                      tRainfall=sum(m2$RAINFALL),
                      evel=m2$elev,
                      slope=m2$SLOPE,
                      aspect=m2$ASPECT,
                      ndays=m2$ndays), '\n'))
cat(head(m2$metout))
