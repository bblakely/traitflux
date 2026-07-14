#Read in ERA5
library(ncdf4)
library(terra)
library(sf)
library(raster)
library(lubridate)

homedir<-"D:/analysis/traitflux"
setwd(homedir)

# Bounding box
bbox<-vect("trait_bbox.geojson"); bbox2<-project(bbox, CRS("+init=epsg:32615") )

# Read in .nc
datadir<-"D:/analysis/traitflux/rsdata/ERA5"
setwd(datadir)
blh<-rast("data_stream-oper_stepType-instant.nc")

#assign time
timestr<-as.numeric(substr(names(blh),start=16, stop=26))
ts<-as.POSIXct(timestr, origin = "1970-01-01", tz = "UTC")
time(blh)<-ts

#reproject and crop
blh.utm<-project(blh, CRS("+init=epsg:32615"))
blh.crop<-crop(blh.utm, bbox2, snap="out")

blh.good<-blh.crop

#blh.res<-resample(blh.crop[[1:100]], evi.good, method="near") #consider doing average or nn to save time

