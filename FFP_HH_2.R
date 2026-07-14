setwd("D:/analysis/traitflux")

## set home directory, get packages#####
library(EBImage)
library(jpeg)
library(fields)
library(zoom)
library(readxl)
library(raster)
library(terra)
library(stringr)
library(ggplot2)
library(esquisse)
library(sf)
library(lubridate)
library(ncdf4)

rt<-"D:/analysis/traitflux"
setwd(rt)


source("FFP_R/calc_footprint_FFP_climatology.R") #footprint function

#####

##read in trait rasters and hyspex LAI + albedo####

setwd(rt)

if(!exists("pheno.trait.ls")){
  
  setwd("rsdata/traitmosaic/0629")
  
  nsc.trait<-terra::rast("NSC.tif")
  lma.trait<-terra::rast("LMA.tif")
  lig.trait<-terra::rast("Lignin.tif")
  n.trait<-terra::rast("Nitrogen.tif")
  pheno.trait<-terra::rast("Phenolics.tif")
  
  setwd(rt)

  
  setwd("rsdata/traitmosaic/0830")
  
  nsc.trait.ls<-terra::rast("NSC.tif")
  lma.trait.ls<-terra::rast("LMA.tif")
  lig.trait.ls<-terra::rast("Lignin.tif")
  n.trait.ls<-terra::rast("Nitrogen.tif")
  pheno.trait.ls<-terra::rast("Phenolics.tif")
  
  setwd(rt)


}

setwd(rt)


if(!exists("lai_hires")){
  
  setwd("rsdata/laimosaic")
  
  lai.jun<-terra::rast("lai_0629.tif"); time(lai.jun)<-180
  lai.aug<-terra::rast("lai_0830.tif"); time(lai.aug)<-242
  lai.aug<-resample(lai.aug, lai.jun)
  
  lai.highres<-c(lai.jun, lai.aug)
  
}

setwd(rt)


if(!exists("alb_hires")){
  
  setwd("rsdata/albmosaic")
  
  
  albhr.jun<-terra::rast("hyspex_albedo_0629.tif"); time(albhr.jun)<-180
  albhr.aug<-terra::rast("hyspex_albedo_0830.tif"); time(albhr.aug)<-242
  albhr.aug<-resample(albhr.aug, albhr.jun)
  
  alb.highres<-c(albhr.jun, albhr.aug)
  
}

setwd(rt)

#####


##Get tower names and metadata + veg indices, LAI, BLH and WiscLand#####

setwd("processed")
sitelist<-substr(list.files(),1,6) #get list of processed sites
setwd(rt)
geom<-read.csv("tower_coords_BB.csv") #from Emily; WGS coords and instrument info for towers
mgmt<-read.csv("Murphy_etal_table2.csv") #from Bailey's paper; stand age for forests and qualitative management info
colnames(mgmt)[1]<-"Site"

#both together in one df, skipping Emily's file name columns
ancil<-merge(mgmt, geom, by.x="Site", by.y="dir_names", all=TRUE)
ancil<-ancil[,c(1:6,8:9, 12:18)]

#Make a points layer for plotting
pts.dat<-ancil[,c(6,8,7)]; colnames(pts.dat)<-c("site", "lon", "lat")
allpts.wgs<-vect(x=pts.dat);crs(allpts.wgs) ="+init=epsg:4326";
allpts<-project(allpts.wgs, CRS("+init=epsg:32615"))


#read in wiscland
if (!exists("wl.fact")){source("D:/analysis/traitflux/wiscland/wiscland2.R")}
#report most common land types


setwd(rt)

#Read in veg indices
if(!exists("evi.good")){source("D:/analysis/traitflux/Read_VI.R")}
#vi.index<-2 #which of the reflectance images to use. 
evirast<-evi.good#[[vi.index]]
ndvirast<-ndvi.good#[[vi.index]]
nirvrast<-nirv.good#[[vi.index]] 

if(!exists("lai.good")){source("D:/analysis/traitflux/Read_LAI.R")}
lairast<-lai.good

if(!exists("blh.good")){source("D:/analysis/traitflux/BLH.R")}
blhrast<-blh.good

if(!exists("alb")){source("D:/analysis/CH19Alb/Read_Albedo.R")
  albrast<-alb
  }

#####


##Set up main loop ####

#make empty items to hold stuff
fprast<-list()
pftholder<-rep(NA, length(sitelist))

#currently hard-coded; consider making flexible at some point

vars<-c("site", "igbp","agemean", "agemax", "agemin", "mgmt")
sitedat<-data.frame(matrix(nrow=length(sitelist), ncol=length(vars))); colnames(sitedat)<-vars


#time range you want to conduct this analysis for
doylist<-c(152:304)#c(180:195)

#####

## Main Site loop: get metadata, calculate footprints and extract FP-weighted values #### 
setwd(rt)
traitdat<-list()
for (f in 1:length(sitelist)){
  
  towername<-sitelist[f]
  
  readname<-paste("processed/", towername, "_processed.csv", sep="")
  data<-read.csv(readname)
  
  readname.meta<-paste("BIF/AMF_",towername,"_BIF_20210730.xlsx", sep="")
  meta<-read_excel(readname.meta)
  
  ### Part 1: get metadata #####
  
  #grab the BIF PFT, coords
  site.pft<-meta$DATAVALUE[meta$VARIABLE=="IGBP"]
  site.lat<-as.numeric(meta$DATAVALUE[meta$VARIABLE=="LOCATION_LAT"])
  site.lon<-as.numeric(meta$DATAVALUE[meta$VARIABLE=="LOCATION_LONG"])
  
  pftholder[f]<-site.pft
  
  #grab Mather-extracted metadata
  site.ht<-geom$z_son[tolower(geom$Tower)==tolower(substr(towername, 4,6))]
  site.can<-geom$veg_h[tolower(geom$Tower)==tolower(substr(towername, 4,6))]
  
  #grab murphy-extracted harvest data
  site.agemean<-ancil$av_age[tolower(ancil$Tower)==tolower(substr(towername, 4,6))]
  site.agemin<-ancil$min_age[tolower(ancil$Tower)==tolower(substr(towername, 4,6))]
  site.agemax<-ancil$max_age[tolower(ancil$Tower)==tolower(substr(towername, 4,6))]
  #should ultimately make a bunch of dummy variables for e.g. has it been harvested
  
  #skip FFP if no instrument/canopy height
  if(length(site.ht)==0|length(site.can)==0){
    print(paste("missing geom for", towername))
    print("##############")
    next
  }
  
  
  ###Part 2: calculate footprints####
  
  #clip timespan to where we actually have data
  bounds<-c(min(which(!is.na(data$FC))),max(which(!is.na(data$FC))))
  dat.clip<-data[bounds[1]:bounds[2],]
  #make a numerical timestamp
  ts<-as.POSIXlt(dat.clip$ts_posix)
  dat.clip$doy<-as.numeric(format(ts, "%j"))
  #flag day and night
  dat.clip$daytime<-1; dat.clip$daytime[dat.clip$SW_IN_f<10]<-0
  
  #get full list of available DOYs
  doys<-unique(dat.clip$doy)
  nlev<-rep(NA, length(doys))
  

  print("########")
  print(paste("starting site", towername))
  print("########")
  
  #Clip to DOYs specified above
  dat<-dat.clip[dat.clip$doy%in%doylist,]
  
  #make a decimal doy timestamp
  dat.ts<-as.POSIXct(substr(dat$TIMESTAMP_END,5, 12), format="%m%d%H%M")
  year(dat.ts)<-2019
  h<-hour(dat.ts); m<-minute(dat.ts); doy<-dat$doy
  dec<-doy+(h/24)+(m/1440)
  dat$decdoy<-dec
  
  #add columns for all the information to be extracted
  newtraitnames<-c("nsc", "lignin","lma", "nitr", "pheno",
                   "nsc_ls", "lignin_ls", "lma_ls", "nitr_ls", "pheno_ls",
                   "evi", "ndvi", "nirv", "nirvp", "lai", "laihr", "alb", "albhr")
  traitcols<-data.frame(matrix(nrow=nrow(dat), ncol=length(newtraitnames), data=NA)); colnames(traitcols)<-newtraitnames
  
  newcompnames<-c("Fir Spruce","Coniferous Forested Wetland", "Pine", "Aspen/Paper Birch" , "Northern Hardwoods","Broad-leaved Deciduous Scrub/Shrub", "Swamp Hardwoods", "Broad-leaved Evergreen Scrub/Shrub", "Open Water", "Red Maple","Mixed Deciduous/Coniferous Forest", "Other")
  compcols<-data.frame(matrix(nrow=nrow(dat), ncol=length(newcompnames), data=0)); colnames(compcols)<-newcompnames
  
  
  dat<-cbind(dat, traitcols, compcols)
  
  if(nrow(dat)==0){
    print(paste("site", towername,"has no data for these DOYs"))
    print("##########")
    next
  }
  
  #Provide FFP parameters####
  
  #Prepare blh:
  #Get tower location (this is repeated below; clean up at some point)
  coord1<-cbind(site.lat, site.lon)
  coord<-data.frame(coord1); colnames(coord)<-c("lat", "lon") #rbind(coord1, coord2))
  loc<-vect(coord, geom=c("lon","lat"),crs="+proj=longlat +datum=WGS84")#SpatialPoints(coords=coord,proj4string=CRS("+proj=longlat"))
  utmloc<-project(loc, CRS("+init=epsg:32615"))
  #extract hourly blh and time
  bl.raw<-unname(extract(blh.good, utmloc, ID=FALSE, raw=TRUE)[1,])
  bl.time<-time(blh.good);bl.time.loc<-with_tz(bl.time,"Etc/GMT+6")
  #Convert bl.time to decimal doy
  h<-hour(bl.time.loc); m<-minute(bl.time.loc); doy<-yday(bl.time.loc)
  dec<-doy+(h/24)+(m/1440); bl.decdoy<-dec
  #Make df
  bl.df<-data.frame(cbind(bl.raw, bl.time, bl.decdoy))
  bl.df$posix<-as.POSIXct(bl.df$bl.time, origin = "1970-01-01", tz = "Etc/GMT+6")
  #clip to same time period
  ind<-which(bl.decdoy%in%dat$decdoy);bl.df<-bl.df[ind,]
  #interpolate to half-hourly
  bl.hh<-approx(y=bl.df$bl.raw,x=bl.df$bl.decdoy, xout=dat$decdoy)
  
  
  #non-varying params
  z<-site.ht #instrument height
  canht<-site.can #canopy height
  d<-0.67*canht #zero plane displacement, estimated
  
  #varying params
  zm<-rep((z-d), nrow(dat))#klujn wants zm to be above displacement height
  z0<-rep((0.15*canht),nrow(dat))#rep((NaN),nrow(dat))#rep((0.15*canht),nrow(dat)) #klujn uses ws OR Z0; Z0 is preferred if both given
  h<-bl.hh$y#h<-rep(1500, nrow(dat)); h[dat$daytime==0]<-300 #should use real BL ht of course, just making it generaically unstable in day, generaically stable at night
  umean<-dat$WS_1_1_1; if (length(umean)<1){umean<-dat$WS}
  ol<-dat$MO_LENGTH_1_1_1
  if(length(ol)!=0){ol[abs(ol)>5000]<-NA}else{ol<-zm/dat$ZL;ol[abs(ol)>5000]<-NA} #clip extreme values, calc from Zl if missing
  sigmav<-dat$V_SIGMA_1_1_1;if(length(sigmav)<1){sigmav<-dat$V_SIGMA};if(all(is.na(sigmav))){source("Calc_sigv.R"); sigmav<-calc.sigv(dat,zi=h)}
  ustar<-dat$USTAR
  wind_dir<-dat$WD_1_1_1
  
  #if((z0[1]*12.5)>=zm[1]){z0<-(zm/12.5)-0.1; print(paste("Adjusting roughness length for",sitename[f]))}
  if((z0[1]*12.5)>=zm[1]){z0<-NA} #if roughness layer doesn't make sense, take it out and let the algo work on WS.
  

  #####

  #Attempt FFP calculation#####
  
  dat.ts<-as.POSIXct(substr(dat$TIMESTAMP_END,5, 12), format="%m%d%H%M")
  year(dat.ts)<-2019

  for(i in (1:nrow(dat))){ #for each half-hourly period

  print(paste("starting", dat.ts[i]))

  fp<-ffp<-fpgr<-NA #set to NA by default so that it's of length 1 when ffp fails

  tryCatch(
    fp<-calc_footprint_FFP_climatology(zm=zm[i],z0=z0[i],umean=umean[i],h=h[i],ol=ol[i],sigmav=sigmav[i],ustar=ustar[i],wind_dir=wind_dir[i], fig=0, r=80, crop=1, pulse=10),
    error = function(e) return(NA)
  )
 
  ######
  
##If FFP worked, proceed to extraction.#####

  if(length(fp)!=1&length(fp)!=4){ #1 happens when fp fails

    #check alignment; these should look similar
    #par(mfrow=c(1,2))
    #image.plot(fp$x_2d[1,], fp$y_2d[,1], fp$fclim_2d) #default from ffp package

    #timg2<-rasterFromXYZ(cbind(rep(2*fp$x_2d[1,], length(fp$y_2d[,1])), rep(2*fp$y_2d[,1], each=length(fp$x_2d[1,])),as.vector(fp$fclim_2d)))
    #spatimg<-rast(timg2)
    #plot(spatimg) #my raster

    #Georeference footprint. Center on tower, add/subtract meters
    coord1<-cbind(site.lat, site.lon)
    coord<-data.frame(coord1); colnames(coord)<-c("lat", "lon") #rbind(coord1, coord2))
    loc<-vect(coord, geom=c("lon","lat"),crs="+proj=longlat +datum=WGS84")#SpatialPoints(coords=coord,proj4string=CRS("+proj=longlat"))
    utmloc<-project(loc, CRS("+init=epsg:32615"))

    #Add x and y coordinates (which are in meters) to tower center points in UTM
    x.new<-2*fp$x_2d[1,]+crds(utmloc)[1]
    y.new<-2*fp$y_2d[,1]+crds(utmloc)[2]

    #Make a georeferenced footprint weight raster from new x, new y, and footprint weights.
    timg3<-rasterFromXYZ(cbind(rep(x.new, length(y.new)), rep(y.new, each=length(x.new)),as.vector(fp$fclim_2d)))
    fpgr<-rast(timg3); crs(fpgr)<-crs(utmloc)
    fprast[[f]]<-fpgr
    
    #plot(evirast); plot(fpgr, add=TRUE, col="white")

    #print sum of FFP; need to address why this is so low
    #print(paste("fp sum:", sum(as.vector(fpgr))))

    
    ##Part 3: Extract information#####
    
    print("Extracting information...")
    
    #extract traits by footprint:
    lma<-lig<-carb<-nitr<-pheno<-nsc<-alum<-boron<-cellulose<-calcium<-fiber<-copper<-NA; #set to NA so that old values aren't copied erroneously
    
    traitextract<-function(rast, ffp=fpgr){

      #crop and resample trait to match footprint weight raster
      trait.crop<-crop(rast, ffp)
      trait.resamp<-resample(trait.crop, ffp)
      #plot(trait.resamp)
      #points(utmloc, col='white', cex=3)

      #do the multiplication
      trait.val<-ffp*trait.resamp

      #might be worth checking to see if this is handling NA properly
      trait<-sum(as.vector(trait.val), na.rm=TRUE)*(1/sum(as.vector(ffp[!is.na(trait.val)])))

      #print(paste("trait value:", trait))

      return(trait)


    }
    
    #Extract traits
    
    #June 29
    nsc<-traitextract(rast=nsc.trait, ffp=fpgr); dat$nsc[i]<-nsc
    lma<-traitextract(rast=lma.trait, ffp=fpgr); dat$lma[i]<-lma
    lig<-traitextract(rast=lig.trait, ffp=fpgr); dat$lignin[i]<-lig
    nitr<-traitextract(rast=n.trait, ffp=fpgr); dat$nitr[i]<-nitr
    pheno<-traitextract(rast=pheno.trait, ffp=fpgr); dat$pheno[i]<-pheno

    
    #August 30
    nsc.ls<-traitextract(rast=nsc.trait.ls, ffp=fpgr); dat$nsc_ls[i]<-nsc.ls
    lma.ls<-traitextract(rast=lma.trait.ls, ffp=fpgr); dat$lma_ls[i]<-lma.ls
    lig.ls<-traitextract(rast=lig.trait.ls, ffp=fpgr); dat$lignin_ls[i]<-lig.ls
    nitr.ls<-traitextract(rast=n.trait.ls, ffp=fpgr); dat$nitr_ls[i]<-nitr.ls
    pheno.ls<-traitextract(rast=pheno.trait.ls, ffp=fpgr); dat$pheno_ls[i]<-pheno.ls


    
    #extract VI's.
    #find closest DOY from vi stack. Depends on all vi having same timestamps (which they do here, because they are all from landsat)
    closest.day.ind<-which(abs(dat$doy[i]-time(evirast))==min(abs(dat$doy[i]-time(evirast))))
    #if equidistant, just pick one
    if(length(closest.day.ind)>1){closest.day.ind<-closest.day.ind[sample(1:length(closest.day.ind), 1)]}
    
    #extract evi by footprint:
    evi<-NA
    #grab the image from nearest day
    evi.fp<-evirast[[closest.day.ind]]
    evi.crop<-crop(evi.fp, fpgr,snap="out")
    evi.resamp<-resample(evi.crop, fpgr)
    evi.val<-fpgr*evi.resamp
    evi<-sum(as.vector(evi.val), na.rm=TRUE)*(1/sum(as.vector(fpgr[!is.na(evi.val)])))
    
    dat$evi[i]<-evi
    
    #extract ndvi by footprint:
    ndvi<-NA
    #grab the image from nearest day
    ndvi.fp<-ndvirast[[closest.day.ind]]
    ndvi.crop<-crop(ndvi.fp, fpgr,snap="out")
    ndvi.resamp<-resample(ndvi.crop, fpgr)
    ndvi.val<-fpgr*ndvi.resamp
    ndvi<-sum(as.vector(ndvi.val), na.rm=TRUE)*(1/sum(as.vector(fpgr[!is.na(ndvi.val)])))
    
    dat$ndvi[i]<-ndvi
    
    #extract nirv by footprint:
    nirv<-NA
    #grab the image from nearest day
    nirv.fp<-nirvrast[[closest.day.ind]]
    nirv.crop<-crop(nirv.fp, fpgr,snap="out")
    nirv.resamp<-resample(nirv.crop, fpgr)
    nirv.val<-fpgr*nirv.resamp
    nirv<-sum(as.vector(nirv.val), na.rm=TRUE)*(1/sum(as.vector(fpgr[!is.na(nirv.val)])))
    
    dat$nirv[i]<-nirv
    
    #calculate nirvp
    nirvp<-NA
    par<-dat$SW_IN[i]*0.48 #yields PAR in Wm-2
    dat$nirvp[i]<-nirv*par
    
    
    #Extract MOSDIS albedo (for radiometer experiment)
    
    alb<-NA
    #find closest DOY from alb stack
    closest.day.ind<-which(abs(dat$doy[i]-time(albrast))==min(abs(dat$doy[i]-time(albrast))))
    #if equidistant, just pick one
    if(length(closest.day.ind)>1){closest.day.ind<-closest.day.ind[sample(1:length(closest.day.ind), 1)]}
    #grab the image from that day
    alb.fp<-albrast[[closest.day.ind]]
    #rest of the steps just like EVI
    alb.crop<-crop(alb.fp, fpgr, snap="out")
    alb.resamp<-resample(alb.crop, fpgr)
    alb.val<-fpgr*alb.resamp
    alb<-sum(as.vector(alb.val), na.rm=TRUE)*(1/sum(as.vector(fpgr[!is.na(alb.val)])))
    
    dat$alb[i]<-alb
    
    
    #extract hyspex albedo
    alb.hr<-NA
    #find closest DOY from alb stack
    closest.day.ind<-which(abs(dat$doy[i]-time(alb.highres))==min(abs(dat$doy[i]-time(alb.highres))))
    #if equidistant, just pick one
    if(length(closest.day.ind)>1){closest.day.ind<-closest.day.ind[sample(1:length(closest.day.ind), 1)]}
    #grab the image from nearest day
    albhr.fp<-alb.highres[[closest.day.ind]]
    albhr.crop<-crop(albhr.fp, fpgr,snap="out")
    albhr.resamp<-resample(albhr.crop, fpgr)
    albhr.val<-fpgr*albhr.resamp
    albhr<-sum(as.vector(albhr.val), na.rm=TRUE)*(1/sum(as.vector(fpgr[!is.na(albhr.val)])))
    
    dat$albhr[i]<-albhr
    
    
    
    
    #extract MODIS lai by footprint
    lai<-NA
    #find closest DOY from lai stack
    closest.day.ind<-which(abs(dat$doy[i]-time(lairast))==min(abs(dat$doy[i]-time(lairast))))
    #if equidistant, just pick one
    if(length(closest.day.ind)>1){closest.day.ind<-closest.day.ind[sample(1:length(closest.day.ind), 1)]}
    #grab the image from that day
    lai.fp<-lairast[[closest.day.ind]]
    #rest of the steps just like EVI
    lai.crop<-crop(lai.fp, fpgr, snap="out")
    lai.resamp<-resample(lai.crop, fpgr)
    lai.val<-fpgr*lai.resamp
    lai<-sum(as.vector(lai.val), na.rm=TRUE)*(1/sum(as.vector(fpgr[!is.na(lai.val)])))
    
    dat$lai[i]<-lai
    
    
    #extract hyspex lai by footprint
    lai.hr<-NA
    #find closest DOY from lai stack
    closest.day.ind<-which(abs(dat$doy[i]-time(lai.highres))==min(abs(dat$doy[i]-time(lai.highres))))
    #if equidistant, just pick one
    if(length(closest.day.ind)>1){closest.day.ind<-closest.day.ind[sample(1:length(closest.day.ind), 1)]}
    #grab the image from nearest day
    laihr.fp<-lai.highres[[closest.day.ind]]
    laihr.crop<-crop(laihr.fp, fpgr,snap="out")
    laihr.resamp<-resample(laihr.crop, fpgr)
    laihr.val<-fpgr*laihr.resamp
    laihr<-sum(as.vector(laihr.val), na.rm=TRUE)*(1/sum(as.vector(fpgr[!is.na(laihr.val)])))
    
    dat$laihr[i]<-laihr
    
    
    #extract composition by footprint
    wl.crop<-crop(wl.fact, fpgr)
    wl.resamp<-resample(wl.crop, fpgr)
    
    #calculate footprint-weighted compositions
    compweights<-zonal(fpgr, wl.resamp, fun="sum"); compweights$pct<-round(compweights$layer/sum(values(fpgr)),3) #gives plausible values,must convince myself I am correct
    #cover types not included in 95% list, put in "other" bin
    other<-sum(compweights$pct[which(!compweights$cover%in%newcompnames)])
    dat$Other[i]<-other
    #checking the results manually####
    # cls$id[cls$cover=="Fir Spruce"]
    # fsmask<-wl.resamp; values(fsmask)<-1; values(fsmask)[values(wl.resamp!=4110)]<-0
    # fsweight<-(ffp*fsmask)
    # sum(values(fsweight), na.rm=TRUE)/sum(values(ffp))
    # 
    # cls$id[cls$cover=="Pine"]
    # pmask<-wl.resamp; values(pmask)<-1; values(pmask)[values(wl.resamp!=4120)]<-0
    # pweight<-(ffp*pmask)
    # sum(values(pweight), na.rm=TRUE)/sum(values(ffp))
    # 
    # cls$id[cls$cover=="Coniferous Forested Wetland"]
    # cfmask<-wl.resamp; values(cfmask)<-1; values(cfmask)[values(wl.resamp!=6410)]<-0
    # cfweight<-(ffp*cfmask)
    # sum(values(cfweight), na.rm=TRUE)/sum(values(ffp))
    #####

    for(c in 1:nrow(compweights)){
      ind<-which(colnames(dat)==compweights$cover[c])
      if(length(ind)>0){dat[i,ind]<-compweights$pct[c]}
    }
    
    
    
  }else{
    print(paste("FFP did not work for", dat.ts[i])); next
  }
  #####
  }

  
  ##Regardless of FFP success, record other metadata####
  
  #input other metadata regardless if FFP success or failure
  sitedat$site[f]<-towername
  sitedat$igbp[f]<-site.pft
  sitedat$agemean[f]<-site.agemean
  sitedat$agemin[f]<-site.agemin
  sitedat$agemax[f]<-site.agemax
  
  traitdat[[f]]<-dat; names(traitdat)[f]<-towername
  
}

traitdat.bk<-traitdat

#assemble into a giant DF

master<-data.frame()
for (i in 1:length(traitdat)){
  
  d<-traitdat[[i]]
  d$sitename<-names(traitdat)[i]
  d$igbp<-sitedat$igbp[i]
  #dat$age<-sitedat$agemean
  
  if(i==1){master<-d}else{master<-merge(master, d, all=TRUE)}
  
}

#explore
library(plotly)
#esquisser()


write.csv(master, "alltower_trait_fullseason3.csv")


