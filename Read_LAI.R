library("sf")
library("terra")
library("raster")
library("dplyr")
library("stringr")

#MODIS

homedir<-"D:/analysis/traitflux"
datadir<-"D:/analysis/traitflux/rsdata/LAI_FPAR"
setwd(homedir)

# bounding box
bbox<-vect("trait_bbox.geojson"); bbox2<-project(bbox, CRS("+init=epsg:32615") )

#QC LUT
qc.lut<-read.csv("MCD15A3H-061-FparLai-QC-lookup.csv")
#GOOD values
vals.cld<-qc.lut$Value[qc.lut$CloudState=="Significant clouds NOT present (clear)" & qc.lut$MODLAND=="Good quality (main algorithm with or without saturation)"]# & cld.lut$Adjacent.to.cloud.shadow=="No"]

#Begin main program,

setwd(datadir)

files<-list.files()


lai.filelist<-str_subset(files, "Lai_500")
qc.filelist<-str_subset(files, "FparLai_QC")

doylist<-usable<-rep(0,length(lai.filelist))
lai.list<-list()

minlai.refdoy<-157 #day of a high-quality raster to use as a minimum plausible LAI. 157 = June 6

for(f in 1:length(lai.filelist)){
  stamp<-substr(lai.filelist[f], 26,32)
  print(paste("DOY", substr(stamp,5,9)))
  
  doy<-as.numeric(substr(stamp,5,9))
  
  lai<-rast(lai.filelist[f])#*0.0001
  
  qcfile<-grepl(stamp,qc.filelist)
  qc<-rast(qc.filelist[qcfile])

  
  lai[!qc%in%vals.cld]<-NA
  
  if(doy==minlai.refdoy){lai.minref<-lai}
  
  if(exists("lai.minref")&doy>minlai.refdoy&doy<241){lai[lai<(lai.minref-2)]<-NA}
  
  lai.proj<-project(lai, CRS("+init=epsg:32615"))
  lai.resamp<-disagg(x=lai.proj,fact=4,method="bilinear")
  
  lai.crop<-crop(lai.resamp, bbox2, snap="out")
  lai<-lai.crop
  
  
  cov<-length(which(is.na(as.vector(lai))))/length(as.vector(lai))
  
  #print(paste(length(which(is.na(as.vector(lai)))), "NA on", length(as.vector(lai)), "pixels;", round(cov*100, 2), "%"))
  
  plot(lai, main=paste("DOY", substr(stamp,5,7), ":", round((1-cov)*100, 1), "% good"), range=c(0,6))
  
  if(doy)
  if(cov<0.25){usable[f]<-1}
  #usable<-1
  
  doylist[f]<-substr(stamp, 5, 7)
  
  
  lai.list[[f]]<-lai
  
  
}

setwd(homedir)

lai.brick<-rast(lai.list)

time(lai.brick)<-as.numeric(doylist)

lai.good<-lai.brick[[which(usable==1)]]; #plot(lai.good, range=c(2,6.5)); lai.bk<-lai.good

lai.fill.es<-approximate(lai.good[[time(lai.good)<=241]]); #plot(lai.fill.es)
lai.fill.ls<-approximate(lai.good[[time(lai.good)>241]]); #plot(lai.fill.ls)
lai.fill<-c(lai.fill.es, lai.fill.ls); plot(lai.fill, range=c(2,6.5))

lai.fill<-focal(lai.fill, w=5, fun="mean", na.policy="only"); plot(lai.fill, range=c(2,6.8))

#lai.fill<-approximate(lai.good); plot(lai.fill)



lai.good<-lai.fill


##for comparison to hyspex lai
# diff<-lai.good[[which(time(lai)==149)]]-lai.good[[which(time(lai)==241)]]
# plot(diff)
# hist(values(diff))
# mean(values(diff))

#less variability, larger average difference


