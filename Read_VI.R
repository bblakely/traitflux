#get EVI data for dim analysis

library("sf")
library("terra")
library("raster")
library("dplyr")
library("stringr")

homedir<-"D:/analysis/traitflux"
datadir<-"D:/analysis/traitflux/rsdata/EVI"
setwd(homedir)


#QC LUT
cld.lut<-read.csv("HLSL30-020-Fmask-lookup.csv")
#GOOD values
vals.cld<-cld.lut$Value[cld.lut$Cloud=="No" & cld.lut$Cloud.shadow=="No"] #& cld.lut$Adjacent.to.cloud.shadow=="No"]

#Begin main program,

setwd(datadir)

files<-list.files()

str_subset(files, "B02")


b2.filelist<-str_subset(files, "B02")
b4.filelist<-str_subset(files, "B04")
b5.filelist<-str_subset(files, "B05")
qc.filelist<-str_subset(files, "Fmask")


#main loop

doylist<-usable<-rep(0,length(b2.filelist))
evi.list<-ndvi.list<-nirv.list<-list()

for(f in 1:length(b2.filelist)){
  
  #announce loop

  stamp<-substr(b2.filelist[f], 19,25)
  print(paste("DOY", substr(stamp,5,9)))
  
  #bands
  
  b2<-rast(b2.filelist[f])#*1000#*0.0001
  
  b4file<-grepl(stamp,b4.filelist)
  b4<-rast(b4.filelist[b4file])#*1000#*0.0001
  
  b5file<-grepl(stamp,b5.filelist)
  b5<-rast(b5.filelist[b5file])#*0.0001
  
  
  #calculations
  
  evi<- 2.5 * (b5 - b4) / ((b5 + 6*b4 - 7.5*b2) + 1); evi[abs(evi)>1]<-NA
  
  ndvi<- (b5 - b4) / (b5 + b4); ndvi[abs(ndvi)>1]<-NA
  
  nirv<-(ndvi-0.08) * b5;
  
  
  
  qcfile<-grepl(stamp,qc.filelist)
  qc<-rast(qc.filelist[qcfile])
  
  
  evi[!qc%in%vals.cld]<-NA
  ndvi[!qc%in%vals.cld]<-NA
  nirv[!qc%in%vals.cld]<- NA
  
  
  cov<-length(which(is.na(as.vector(evi))))/length(as.vector(evi))

  #print(paste(length(which(is.na(as.vector(evi2)))), "NA on", length(as.vector(evi2)), "pixels;", round(cov*100, 2), "%"))
  
  par(mfrow=c(2,2))
  plot(evi, main=paste("EVI DOY", substr(stamp,5,7), ":", round((1-cov)*100, 1), "% good"))
  plot(ndvi, main=paste("NDVI DOY", substr(stamp,5,7), ":", round((1-cov)*100, 1), "% good"))
  plot(nirv, main=paste("NIRV DOY", substr(stamp,5,7), ":", round((1-cov)*100, 1), "% good"))
  
  
  if(cov<0.95){usable[f]<-1}
  #usable<-1
  
  doylist[f]<-substr(stamp, 5, 7)
  
  evi.list[[f]]<-evi
  ndvi.list[[f]]<-ndvi
  nirv.list[[f]]<-nirv
  
}

evi.good<-rast(evi.list[which(usable==1)])
evi.good<-approximate(evi.good, method="linear");time(evi.good)<-as.numeric(doylist[which(usable==1)])



ndvi.good<-rast(ndvi.list[which(usable==1)])
ndvi.good<-approximate(ndvi.good, method="linear");time(ndvi.good)<-as.numeric(doylist[which(usable==1)])



nirv.good<-rast(nirv.list[which(usable==1)])
nirv.good<-approximate(nirv.good, method="linear");time(nirv.good)<-as.numeric(doylist[which(usable==1)])




doys<-substr(varnames(evi.good), 23,25)


