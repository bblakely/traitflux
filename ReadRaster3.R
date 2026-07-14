#New raster reader for all local + august traits

library(stringr)
library(raster)
library(rgl)
library(caTools)
library(terra)
library(EBImage)
library(jpeg)
library(fields)
library(zoom)
library(readxl)
library(terra)

setwd("D:/Analysis/traitflux")

#June traits, on this machine
drivedir<-"D:/traitdata/0629"
traitlist<-c("Lignin", "LMA", "Phenolics", "Nitrogen", "NSC")

for(t in c(1:length(traitlist))){
  
  setwd(drivedir)
  files<-list.files()
  
  trait<-traitlist[t]
  traitfiles<-files[str_detect(files, traitlist[t])]
  hdrs<-traitfiles[str_detect(traitfiles,".hdr")]
  xmls<-traitfiles[str_detect(traitfiles,".aux.xml")]
  files<-setdiff(traitfiles, hdrs); files<-setdiff(files, xmls)
  
  rm(mo) #clear out variable name so we don't get repeats
  
  print(paste("starting processing for", traitlist[t]))
  
  print("reading in files")
  allrast<-list()
  for (i in 1:(length(files))){ 
    #print(paste("reading", files[i]))
    readin<-try(rst<-terra::rast(files[i]))
    if(class(readin)=="try-error"){print(paste("could not read", files[i]))}else{allrast[[i]]<-rst}
    
  }
  
  
  
  #back to regular directory so I'm not fucking about in drive
  setwd("D:/Analysis/traitflux")
  
  #reduce list of rasters to ones that worked
  notempty<-sapply(allrast, function(x) !is.null(x))
  allrast<-allrast[notempty]
  
  print("conducting QC")
  
  #clean them up
  cleanrast<-list()
  for(r in 1:length(allrast)){
    print(paste("file",r))
    rast<-allrast[[r]]
    traitrast<-rast[[1]]
    rg<-rast$range_mask; ndi<-rast$ndi; cld<-rast$cloud; edg<-rast$neon_edge
    traitrast[rg<1|ndi<1|cld<1|edg<1]<-NA
    cleanrast[r]<-traitrast
  }
  
  print("beginning mosaicing")
  
  rastlist<-sprc(cleanrast)
  
  mo<-mosaic(rastlist, fun="mean")
  
  #Mosaic.Slow af
  # babymo1<-mosaic(cleanrast[[1]], cleanrast[[2]], cleanrast[[3]],cleanrast[[4]],cleanrast[[5]],cleanrast[[6]])
  # print("first submosaic complete")
  # babymo2<-mosaic(cleanrast[[7]], cleanrast[[8]], cleanrast[[9]],cleanrast[[10]],cleanrast[[11]],cleanrast[[12]])
  # print("second submosaic complete")
  # babymo3<-mosaic(cleanrast[[13]], cleanrast[[14]], cleanrast[[15]])#,cleanrast[[16]],cleanrast[[17]],cleanrast[[18]])
  #print("third submosaic complete, final mosaicing...")
  #mo<-mosaic(babymo1,babymo2, babymo3)
  
  plot(mo)
  
  writeRaster(mo, paste("D:/Analysis/traitflux/rsdata/traitmosaic/0629/",trait,".tif", sep=""), overwrite=TRUE)
  print(paste("raster has been written for", traitlist[t]))
  
}


#Same thing again for August
drivedir<-"D:/traitdata/0830"
traitlist<-c("Lignin", "LMA", "Phenolics", "Nitrogen", "NSC")

for(t in c(2:3)){#c(1:length(traitlist))){
  
  setwd(drivedir)
  files<-list.files()
  
  trait<-traitlist[t]
  traitfiles<-files[str_detect(files, traitlist[t])]
  hdrs<-traitfiles[str_detect(traitfiles,".hdr")]
  xmls<-traitfiles[str_detect(traitfiles,".aux.xml")]
  files<-setdiff(traitfiles, hdrs); files<-setdiff(files, xmls)
  
  rm(mo) #clear out variable name so we don't get repeats
  
  print(paste("starting processing for", traitlist[t]))
  
  print("reading in files")
  allrast<-list()
  for (i in 1:(length(files))){ 
    #print(paste("reading", files[i]))
    readin<-try(rst<-terra::rast(files[i]))
    if(class(readin)=="try-error"){print(paste("could not read", files[i]))}else{allrast[[i]]<-rst}
    
  }
  
  
  
  #back to regular directory so I'm not fucking about in drive
  setwd("D:/Analysis/traitflux")
  
  #reduce list of rasters to ones that worked
  notempty<-sapply(allrast, function(x) !is.null(x))
  allrast<-allrast[notempty]
  
  print("conducting QC")
  
  #clean them up
  cleanrast<-list()
  for(r in 1:length(allrast)){
    print(paste("file",r))
    rast<-allrast[[r]]
    traitrast<-rast[[1]]
    rg<-rast$range_mask; ndi<-rast$ndi; cld<-rast$cloud; edg<-rast$neon_edge
    traitrast[rg<1|ndi<1|cld<1|edg<1]<-NA
    cleanrast[r]<-traitrast
  }
  
  print("beginning mosaicing")
  
  rastlist<-sprc(cleanrast)
  
  mo<-mosaic(rastlist, fun="mean")
  
  
  plot(mo)
  
  writeRaster(mo, paste("D:/Analysis/traitflux/rsdata/traitmosaic/0830/",trait,".tif", sep=""), overwrite=TRUE)
  print(paste("raster has been written for", traitlist[t]))
  
}
