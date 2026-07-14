#hyspex lai processing


drivedir<-"D:/Analysis/traitflux/rsdata/hyspex_lai"
datelist<-c("0629", "0830")

for(t in c(1:length(datelist))){

  subdir<-paste(drivedir,datelist[t], sep="/")
  
  setwd(subdir)
  files<-list.files()
  
  #trait<-traitlist[t]
  laifiles<-files#[str_detect(files, traitlist[t])]
  hdrs<-laifiles[str_detect(laifiles,".hdr")]
  files<-setdiff(laifiles, hdrs)
  
  rm(mo) #clear out variable name so we don't get repeats
  
  print(paste("starting processing for lai",  datelist[t]))
  
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
    lairast<-rast[[1]]
    rg<-rast$range_mask; ndi<-rast$ndi; edg<-rast$neon_edge #cld<-rast$cloud; edg<-rast$neon_edge
    lairast[rg<1|ndi<1|edg<1]<-NA
    cleanrast[r]<-lairast
  }
  
  print("beginning mosaicing")
  
  rastlist<-sprc(cleanrast)
  
  mo<-mosaic(rastlist, fun="mean")
  
  
  plot(mo)
  
  writeRaster(mo, paste("D:/Analysis/traitflux/rsdata/laimosaic/lai_",datelist[t],".tif", sep=""), overwrite=TRUE)
  print(paste("raster has been written for", datelist[t]))
  
}

#exploratory plotting and such (to be deleted or commented out)
# 
# lai_0629<-rast("D:/Analysis/traitflux/rsdata/laimosaic/lai_0629.tif")
# lai_0830<-rast("D:/Analysis/traitflux/rsdata/laimosaic/lai_0830.tif")
# 
# lai_aug<-resample(lai_0830, lai_0629) #slow (~8sec)
# 
# 
# lai.chg<-lai_0629-lai_aug; 
# plot(lai.chg, range=c(-4,3))
# 
# hist(sample(values(lai.chg), 3000), breaks=30)
# mean(sample(values(lai.chg), 3000), na.rm=TRUE)

#lai an average of 0.2 higher in Aug, but variability much larger than mean change

