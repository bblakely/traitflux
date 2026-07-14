#AF Cheesehead

library(stringr)
library(REddyProc)
library(dplyr)
library(readxl)
library(bigleaf)


setwd("D:/Analysis/traitflux")
wd<-getwd()

#read in pheno data
#pheno<-read_excel("Cheesehead_Tree_Phenology_Data_full.xlsx")

setwd("CHF")
files<-filelist<-list.files()
print(paste("There are", length(files), "raw files"))
setwd(wd)

year<-2019
plots<-TRUE #set to true to generate plots

out.er<-out.gpp<-out.carb<-data.frame(matrix(nrow=14592, ncol=length(filelist), data=NA)); colnames(out.carb)<-substr(filelist, 5,10)

for(i in 1:length(files)){
  
#Read in data#####
  
  print("###########")
  print(paste("starting", filelist[i]))
  print(paste("( site", i, ")"))
  
  setwd("CHF")
  dat.raw<-read.csv(filelist[i], skip=2, header=TRUE) #choose which file, will loop eventually
  dat<-dat.raw; dat[dat==-9999]<-NA
  setwd(wd)
  
  failcount<-0
  #get timestamp into posix form
  
  ts.raw<-dat$TIMESTAMP_START
  ts.posix<-as.POSIXct(as.character(ts.raw), #start time based posix ts
                       format="%Y%m%d%H%M", tz="UTC")
  ts.posix.e<-as.POSIXct(as.character(dat$TIMESTAMP_END), #end time based posx ts (reddyproc wants this one)
                         format="%Y%m%d%H%M", tz="UTC")
  
  #put the posix end time in the flux dataset
  dat$ts_posix<-ts.posix.e

  #####
  
  #Check variables and rename where necessary#####  
  
  #Check variables and adjust names accordingly 
  #(i.e. use 1_1_1 if no basic vars available, or _PI_F if there are no raw vars avaiulbale)  
  
  #list of varaibles necessary for flux processing    
  vars<-c("FC", "LE", "H", "SW_IN", "USTAR", "TA","RH") #originally had netrad, ts, swc
  
  datvars<-colnames(dat)
  
  if(all(vars%in%datvars)==FALSE){
    
    missing<-(vars[which(vars%in%datvars==FALSE)])
    print("missing variables named")
    print(missing)
    
    
    for (v in 1:length(missing)){
      
      #a few weird ones...
      if(missing[v]=="TA" & "T_SONIC"%in%colnames(dat) & !"TA_1_1_1"%in%colnames(dat)){
        print("replacing name T_SONIC with name TA")
        colnames(dat)[colnames(dat)=="T_SONIC"]<-"TA"
        next}
      
    
      
      opt<-which(substr(datvars, 0, nchar(missing[v]))==missing[v])
      
      if(length(opt)==0){ #ie there are no options for a necessary variable
        print(paste("variable", missing[v], "is entirely absent from dataset", filelist[i]))
        failcount<-1
        next
      }
      
      #print(datvars[opt])
      
      if(length(datvars[opt])>1){ #if there's more than one option...
        
        #check for 1_1_1 positional qualifier
        
        list<-datvars[opt]
        var.syn<-paste(missing[v],"_1_1_1", sep="")  #prioritize variables with positional qualifiers
        
        if(var.syn%in%list){
          
          ind<-which(colnames(dat)==var.syn)
          print(paste("replacing name", colnames(dat)[ind], "with name", missing[v] ))
          colnames(dat)[ind]<-missing[v]
          
        }else{ #if no version with a positional qualifier, look for a PI-filledvariable
          
          var.syn.f<-paste(missing[v],"_PI_F", sep="")
          
          if(var.syn.f%in%list){
            
            ind<-which(colnames(dat)==var.syn.f)
            print(paste("replacing name", colnames(dat)[ind], "with name", missing[v] ))
            
            colnames(dat)[ind]<-missing[v]
            
            
          }else{ #if no pi_f version, try for just a regular pi version
            
            var.syn.pi<-paste(missing[v],"_PI", sep="")
            
            if(var.syn.pi%in%list){
              
              ind<-which(colnames(dat)==var.syn.pi)
              print(paste("replacing name", colnames(dat)[ind], "with name", missing[v] ))
              colnames(dat)[ind]<-missing[v]
              
            } else{print(paste("no alternative variable could be found for", missing[v]))
              failcount<-1}
            
          }
          
        }
        
        
        
        
        
      }else{
        
        if(length(opt)==1){
          print(paste("replacing name", colnames(dat)[opt], "with name",missing[v] ))
          colnames(dat)[opt]<-missing[v]
        }
        
        #should add some kind of flag here
        
      }
      
    }
    
    if(failcount==1){
      print(paste("dataset", filelist[i], "is missing an essential variable; skipping to next dataset"))
      next } #if an essential variable was missing, stop further processing
    
    
  }else{print(paste("dataset", filelist[i], "has all necessary variables"))}
  
  #test for all-NA columns
  dat.test<-dat[colnames(dat)%in%vars]
  means<-colMeans(dat.test, na.rm=TRUE)
  if(any(is.na(means))){
    print(paste("dataset", filelist[i], "has all NA for a necessary variable; skipping to next dataset"))
    next
  }
  
  #some basic plots for general idea of data continuity
  # if (plots==TRUE){
  #   par(mfrow=c(2,2))
  #   plot(dat$FC~ts.posix.e, main="carbon flux")
  #   plot(dat$LE~ts.posix.e, main="LE")
  #   plot (dat$SW_IN~ts.posix.e, main="shortwave in")
  #   plot (dat$TA~ts.posix.e, main="air temp")
  #   
  # }
  
#} #comment in or out to do filtering without processing
  
#####
  
#Process fluxes
  #Calculate VPD:
  dat$VPD_C<-fCalcVPDfromRHandTair(rH = dat$RH, Tair = dat$TA); dat$VPD_C[dat$VPD_C<0]<-0 #for error at high humidity
  
  
  #Make the instance of processing
  ep<-sEddyProc$new(ID="TST", Data=dat, c("FC", "LE", "H", "SW_IN", "USTAR", "TA","RH", "VPD_C"), ColPOSIXTime="ts_posix", DTS=48) #had netrad, swc, ts
  
  
  #reddyproc plots
  par(mfrow=c(1,1))
  ep$sPlotFingerprintY("FC", Year=year)
  #ep$sPlotDiurnalCycle("FC") #this one saves plots out to file
  #ep$sPlotDiurnalCycle("SW_IN")
  #ep$sPlotDiurnalCycle("TA")
  
  #ustar estimation
  #seasonsplit<-as.factor(c(2017001, 2017004, 2017006, 2017009, 2017012)) #doesn't currently work, but a more reasonable split for arctic ecosystems
  try(ep$sEstimateUstarScenarios(UstarColName = "USTAR", NEEColName = "FC",TempColName = "TA", RgColName = "SW_IN",
                             nSample = 100L, probs = c(0.5), ctrlUstarSub = usControlUstarSubsetting(taClasses=3, minRecordsWithinTemp = 50,  minRecordsWithinSeason =80)))
  
  if(is.na(all(ep$sGetUstarScenarios()$uStar))){ #if ustar didn't work, remake the instance and calcualte a single u*
    ep<-sEddyProc$new(ID="TST", Data=dat, c("FC", "LE", "H", "SW_IN", "USTAR", "TA","RH", "VPD_C"), ColPOSIXTime="ts_posix", DTS=48) #had netrad, swc, ts
    ep$sEstUstarThold(UstarColName = "USTAR", NEEColName = "FC",TempColName = "TA", RgColName = "SW_IN")}
                    
  #Prints out ustar threshholds, if inspection is desired
  ep$sGetEstimatedUstarThresholdDistribution() 
  ep$sGetUstarScenarios()
  
  #Gapfilling
  # The "ustarscens" just means it's doing gapfilling for multiple levels of ustar.
  
  try(ep$sMDSGapFillUStarScens(fluxVar='FC', uStarVar="USTAR", RgColName="SW_IN",swThr=50,
                           V1="SW_IN", T1=50, V2="RH", T2=20, V3="TA", T3=2.5))
  
  
  
  try(ep$sMDSGapFillUStarScens(fluxVar='LE', uStarVar="USTAR", RgColName="SW_IN",swThr=50,
                           V1="SW_IN", T1=50, V2="RH", T2=20, V3="TA", T3=2.5))
  
  
  try(ep$sMDSGapFillUStarScens(fluxVar='H', uStarVar="USTAR", RgColName="SW_IN",swThr=50,
                           V1="SW_IN", T1=50, V2="RH", T2=20, V3="TA", T3=2.5))
  
  
  fcplot<-try(ep$sPlotFingerprintY('FC_U50_f', Year = year))
  leplot<-try(ep$sPlotFingerprintY('LE_U50_f', Year = year))
  
  
  if(class(leplot)=="try-error"){
    leplot<-try(ep$sPlotFingerprintY('LE_uStar_f', Year = year))
    }
  
  if(class(leplot)=="try-error"){
    print(paste("Eddypro processing failed for",sitename))
    next}
  

  #ep$sPlotDiurnalCycle("FC_U50_f")
  
  
  #Prepare filled met data for partitioning
  #read the badm for location info
  
  # setwd("/Users/bethanyblakely/Desktop/Analysis/ABOVE/Sites1")
  # 
  # sitename<-substr(filelist[i],1,15)
  # badmdir<-(list.files(pattern=sitename))
  # setwd(badmdir)
  # badm<-read_excel(list.files())
  # lat<-as.numeric(badm$DATAVALUE[badm$VARIABLE=="LOCATION_LAT"])
  # lon<-as.numeric(badm$DATAVALUE[badm$VARIABLE=="LOCATION_LONG"])
  # 
  # setwd(paste("/Users/bethanyblakely/Desktop/Analysis/ABOVE/writeout/",year, sep=""))
  
  #Shortcut: just pick a cheesehead adjacent spot. All towers close together
  
  lat<-45.89; lon<-(-90.21)
  
  #continue with processing
  ep$sSetLocationInfo(LatDeg = 	lat, LongDeg = lon, TimeZoneHour = 0)  
  ep$sMDSGapFill('TA', FillAll = FALSE,  minNWarnRunLength = NA)     
  ep$sMDSGapFill('RH', FillAll = FALSE,  minNWarnRunLength = NA) 
  ep$sMDSGapFill('SW_IN', FillAll = FALSE,  minNWarnRunLength = NA) 
  ep$sMDSGapFill('VPD_C', FillAll = FALSE,  minNWarnRunLength = NA) 
  

  
  
  partition<-try(ep$sc(FluxVar="FC_U50_f", QFFluxVar="FC_U50_fqc",TempVar="TA_f",QFTempVar="TA_fqc", RadVar="SW_IN"))
  #ep$sGLFluxPartition(NEEVar="FC_U50_f", QFNEEVar="FC_U50_fqc",NEESdVar.s= "FC_U50_fsd", 
  #                    TempVar="TA_f",QFTempVar="TA_fqc", VPDVar = "VPD_C_f", QFVPDVar= "VPD_C_fqc", RadVar="SW_IN",QFRadVar= "SW_IN_fqc")
  
  if(class(partition)=="try-error"){ #if you had to do single ustar, _U50 aren't present. Use _uStar. 
    partition<-try(ep$sMRFluxPartitionUStarScens(FluxVar="FC_uStar_f", QFFluxVar="FC_uStar_fqc",TempVar="TA_f",QFTempVar="TA_fqc", RadVar="SW_IN"))
    
  }
  
  if(class(partition)=="try-error"){print(paste("Partitioning has failed for",sitename))}
  
  #ep$sPlotFingerprintY('GPP_U50_f', Year = year)
  #ep$sPlotFingerprintY('Reco_U50', Year = year)
  
  dat.derived<-ep$sExportResults()
  dat.out<-cbind(dat, dat.derived)
  sitename<-substr(files[i], 5,10)
  
  #Check fill
  
  par(mfrow=c(2,2), mar=c(3,3,1,1))
  dat.clip<-dat.derived[8000:13800,]
  
  plot(dat.clip$FC_U50_f, col='forest green', type='l')
  points(dat.clip$FC_U50_orig, pch="+")

  plot(dat.clip$TA_f, col='darkred', type='l')
  points(dat.clip$TA_orig, pch="+")
  text(1000,0, sitename)
  
  plot(dat.clip$RH_f, col='darkblue', type='l')
  points(dat.clip$RH_orig, pch="+")

  plot(dat.clip$SW_IN_f, col='orange', type='l')
  points(dat.clip$SW_IN_orig, pch="+")
  
  
  #put the carbon into the dataframe for basic plotting after the fact
  out.carb[,i]<-dat.out$FC_U50_f
  out.gpp[,i]<-dat.out$GPP_U50_f
  out.er[,i]<-dat.out$Reco_U50
  
  write.csv(dat.out, file=paste("processed/",sitename, "_processed.csv", sep=""))
}

