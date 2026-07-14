#Gradient boosting regression
  #This one does many model runs and extracts summaries of improvements to R2 
  #and extracts beeswarm and importance plots for each run
  
  #Jan 22 - color coded importance plots not working for ER or an VI runs
  setwd("D:/Analysis/traitflux")

#library(gbm)
library(shapviz)
library(xgboost)
library(dplyr)
library(SHAPforxgboost)
library(caret)
library(gridExtra)

#Data prep: read in, remove bad values ####

master<-read.csv("alltower_trait_fullseason3.csv")
master<-master[,which(colnames(master)!="X" & colnames(master)!="X.1")]
master.bk<-master

master$SG_1_1_1[abs(master$SG_1_1_1)>300]<-NA
master$SWC_1_1_1[abs(master$SWC_1_1_1)<20]<-NA
master$LE[master$LE<(-100)|master$LE>600]<-NA

#calculate canopy surface temp, assuming 99% emissivity

master$TC<-((master$LW_OUT_1_1_1/(5.67E-8*0.99))^0.25)-273.15

#####


#Choose days and land cover types ####

#master<-master[master$doy>=182 & master$doy<=258 ,] #days
master<-master[master$igbp!="WAT"&master$igbp!="WET"&master$igbp!="GRA",] #cover types

doylist<-c(182:243) #182:243 splits the distance between flights
master<-master[master$doy%in%doylist,]
#####

#####


#Data prep: cleaning and variable selection #####

#half hours with no footprint - this eliminates about half the data - 40% of daytime and 60% of nighttime
compind<-(103:114); traitind<-c(85:94)
#vegind<-which(colnames(master)%in%c("evi", "ndvi", "lai", "nirv", "nirvp"))
#metind<-which(colnames(master)%in%c("SW_IN","TA","VPD_C", "TS_1_1_1","SWC_1_1_1","WS_1_1_1","TC"))  


colnames(master[,compind]) #columns with composition (check)
colnames(master[,traitind]) #columns with traits (check)
nocomp<-which(rowSums(master[,compind])==0); notrait<-which(rowSums(master[,traitind], na.rm=TRUE)==0)
nolai<-which(is.na(master$lai))
novi<-which(is.na(master$evi))
missingvals<-unique(c(nocomp, notrait, nolai, novi))
master<-master[-missingvals,]



###Choose and apply trait interpolation policy#####

#set trait.policy as follows
# - "early" : June 29
# - "late"  : August 30
# - "nearest" : choose closest date
# - "interp" : interpolate between values based on date (all past Aug 30 get Aug 30 values)

trait.policy<-"early"

#get early season traits relabeled as such
es.traitind<-c(85:89); ls.traitind<-c(90:94)
colnames(master)[es.traitind]
colnames(master)[es.traitind]<-paste(colnames(master)[es.traitind], "es", sep="_")

early.labs<-colnames(master)[es.traitind]; early.labs
late.labs<-colnames(master)[ls.traitind]; late.labs

#apply policy
trait.df<-data.frame(matrix(nrow=nrow(master), ncol=length(es.traitind)))
colnames(trait.df)<-substr(early.labs,1,nchar(early.labs)-3)

trait.es<-master[,es.traitind]; trait.ls<-master[,ls.traitind]

if(trait.policy=="early"){trait.df<-master[colnames(master)%in%early.labs]}
if(trait.policy=="late"){trait.df<-master[colnames(master)%in%late.labs]}

if(trait.policy=="nearest"){
  print("selecting trait values")
  for(i in 1:nrow(master)){
    if(i%%1000==0){print(i)}
    doy<-master$doy[i]
    diff<-abs(doy-c(180,242)); if(diff[1]-diff[2]==0){diff[1]<-diff[1]+1} #for the day right in the middle, just give it late season
    if(which(diff==min(diff))==1){trait.df[i,]<-trait.es[i,]}else{trait.df[i,]<-trait.ls[i,]}
  }
}

if(trait.policy=="interp"){
  print("interpolating traits")
  for(i in 1:nrow(master)){
    
    if(i%%1000==0){print(i)}
    
    doy<-master$doy[i]
    ls.wt<-(doy-180)/(242-180); es.wt<-(1-ls.wt)
    if(ls.wt>1){ls.wt<-1; es.wt<-0}
    
    trait.df[i,]<-(trait.es[i,]*es.wt)+(trait.ls[i,]*ls.wt)
  }
}

colnames(trait.df)<-substr(early.labs,1,nchar(early.labs)-3)


master<-cbind(master, trait.df)


colnames(master)[242:246]



traitind<-which(colnames(master)%in%c("nsc","lignin","lma","nitr","pheno"))
vegind<-which(colnames(master)%in%c("evi", "ndvi", "lai", "nirv"))
metind<-which(colnames(master)%in%c("SW_IN","TA","VPD_C", "TS_1_1_1","SWC_1_1_1","WS_1_1_1","TC"))  


length(which(is.na(master$evi)))

#corrplot if ya want
library(corrplot); par(mfrow=c(1,1))
corrdat<-master[,c(compind, traitind, vegind)]
corrdat<-cor(corrdat)
corrplot(corrdat, type="upper", method="color", insig="pch") #nirvp drops because it has some NAs due to missing shortwave



#Create color coding by variable type (for later plotting)
comp.lab<-rep("forestgreen", length(compind)); names(comp.lab)<-colnames(master[,compind])
trait.lab<-rep("red4", length(traitind)); names(trait.lab)<-colnames(master[,traitind])
vi.lab<-rep("blue3", length(vegind)); names(vi.lab)<-colnames(master[,vegind])
met.lab<-rep("purple4", length(metind));names(met.lab)<-colnames(master[,metind])
var.lab<-c(comp.lab, trait.lab, vi.lab, met.lab, "doy"="darkgray", "nirvp"="darkblue")



#remove variable categories definitely not useful - seasonal, processing, etc. Also taking out silly traits here.####
patterns<-c("_qc", "_fall","_fmeth", "_fnum", "_fwin", "X", "_fqc", "fsd", "Thres", "SSITC", "SIGMA", "SPEC", "PRI","calcium","copper","alum","carbon", "fiber", "cellulose", "boron");patterns2<-c("PotRad", "_ref","daytime","TIMESTAMP_","posix", "FP_","sitename", "igbp", "season")
patterns<-c(patterns, patterns2) #option to not include the second list which I could see being useful at some point
nogood<-grep(paste(patterns,collapse="|"), colnames(master))

labels<-grep(paste(patterns2,collapse="|"), colnames(master)) #keep useful columns
data.lab<-master[,labels]

data.int<-master[,-nogood]; #remove all non-numeric ones from the analysis set

#####

#Choose fill value retention policy####

#rename gpp to not have f in it, and remove spots likey to be fill
data.int$GPP_BB<-data.int$GPP_uStar_f; data.int$GPP_BB[is.na(data.int$Reco_U50)]<-NA

orig.int<-c("_f", "_NEW") #to remove filled variables, i.e. use actual obs
fill.int<-c("_orig", "_NEW") #to remove orig and use fill
all.int<-c(orig.int, fill.int) #to remove both and only use raw data

rep<-grep(paste(all.int,collapse="|"), colnames(data.int))

data.unique<-data.int[,-rep]

#####

#Data prep: remove so-correlated-as-to-be-redundant variables:####

red.ind<-c("TS_1_2_1","TS_1_3_1","TS_1_4_1","TAU_1_1_1", "ZL_1_1_1","USTAR","RH_1_2_1","RH_1_3_1","RH",
           "T_SONIC_1_1_1","TA_1_2_1", "TA_1_3_1", "SW_OUT_1_1_1", "NETRAD", "SG_1_1_1", "WD_1_1_1","decdoy", "PA_1_1_1") #"LW_OUT_1_1_1", "decdoy"
red<-grep(paste(red.ind,collapse="|"), colnames(data.unique))
data.unique<-data.unique[,-red]

#now remove all columns with no data (leftover from merge with ultimately-discarded towers)

nodat<-which(colSums(data.unique, na.rm=TRUE)==0)
data.unique<-data.unique[,-nodat]

#####

#further variables, usually all false####
use.alb<-FALSE #albedo
use.laihr<-FALSE #high spatial res LAI (2 timepoints only)
use.seas<-FALSE #early or late season trait values

if(use.alb==FALSE){
  alb.int<-c("alb", "albhr")
  alb<-grep(paste(alb.int,collapse="|"), colnames(data.unique))
  data.unique<-data.unique[,-alb]
}

#if you do use the high res lai, don't use the low res one
if(use.laihr==FALSE){
  lai.int<-c("laihr")
  laihr<-grep(paste(lai.int,collapse="|"), colnames(data.unique))
  data.unique<-data.unique[,-laihr]
}else{
  lai.int<-c("lai")
  laihr<-which(colnames(data.unique)==lai.int)
  data.unique<-data.unique[,-laihr]
}

if(use.seas==FALSE){
  seas.int<-c("_es", "_ls")
  seas<-grep(paste(seas.int,collapse="|"), colnames(data.unique))
  data.unique<-data.unique[,-seas]
}

#####


# Model setup ####
results<-data.frame() #empty var for results
n.iter<-50 #number of iterations of model runs
nrounds<-500 #number of trees to do per model
plots<-FALSE
incl.nobio<-FALSE

#plot holders
er.bees<-gpp.bees<-le.bees<-h.bees<-list()
er.imp<-gpp.imp<-le.imp<-h.imp<-list()

#parameter list
xgb_params <- list(
  eta = 0.07,
  max_depth = 5,
  subsample = 0.6,
  colsample_bytree = 1,
  min_child_weight = 10, 
  gamma=0
  
)

traituse<-c(TRUE, FALSE); compuse<-c(TRUE, FALSE); eviuse<-c(TRUE, FALSE)
sets<-expand.grid(traituse, compuse, eviuse); colnames(sets)<-c("trait", "comp", "vi")
if(incl.nobio==FALSE){sets<-expand.grid(traituse, compuse);colnames(sets)<-c("trait", "comp")}



#do you want to hold out whole towers (but use every day of the record) 
#or whole days (but use every tower)?
split.type<-"tower" #doy 



data.unique.bk<-data.unique

#####

setlist.er<-list()

for(s in 1:nrow(sets)){

  #s<-1 #originally parameterized on full model
  
  data.unique<-data.unique.bk #reset input data columns
  
  gridrow<-sets[s,]
  
  #do you want traits and/or composition in the predictors?
  use.trait<-sets$trait[s]
  use.comp<- sets$comp[s]
  use.evi<-sets$vi[s]; if(incl.nobio==FALSE){use.evi=TRUE}
  ####
  
  
  #####
  
  #Apply choice of model driver sets####
  
  if(use.comp==FALSE){
    comp.int<-c("Fir.Spruce","Coniferous.Forested.Wetland", "Pine", "Aspen.Paper.Birch" , "Northern.Hardwoods","Broad.leaved.Deciduous.Scrub.Shrub", "Swamp.Hardwoods", "Broad.leaved.Evergreen.Scrub.Shrub", "Open.Water", "Red.Maple","Mixed.Deciduous.Coniferous.Forest", "Other")
    comp<-grep(paste(comp.int,collapse="|"), colnames(data.unique))
    data.unique<-data.unique[,-comp]
  }
  
  if(use.trait==FALSE){
    trait.int<-c("nsc","lignin","lma","nitr","pheno" )
    trait<-grep(paste(trait.int,collapse="|"), colnames(data.unique))
    data.unique<-data.unique[,-trait]
    
  }
  
  if(use.evi==FALSE){
    evi.int<-c("evi", "nirv", "nirvp", "ndvi")
    evi<-grep(paste(evi.int,collapse="|"), colnames(data.unique))
    data.unique<-data.unique[,-evi]
    
  }
  
  
  #label run type
  if(use.trait==TRUE & use.comp==FALSE & use.evi=="TRUE"){setlab<-"trait"}
  if(use.trait==FALSE & use.comp==TRUE & use.evi=="TRUE"){setlab<-"comp"}
  if(use.trait==TRUE & use.comp==TRUE & use.evi=="TRUE"){setlab<-"both"}
  if(use.trait==FALSE & use.comp==FALSE & use.evi=="TRUE"){setlab<-"vi"}
  
  if(use.trait==TRUE & use.comp==FALSE & use.evi=="FALSE"){setlab<-"trait.only"}
  if(use.trait==FALSE & use.comp==TRUE & use.evi=="FALSE"){setlab<-"comp.only"}
  if(use.trait==TRUE & use.comp==TRUE & use.evi=="FALSE"){setlab<-"both.only"}
  if(use.trait==FALSE & use.comp==FALSE & use.evi=="FALSE"){setlab<-"no.bio"}
  
  
  print(paste("starting runs for model type", setlab))
  
  #####
  
  
  
  #Function to remove NANs####
  kill.nan<-function(dat, var){
    
    var.ind<-which(colnames(dat)==var)
    
    nona.ind<-which(!is.na(dat.in[var.ind]))
    dat.out<-dat.in[nona.ind,]
    
  }
  
  
  #####
  
  
  #Modeling:
  
  

  
  #ER####
  
  print("starting ER")
  #Setup####
  #Flux-specific variable removal
  circ<-c("FC", "E_0", "resid", "GPP_", "LE","H", "NEE", "SC_", "Reco_uStar", "CO2","LW_IN_1_1_1", "LW_OUT_1_1_1", "G_1_1_1", "nirvp","nirv")
  circ.ind<-grep(paste(circ,collapse="|"), colnames(data.unique))
  dat.in<-data.unique[,-circ.ind]
  colnames(dat.in)
  
  
  dat.er<-kill.nan(dat.in, var="Reco_U50") #unsplit dataset
  dat.lab.er<-data.lab; dat.lab.er<-dat.lab.er[!is.na(dat.in$Reco_U50),]
  
  #take out weird towers
  norm<-which(!dat.lab.er$sitename%in%c("US-PFj", "US-PFh", "US-PFL"))
  dat.er<-dat.er[norm,]; dat.lab.er<-dat.lab.er[norm,]
  
  
  avail.twr<-unique(dat.lab.er$sitename)
  leaveout.count<-floor(length(avail.twr)*0.3)
  pairs<-combn(avail.twr, leaveout.count)
  #optionally, sample to the number of pairs for gpp
  #pairs<-pairs[,sample(1:ncol(pairs), min.pairs)]
  
  r2dat.er<-data.frame(matrix(nrow=ncol(pairs), ncol=2)); colnames(r2dat.er)<-c("train", "test")
  
  
#gpp grid search modified
  
  
  gridresult<-list()
  testresult<-rep(0, ncol(pairs))
  #result.grid<-data.frame(matrix(nrow=9), ncol=ncol(pairs))
  
  for(i in 1:ncol(pairs)){ #for each set of two left out towers
    print(paste("leaving out", pairs[1,i], "and", pairs[2,i], "; pair", i, "of", ncol(pairs)))
    leaveout<-c(pairs[,i])
    train.ind<-which(!dat.lab.er$sitename%in%leaveout) 
    
    
    y_all <- dat.er$Reco_U50
    X_all <- dat.er %>% select(-Reco_U50)
    
    xgb_all<-xgb.DMatrix(data = as.matrix(X_all), label = y_all)
    
    dat.train<-dat.er[train.ind,] #training data
    dat.test<-dat.er[-train.ind,] #testing data
    
    #xgboost
    y_train <- dat.train$Reco_U50
    y_test <- dat.test$Reco_U50
    X_train <- dat.train %>% select(-Reco_U50)
    X_test <- dat.test %>% select(-Reco_U50)
    
    xgb_train <- xgb.DMatrix(data = as.matrix(X_train), label = y_train)
    xgb_test <- xgb.DMatrix(data = as.matrix(X_test), label = y_test)
    
    
    ##Grid search for params####
    param.grid<-expand.grid(nrounds = c(1000),
                            max_depth = c(8, 10),
                            eta = c(0.01, 0.03),
                            gamma = c(0,10),
                            subsample = c(0.7),
                            min_child_weight = c(5,10,15),
                            colsample_bytree = c(0.7))
    
    train_control = trainControl(method = "cv", number = 5, search = "grid")
    
    start<-(Sys.time())
    # training a XGboost Regression tree model while tuning parameters
    xgb2 = train(Reco_U50~., data = dat.train, metric="Rsquared", method = "xgbTree", trControl = train_control,tuneGrid = param.grid, na.action = na.omit, verbosity=0)#tuneGrid = param.grid,
    print(Sys.time()-start)
    
#summarising the results
    print(xgb2)
    
    tune<-xgb2$bestTune; tune$leaveout1<-pairs[1,i]; tune$leaveout2<-pairs[2,i]
    
    gridresult[[i]]<-xgb2
    
    if(i==1){result.grid<-tune}else{result.grid<-rbind(result.grid, tune)}
    
    xgb.testmod<-xgb.train(
      data = xgb_train,
      params=as.list(gridresult[[i]]$bestTune),
      nrounds=1000,
      verbosity = 1
    )
    
    
    xgb_preds <- predict(xgb.testmod, as.matrix(X_test))
    rmse <- sqrt(mean((xgb_preds - y_test)^2)); #rmse
    r2<-cor(xgb_preds, y_test)^2; r2
    
    xgb_preds_tr<-predict(xgb.testmod, as.matrix(X_train))
    rmse.tr <- sqrt(mean((xgb_preds_tr - y_train)^2)); #rmse.tr
    r2.tr<-cor(xgb_preds_tr, y_train)^2; r2.tr
    
    testresult[i]<-r2
    
    
    #####
    
  }#comment out to do more than parameterize
  
  out.set<-cbind(result.grid, testresult)
  setlist.er[[s]]<-out.set
  
  }#to parameterize across model sets
  

#combine and write out results
results.both<-data.frame(setlist.er[[1]]$testresult);colnames(results.both)<-"test";results.both$lab<-"both"
results.both$leaveout1<-setlist.er[[1]]$leaveout1;results.both$leaveout2<-setlist.er[[1]]$leaveout2;

results.comp<-data.frame(setlist.er[[2]]$testresult);colnames(results.comp)<-"test"; results.comp$lab<-"comp"
results.comp$leaveout1<-setlist.er[[2]]$leaveout1;results.comp$leaveout2<-setlist.er[[2]]$leaveout2;

results.trait<-data.frame(setlist.er[[3]]$testresult);colnames(results.trait)<-"test"; results.trait$lab<-"trait"
results.trait$leaveout1<-setlist.er[[3]]$leaveout1;results.trait$leaveout2<-setlist.er[[3]]$leaveout2;

results.vi<-data.frame(setlist.er[[4]]$testresult);colnames(results.vi)<-"test"; results.vi$lab<-"vi"
results.vi$leaveout1<-setlist.er[[4]]$leaveout1;results.vi$leaveout2<-setlist.er[[4]]$leaveout2;

results.er.opt<-rbind(results.both, results.comp, results.trait, results.vi)

write.csv(results.er.opt, "optimized_er_results.csv", row.names=FALSE)

ggplot(results.er.opt) +
  aes(x = lab, y = test,colour = lab) +
  stat_summary(position=position_dodge(width=0.7), size=1.3,fun="median", fun.min="min", fun.max="max") +
  scale_fill_hue(direction = 1) +
  ylab("testing R2")+
  ylim(0, 1)+
  theme_minimal()+
  labs(color="model")+
  ylab(bquote("testing"~R^2))+
  theme(axis.title = element_text(size=22), 
        axis.text=element_text(size=16), 
        legend.title=element_text(size=20),
        legend.text=element_text(size=16))

 