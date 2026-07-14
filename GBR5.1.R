#Gradient boosting regression
#single run through of each model to get plots quickly; use in conjunction with GBR_summaries

setwd("D:/Analysis/traitflux")

#library(gbm)
library(shapviz)
library(xgboost)
library(SHAPforxgboost)
library(caret)
library(pals)
library(gridExtra)
library(dplyr)

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


#Data prep: cleaning and variable selection #####

#half hours with no footprint - this eliminates about half the data - 40% of daytime and 60% of nighttime
compind<-(103:114); traitind<-c(85:94)

colnames(master[,compind]) #columns with composition (check)
colnames(master[,traitind]) #columns with traits (check)
nocomp<-which(rowSums(master[,compind])==0); notrait<-which(rowSums(master[,traitind], na.rm=TRUE)==0)
nolai<-which(is.na(master$lai))
novi<-which(is.na(master$ndvi))
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


#####

traitind<-which(colnames(master)%in%c("nsc","lignin","lma","nitr","pheno"))
vegind<-which(colnames(master)%in%c("evi", "ndvi", "lai", "nirv"))

length(which(is.na(master$evi)))



#corrplot if ya want
library(corrplot); par(mfrow=c(1,1))
corrdat<-data.frame(master[,c(compind, traitind, vegind)])
corrdat<-cor(corrdat)
corrplot(corrdat, type="upper", method="color", insig="pch")

 corrdat2<-master[,c(compind, traitind, vegind)]
 #corrdat2<-corrdat2[which(!colnames(corrdat2)%in%c("boron", "carbon", "calcium","copper","alum", "fiber", "cellulose"))]
 corrdat2<-cor(corrdat2)# corrplot(corrdat2, type="upper", method="number", number.cex=0.5)
 corrplot(corrdat2, type="upper", method="color", insig="pch", tl.col="black")
 
#Data preparaton - additional#### 
 
#remove variable categories definitely not useful - seasonal, processing, etc. Also taking out silly traits here.
patterns<-c("_qc", "_fall","_fmeth", "_fnum", "_fwin", "X", "_fqc", "fsd", "Thres", "SSITC", "SIGMA", "SPEC", "PRI","calcium","copper","alum","fiber", "cellulose", "boron");patterns2<-c("PotRad", "_ref","daytime","TIMESTAMP_","posix", "FP_","sitename", "igbp", "season")
patterns<-c(patterns, patterns2) #option to not include the second list which I could see being useful at some point
nogood<-grep(paste(patterns,collapse="|"), colnames(master))

labels<-grep(paste(patterns2,collapse="|"), colnames(master)) #keep useful columns
data.lab<-master[,c(labels, which(colnames(master)=="nsc"), which(colnames(master)=="evi"))] #nsc has all unique values; use something like an ID code later for merging

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
           "T_SONIC_1_1_1","TA_1_2_1", "TA_1_3_1", "SW_OUT_1_1_1", "NETRAD", "SG_1_1_1","decdoy", "PA_1_1_1","WD_1_1_1") #"LW_OUT_1_1_1", "decdoy" "WD_1_1_1"
red<-grep(paste(red.ind,collapse="|"), colnames(data.unique))
data.unique<-data.unique[,-red]

#now remove all columns with no data (leftover from merge with ultimately-discarded towers)

nodat<-which(colSums(data.unique, na.rm=TRUE)==0)
data.unique<-data.unique[,-nodat]

data.unique.bk<-data.unique
#####


# Model setup ####

results<-data.frame() #empty var for results
n.iter<-50 #number of iterations of model runs
nrounds<-500 #number of trees to do per model
plots<-FALSE

#parameter list
xgb_params <- list(
  eta = 0.07,
  max_depth = 5,
  subsample = 0.6,
  colsample_bytree = 1,
  min_child_weight = 10, 
  gamma=0
  
)

use.trait<-TRUE
use.comp<-FALSE
use.evi<-TRUE

#do you want to hold out whole towers (but use every day of the record) or whole days (but use every tower)?
split.type<-"doy" #doy 


#further variables, usually all false
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
  evi.int<-c("evi", "nirv", "nirvp", "ndvi" )
  evi<-grep(paste(evi.int,collapse="|"), colnames(data.unique))
  data.unique<-data.unique[,-evi]
  
}


#label run type
if(use.trait==TRUE & use.comp==FALSE & use.evi=="TRUE"){setlab<-"trait"}
if(use.trait==FALSE & use.comp==TRUE & use.evi=="TRUE"){setlab<-"comp"}
if(use.trait==TRUE & use.comp==TRUE & use.evi=="TRUE"){setlab<-"both"}
if(use.trait==FALSE & use.comp==FALSE & use.evi=="TRUE"){setlab<-"evi"}

if(use.trait==TRUE & use.comp==FALSE & use.evi=="FALSE"){setlab<-"trait.only"}
if(use.trait==FALSE & use.comp==TRUE & use.evi=="FALSE"){setlab<-"comp.only"}
if(use.trait==TRUE & use.comp==TRUE & use.evi=="FALSE"){setlab<-"both.only"}
if(use.trait==FALSE & use.comp==FALSE & use.evi=="FALSE"){setlab<-"no.bio"}

data.plotexp<-cbind(data.unique, data.lab$sitename, data.lab$igbp); colnames(data.plotexp)[c(ncol(data.plotexp)-1, ncol(data.plotexp))]<-c("sitename", "igbp")

print(paste("starting runs for model type", setlab))

#####


#Function to remove NANs####
kill.nan<-function(dat, var){

  var.ind<-which(colnames(dat)==var)

  nona.ind<-which(!is.na(dat.in[var.ind]))
  dat.out<-dat.in[nona.ind,]

}

#####


#testing area#

# dat.lab.2<-data.lab; colnames(dat.lab.2)[19:20]<-c("nsc_lab","evi_lab")
# dat<-cbind(dat.lab.2, data.unique)
# 
# dat$ts<-as.POSIXct(dat$ts_posix, format="%Y-%m-%d %H:%M:%S")




#Modeling:


#GPP#####

circ<-c("FC", "E_0", "resid", "GPP_uStar_f", "LE","H", "NEE", "SC_", "Reco","CO2", "nirvp", "G_1_1_1", "LW_IN_1_1_1", "LW_OUT_1_1_1", "nirv")
circ.ind<-grep(paste(circ,collapse="|"), colnames(data.unique))
dat.in<-data.unique[,-circ.ind]
colnames(dat.in)


dat.gpp<-kill.nan(dat=dat.in, var="GPP_BB"); 
dat.lab.gpp<-data.lab; dat.lab.gpp<-dat.lab.gpp[!is.na(dat.in$GPP_BB),]


#withhold full days for testing
if(split.type=="doy"){
doy.train<-sample(unique(dat.gpp$doy), round(0.7*length(unique(dat.gpp$doy)))) #use 70% of days (i.e. withhold whole days) to make it harder for the model
train.ind<-which(dat.gpp$doy%in%doy.train) # sample(1:nrow(dat.gpp), 0.7*nrow(dat.gpp))
}

#withhold full towers for testing
if(split.type=="tower"){
twr.train<-sample(unique(dat.lab.gpp$sitename), round(0.8*length(unique(dat.lab.gpp$sitename))))
train.ind<-which(dat.lab.gpp$sitename%in%twr.train) # sample(1:nrow(dat.gpp), 0.7*nrow(dat.gpp))
}


y_all <- dat.gpp$GPP_BB
X_all <- dat.gpp %>% dplyr::select(-GPP_BB)


xgb_all<-xgb.DMatrix(data = as.matrix(X_all), label = y_all)


dat.train<-dat.gpp[train.ind,] #training data
dat.test<-dat.gpp[-train.ind,] #testing data

print(unique(dat.lab.gpp$sitename[-train.ind]))

#xgboost
y_train <- dat.train$GPP_BB
y_test <- dat.test$GPP_BB
X_train <- dat.train %>%  dplyr::select(-GPP_BB)
X_test <- dat.test %>%  dplyr::select(-GPP_BB)

xgb_train <- xgb.DMatrix(data = as.matrix(X_train), label = y_train)
xgb_test <- xgb.DMatrix(data = as.matrix(X_test), label = y_test)


#parameter list
xgb_params <- list(
  eta = 0.05,
  max_depth = 8,
  subsample = 0.6,
  colsample_bytree = 0.6,
  min_child_weight = 10, 
  gamma=10
  
)


######


xgb_model <- xgb.train( #params = xgb_params,
  data = xgb_train,
  params=xgb_params,
  verbosity = 1,
  nround=500
  
)

xgb_preds <- predict(xgb_model, as.matrix(X_test))
rmse <- sqrt(mean((xgb_preds - y_test)^2)); #rmse
r2<-cor(xgb_preds, y_test)^2; r2

xgb_preds_tr<-predict(xgb_model, as.matrix(X_train))
rmse.tr <- sqrt(mean((xgb_preds_tr - y_train)^2)); #rmse.tr
r2.tr<-cor(xgb_preds_tr, y_train)^2; r2.tr

r2rec<-c(r2.tr, r2)

#xgb.importance(xgb_model)
#Plots####
#if(plots==TRUE){


viz<-shapviz(xgb_model, X_pred=as.matrix(X_train))#, interactions=TRUE)
#imp<-sv_importance(viz, kind="beeswarm"); imp

imp2<-sv_importance(viz, kind="bar", show_numbers = TRUE); imp2+theme_minimal()
#imp2<-sv_waterfall(viz,row_id = 3); imp2

  sv_dependence(viz, "SW_IN", "doy")+
    ylab("GPP SHAP value")+
    theme_minimal()+
    theme(axis.title = element_text(size=24), 
          axis.text=element_text(size=22), 
          legend.position="none")
  
    
  dev.copy(png, "D:/Analysis/traitflux/plots/sw_gpp.png", height=600, width=400)
  dev.off()

 
 if(setlab=="trait"){viz$X<-merge(viz$X, data.lab, by=c("nsc", "evi"), sort=FALSE, no.dups=FALSE)}
 if(setlab=="comp"){viz$X<-merge(viz$X, data.lab, by="evi", sort=FALSE)}
 
 viz$X$sitecode<-as.factor(viz$X$sitename)
 mycol<-glasbey(n=length(levels(viz$X$sitecode)))
 names(mycol)<-levels(viz$X$sitecode)
 
 viz$X$pftcode<-as.factor(viz$X$igbp)
 igbpcol<-glasbey(n=length(levels(viz$X$pftcode)))
 names(igbpcol)<-levels(viz$X$pftcode)
 
 samp<-sample(1:nrow(viz$X), 3000)
 
 #sv_interaction(viz, "lai", kind="bar", max_display=10)
 
 #trait plots####
 
 
 gpp.lma<-sv_dependence(viz[samp],"lma")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
   scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+
   scale_color_manual(values=mycol, name="site")+ylim(-2,3)+ylab("GPP SHAP value")+
   theme_minimal()+
   theme(axis.title = element_text(size=20), 
         axis.text=element_text(size=18),
         legend.position="none")
 
 gpp.nitr<-sv_dependence(viz[samp],"nitr")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
   scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+
   scale_color_manual(values=mycol, name="site")+ylim(-2,3)+xlab("nitrogen")+ylab("GPP SHAP value")+
   theme_minimal()+
   theme(axis.title = element_text(size=20), 
         axis.text=element_text(size=18),
         legend.position="none")
 

 grid.arrange(gpp.lma, gpp.nitr)
 
 dev.copy(png, "D:/Analysis/traitflux/plots/gpp_les.png", height=550, width=400)
 dev.off()
 
 # sv_dependence(viz[samp],"nitr")+aes(shape=viz$X$pftcode[samp], color=viz$X$pftcode[samp])+
 #   scale_shape_manual(values=1:nlevels(viz$X$pftcode[samp]), name="pft")+
 #   scale_color_manual(values=igbpcol, name="pft")+
 #   theme_minimal()
 # 
 
 gpp.nsc<-sv_dependence(viz[samp],"nsc")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
   scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+
   scale_color_manual(values=mycol, name="site")+ylab("GPP SHAP value")+ylim(-3,6)+
   theme_minimal()+
   theme(axis.title = element_text(size=20), 
         axis.text=element_text(size=18),
         legend.position="none")
 
 gpp.pheno<-sv_dependence(viz[samp],"pheno")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
   scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+ylim(-3,6)+
   scale_color_manual(values=mycol, name="site")+ylab("GPP SHAP value")+xlab("phenolics")+
   theme_minimal()+
   theme(axis.title = element_text(size=20), 
         axis.text=element_text(size=18),
         legend.position="none")
 
 gpp.lignin<-sv_dependence(viz[samp],"lignin")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
   scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+ylim(-3,6)+
   scale_color_manual(values=mycol, name="site")+ylab("GPP SHAP value")+
   theme_minimal()+
   theme(axis.title = element_text(size=20), 
         axis.text=element_text(size=18),
         legend.position="none")
 
 ggarrange(gpp.lignin,gpp.nsc, gpp.pheno)
 
 dev.copy(png, "D:/Analysis/traitflux/plots/gpp_alloc.png", height=600, width=350)
 dev.off()
 
 #LAI plots#####
 
 gpp.lai<-sv_dependence(viz[samp],"lai")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
   scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+
   scale_color_manual(values=mycol, name="site")+
   theme_minimal(); gpp.lai
 
 gpp.lai.pft<-sv_dependence(viz[samp],"lai")+aes(shape=viz$X$pftcode[samp], color=viz$X$pftcode[samp])+
   scale_shape_manual(values=1:nlevels(viz$X$pftcode[samp]), name="pft")+
   scale_color_manual(values=igbpcol, name="pft")+
   theme_minimal(); gpp.lai.pft
 
 grid.arrange(gpp.lai, gpp.lai.pft)
 
 dev.copy(png, "D:/Analysis/traitflux/plots/gpp_lai.png", height=550, width=500)
 dev.off()

 #####

 #comp plots####
 
 sv_dependence(viz[samp],"Pine")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
   scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+
   scale_color_manual(values=mycol, name="site")+ylim(-2,2)+
   theme_minimal()
 
 
 sv_dependence(viz[samp],"Fir.Spruce")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
   scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+
   scale_color_manual(values=mycol, name="site")+ylim(-2,2)+
   theme_minimal()
 
 sv_dependence(viz[samp],"Aspen.Paper.Birch")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
   scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+
   scale_color_manual(values=mycol, name="site")+
   theme_minimal()
 
 sv_dependence(viz[samp],"Broad.leaved.Deciduous.Scrub.Shrub")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
   scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+
   scale_color_manual(values=mycol, name="site")+
   theme_minimal()
 
 sv_dependence(viz[samp],"Other")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
   scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+
   scale_color_manual(values=mycol, name="site")+
   theme_minimal()
 
 sv_dependence(viz[samp],"Coniferous.Forested.Wetland")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
   scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+
   scale_color_manual(values=mycol, name="site")+
   theme_minimal()
 #####

rm("dat.in")

#####



#ER####

print("starting ER")

#Flux-specific variable removal
circ<-c("FC", "E_0", "resid", "GPP_", "LE","H", "NEE", "SC_", "Reco_uStar", "CO2","LW_IN_1_1_1", "LW_OUT_1_1_1", "G_1_1_1", "nirvp","nirv")
circ.ind<-grep(paste(circ,collapse="|"), colnames(data.unique))
dat.in<-data.unique[,-circ.ind]
colnames(dat.in)


dat.er<-kill.nan(dat.in, var="Reco_U50") #unsplit dataset
dat.lab.er<-data.lab; dat.lab.er<-dat.lab.er[!is.na(dat.in$Reco_U50),]


#try taking out pfj, pfh
#dat.er<-dat.er[dat.lab.er$sitename!="US-PFj" & dat.lab.er$sitename!="US-PFh",]
#dat.lab.er<-dat.lab.er[dat.lab.er$sitename!="US-PFj"& dat.lab.er$sitename!="US-PFh",]


#withhold full days for testing
if(split.type=="doy"){
doy.train<-sample(unique(dat.er$doy), round(0.7*length(unique(dat.er$doy)))) #use 70% of days (i.e. withhold whole days) to make it harder for the model
train.ind<-which(dat.er$doy%in%doy.train) # sample(1:nrow(dat.er), 0.7*nrow(dat.er))
}

#withhold full towers for testing
if(split.type=="tower"){
  twr.train<-sample(unique(dat.lab.er$sitename), round(0.8*length(unique(dat.lab.er$sitename))))
  train.ind<-which(dat.lab.er$sitename%in%twr.train)
}


y_all <- dat.er$Reco_U50
X_all <- dat.er %>% dplyr::select(-Reco_U50)

xgb_all<-xgb.DMatrix(data = as.matrix(X_all), label = y_all)


dat.train<-dat.er[train.ind,] #training data
dat.test<-dat.er[-train.ind,] #testing data

if(split.type=="tower"){print(unique(dat.lab.er$sitename[-train.ind]))}

#xgboost
y_train <- dat.train$Reco_U50
y_test <- dat.test$Reco_U50
X_train <- dat.train %>% dplyr::select(-Reco_U50)
X_test <- dat.test %>% dplyr::select(-Reco_U50)

xgb_train <- xgb.DMatrix(data = as.matrix(X_train), label = y_train)
xgb_test <- xgb.DMatrix(data = as.matrix(X_test), label = y_test)

#parameter list
xgb_params <- list(
  eta = 0.07,
  max_depth = 5,
  subsample = 0.6,
  colsample_bytree = 1,
  min_child_weight = 10, 
  gamma=0
  
)

xgb_model <- xgb.train( #params = xgb_params,
  data = xgb_train,
  params=xgb_params,
  nrounds=500,
  verbosity = 1,

)

xgb_preds <- predict(xgb_model, as.matrix(X_test))
rmse <- sqrt(mean((xgb_preds - y_test)^2)); #rmse
r2<-cor(xgb_preds, y_test)^2; r2

xgb_preds_tr<-predict(xgb_model, as.matrix(X_train))
rmse.tr <- sqrt(mean((xgb_preds_tr - y_train)^2)); #rmse.tr
r2.tr<-cor(xgb_preds_tr, y_train)^2; r2.tr

r2rec<-c(r2.tr, r2)


#####
#Plots####

viz<-shapviz(xgb_model, X_pred=as.matrix(X_train))
#imp<-sv_importance(viz, kind="beeswarm"); imp

imp2<-sv_importance(viz, kind="bar", show_numbers = TRUE); imp2+theme_minimal()

sv_dependence(viz, "TA", "doy")+
  ylab("ER SHAP value")+
  theme_minimal()+
  theme(axis.title = element_text(size=24), 
        axis.text=element_text(size=22), 
        legend.position="none")

dev.copy(png, "D:/Analysis/traitflux/plots/ta_er.png", height=600, width=400)
dev.off()

sv_dependence(viz, "TC", "doy")+
  theme_minimal()+
  theme(axis.title = element_text(size=24), 
        axis.text=element_text(size=22), 
        legend.position="none")

dev.copy(png, "D:/Analysis/traitflux/plots/tc_er.png", height=600, width=400)
dev.off()

#imp2<-sv_waterfall(viz,row_id = 3); imp2

# #sv_importance(viz)
  sv_dependence(viz, "nsc", "TC", alpha=0.5)
#  sv_dependence(viz, "TC", "doy")
# sv_dependence(viz, "lma", "TA")
#  sv_dependence(viz, "lai", "doy")

 
 #Fancy dependences (Slow!)

  if(setlab=="trait"){viz$X<-merge(viz$X, data.lab, by=c("nsc", "evi"), sort=FALSE)}

  if(setlab=="comp"){viz$X<-merge(viz$X, data.lab, by="evi", sort=FALSE)}
 
 viz$X$sitecode<-as.factor(viz$X$sitename)
 mycol<-glasbey(n=length(levels(viz$X$sitecode)))
 names(mycol)<-levels(viz$X$sitecode)
 
 viz$X$pftcode<-as.factor(viz$X$igbp)
 igbpcol<-glasbey(n=length(levels(viz$X$pftcode)))
 names(igbpcol)<-levels(viz$X$pftcode)
 
 samp<-sample(1:nrow(viz$X), 3000)
 
 ##trait plots####
 
 er.lma<-sv_dependence(viz[samp],"lma")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
   scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+ylim(-1.5, 2)+
   scale_color_manual(name="site", values=mycol)+
   ylab("ER SHAP value")+
   theme_minimal()+
   theme(axis.title = element_text(size=20), 
         axis.text=element_text(size=18),
         legend.title=element_text(size=20),
         legend.text=element_text(size=18))
 
 
 
 er.nitr<-sv_dependence(viz[samp],"nitr")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
   scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+ylim(-1.5, 2)+
   scale_color_manual(values=mycol, name="site")+xlab("nitrogen")+
   ylab("ER SHAP value")+
   theme_minimal()+
   theme(axis.title = element_text(size=20), 
         axis.text=element_text(size=18),
         legend.position = "none")
 
 
 ggarrange(er.lma, er.nitr)
 
 dev.copy(png, "D:/Analysis/traitflux/plots/er_les.png", height=550, width=500)
 dev.off()

 
 ylim<-c(-2.5, 3)#c(-1,1) #c(-2.5, 3)
 

 er.nsc<-sv_dependence(viz[samp],"nsc")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
   scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+ylim(ylim)+
   scale_color_manual(values=mycol, name="site")+
   ylab("ER SHAP value")+
   theme_minimal()+
   theme(axis.title = element_text(size=20), 
         axis.text=element_text(size=18),
         legend.title=element_text(size=20),
         legend.text=element_text(size=18))
 
 er.pheno<-sv_dependence(viz[samp],"pheno")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
   scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+ylim(ylim)+
   scale_color_manual(values=mycol, name="site")+
   ylab("ER SHAP value")+xlab("phenolics")+
   theme_minimal()+
   theme(axis.title = element_text(size=20), 
         axis.text=element_text(size=18), 
         legend.position="none")
 
 er.lignin<-sv_dependence(viz[samp],"lignin")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
   scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+ylim(ylim)+
   scale_color_manual(values=mycol, name="site")+
   ylab("ER SHAP value")+
   theme_minimal()+
   theme(axis.title = element_text(size=20), 
         axis.text=element_text(size=18),
         legend.position="none")
 
 ggarrange(er.lignin, er.nsc, er.pheno)
 
 
 dev.copy(png, "D:/Analysis/traitflux/plots/er_alloc.png", height=600, width=450)
 dev.off()

 
 sv_dependence(viz[samp],"lignin")+aes(shape=viz$X$igbp[samp], color=viz$X$igbp[samp])+
   scale_shape_manual(values=1:nlevels(viz$X$igbp[samp]), name="site")+
   scale_color_manual(values=igbpcol, name="site")+
   ylab("Lignin - ER SHAP")+
   theme_minimal()+
   theme(axis.title = element_text(size=24), 
        axis.text=element_text(size=22))
 
 er.lai<-sv_dependence(viz[samp],"lai")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
   scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+
   scale_color_manual(values=mycol, name="site")+
   ylab("LAI - ER SHAP")+
   theme_minimal()+
   theme(axis.title = element_text(size=24), 
         axis.text=element_text(size=22))
 
 er.lai.pft<-sv_dependence(viz[samp],"lai")+aes(shape=viz$X$pftcode[samp], color=viz$X$pftcode[samp])+
   scale_shape_manual(values=1:nlevels(viz$X$pftcode[samp]), name="igbp")+
   scale_color_manual(values=igbpcol, name="igbp")+
   ylab("LAI- ER SHAP")+
   theme_minimal()+
   theme(axis.title = element_text(size=24), 
         axis.text=element_text(size=22))
 
 grid.arrange(er.lai, er.lai.pft)
 
 dev.copy(png, "D:/Analysis/traitflux/plots/er_lai.png", height=550, width=500)
 dev.off()
 
 # sv_dependence(viz[samp],"doy")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
 #   scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+
 #   scale_color_manual(values=mycol, name="site")+
 #   ylab("LAI - ER SHAP")+
 #   theme_minimal()
 # 
 # sv_dependence(viz[samp],"doy")+aes(shape=viz$X$pftcode[samp], color=viz$X$pftcode[samp])+
 #   scale_shape_manual(values=1:nlevels(viz$X$pftcode[samp]), name="igbp")+
 #   scale_color_manual(values=igbpcol, name="igbp")+
 #   ylab("LAI- ER SHAP")+
 #   theme_minimal()
 
 #####
 
 #comp plots####
 
 sv_dependence(viz[samp],"Pine")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
   scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+
   scale_color_manual(values=mycol, name="site")+
   theme_minimal()
 
 sv_dependence(viz[samp],"Aspen.Paper.Birch")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
   scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+
   scale_color_manual(values=mycol, name="site")+
   theme_minimal()
 
 sv_dependence(viz[samp],"Broad.leaved.Deciduous.Scrub.Shrub")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
   scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+
   scale_color_manual(values=mycol, name="site")+
   theme_minimal()
 
 sv_dependence(viz[samp],"Other")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
   scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+
   scale_color_manual(values=mycol, name="site")+
   theme_minimal()
 
 sv_dependence(viz[samp],"Fir.Spruce")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
   scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+
   scale_color_manual(values=mycol, name="site")+ylim(-2,2)+
   theme_minimal()
 
 sv_dependence(viz[samp],"Coniferous.Forested.Wetland")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
   scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+
   scale_color_manual(values=mycol, name="site")+ylim(-2,2)+
   theme_minimal()
 #####
 


#####
 rm("dat.in", "xgb_model", "viz")

 
 
 
 ##Water####

print("starting LE")

circ<-c("FC", "E_0", "resid", "GPP_","H", "NEE", "SC_", "Reco_uStar", "CO2", "LE_U50_orig", "LE_uStar", "Reco_U50", "SLE_1_1_1", "nirvp","LW_IN_1_1_1", "LW_OUT_1_1_1", "G_1_1_1","nirv", "ndvi") #, "nirv", "ndvi"
circ.ind<-grep(paste(circ,collapse="|"), colnames(data.unique))
dat.in<-data.unique[,-circ.ind]
colnames(dat.in)


dat.le<-kill.nan(dat.in, var="LE")
dat.lab.le<-data.lab; dat.lab.le<-dat.lab.le[!is.na(dat.in$LE),]


#withhold full days for testing
if(split.type=="doy"){
  doy.train<-sample(unique(dat.le$doy), round(0.7*length(unique(dat.le$doy)))) #use 70% of days (i.e. withhold whole days) to make it harder for the model
  train.ind<-which(dat.le$doy%in%doy.train) # sample(1:nrow(dat.er), 0.7*nrow(dat.er))
}

#withhold full towers for testing
if(split.type=="tower"){
  twr.train<-sample(unique(dat.lab.le$sitename), round(0.8*length(unique(dat.lab.le$sitename))))
  train.ind<-which(dat.lab.le$sitename%in%twr.train)
}


if(split.type=="tower"){print(unique(dat.lab.le$sitename[-train.ind]))}


#train.ind<-sample(1:nrow(dat.le), 0.7*nrow(dat.le))
dat.train<-dat.le[train.ind,]
dat.test<-dat.le[-train.ind,]


#xgboost
y_train <- dat.train$LE
y_test <- dat.test$LE
X_train <- dat.train %>% dplyr::select(-LE)
X_test <- dat.test %>% dplyr::select(-LE)

xgb_train <- xgb.DMatrix(data = as.matrix(X_train), label = y_train)
xgb_test <- xgb.DMatrix(data = as.matrix(X_test), label = y_test)

xgb_params <- list(
  eta = 0.05,
  max_depth = 8,
  gamma = 10,
  subsample = 0.6,
  colsample_bytree = 0.6,
  min_child_weight=10
)


xgb_model <- xgb.train(
  data = xgb_train,
  params=xgb_params,
  nrounds = 500,
  verbose = 1
)



xgb_preds <- predict(xgb_model, as.matrix(X_test))
rmse <- sqrt(mean((xgb_preds- y_test)^2)); rmse
r2<-cor(xgb_preds, y_test)^2; r2

xgb_preds_tr<-predict(xgb_model, as.matrix(X_train))
r2.tr<-cor(xgb_preds_tr, y_train)^2; r2.tr


r2rec<-c(r2.tr, r2)



#plots####
#if(plots==TRUE){

viz<-shapviz(xgb_model, X_pred=as.matrix(X_train))
#imp<-sv_importance(viz, kind="beeswarm"); imp

#imp1<-sv_waterfall(viz); imp1
imp2<-sv_importance(viz,show_numbers = TRUE); imp2


sv_dependence(viz, "VPD_C", "doy")+
  ylab("LE SHAP value")+
  theme_minimal()+
  theme(axis.title = element_text(size=24), 
      axis.text=element_text(size=22), 
      legend.title=element_text(size=24),
      legend.text=element_text(size=22))
dev.copy(png, "D:/Analysis/traitflux/plots/vpd_le.png", height=600, width=500)
dev.off()

# sv_dependence(viz, "TC")
# sv_dependence(viz, "doy")
# sv_dependence(viz, "SW_IN", "doy")
 sv_dependence(viz,  "VPD_C","doy", alpha=0.5)+ylab("SHAP value for LE")+theme_minimal()
 
 #Fancy dependences (Slow!)

 if(setlab=="trait"){viz$X<-merge(viz$X, data.lab, by=c("nsc", "evi"), sort=FALSE)}

 if(setlab=="comp"){viz$X<-merge(viz$X, data.lab, by="evi", sort=FALSE)}

 viz$X$sitecode<-as.factor(viz$X$sitename)
 mycol<-glasbey(n=length(levels(viz$X$sitecode)))
 names(mycol)<-levels(viz$X$sitecode)
 
 viz$X$pftcode<-as.factor(viz$X$igbp)
 igbpcol<-glasbey(n=length(levels(viz$X$pftcode)))
 names(igbpcol)<-levels(viz$X$pftcode)
 
 
 samp<-sample(1:nrow(viz$X), 3000)
 
 #trait plots####
 
 le.lma<-sv_dependence(viz[samp],"lma")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
   scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+ylim(-15,50)+
   scale_color_manual(values=mycol, name="site")+ylab("LE SHAP value")+
   theme_minimal()+
   theme(axis.title = element_text(size=20), 
         axis.text=element_text(size=18),
         legend.position="none")
 
 le.nitr<-sv_dependence(viz[samp],"nitr")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
   scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+ylim(-15,50)+
   scale_color_manual(values=mycol, name="site")+ylab("LE SHAP value")+xlab("nitrogen")+
   theme_minimal()+
   theme(axis.title = element_text(size=20), 
         axis.text=element_text(size=18),
         legend.position="none")
 
 ggarrange(le.lma, le.nitr)
 
 
 dev.copy(png, "D:/Analysis/traitflux/plots/le_les.png", height=550, width=400)
 dev.off()
 
 le.nsc<-sv_dependence(viz[samp],"nsc")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
   scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+ylim(-30, 40)+
   scale_color_manual(values=mycol, name="site")+ylab("LE SHAP value")+
   theme_minimal()+
   theme(axis.title = element_text(size=20), 
         axis.text=element_text(size=18),
         legend.position="none")
 
 le.pheno<-sv_dependence(viz[samp],"pheno")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
   scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+ylim(-30, 40)+
   scale_color_manual(values=mycol, name="site")+xlab("phenolics")+ylab("LE SHAP value")+
   theme_minimal()+
   theme(axis.title = element_text(size=20), 
         axis.text=element_text(size=18),
         legend.position="none")
 
 le.lig<-sv_dependence(viz[samp],"lignin")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
   scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+ylim(-30, 40)+
   scale_color_manual(values=mycol, name="site")+ylab("LE SHAP value")+
   theme_minimal()+
   theme(axis.title = element_text(size=20), 
         axis.text=element_text(size=18),
         legend.position="none")
 
 ggarrange(le.lig, le.nsc, le.pheno)
 
 
 dev.copy(png, "D:/Analysis/traitflux/plots/le_alloc.png", height=600, width=400)
 dev.off()
 
 
 
 le.lai<-sv_dependence(viz[samp],"lai")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
   scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+
   scale_color_manual(values=mycol, name="site")+
   theme_minimal()
 
 le.lai.pft<-sv_dependence(viz[samp],"lai")+aes(shape=viz$X$pftcode[samp], color=viz$X$pftcode[samp])+
   scale_shape_manual(values=1:nlevels(viz$X$pftcode[samp]), name="pft")+
   scale_color_manual(values=igbpcol, name="pft")+
   theme_minimal()
 

 grid.arrange(le.lai, le.lai.pft)
 
 
 dev.copy(png, "D:/Analysis/traitflux/plots/le_lai.png", height=550, width=500)
 dev.off()
 
#####

 #comp plots####
 
 sv_dependence(viz[samp],"Pine")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
   scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+
   scale_color_manual(values=mycol, name="site")+
   theme_minimal()
 
 sv_dependence(viz[samp],"Aspen.Paper.Birch")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
   scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+
   scale_color_manual(values=mycol, name="site")+
   theme_minimal()
 
 sv_dependence(viz[samp],"Broad.leaved.Deciduous.Scrub.Shrub")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
   scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+
   scale_color_manual(values=mycol, name="site")+
   theme_minimal()
 
 sv_dependence(viz[samp],"Fir.Spruce")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
   scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+
   scale_color_manual(values=mycol, name="site")+
   theme_minimal()
 
 sv_dependence(viz[samp],"Other")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
   scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+
   scale_color_manual(values=mycol, name="site")+
   theme_minimal()
 
 sv_dependence(viz[samp],"Coniferous.Forested.Wetland")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
   scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+
   scale_color_manual(values=mycol, name="site")+
   theme_minimal()
 #####

rm("dat.in", "xgb_model", "viz")
#####

 
 
 

##Sensible heat####
print("starting H")

circ<-c("FC", "E_0", "resid", "GPP_","LE", "NEE", "SC_", "Reco_uStar", "CO2", "LE_U50_orig", "LE_uStar", "Reco_U50", "SLE_1_1_1", "SH_1_1_1", "nirvp","LW_IN_1_1_1", "LW_OUT_1_1_1", "G_1_1_1", "nirv")
circ.ind<-grep(paste(circ,collapse="|"), colnames(data.unique))
dat.in<-data.unique[,-circ.ind]
colnames(dat.in)


dat.h<-kill.nan(dat.in, var="H")
dat.lab.h<-data.lab; dat.lab.h<-dat.lab.h[!is.na(dat.in$H),]



#withhold full days for testing
if(split.type=="doy"){
  doy.train<-sample(unique(dat.h$doy), round(0.7*length(unique(dat.h$doy)))) #use 70% of days (i.e. withhold whole days) to make it harder for the model
  train.ind<-which(dat.h$doy%in%doy.train) # sample(1:nrow(dat.er), 0.7*nrow(dat.er))
}

#withhold full towers for testing
if(split.type=="tower"){
  twr.train<-sample(unique(dat.lab.h$sitename), round(0.8*length(unique(dat.lab.h$sitename))))
  train.ind<-which(dat.lab.h$sitename%in%twr.train)
}


if(split.type=="tower"){print(unique(dat.lab.h$sitename[-train.ind]))}


#train.ind<-sample(1:nrow(dat.h), 0.7*nrow(dat.h))
dat.train<-dat.h[train.ind,]
dat.test<-dat.h[-train.ind,]


#xgboost
y_train <- dat.train$H
y_test <- dat.test$H
X_train <- dat.train %>% dplyr::select(-H)
X_test <- dat.test %>% dplyr::select(-H)

xgb_train <- xgb.DMatrix(data = as.matrix(X_train), label = y_train)
xgb_test <- xgb.DMatrix(data = as.matrix(X_test), label = y_test)




xgb_params <- list(
  eta = 0.05,
  max_depth = 8,
  subsample = 0.6,
  colsample_bytree = 0.6,
  min_child_weight = 1,
  gamma=0

)

xgb_model <- xgb.train( #params = xgb_params,
  data = xgb_train,
  params=xgb_params,
  nrounds = 500,
  verbose = 1
)



xgb_preds <- predict(xgb_model, as.matrix(X_test))
rmse <- sqrt(mean((xgb_preds- y_test)^2)); rmse
r2<-cor(xgb_preds, y_test)^2; r2

xgb_preds_tr<-predict(xgb_model, as.matrix(X_train))
r2.tr<-cor(xgb_preds_tr, y_train)^2; r2.tr

r2rec<-c(r2.tr, r2)

#plots####
#if(plots==TRUE){

viz<-shapviz(xgb_model, X_pred=as.matrix(X_train))
#imp<-sv_importance(viz, kind="beeswarm"); imp

#imp1<-sv_waterfall(viz); imp1
imp2<-sv_importance(viz,show_numbers = TRUE); imp2



if(setlab=="trait"){viz$X<-merge(viz$X, data.lab, by=c("nsc", "evi"), sort=FALSE)}

if(setlab=="comp"){viz$X<-merge(viz$X, data.lab, by="evi", sort=FALSE)}

viz$X$sitecode<-as.factor(viz$X$sitename)
mycol<-glasbey(n=length(levels(viz$X$sitecode)))
names(mycol)<-levels(viz$X$sitecode)

viz$X$pftcode<-as.factor(viz$X$igbp)
igbpcol<-glasbey(n=length(levels(viz$X$pftcode)))
names(igbpcol)<-levels(viz$X$pftcode)


samp<-sample(1:nrow(viz$X), 3000)


#trait plots####

h.lma<-sv_dependence(viz[samp],"lma")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
  scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+ylim(-10,20)+
  scale_color_manual(values=mycol, name="site")+ylab("H SHAP value")+
  theme_minimal()+
  theme(axis.title = element_text(size=20), 
        axis.text=element_text(size=18),
        legend.title=element_text(size=20),
        legend.text=element_text(size=18))


h.nitr<-sv_dependence(viz[samp],"nitr")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
  scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+ylim(-10,20)+
  scale_color_manual(values=mycol, name="site")+xlab("nitrogen")+ylab("H SHAP value")+
  theme_minimal()+
  theme(axis.title = element_text(size=20), 
        axis.text=element_text(size=18),
        legend.position="none")

ggarrange(h.lma, h.nitr)

dev.copy(png, "D:/Analysis/traitflux/plots/h_les.png", height=550, width=500)
dev.off()

# sv_dependence(viz[samp],"nitr")+aes(shape=viz$X$pftcode[samp], color=viz$X$pftcode[samp])+
#   scale_shape_manual(values=1:nlevels(viz$X$pftcode[samp]), name="pft")+
#   scale_color_manual(values=igbpcol, name="pft")+
#   theme_minimal()

h.nsc<-sv_dependence(viz[samp],"nsc")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
  scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+ylab("H SHAP value")+
  scale_color_manual(values=mycol, name="site")+ylim(-10,20)+
  theme_minimal()+
  theme(axis.title = element_text(size=20), 
        axis.text=element_text(size=18),
        legend.title=element_text(size=20),
        legend.text=element_text(size=18))


h.pheno<-sv_dependence(viz[samp],"pheno")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
  scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+ylab("H SHAP value")+xlab("phenolics")+
  scale_color_manual(values=mycol, name="site")+ylim(-10,20)+
  theme_minimal()+
  theme(axis.title = element_text(size=20), 
        axis.text=element_text(size=18),
        legend.position="none")

h.lig<-sv_dependence(viz[samp],"lignin")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
  scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+ylab("H SHAP value")+
  scale_color_manual(values=mycol, name="site")+ylim(-10,20)+
  theme_minimal()+
  theme(axis.title = element_text(size=20), 
        axis.text=element_text(size=18),
        legend.position="none")
  
ggarrange(h.lig, h.nsc, h.pheno)

dev.copy(png, "D:/Analysis/traitflux/plots/h_alloc.png", height=600, width=500)
dev.off()


h.lai<-sv_dependence(viz[samp],"lai")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
  scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+ylim(-15,15)+
  scale_color_manual(values=mycol, name="site")+
  theme_minimal()

h.lai.pft<-sv_dependence(viz[samp],"lai")+aes(shape=viz$X$pftcode[samp], color=viz$X$pftcode[samp])+
  scale_shape_manual(values=1:nlevels(viz$X$pftcode[samp]), name="site")+ylim(-15,15)+
  scale_color_manual(values=igbpcol, name="site")+
  theme_minimal()


grid.arrange(h.lai, h.lai.pft)

dev.copy(png, "D:/Analysis/traitflux/plots/h_lai.png", height=550, width=500)
dev.off()



sv_dependence(viz[samp],"nirv")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
  scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+
  scale_color_manual(values=mycol, name="site")+
  theme_minimal()

sv_dependence(viz[samp],"nirv")+aes(shape=viz$X$pftcode[samp], color=viz$X$pftcode[samp])+
  scale_shape_manual(values=1:nlevels(viz$X$pftcode[samp]), name="site")+
  scale_color_manual(values=igbpcol, name="site")+
  theme_minimal()

#####
#comp plots####

sv_dependence(viz[samp],"Pine")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
  scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+
  scale_color_manual(values=mycol, name="site")+
  theme_minimal()

sv_dependence(viz[samp],"Aspen.Paper.Birch")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
  scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+
  scale_color_manual(values=mycol, name="site")+
  theme_minimal()

sv_dependence(viz[samp],"Broad.leaved.Deciduous.Scrub.Shrub")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
  scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+
  scale_color_manual(values=mycol, name="site")+
  theme_minimal()

sv_dependence(viz[samp],"Fir.Spruce")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
  scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+
  scale_color_manual(values=mycol, name="site")+
  theme_minimal()

sv_dependence(viz[samp],"Other")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
  scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+
  scale_color_manual(values=mycol, name="site")+
  theme_minimal()

sv_dependence(viz[samp],"Coniferous.Forested.Wetland")+aes(shape=viz$X$sitecode[samp], color=viz$X$sitecode[samp])+
  scale_shape_manual(values=1:nlevels(viz$X$sitecode[samp]), name="site")+
  scale_color_manual(values=mycol, name="site")+
  theme_minimal()
#####

# }

#####

rm(dat.in)

#plot each time a new model set added

ggplot(results) +
  aes(x = flux, y = test, fill = set) +
  geom_boxplot() +
  scale_fill_hue(direction = 1) +
  theme_minimal()

#}

#closes model loop

ggplot(results) +
  aes(x = flux, y = test, fill = set) +
  geom_boxplot() +
  scale_fill_hue(direction = 1) +
  ylab("testing R2")+
  ylim(0.2, 1)+
  theme_minimal()

results.clear<-results[results$set%in%c("both", "trait", "evi", "comp"),]

ggplot(results.clear) +
  aes(x = flux, y = test, fill = set) +
  geom_boxplot() +
  scale_fill_hue(direction = 1) +
  ylab("testing R2")+
  ylim(0.2, 1)+
  theme_minimal()

library(ggplot2)

ggplot(master) +
  aes(x = lignin, y = pheno, colour = sitename) +
  geom_point() +
  scale_color_hue(direction = 1) +
  theme_minimal() +
  facet_wrap(vars(igbp))

ggplot(master) +
  aes(x = nitr, y = lma, colour = sitename) +
  geom_point() +
  scale_color_hue(direction = 1) +
  theme_minimal() +
  facet_wrap(vars(igbp))

