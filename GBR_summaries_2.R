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
plots<-TRUE
incl.nobio<-FALSE

#plot holders
er.bees<-gpp.bees<-le.bees<-h.bees<-list()
er.imp<-gpp.imp<-le.imp<-h.imp<-list()

#parameter list; get set by individual model now
# xgb_params <- list(
#   eta = 0.07,
#   max_depth = 5,
#   subsample = 0.6,
#   colsample_bytree = 1,
#   min_child_weight = 10, 
#   gamma=0
#   
# )

traituse<-c(TRUE, FALSE); compuse<-c(TRUE, FALSE); eviuse<-c(TRUE, FALSE)
sets<-expand.grid(traituse, compuse, eviuse); colnames(sets)<-c("trait", "comp", "vi")
if(incl.nobio==FALSE){sets<-expand.grid(traituse, compuse);colnames(sets)<-c("trait", "comp")}

data.unique.bk<-data.unique

#####

for(s in 1:nrow(sets)){

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

#vars to leave out to avoid moderate correlations:
# ,"lma", "pheno", "ndvi"

#GPP#####

print("starting GPP")

circ<-c("FC", "E_0", "resid", "GPP_uStar_f", "LE","H", "NEE", "SC_", "Reco","CO2", "nirvp", "G_1_1_1", "LW_IN_1_1_1", "LW_OUT_1_1_1", "nirv")
circ.ind<-grep(paste(circ,collapse="|"), colnames(data.unique))
dat.in<-data.unique[,-circ.ind]
colnames(dat.in)


dat.gpp<-kill.nan(dat=dat.in, var="GPP_BB"); 
dat.lab.gpp<-data.lab; dat.lab.gpp<-dat.lab.gpp[!is.na(dat.in$GPP_BB),]


r2dat.gpp<-data.frame(matrix(nrow=n.iter, ncol=2)); colnames(r2dat.gpp)<-c("train", "test")
n.iter<-n.iter

for(i in 1:n.iter){
  
  print(i)
  
#gpp.model<-function(dat.gpp){

#withhold full days for testing
doy.train<-sample(unique(dat.gpp$doy), round(0.7*length(unique(dat.gpp$doy)))) #use 70% of days (i.e. withhold whole days) to make it harder for the model
train.ind<-which(dat.gpp$doy%in%doy.train) # sample(1:nrow(dat.gpp), 0.7*nrow(dat.gpp))

y_all <- dat.gpp$GPP_BB
X_all <- dat.gpp %>% dplyr::select(-GPP_BB)

xgb_all<-xgb.DMatrix(data = as.matrix(X_all), label = y_all)


dat.train<-dat.gpp[train.ind,] #training data
dat.test<-dat.gpp[-train.ind,] #testing data

#xgboost
y_train <- dat.train$GPP_BB
y_test <- dat.test$GPP_BB
X_train <- dat.train %>% select(-GPP_BB)
X_test <- dat.test %>% select(-GPP_BB)

xgb_train <- xgb.DMatrix(data = as.matrix(X_train), label = y_train)
xgb_test <- xgb.DMatrix(data = as.matrix(X_test), label = y_test)

#parameter list
xgb_params <- list(
  eta = 0.05,
  max_depth = 8,
  subsample = 0.6,
  colsample_bytree = 0.6,
  min_child_weight = 15, 
  gamma=10
  
)


######


xgb_model <- xgb.train( #params = xgb_params,
  data = xgb_train,
  params=xgb_params,
  verbosity = 1,
  nrounds=500
  
)

xgb_preds <- predict(xgb_model, as.matrix(X_test))
rmse <- sqrt(mean((xgb_preds - y_test)^2)); #rmse
r2<-cor(xgb_preds, y_test)^2; r2

xgb_preds_tr<-predict(xgb_model, as.matrix(X_train))
rmse.tr <- sqrt(mean((xgb_preds_tr - y_train)^2)); #rmse.tr
r2.tr<-cor(xgb_preds_tr, y_train)^2; r2.tr

r2rec<-c(r2.tr, r2)


r2rec<-c(r2.tr, r2)
r2dat.gpp[i,]<-r2rec

##save plots
runlab<-paste("gpp", setlab)

if(i==1){  
  print("saving initial plots")
  viz<-shapviz(xgb_model, X_pred=as.matrix(X_train))
  imp<-sv_importance(viz, kind="beeswarm")
  gpp.bees[[s]]<-imp+xlab(paste("SHAP values for", runlab)); names(gpp.bees)[[s]]<-runlab
  
  
  impvar<-names(sort(colMeans(abs(viz$S), na.rm=TRUE), decreasing=TRUE))[1:9] #12 bars used by default importance plot
  sublab<-var.lab[names(var.lab)%in%impvar]; sublab<-sublab[impvar]
  imp2<-sv_importance(viz, kind="bar", show_numbers = TRUE, fill=sublab, max_display=9L); gpp.imp[[s]]<-imp2
  gpp.imp[[s]]<-imp2+xlab(paste("mean |SHAP values| for", runlab)); names(gpp.imp)[[s]]<-runlab
  
}

if(i>1 & r2>=max(r2dat.gpp$test, na.rm=TRUE)){
  
  print("best r2 so far, replacing plots...")
  
  viz<-shapviz(xgb_model, X_pred=as.matrix(X_train))
  imp<-sv_importance(viz, kind="beeswarm")
  gpp.bees[[s]]<-imp+xlab(paste("SHAP values for", runlab)); names(gpp.bees)[[s]]<-runlab
  
  
  impvar<-names(sort(colMeans(abs(viz$S), na.rm=TRUE), decreasing=TRUE))[1:9] #12 bars used by default importance plot
  sublab<-var.lab[names(var.lab)%in%impvar]; sublab<-sublab[impvar]
  imp2<-sv_importance(viz, kind="bar", show_numbers = TRUE, fill=sublab, max_display=9L); gpp.imp[[s]]<-imp2
  gpp.imp[[s]]<-imp2+xlab(paste("mean |SHAP values| for", runlab)); names(gpp.imp)[[s]]<-runlab
  
}


}


r2dat.gpp$diff<-r2dat.gpp$train-r2dat.gpp$test
r2dat.gpp$flux<-"GPP"; r2dat.gpp$set<-setlab
r2dat.gpp

results<-rbind(results, r2dat.gpp)

#Plots####
if(plots==TRUE){


viz<-shapviz(xgb_model, X_pred=as.matrix(X_train))
#imp<-sv_importance(viz, kind="beeswarm"); imp
imp2<-sv_importance(viz, kind="bar", show_numbers = TRUE); imp2
#imp2<-sv_waterfall(viz,row_id = 3); imp2

}

rm("dat.in")

#####


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
#Main loop#####

r2dat.er<-data.frame(matrix(nrow=n.iter, ncol=2)); colnames(r2dat.er)<-c("train", "test")
n.iter<-n.iter


for(i in 1:n.iter){

print(i)


#withhold full days for testing
doy.train<-sample(unique(dat.er$doy), round(0.7*length(unique(dat.er$doy)))) #use 70% of days (i.e. withhold whole days) to make it harder for the model
train.ind<-which(dat.er$doy%in%doy.train) # sample(1:nrow(dat.er), 0.7*nrow(dat.er))

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

##Grid search#### 


# param.grid<-expand.grid(nrounds = c(50, 100, 500),
#                         max_depth = c(3, 5, 8),
#                         eta = c(0.05, 0.07),
#                         gamma = c(0,10),
#                         subsample = c(0.6,1),
#                         min_child_weight = c(1,5,10),
#                         colsample_bytree = c(0.6,1))
# 
# 
# param.grid.rate<-expand.grid(nrounds = 2000,
#                              max_depth = 3,
#                              eta = c(0.01, 0.03, 0.05, 0.07, 0.1),
#                              gamma = 0,
#                              subsample = 1,
#                              min_child_weight = 1,
#                              colsample_bytree = 0.6)
# 
# 
# 
# option 1###
# 
# output<-data.frame(matrix(ncol=5,nrow=nrow(param.grid))); colnames(output)<-c("train", "test", "eta", "depth", "nrounds")
#
# for(i in 1:nrow(param.grid.rate)){
#
#   print(i)
#
#   xgb_params <- list(
#   #booster = "gbtree",
#   eta = param.grid.rate$eta[i],
#   max_depth = param.grid.rate$max_depth[i]
#   #gamma = 4,
#   #subsample = 0.75,
#   #colsample_bytree = 1,
#   #objective = "reg:squarederror"
# )
#
# xgb_model <- xgb.cv(
#   data = xgb_all,
#   params=xgb_params,
#   nrounds = param.grid.rate$nrounds[i],
#   verbose = 1,
#   nfold=3,
#   #prediction = TRUE,
#   early_stopping_rounds = 10
# )
#
# xvalidationScores <- as.data.frame(xgb_model$evaluation_log)
#
# rmse <- tail(xvalidationScores$test_rmse_mean, 1)
# trmse <- tail(xvalidationScores$train_rmse_mean,1)
#
# output[i,] <- c(rmse, trmse, param.grid$eta[i], param.grid$max_depth[i], param.grid$nrounds[i])
#
# plot(xvalidationScores$train_rmse_mean, ylab="RMSE", xlab="iteration")
# points(xvalidationScores$test_rmse_mean, col='red')
# }
#


#Option 2##


#train_control = trainControl(method = "cv", number = 3, search = "grid")

# # training a XGboost Regression tree model while tuning parameters
# xgb2 = train(Reco_U50~., data = dat.train, metric="Rsquared", method = "xgbTree", trControl = train_control, tuneGrid = param.grid, na.action = na.omit)
# 
# # summarising the results
# print(xgb2)
# 
# ####

#parameter list
xgb_params <- list(
  eta = 0.03, #0.07
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
r2dat.er[i,]<-r2rec

##save plots
runlab<-paste("ER", setlab)

if(i==1){  
print("saving initial plots")
viz<-shapviz(xgb_model, X_pred=as.matrix(X_train))
imp<-sv_importance(viz, kind="beeswarm")
er.bees[[s]]<-imp+xlab(paste("SHAP values for", runlab)); names(er.bees)[[s]]<-runlab


impvar<-names(sort(colMeans(abs(viz$S), na.rm=TRUE), decreasing=TRUE))[1:9] #12 bars used by default importance plot
sublab<-var.lab[names(var.lab)%in%impvar]; sublab<-sublab[impvar]
imp2<-sv_importance(viz, kind="bar", show_numbers = TRUE, fill=sublab, max_display=9L); er.imp[[s]]<-imp2
er.imp[[s]]<-imp2+xlab(paste("mean |SHAP values| for", runlab)); names(er.imp)[[s]]<-runlab

}

if(i>1 & r2>=max(r2dat.er$test, na.rm=TRUE)){
  
  print("best r2 so far, replacing plots...")
  
  viz<-shapviz(xgb_model, X_pred=as.matrix(X_train))
  imp<-sv_importance(viz, kind="beeswarm")
  er.bees[[s]]<-imp+xlab(paste("SHAP values for", runlab)); names(er.bees)[[s]]<-runlab
  
  
  impvar<-names(sort(colMeans(abs(viz$S), na.rm=TRUE), decreasing=TRUE))[1:9] #12 bars used by default importance plot
  sublab<-var.lab[names(var.lab)%in%impvar]; sublab<-sublab[impvar]
  imp2<-sv_importance(viz, kind="bar", show_numbers = TRUE, fill=sublab, max_display=9L); er.imp[[s]]<-imp2
  er.imp[[s]]<-imp2+xlab(paste("mean |SHAP values| for", runlab)); names(er.imp)[[s]]<-runlab
  
}



}


r2dat.er$diff<-r2dat.er$train-r2dat.er$test
r2dat.er$flux<-"ER"; r2dat.er$set<-setlab
r2dat.er

results<-rbind(results, r2dat.er)


#####
#Plots####
if(plots==TRUE){

viz<-shapviz(xgb_model, X_pred=as.matrix(X_train))
#imp<-sv_importance(viz, kind="beeswarm"); imp

imp2<-sv_importance(viz, kind="bar", show_numbers = TRUE); 

imp2+
xlab("mean |SHAP value| for ER")+
theme_minimal()


}

#####

rm("dat.in")

# }  #loop cutoff to run for only GPP and ER



##Water####

print("starting LE")

circ<-c("FC", "E_0", "resid", "GPP_","H", "NEE", "SC_", "Reco_uStar", "CO2", "LE_U50_orig", "LE_uStar", "Reco_U50", "SLE_1_1_1", "nirvp","LW_IN_1_1_1", "LW_OUT_1_1_1", "G_1_1_1","nirv")
circ.ind<-grep(paste(circ,collapse="|"), colnames(data.unique))
dat.in<-data.unique[,-circ.ind]
colnames(dat.in)


dat.le<-kill.nan(dat.in, var="LE")

r2dat.le<-data.frame(matrix(nrow=n.iter, ncol=2)); colnames(r2dat.le)<-c("train", "test")
n.iter<-n.iter


for(i in 1:n.iter){

  #le.model<-function(dat.le){

doy.train<-sample(unique(dat.le$doy), round(0.7*length(unique(dat.le$doy)))) #use 70% of days (i.e. withhold whole days) to make it harder for the model
train.ind<-which(dat.le$doy%in%doy.train) # sample(1:nrow(dat.er), 0.7*nrow(dat.er))

#train.ind<-sample(1:nrow(dat.le), 0.7*nrow(dat.le))
dat.train<-dat.le[train.ind,]
dat.test<-dat.le[-train.ind,]


#xgboost
y_train <- dat.train$LE
y_test <- dat.test$LE
X_train <- dat.train %>% select(-LE)
X_test <- dat.test %>% select(-LE)

xgb_train <- xgb.DMatrix(data = as.matrix(X_train), label = y_train)
xgb_test <- xgb.DMatrix(data = as.matrix(X_test), label = y_test)

# ##Grid search for params####
# 
# param.grid<-expand.grid(nrounds = c(50, 100, 500),
#                         max_depth = c(3, 5, 8),
#                         eta = c(0.05, 0.07),
#                         gamma = c(0,10),
#                         subsample = c(0.6,1),
#                         min_child_weight = c(1,5,10),
#                         colsample_bytree = c(0.6,1))
# 
# train_control = trainControl(method = "cv", number = 3, search = "grid")
# 
# # training a XGboost Regression tree model while tuning parameters
#  xgb2 = train(LE~., data = dat.train, metric="Rsquared", method = "xgbTree", trControl = train_control, tuneGrid = param.grid, na.action = na.omit)
# 
# #summarising the results
# print(xgb2)
# gridresult<-xgb2$results
# #####

xgb_params <- list(
  eta = 0.05, #0.05
  max_depth = 8,
  gamma = 10,
  subsample = 0.6,
  colsample_bytree = 0.6,
  min_child_weight=10
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
r2dat.le[i,]<-r2rec


##save plots
runlab<-paste("LE", setlab)

if(i==1){  
  print("saving initial plots")
  viz<-shapviz(xgb_model, X_pred=as.matrix(X_train))
  imp<-sv_importance(viz, kind="beeswarm")
  le.bees[[s]]<-imp+xlab(paste("SHAP values for", runlab)); names(le.bees)[[s]]<-runlab
  
  
  
  
  impvar<-names(sort(colMeans(abs(viz$S), na.rm=TRUE), decreasing=TRUE))[1:9] #12 bars used by default importance plot
  sublab<-var.lab[names(var.lab)%in%impvar]; sublab<-sublab[impvar]
  imp2<-sv_importance(viz, kind="bar", show_numbers = TRUE, fill=sublab, max_display=9L); le.imp[[s]]<-imp2
  le.imp[[s]]<-imp2+xlab(paste("mean |SHAP values| for", runlab)); names(le.imp)[[s]]<-runlab
  
}

if(i>1 & r2>=max(r2dat.le$test, na.rm=TRUE)){
  
  print("best r2 so far, replacing plots...")
  
  viz<-shapviz(xgb_model, X_pred=as.matrix(X_train))
  imp<-sv_importance(viz, kind="beeswarm")
  le.bees[[s]]<-imp+xlab(paste("SHAP values for", runlab)); names(le.bees)[[s]]<-runlab
  
  
  impvar<-names(sort(colMeans(abs(viz$S), na.rm=TRUE), decreasing=TRUE))[1:9] #12 bars used by default importance plot; 9 is minimum variables used across all models
  sublab<-var.lab[names(var.lab)%in%impvar]; sublab<-sublab[impvar]
  imp2<-sv_importance(viz, kind="bar", show_numbers = TRUE, fill=sublab, max_display=9L); le.imp[[s]]<-imp2
  le.imp[[s]]<-imp2+xlab(paste("mean |SHAP values| for", runlab)); names(le.imp)[[s]]<-runlab
  
}


}


r2dat.le$diff<-r2dat.le$train-r2dat.le$test
r2dat.le$flux<-"LE"; r2dat.le$set<-setlab
r2dat.le

results<-rbind(results, r2dat.le)

#plots####
if(plots==TRUE){

viz<-shapviz(xgb_model, X_pred=as.matrix(X_train))
imp<-sv_importance(viz, kind="beeswarm"); imp

#imp1<-sv_waterfall(viz); imp1
imp2<-sv_importance(viz,show_numbers = TRUE); imp2

}

rm("dat.in")
#####


##Sensible heat####
print("starting H")

circ<-c("FC", "E_0", "resid", "GPP_","LE", "NEE", "SC_", "Reco_uStar", "CO2", "LE_U50_orig", "LE_uStar", "Reco_U50", "SLE_1_1_1", "SH_1_1_1", "nirvp","LW_IN_1_1_1", "LW_OUT_1_1_1", "G_1_1_1", "nirv")
circ.ind<-grep(paste(circ,collapse="|"), colnames(data.unique))
dat.in<-data.unique[,-circ.ind]
colnames(dat.in)


dat.h<-kill.nan(dat.in, var="H")

r2dat.h<-data.frame(matrix(nrow=n.iter, ncol=2)); colnames(r2dat.h)<-c("train", "test")
n.iter<-n.iter

for(i in 1:n.iter){
  

#h.model<-function(dat.h){

doy.train<-sample(unique(dat.h$doy), round(0.7*length(unique(dat.h$doy)))) #use 70% of days (i.e. withhold whole days) to make it harder for the model
train.ind<-which(dat.h$doy%in%doy.train) # sample(1:nrow(dat.er), 0.7*nrow(dat.er))

#train.ind<-sample(1:nrow(dat.h), 0.7*nrow(dat.h))
dat.train<-dat.h[train.ind,]
dat.test<-dat.h[-train.ind,]


#xgboost
y_train <- dat.train$H
y_test <- dat.test$H
X_train <- dat.train %>% select(-H)
X_test <- dat.test %>% select(-H)

xgb_train <- xgb.DMatrix(data = as.matrix(X_train), label = y_train)
xgb_test <- xgb.DMatrix(data = as.matrix(X_test), label = y_test)


# ##Grid search for params####
# 
# param.grid<-expand.grid(nrounds = c(50, 100, 500),
#                         max_depth = c(3, 5, 8),
#                         eta = c(0.05, 0.07),
#                         gamma = c(0,10),
#                         subsample = c(0.6,1),
#                         min_child_weight = c(1,5,10),
#                         colsample_bytree = c(0.6,1))
# 
# train_control = trainControl(method = "cv", number = 3, search = "grid")
# 
# # training a XGboost Regression tree model while tuning parameters
#  xgb2 = train(H~., data = dat.train, metric="Rsquared", method = "xgbTree", trControl = train_control, tuneGrid = param.grid, na.action = na.omit)
# 
# #summarising the results
# print(xgb2)
# gridresult<-xgb2$results
# #####

xgb_params <- list(
  eta = 0.05, #0.05
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
r2dat.h[i,]<-r2rec

##save plots
runlab<-paste("H", setlab)

if(i==1){  
  print("saving initial plots")
  viz<-shapviz(xgb_model, X_pred=as.matrix(X_train))
  imp<-sv_importance(viz, kind="beeswarm")
  h.bees[[s]]<-imp+xlab(paste("SHAP values for", runlab)); names(h.bees)[[s]]<-runlab
  
  #Improtance plot coded by variable type
  impvar<-names(sort(colMeans(abs(viz$S), na.rm=TRUE), decreasing=TRUE))[1:9] #12 bars used by default importance plot
  sublab<-var.lab[names(var.lab)%in%impvar]; sublab<-sublab[impvar]
  imp2<-sv_importance(viz, kind="bar", show_numbers = TRUE,fill=sublab, max_display=9L); h.imp[[s]]<-imp2
  h.imp[[s]]<-imp2+xlab(paste("mean |SHAP values| for", runlab)); names(h.imp)[[s]]<-runlab
  
}

if(i>1 & r2>=max(r2dat.h$test, na.rm=TRUE)){
  
  print("best r2 so far, replacing plots...")
  
  viz<-shapviz(xgb_model, X_pred=as.matrix(X_train))
  imp<-sv_importance(viz, kind="beeswarm")
  h.bees[[s]]<-imp+xlab(paste("SHAP values for", runlab)); names(h.bees)[[s]]<-runlab
  
  
  impvar<-names(sort(colMeans(abs(viz$S), na.rm=TRUE), decreasing=TRUE))[1:9] #12 bars used by default importance plot
  sublab<-var.lab[names(var.lab)%in%impvar]; sublab<-sublab[impvar]
  imp2<-sv_importance(viz, kind="bar", show_numbers = TRUE, fill=sublab, max_display=9L); h.imp[[s]]<-imp2
  h.imp[[s]]<-imp2+xlab(paste("mean |SHAP values| for", runlab)); names(h.imp)[[s]]<-runlab
  
}

}


r2dat.h$diff<-r2dat.h$train-r2dat.h$test
r2dat.h$flux<-"H"; r2dat.h$set<-setlab
r2dat.h

results<-rbind(results, r2dat.h)

#plots####
if(plots==TRUE){

viz<-shapviz(xgb_model, X_pred=as.matrix(X_train))
imp<-sv_importance(viz, kind="beeswarm"); imp
#imp1<-sv_waterfall(viz); imp1

imp2<-sv_importance(viz,show_numbers = TRUE); imp2
 
}

#####

rm(dat.in)

# ggplot(results) +
#   aes(x = flux, y = test, fill = set) +
#   geom_boxplot() +
#   scale_fill_hue(direction = 1) +
#   theme_minimal()

} #end model type loop



#Summary Plots#####

library(Hmisc)

write.csv(results,"plots/full_results_nomodcor.csv")


#rename "vi" to "base" to not exclude LAI
results$set[results$set=="vi"]<-"base"


ggplot(results) +
  aes(x = flux, y = test,colour = set) +
  stat_summary(position=position_dodge(width=0.5), fun="median", fun.min="min", fun.max="max") +
  scale_fill_hue(direction = 1) +
  ylab("testing R2")+
  ylim(0.6, 1)+
  theme_minimal()


ggplot(results) +
  aes(x = flux, y = test,colour = set) +
  stat_summary(position=position_dodge(width=0.7), size=1.3,fun="median", fun.min="min", fun.max="max") +
  scale_fill_hue(direction = 1) +
  ylab("testing R2")+
  ylim(0.6, 1)+
  theme_minimal()+
  labs(color="model")+
  ylab(bquote("testing"~R^2))+
  theme(axis.title = element_text(size=22), 
        axis.text=element_text(size=16), 
        legend.title=element_text(size=20),
        legend.text=element_text(size=16))

dev.copy(png, "D:/Analysis/traitflux/plots/performance_all_ogparam.png", height=525, width=750)
dev.off()


results.clear<-results[results$set%in%c("both", "trait", "base", "comp"),]


ggplot(results.clear) +
  aes(x = flux, y = test, colour = set) +
  #stat_summary(position=position_dodge(width=0.5)) +
  stat_summary(position=position_dodge(width=0.5), size=1.3, fun="median", fun.min="min", fun.max="max") +
  scale_fill_hue(direction = 1) +
  ylim(0.6, 1)+
  theme_minimal()+
  labs(color="model")+
  ylab(bquote("testing"~R^2))+
  theme(axis.title = element_text(size=20), 
        axis.text=element_text(size=14), 
        legend.title=element_text(size=18),
        legend.text=element_text(size=14))


dev.copy(png, "D:/Analysis/traitflux/plots/performance_clear_ogparam.png", height=350, width=500)
dev.off()


# results.simp<-results[results$set%in%c("trait", "base"),]
# 
# ggplot(results.simp) +
#   aes(x = flux, y = test, color = set) +
#   stat_summary(position=position_dodge(width=0.5), size=1.3, fun="median", fun.min="min", fun.max="max") +
#   scale_fill_hue(direction = 1) +
#   ylim(0.6, 1)+
#   theme_minimal()+
#   labs(color="model")+
#   ylab(bquote("testing"~R^2))+
#   theme(axis.title = element_text(size=20), 
#         axis.text=element_text(size=14), 
#         legend.title=element_text(size=18),
#         legend.text=element_text(size=14))
# 
# dev.copy(png, "D:/Analysis/traitflux/plots/performance_simp.png", height=350, width=500)
# dev.off()

results.out<-aggregate(cbind(results$train, results$test, results$diff)~results$flux+results$set, FUN="median", na.rm=TRUE)
colnames(results.out)<-c("flux", "set", "training R2", "testing R2", "trainR2 - test R2")
write.csv(results.out,"plots/results_summary_median.csv")

#####

# write.csv(results, "results.csv")
# implist<-list(gpp.imp, er.imp, h.imp,le.imp)
# save(implist, file="summer_importance.rdata")
# beelist<-list(gpp.bees, er.bees, h.bees,le.bees)
# save(beelist, file="summer_beeswarm.rdata")

library(egg)

ggarrange(
gpp.imp[[4]]+theme_minimal(),
er.imp[[4]]+theme_minimal(),
le.imp[[4]]+theme_minimal(),
h.imp[[4]]+theme_minimal(),
ncol=2, nrow=2
)

dev.copy(png, "D:/Analysis/traitflux/plots/trait_imp.png", height=600, width=750)
dev.off()


ggarrange(
  gpp.imp[[3]]+theme_minimal(),
  er.imp[[3]]+theme_minimal(),
  le.imp[[3]]+theme_minimal(),
  h.imp[[3]]+theme_minimal(),
  ncol=2, nrow=2
)

dev.copy(png, "D:/Analysis/traitflux/plots/vi_imp.png", height=600, width=750)
dev.off()



ggarrange(
  gpp.imp[[1]]+theme_minimal(),
  er.imp[[1]]+theme_minimal(),
  le.imp[[1]]+theme_minimal(),
  h.imp[[1]]+theme_minimal(),
  ncol=2, nrow=2
)

ggarrange(
  gpp.imp[[2]]+theme_minimal(),
  er.imp[[2]]+theme_minimal(),
  le.imp[[2]]+theme_minimal(),
  h.imp[[2]]+theme_minimal(),
  ncol=2, nrow=2
)
