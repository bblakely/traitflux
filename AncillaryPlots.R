#Ancillary plots for traits paper


library(data.table)
library(esquisse)
library(ggplot2)
library(egg)


setwd("D:/analysis/traitflux")

#Read in data and get it to match what was used for modeling
master<-read.csv("alltower_trait_fullseason3.csv")
alldat<-master
alldat$ts<-as.POSIXct(alldat$ts_posix, format="%Y-%m-%d %H:%M:%S")



#Data prep: read in, remove bad values ####

#master<-read.csv("alltower_trait_fullseason2.csv")
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







##Composition (replaces "CompositonPlot.R")####

sites<-unique(master$sitename)


pft.dat<-data.frame()
for(i in 1:length(sites)){
  site<-sites[i]
  pft<-unique(master$igbp[master$sitename==site])
  pft.dat[i,1:2]<-c(site,pft)
  
}

colnames(pft.dat)<-c("site", "ptf")

comps<-c("Fir.Spruce", "Coniferous.Forested.Wetland","Pine","Aspen.Paper.Birch",
         "Northern.Hardwoods","Broad.leaved.Deciduous.Scrub.Shrub","Swamp.Hardwoods", 
         "Broad.leaved.Evergreen.Scrub.Shrub","Open.Water","Red.Maple",
         "Mixed.Deciduous.Coniferous.Forest","Other")

compcols<-which(colnames(master)%in%c("Fir.Spruce", "Coniferous.Forested.Wetland","Pine","Aspen.Paper.Birch",
                                      "Northern.Hardwoods","Broad.leaved.Deciduous.Scrub.Shrub","Swamp.Hardwoods", 
                                      "Broad.leaved.Evergreen.Scrub.Shrub","Open.Water","Red.Maple",
                                      "Mixed.Deciduous.Coniferous.Forest","Other"))

meanpft<-aggregate(master[,compcols], by=list(master$sitename), FUN="mean")



long <- melt(setDT(meanpft), id.vars = "Group.1")

ggplot(long) +
  aes(x = Group.1, y = value, fill = variable) +
  geom_col(position = "fill") +
  scale_fill_brewer(palette = "Paired", direction = 1) +
  xlab("CH19 tower")+
  ylab("")+
  labs(fill="vegetation composition")+
  theme_minimal()




##Trait levels by tower####


library(ggplot2)

lig.site<-ggplot(master) +
 aes(x = sitename, y = lignin, fill = igbp) + #fill=sitename
 geom_violin(adjust = 1L) +
 scale_fill_hue(direction = 1) +
 xlab("")+
 theme_minimal()

nsc.site<-ggplot(master) +
  aes(x = sitename, y = nsc, fill = igbp) + #fill=sitename
  geom_violin(adjust = 1L) +
  scale_fill_hue(direction = 1) +
  xlab("")+
  theme_minimal()

pheno.site<-ggplot(master) +
  aes(x = sitename, y = pheno, fill = igbp) + #fill=sitename
  geom_violin(adjust = 1L) +
  scale_fill_hue(direction = 1) +
  xlab("")+
  theme_minimal()

ggarrange(lig.site, nsc.site, pheno.site, nrow=3, ncol=1)


lma.site<-ggplot(master) +
  aes(x = sitename, y = lma, fill = igbp) + #fill=sitename
  geom_violin(adjust = 1L) +
  scale_fill_hue(direction = 1) +
  xlab("")+
  theme_minimal()

nitr.site<-ggplot(master) +
  aes(x = sitename, y = nitr, fill = igbp) + #fill=sitename
  geom_violin(adjust = 1L) +
  scale_fill_hue(direction = 1) +
  xlab("")+
  theme_minimal()

ggarrange(nitr.site,lma.site, nrow=2, ncol=1)

#####


#shows ER much more site-specific than GPP
library(ggplot2)

ggplot(alldat) +
 aes(x = ts, y = Reco_U50, colour = sitename) +
 geom_point() +
 scale_color_hue(direction = 1) +
 theme_minimal()


#comparative


gpp_trace<-ggplot(alldat[!is.na(alldat$FC)&alldat$ts_posix<("2019-07-15 00:00:00"),]) +
 aes(x = ts, y = GPP_U50_f, colour = sitename) +
 geom_point()+
 ylim(-12, 80)+
 scale_color_hue(direction = 1, na.translate=FALSE) +
  xlab("")+
  ylab("GPP")+
  labs(colour="site")+
  guides(colour=guide_legend(ncol=2))+
 theme_minimal()


le_trace<-ggplot(alldat[alldat$ts_posix<("2019-07-15 00:00:00"),]) +
  aes(x = ts, y = LE, colour = sitename) +
  geom_point()+
  ylim(-10,650)+
  scale_color_hue(direction = 1) +
  xlab("")+
  labs(colour="site")+
  theme_minimal()+
  theme(legend.position = "none")

h_trace<-ggplot(alldat[alldat$ts_posix<("2019-07-15 00:00:00"),]) +
  aes(x = ts, y = H, colour = sitename) +
  geom_point()+
  ylim(-10,650)+
  scale_color_hue(direction = 1) +
  xlab("")+
  labs(colour="site")+
  guides(colour=guide_legend(ncol=2))+
  theme_minimal()+
  theme(legend.position = "none")


er_trace<-ggplot(alldat[alldat$ts_posix<("2019-07-15 00:00:00"),]) +
  aes(x = ts, y = Reco_U50, colour = sitename) +
  geom_point()+
  ylim(-12, 80)+
  scale_color_hue(direction = 1) +
  xlab("")+
  ylab("ER")+
  labs(colour="site")+
  theme_minimal()+
  theme(legend.position = "none")

ggarrange(er_trace, gpp_trace, le_trace, h_trace, ncol=2, nrow=2)

dev.copy(png, "D:/Analysis/traitflux/plots/compare_traces.png", height=400, width=800)
dev.off()


library(ggplot2)

ggplot(alldat) +
 aes(x = ts, y = lai, colour = sitename) +
 geom_point() +
 scale_color_hue(direction = 1) +
 theme_minimal() +
 facet_wrap(vars(sitename))

ggplot(alldat) +
 aes(x = ts, y = evi, colour = sitename) +
 geom_point() +
 ylim(0.2,0.8)+
 scale_color_hue(direction = 1) +
 theme_minimal() +
 facet_wrap(vars(sitename))


