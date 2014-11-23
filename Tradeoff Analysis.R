################## SNAP DPSA Tradeoff Analayis Framework #######################

# This wrapper enables the use of the Master simulation framework developed by Cody to test a variety of management interventions. 
#Version 1.0 will be for one species under multiple management options, with potential for life history uncertainty 

rm(list=ls())
sapply(list.files(pattern="[.]R$", path="Functions", full.names=TRUE), source)
source('MasterForTradeoff.R')
library(animation)
library(caTools)
library(plyr)
library(ggplot2)
library(lattice)
library(ggthemes)
source("lenwei.R")
source("VisualizeMovement.R")
source("movArray.R")
source("InitialPop.R")
source("Recruitment.R")
source("samplingFunc.R")


# Storage Settings --------------------------------------------------------


Site<- 'Belize'

Species<- 'Lobster'

RunName<- 'Test Run'


Fishery<- paste(Site,Species)

FisheryPlace<- paste(Site,Species,sep='/')

RunName<- paste(FisheryPlace,RunName,sep='/')

dir.create(RunName)

FigureFolder<- paste(RunName,'/Figures',sep='')

ResultFolder<- paste(RunName,'/Results',sep='')

dir.create(FigureFolder)

dir.create(ResultFolder)


# Load Fishery ------------------------------------------------------------

OriginalWorkingDir<- getwd()

setwd(FisheryPlace)

Life<-read.csv("LifeHistory.csv")                 # life history characteristics
SimCTL<-read.csv("GrandSimCtlBLZ.csv",header=F)               # simulation controls
Fleets<-read.csv("Fleets.csv",header=F)                   # fleet characteristics  
season<-read.csv("season.csv",header=F)           # fishing seasons by fleet
Samp <- read.csv("SamplingParams.csv")            # sampling controls for management
NoTakeZoneNULL<-read.csv("notakezoneBLZNULL.csv",header=F)   # marine protected areas (0=open access, 1=MPA, 2=TURF?)
NoTakeZoneImp<-read.csv("notakezoneBLZ.csv",header=F)   # marine protected areas (0=open access, 1=MPA, 2=TURF?)
habitat<-read.csv("habitatBLZ.csv",header=F)       
ManageStrats<- read.csv('ManagementStrategies.csv')
setwd(OriginalWorkingDir)

# Set Management ----------------------------------------------------------

Management<- NULL

Management$SizeLimit<- 181

Management$NTZ<- NoTakeZoneImp

Management$Season<- season

Management$Quota<- 10000

Management$Effort<- 5

Management$Gear<- 0.5

Management$Tax<- 1.1


# Run Tradeoffs -----------------------------------------------------------

ManageSims<- list()

ManageResults<-as.data.frame(matrix(NA,nrow=dim(ManageStrats)[1],ncol=6))

colnames(ManageResults) <- c('ManagementPlan','Catch','FishingCost','FishingProfit','ManagementCost','SpawningBiomass')

for (i in 1:dim(ManageStrats)[1]) #Can replace this with mclapply later if this takes too long, easier to debug this way
{
  show(paste(round(100*(i/dim(ManageStrats)[1])),'% done with management iterations',sep=''))
  
  ManageSims[[i]]<-Master(Life,SimCTL,Fleets,season,Samp,ManageStrats[i,],Management,NoTakeZoneNULL,NoTakeZoneImp,habitat,Graphs=F,GraphsFish=F,PrintLifeHistory=F)
  
  ManageResults$ManagementPlan[i]<- as.character(ManageStrats$Plan[i] )
  
  ManageResults$Catch[i] <- sum(ManageSims[[i]]$CatchByFisher,na.rm=T)
  
  ManageResults$FishingCost[i] <- sum(ManageSims[[i]]$CostByFisher,na.rm=T)
  
  ManageResults$FishingProfit[i] <- sum(ManageSims[[i]]$ProfitByFisher,na.rm=T)
  
  ManageResults$ManagementCost[i] <- sum(ManageSims[[i]]$CostOfManagement,na.rm=T)
  
  BioPath<- (ManageSims[[i]]$SpawningBiomass)
  
  ManageResults$SpawningBiomass[i] <- BioPath[length(BioPath)]
  
}

m<- ManageSims[[2]]

pdf(file=paste(FigureFolder,'Profit and Biomass Tradeoff.pdf',sep='/'),width=8,height=6)
print(ggplot(ManageResults,aes(FishingProfit,SpawningBiomass,size=ManagementCost)) +
        #   geom_point(aes(shape=ManagementPlan,color=ManagementCost),size=5)+
        geom_point(aes(color=ManagementPlan,size=ManagementCost),size=5)+
        xlab('Cumulative Fishing Profits')+
        ylab('Final Spawning Biomass'))
#   theme_economist()+scale_colour_economist())
dev.off()


# pdf(file=paste(FigureFolder,'Profit and Biomass Tradeoff Panel.pdf',sep='/'),width=8,height=6)
# print(ggplot(ManageResults,aes(FishingProfit,SpawningBiomass)) +
#         geom_point(aes(color=ManagementCost),size=10)+
#        facet_wrap(~ManagementPlan)+
#         xlab('Cumulative Fishing Profits')+
#         ylab('Final Spawning Biomass'))
# dev.off()


# ggplot(nmmaps, aes(date,temp))+geom_point(color="chartreuse4")+ facet_wrap(~year, ncol=2)

