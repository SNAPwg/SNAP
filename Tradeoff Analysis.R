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

Management$SizeLimit<- 90

Management$NTZ<- NoTakeZoneImp

Management$Season<- season

Management$Quota<- 10000

Management$Effort<- 1

Management$Gear<- 0.5

Management$Tax<- 1.1


# Run Tradeoffs -----------------------------------------------------------
DefaultFish<- Life
DefaultFleet<- Fleets

ManageSims<- list()

for (i in 1:dim(ManageStrats)[1]) #Can replace this with mclapply later if this takes too long, easier to debug this way
{

  
ManageSims[[i]]<-Master(Life,SimCTL,Fleets,season,Samp,ManageStrats[i,],Management,NoTakeZoneNULL,NoTakeZoneImp,habitat,Graphs=F,GraphsFish=F,PrintLifeHistory=T)

ManageSims[[i]]$CostOfManagement

OAtotCatch<-apply(ManageSims[[i]]$CatchByFisher,2,sum,na.rm=T)

Life<- DefaultFish

Fleets<- DefaultFleet

}

# Define life history uncertainty


