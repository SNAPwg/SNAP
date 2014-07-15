#==specify options (e.g. what management tactic is going to be searched over)
#==make storage facility for key output (e.g. total cost vs. total profit; depletion, spr)
#==create directories
#==write csv files
#==scenarios: changing MPA size, changing size limit, changing season length, changing effort input through fishers, etc.

require(animation)
require(caTools)
source("Master.R")
source("lenwei.R")
source("VisualizeMovement.R")
source("movArray.R")
source("InitialPop.R")
source("Recruitment.R")
source("samplingFunc.R")
#rm(list=ls())

Graphs<-F
GraphsFish<-F
PrintLifeHistory<-F

#==open access===========================

Life<-read.csv("LifeHistory.csv")                 # life history characteristics
SimCTL<-read.csv("GrandSimCtl.csv")               # simulation controls
Fleets<-read.csv("Fleets2.csv")                   # fleet characteristics  
season<-read.csv("SeasonNULL.csv",header=F)           # fishing seasons by fleet
Samp <- read.csv("SamplingParams.csv")            # sampling controls for management
NoTakeZone<-read.csv("notakezoneNULL.csv",header=F)   # marine protected areas (0=open access, 1=MPA, 2=TURF?)
habitat<-read.csv("habitat.csv",header=F)         # habitat quality (recruitment suitability)

OpenAccess<-Master(Life,SimCTL,Fleets,season,Samp,NoTakeZone,habitat,Graphs=F,GraphsFish=F,PrintLifeHistory=F)

#==MPA====

Life<-read.csv("LifeHistory.csv")                 # life history characteristics
SimCTL<-read.csv("GrandSimCtl.csv")               # simulation controls
Fleets<-read.csv("Fleets2.csv")                   # fleet characteristics  
season<-read.csv("Season.csv",header=F)           # fishing seasons by fleet
Samp <- read.csv("SamplingParams.csv")            # sampling controls for management
NoTakeZone<-read.csv("notakezone.csv",header=F)   # marine protected areas (0=open access, 1=MPA, 2=TURF?)
habitat<-read.csv("habitat.csv",header=F)         # habitat quality (recruitment suitability)

#==Seasons

#==different size selection

#==do four different runs

#==collate the results


#==make a graph




outs<-Master(Life,SimCTL,Fleets,season,Samp,NoTakeZone,habitat,Graphs=F,GraphsFish=F,PrintLifeHistory=F)

totCatch<-apply(CatchByFisher,2,sum,na.rm=T)
totCost<-apply(CostByFisher,2,sum,na.rm=T)
totProfit<-apply(ProfitByFisher,2,sum,na.rm=T)

par(mfrow=c(4,1),mar=c(.1,4,.1,.1))
plot(totCatch,type="b",xaxt='n',las=2)
plot(totCost,lty=2,type="b",xaxt='n',las=2)
plot(totProfit,lty=2,type="b",xaxt='n',las=2)
lines(CostOfManagement,lty=2,type="b",col=2)
plot(SpawningBiomass[burn:simTime],type="b",xaxt='n',las=2,ylab="SpawningBio",ylim=c(0,max(SpawningBiomass[burn:simTime],na.rm=T)))




#ani.options(interval=.15)  
#ani.replay()
#saveHTML(ani.replay(), img.name = "record_plot_oldf",outdir = getwd(),interval=0.05)
