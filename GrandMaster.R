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
SimCTL<-read.csv("GrandSimCtl.csv",header=F)               # simulation controls
Fleets<-read.csv("Fleets.csv",header=F)                   # fleet characteristics  
season<-read.csv("seasonNULL.csv",header=F)           # fishing seasons by fleet
Samp <- read.csv("SamplingParams.csv")            # sampling controls for management
NoTakeZone<-read.csv("notakezoneNULL.csv",header=F)   # marine protected areas (0=open access, 1=MPA, 2=TURF?)
habitat<-read.csv("habitatNULL.csv",header=F)         # habitat quality (recruitment suitability)

OpenAccess<-Master(Life,SimCTL,Fleets,season,Samp,NoTakeZone,habitat,Graphs=F,GraphsFish=F,PrintLifeHistory=F)
OAtotCatch<-apply(OpenAccess$CatchByFisher,2,sum,na.rm=T)
OAtotCost<-apply(OpenAccess$CostByFisher,2,sum,na.rm=T)
OAtotProfit<-apply(OpenAccess$ProfitByFisher,2,sum,na.rm=T)

#==MPA====
Life<-read.csv("LifeHistory.csv")                 # life history characteristics
SimCTL<-read.csv("GrandSimCtl.csv",header=F)               # simulation controls
Fleets<-read.csv("Fleets.csv",header=F)                   # fleet characteristics  
season<-read.csv("seasonNULL.csv",header=F)           # fishing seasons by fleet
Samp <- read.csv("SamplingParams.csv")            # sampling controls for management
NoTakeZone<-read.csv("notakezone.csv",header=F)   # marine protected areas (0=open access, 1=MPA, 2=TURF?)
habitat<-read.csv("habitatNULL.csv",header=F)         # habitat quality (recruitment suitability)

halfMPA<-Master(Life,SimCTL,Fleets,season,Samp,NoTakeZone,habitat,Graphs=F,GraphsFish=F,PrintLifeHistory=F)
MPAtotCatch<-apply(halfMPA$CatchByFisher,2,sum,na.rm=T)
MPAtotCost<-apply(halfMPA$CostByFisher,2,sum,na.rm=T)
MPAtotProfit<-apply(halfMPA$ProfitByFisher,2,sum,na.rm=T)

#==Seasons
Life<-read.csv("LifeHistory.csv")                 # life history characteristics
SimCTL<-read.csv("GrandSimCtl.csv",header=F)               # simulation controls
Fleets<-read.csv("Fleets.csv",header=F)                   # fleet characteristics  
season<-read.csv("season.csv",header=F)           # fishing seasons by fleet
Samp <- read.csv("SamplingParams.csv")            # sampling controls for management
NoTakeZone<-read.csv("notakezoneNULL.csv",header=F)   # marine protected areas (0=open access, 1=MPA, 2=TURF?)
habitat<-read.csv("habitatNULL.csv",header=F)         # habitat quality (recruitment suitability)

Seas<-Master(Life,SimCTL,Fleets,season,Samp,NoTakeZone,habitat,Graphs=F,GraphsFish=F,PrintLifeHistory=F)
StotCatch<-apply(Seas$CatchByFisher,2,sum,na.rm=T)
StotCost<-apply(Seas$CostByFisher,2,sum,na.rm=T)
StotProfit<-apply(Seas$ProfitByFisher,2,sum,na.rm=T)

#==plot it all
burnIn    <-SimCTL[grep('burn',SimCTL[,2]),1]
simTimePlt <-SimCTL[grep('simTime',SimCTL[,2]),1]
par(mfrow=c(5,1),mar=c(.1,4,.1,.1))
plot(OAtotCatch,type="b",xaxt='n',las=2,ylim=c(0,max(OAtotCatch,MPAtotCatch,StotCatch,na.rm=T)))
 lines(MPAtotCatch,type="b",col=2)
 lines(StotCatch,type="b",col=3)
plot(OAtotCost,lty=2,type="b",xaxt='n',las=2,ylim=c(0,max(OAtotCost,MPAtotCost,StotCost,na.rm=T)))
lines(MPAtotCost,type="b",col=2)
lines(StotCost,type="b",col=3)
plot(OAtotProfit,lty=2,type="b",xaxt='n',las=2,ylim=c(0,max(OAtotProfit,MPAtotProfit,StotProfit,na.rm=T)))
lines(MPAtotProfit,type="b",col=2)
lines(StotProfit,type="b",col=3)
plot(OpenAccess$CostOfManagement,type="b",xaxt='n',las=2,ylim=c(0,max(OpenAccess$CostOfManagement,halfMPA$CostOfManagement,Seas$CostOfManagement,na.rm=T)))
lines(halfMPA$CostOfManagement,type="b",col=2)
lines(Seas$CostOfManagement,type="b",col=3)
plot(OpenAccess$SpawningBiomass[burnIn:simTimePlt],type="b",xaxt='n',las=2,ylab="SpawningBio",ylim=c(0,max(OpenAccess$SpawningBiomass[burnIn:simTimePlt],na.rm=T)))
lines(halfMPA$SpawningBiomass[burnIn:simTimePlt],type="b",col=2)
lines(Seas$SpawningBiomass[burnIn:simTimePlt],type="b",col=3)



#saveHTML(ani.replay(), img.name = "record_plot_oldf",outdir = getwd(),interval=0.05)
