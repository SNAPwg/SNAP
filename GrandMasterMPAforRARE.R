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
PrintLifeHistory<-T

#==open access===========================

Life<-read.csv("LifeHistoryRARE.csv")                 # life history characteristics
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

tempMat<-matrix(0,ncol=10,nrow=10)
for(x in 1:275)
  tempMat<-tempMat+OpenAccess$SpaceEffort[x,,,]
sum(tempMat)
dev.off()
filled.contour(tempMat)

#==MPA====
Life<-read.csv("LifeHistoryRARE.csv")                 # life history characteristics
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

tempMat<-matrix(0,ncol=10,nrow=10)
for(x in 1:275)
  tempMat<-tempMat+halfMPA$SpaceEffort[x,,,]


filled.contour(tempMat)
sum(tempMat)

#==plot it all
burnIn    <-SimCTL[grep('burn',SimCTL[,2]),1]
simTimePlt <-SimCTL[grep('simTime',SimCTL[,2]),1]

par(mfrow=c(6,1),mar=c(.1,4,.1,.1))
plot(OAtotCatch,type="b",xaxt='n',las=2,ylim=c(0,max(OAtotCatch,MPAtotCatch,na.rm=T)),ylab="Total Catch",pch=16)
 lines(MPAtotCatch,type="b",col=2,pch=16)

plot(OAtotCost,lty=2,type="b",pch=16,xaxt='n',las=2,ylim=c(0,max(OAtotCost,MPAtotCost,na.rm=T)),ylab="Total cost of fishing")
lines(MPAtotCost,type="b",col=2,pch=16)

# plot(OAtotProfit,lty=2,type="b",pch=16,xaxt='n',las=2,ylim=c(0,max(OAtotProfit,MPAtotProfit,na.rm=T)),ylab="Total profit of fishing")
plot(OAtotProfit,lty=2,type="b",pch=16,xaxt='n',las=2,ylim=c(0,200),ylab="Total profit of fishing")
lines(MPAtotProfit,type="b",col=2,pch=16)
# lines(MPAtotProfit,type="b",col=2,pch=16)

plot(OpenAccess$CostOfManagement,pch=16,type="b",xaxt='n',las=2,ylim=c(0,max(OpenAccess$CostOfManagement,halfMPA$CostOfManagement,na.rm=T)),
     ylab="Cost of MPA")
lines(halfMPA$CostOfManagement,pch=16,type="b",col=2)

plot(OpenAccess$SpawningBiomass[burnIn:simTimePlt],pch=16,type="b",xaxt='n',las=2,
#      ylab="SpawningBio",ylim=c(0,max(OpenAccess$SpawningBiomass[burnIn:simTimePlt],na.rm=T)))
     ylab="SpawningBio",ylim=c(0,2100))
lines(halfMPA$SpawningBiomass[burnIn:simTimePlt],pch=16,type="b",col=2)

plot(halfMPA$InsideMPAspbio[burnIn:simTimePlt])
lines(halfMPA$OutsideMPAspbio[burnIn:simTimePlt])

legend("topright",col=c(1,2),pch=16,legend=c("Open access","MPA"),bty='n')



DiffMPA<-cbind((OAtotCatch-MPAtotCatch)/OAtotCatch,(OAtotProfit-MPAtotProfit)/OAtotProfit,
 (OpenAccess$SpawningBiomass[burnIn:simTimePlt]-halfMPA$SpawningBiomass[burnIn:simTimePlt])/OpenAccess$SpawningBiomass[burnIn:simTimePlt]
step1<-251
step2<-281
sum(OAtotCatch[step1:step2])
sum(MPAtotCatch[step1:step2])

sum(OAtotProfit[step1:step2])
sum(MPAtotProfit[step1:step2])

sum(OpenAccess$SpawningBiomass[step1:step2])
sum(halfMPA$SpawningBiomass[step1:step2])
sum((OpenAccess$SpawningBiomass[step1:step2]-halfMPA$SpawningBiomass[step1:step2]))/(OpenAccess$SpawningBiomass[step1:step2]))

dev.off()
boxplot(DiffMPA,ylim=c(-2,2))

#saveHTML(ani.replay(), img.name = "record_plot_oldf",outdir = getwd(),interval=0.05)
