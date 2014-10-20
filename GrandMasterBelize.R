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
SimCTL<-read.csv("GrandSimCtlBLZ.csv",header=F)               # simulation controls
Fleets<-read.csv("Fleets.csv",header=F)                   # fleet characteristics  
season<-read.csv("season.csv",header=F)           # fishing seasons by fleet
Samp <- read.csv("SamplingParams.csv")            # sampling controls for management
NoTakeZoneNULL<-read.csv("notakezoneBLZNULL.csv",header=F)   # marine protected areas (0=open access, 1=MPA, 2=TURF?)
NoTakeZoneImp<-read.csv("notakezoneBLZ.csv",header=F)   # marine protected areas (0=open access, 1=MPA, 2=TURF?)
habitat<-read.csv("habitatBLZ.csv",header=F)         # habitat quality (recruitment suitability)
filled.contour(z=habitat)
#============================================================
#==find Rzero
#=run with specific R0 once
# OpenAccess<-Master(Life,SimCTL,Fleets,season,Samp,NoTakeZoneNULL,NoTakeZoneImp,habitat,Graphs,GraphsFish,PrintLifeHistory)
# OAtotCatch<-apply(OpenAccess$CatchByFisher,2,sum,na.rm=T)
# compCatch<-sum(OAtotCatch[502:514])
# #==set catch in last year
# EqCatch<-5071
# highR0<-2*Life[18,1]
# lowR0<-1
# RzeroTest<-Life[18,1]
# RzeroFin<-12700.73
# for(x in 1:20)
# {
# #==change Rzero in "life.csv"
# if(compCatch<EqCatch)
#   lowR0<-RzeroTest
# if(compCatch>=EqCatch)
#   highR0<-RzeroTest
# 
#  Life[18,1]<-(highR0+lowR0)/2
#  RzeroTest<-Life[18,1]
#  OpenAccess<-Master(Life,SimCTL,Fleets,season,Samp,NoTakeZoneNULL,NoTakeZoneImp,habitat,Graphs,GraphsFish,PrintLifeHistory)
#  OAtotCatch<-apply(OpenAccess$CatchByFisher,2,sum,na.rm=T)
#  compCatch<-sum(OAtotCatch[502:514])
#  print(compCatch)
#  print(EqCatch)
#  print(Life[18,1])
# }
# 
# #===================================================================
# #==find costs given depletion
# 





OpenAccess<-Master(Life,SimCTL,Fleets,season,Samp,NoTakeZoneNULL,NoTakeZoneImp,habitat,Graphs,GraphsFish,PrintLifeHistory)
OAtotCatch<-apply(OpenAccess$CatchByFisher,2,sum,na.rm=T)
OAtotCost<-apply(OpenAccess$CostByFisher,2,sum,na.rm=T)
OAtotProfit<-apply(OpenAccess$ProfitByFisher,2,sum,na.rm=T)

burnInt   <-SimCTL[grep('burn',SimCTL[,2]),1]
simTimePlt <-SimCTL[grep('simTime',SimCTL[,2]),1]
initManage<-SimCTL[grep('initManage',SimCTL[,2]),1]   # year in which to initiate management
yearMark2  <-SimCTL[grep('yearMark',SimCTL[,2]),1]    	# number of time steps in a year
SpcRow    <-SimCTL[grep('SpaceR',SimCTL[,2]),1]	    	# Rows in the grid space
SpcCol	  <-SimCTL[grep('SpaceC',SimCTL[,2]),1]  			# cols in the grid spcae

tempMat<-matrix(0,ncol=SpcCol,nrow=SpcRow)
for(x in (initManage-burnInt):(simTimePlt-burnInt))
  tempMat<-tempMat+OpenAccess$SpaceEffort[x,,,]

tempMat2<-matrix(0,ncol=SpcCol,nrow=SpcRow)
for(x in 1:(initManage-burnInt))
  tempMat2<-tempMat2+OpenAccess$SpaceEffort[x,,,]

tempMat3<-matrix(0,ncol=SpcCol,nrow=SpcRow)
for(x in (simTimePlt-burnInt-36):(simTimePlt-burnInt))
  tempMat3<-tempMat3+OpenAccess$SpaceEffort[x,,,]


# dev.off()
# filled.contour(tempMat)
# filled.contour(tempMat2)
# filled.contour(tempMat3)


par(mfrow=c(7,1),mar=c(.1,4,.1,.1))
plot(OAtotCatch,type="b",xaxt='n',las=2,ylim=c(0,max(OAtotCatch,na.rm=T)),ylab="Total Catch",pch=16)
abline(v=initManage-burnInt,col=2,lty=2)
legend("topright",col=2,lty=2,"Management implemented",bty='n')
plot(OAtotCost,lty=2,type="b",pch=16,xaxt='n',las=2,ylim=c(0,max(OAtotCost,na.rm=T)),ylab="Total cost of fishing")
abline(v=initManage-burnInt,col=2,lty=2)

plot(OAtotProfit,lty=2,type="b",pch=16,xaxt='n',las=2,ylim=c(0,max(OAtotProfit,na.rm=T)),ylab="Total profit of fishing")
abline(v=initManage-burnInt,col=2,lty=2)

plot(OpenAccess$CostOfManagement[burnInt:(simTimePlt)],pch=16,type="b",xaxt='n',las=2,ylim=c(0,max(OpenAccess$CostOfManagement,na.rm=T)),
     ylab="Cost of MPA")
abline(v=initManage-burnInt,col=2,lty=2)

plot(OpenAccess$SpawningBiomass[burnInt:simTimePlt],pch=16,type="b",xaxt='n',las=2,
     ylab="SpawningBio",ylim=c(0,max(OpenAccess$SpawningBiomass[burnInt:simTimePlt],na.rm=T)))
abline(v=initManage-burnInt,col=2,lty=2)

plot(OpenAccess$ExploitableNumbers[burnInt:simTimePlt],pch=16,type="b",xaxt='n',las=2,
     ylab="Exploitable Numbers",ylim=c(0,max(OpenAccess$ExploitableNumbers[burnInt:simTimePlt],na.rm=T)))

plot(OpenAccess$OutsideMPAspbio[burnInt:simTimePlt],ylim=c(0,max(OpenAccess$OutsideMPAspbio[burnInt:simTimePlt],OpenAccess$InsideMPAspbio[burnInt:simTimePlt])),type='b',pch=16)
lines(OpenAccess$InsideMPAspbio[burnInt:simTimePlt],type='b',pch=16,col=2)
abline(v=initManage-burnInt,col=2,lty=2)

legend("topright",col=c(1,2),pch=16,legend=c("Outside MPA","Inside MPA"),bty='n')
sum(OAtotCatch[191:201])

#== calculate statistics for before and after managmeent implementation
#==find a range of years in each management regime for comparison
#==when recruitment happens once a year, they must both start at the same time in the reproductive cycle

timeBack<-60

timeInd11<-initManage-timeBack-burnInt
timeInd12<-initManage-burnInt

timeInd21<-simTimePlt-burnInt-timeBack
timeInd22<-simTimePlt-burnInt

timeInd21<-initManage+timeBack-burnInt
timeInd22<-initManage-burnInt+1

plot(OAtotCatch,type="l")
points(y=100,x=timeInd11)
points(y=100,x=timeInd12)
points(y=100,x=timeInd21)
points(y=100,x=timeInd22)
# timeInd21<-simTimePlt-burnInt-timeBack
# timeInd22<-simTimePlt-burnInt
sum(OAtotCatch[(timeInd11):(timeInd12)])
sum(OAtotCatch[(timeInd21):(timeInd22)])

sum(OAtotCost[(timeInd11):(timeInd12)])
sum(OAtotCost[(timeInd21):(timeInd22)])

sum(OAtotProfit[(timeInd11):(timeInd12)])
sum(OAtotProfit[(timeInd21):(timeInd22)])

sum(OpenAccess$SpawningBiomass[(timeInd11):(timeInd12)])
sum(OpenAccess$SpawningBiomass[(timeInd21):(timeInd22)])

sum(OpenAccess$OutsideMPAspbio[(timeInd11):(timeInd12)])
sum(OpenAccess$OutsideMPAspbio[(timeInd21):(timeInd22)])

sum(OpenAccess$InsideMPAspbio[(timeInd11):(timeInd12)])
sum(OpenAccess$InsideMPAspbio[(timeInd21):(timeInd22)])

Cat1<-mean(OAtotCatch[(timeInd11):(timeInd12)])
Cat2<-mean(OAtotCatch[(timeInd21):(timeInd22)])
CatchChange<-(Cat2-Cat1)/Cat1

Cost1<-mean(OAtotCost[(timeInd11):(timeInd12)])
Cost2<-mean(OAtotCost[(timeInd21):(timeInd22)])
CostChange<-(Cost2-Cost1)/Cost1

Prof1<-mean(OAtotProfit[(timeInd11):(timeInd12)])
Prof2<-mean(OAtotProfit[(timeInd21):(timeInd22)])
ProfitChange<-(Prof2-Prof1)/Prof1

SpB1<-mean(OpenAccess$SpawningBiomass[(timeInd11):(timeInd12)])
SpB2<-mean(OpenAccess$SpawningBiomass[(timeInd21):(timeInd22)])
SpBChange<-(SpB2-SpB1)/SpB1

MPAout1<-mean(OpenAccess$OutsideMPAspbio[(timeInd11):(timeInd12)])
MPAout2<-mean(OpenAccess$OutsideMPAspbio[(timeInd21):(timeInd22)])
MPAoutChg<-(MPAout2-MPAout1)/MPAout1

MPAin1<-mean(OpenAccess$InsideMPAspbio[(timeInd11):(timeInd12)])
MPAin2<-mean(OpenAccess$InsideMPAspbio[(timeInd21):(timeInd22)])
MPAinChg<-(MPAin2-MPAin1)/MPAin1

round(rbind(CatchChange,CostChange,ProfitChange,SpBChange,MPAoutChg,MPAinChg),2)
#saveHTML(ani.replay(), img.name = "record_plot_oldf",outdir = getwd(),interval=0.05)
