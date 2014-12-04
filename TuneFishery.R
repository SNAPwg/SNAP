################## SNAP DPSA Model Tuning Framework #######################

# This wrapper enables the use of the Master simulation framework developed by Cody to test a variety of management interventions. 
#Version 1.0 will be for one species under multiple management options, with potential for life history uncertainty 

rm(list=ls())
sapply(list.files(pattern="[.]R$", path="Functions", full.names=TRUE), source)
library(animation)
library(caTools)
library(plyr)
library(ggplot2)
library(gridExtra)


library(lattice)
library(ggthemes)
source("lenwei.R")
source("VisualizeMovement.R")
source("movArray.R")
source("InitialPop.R")
source("Recruitment.R")
source("samplingFunc.R")


# Storage Settings --------------------------------------------------------


Site<- 'Peru'

Species<- 'Lorna Drum'

RunName<- 'Test'


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
SimCTL<-read.csv("GrandSimCtlPeru.csv",header=F)               # simulation controls
Fleets<-read.csv("Fleets.csv",header=F)                   # fleet characteristics  
season<-read.csv("season.csv",header=F)           # fishing seasons by fleet
Samp <- read.csv("SamplingParams.csv")            # sampling controls for management
NoTakeZoneNULL<-read.csv("notakezonePeruNULL.csv",header=F)   # marine protected areas (0=open access, 1=MPA, 2=TURF?)
NoTakeZoneImp<-read.csv("notakezonePeru.csv",header=F)   # marine protected areas (0=open access, 1=MPA, 2=TURF?)
habitat<-read.csv("habitatPeru.csv",header=F)       
ManageStrats<- read.csv('ManagementStrategies.csv')
setwd(OriginalWorkingDir)

# Set Management ----------------------------------------------------------

Management<- NULL

Management$SizeLimit<- 320

Management$NTZ<- NoTakeZoneImp

shortseason<- season

shortseason[,2]<- NA

shortseason[1:4,2]<- 1 

Management$Season<- shortseason

Management$Quota<- 10000

Management$Effort<- 12

Management$Gear<- 1.5

Management$Tax<- 1.3


# Run Tradeoffs -----------------------------------------------------------


ManageResults<-as.data.frame(matrix(NA,nrow=dim(ManageStrats)[1],ncol=6))

colnames(ManageResults) <- c('ManagementPlan','Catch','FishingCost','FishingProfit','ManagementCost','SpawningBiomass')

#Turn Off Management 

ManageSims<-Master(Life,SimCTL,Fleets,season,Samp,ManageStrats[1,],Management,NoTakeZoneNULL,NoTakeZoneImp,habitat,Graphs=F,GraphsFish=F,PrintLifeHistory=T)

TimeLength<- dim(ManageSims[[1]]$CatchByFisher)[2] #time run for fishing scenarios

SimLength<- SimCTL[grep('simTime',SimCTL[,2]),1] #Complete time with burn in

BurnIn<- SimCTL[grep('burn',SimCTL[,2]),1] #Burn in period


MSE <- ldply(ManageSims, function(mod2) 
{
  
  data.frame(rep(as.character(mod2$Fishery$ManagementPlan[1]),TimeLength),
             rep(mod2$Fishery$Iteration,TimeLength),
             1:TimeLength,
             colSums(mod2$CatchByFisher,na.rm=T),
             colSums(mod2$ProfitByFisher,na.rm=T),
             colSums(mod2$CostByFisher,na.rm=T),
             (mod2$CostOfManagement[BurnIn:SimLength]),
             mod2$SpawningBiomass[BurnIn:SimLength],
             rep(mod2$Fishery$mat50,TimeLength),
             rep(mod2$Fishery$VirSpBio,TimeLength),
             rep(sum(mod2$CostOfManagement[BurnIn:SimLength],na.rm=T),TimeLength)
  )
}
)


colnames(MSE)<- c('ManagementPlan','Iteration','TimeStep','Catch','Profits','FishingCosts',
                  'ManagementCosts','SpawningBiomass','mat50','VirginSpawningBiomass','TotalManagementCosts')

MSE$ManagementCosts[is.na(MSE$ManagementCosts)]<- 0

MSE$Year<- floor(MSE$TimeStep/12)


MSE$RunName<- paste(MSE$ManagementPlan,MSE$Iteration,sep='-')

MSE_Totals<- ddply(MSE,c('ManagementPlan','Iteration'),summarize,TotalManagementCosts=sum(ManagementCosts,na.rm=T),TotalProfits=sum(Profits,na.rm=T),
                   NPV=sum(Profits*(1+0.05)^-TimeStep),FinalSSB=mean(SpawningBiomass[Year==max(Year)]),FinalProfits=Profits[TimeStep==max(TimeStep)])


MSE_ByMonth<- ddply(MSE,c('ManagementPlan','TimeStep'),summarize,Costs=mean(ManagementCosts,na.rm=T),Profits=mean(Profits,na.rm=T),
                    NPV=mean(Profits*(1+0.05)^-TimeStep),SSB=mean(SpawningBiomass))

MSE_ByYear<- ddply(MSE,c('ManagementPlan','Iteration','Year'),summarize,ManagementCosts=sum(ManagementCosts,na.rm=T),Profits=sum(Profits,na.rm=T),
                   PV=sum(Profits*(1+0.05)^-TimeStep),SSB=mean(SpawningBiomass),Depletion=mean(SpawningBiomass)/mean(VirginSpawningBiomass),Catch=sum(Catch,na.rm=T))


MSE_ByYear<- ddply(MSE_ByYear,c('ManagementPlan','Year'),summarize,ManagementCosts=mean(ManagementCosts,na.rm=T),Profits=mean(Profits,na.rm=T),
                   PV=mean(PV),SSB=mean(SSB),Depletion=mean(Depletion),Catch=mean(Catch))


pdf(file=paste(FigureFolder,'Final Profit and Final Biomass Tradeoff No Cost.pdf',sep='/'),width=8,height=6)
print(ggplot(subset(MSE,TimeStep==max(TimeStep)),aes(Profits,SpawningBiomass))+
        geom_point(aes(color=ManagementPlan),size=10,alpha=0.6)+
        scale_size_continuous(range = c(5, 15))+
        xlab('Final Fishing Profits')+
        ylab('Final Spawning Biomass'))
dev.off()


pdf(file=paste(FigureFolder,'Final Profit and Final Biomass Tradeoff.pdf',sep='/'),width=8,height=6)
print(ggplot(subset(MSE,TimeStep==max(TimeStep)),aes(Profits,SpawningBiomass))+
        geom_point(aes(color=ManagementPlan,size=TotalManagementCosts),alpha=0.6)+
        scale_size_continuous(range = c(5, 15))+
        xlab('Final Fishing Profits')+
        ylab('Final Spawning Biomass'))
dev.off()

pdf(file=paste(FigureFolder,'Cumulative Tradeoffs.pdf',sep='/'),width=8,height=6)
print(ggplot(MSE_Totals,aes(TotalProfits,FinalSSB))+
        geom_point(aes(color=ManagementPlan),size=10,alpha=0.6)+
        xlab('Cumulative Fishing Profits')+
        ylab('Final Spawning Biomass'))
dev.off()


pdf(file=paste(FigureFolder,'Cumulative Tradeoffs 2.pdf',sep='/'),width=8,height=6)
print(ggplot(MSE_Totals,aes(TotalProfits,FinalSSB))+
        geom_point(aes(color=ManagementPlan,size=TotalManagementCosts),alpha=0.6)+
        xlab('Cumulative Fishing Profits')+
        ylab('Final Spawning Biomass')+
        scale_size_continuous(range = c(5, 15)))


dev.off()


pdf(file=paste(FigureFolder,'Cumulative Tradeoffs 3.pdf',sep='/'),width=8,height=6)
print(ggplot(MSE_Totals,aes(NPV,FinalSSB))+
        geom_point(aes(color=ManagementPlan,size=TotalManagementCosts),alpha=0.6)+
        xlab('Net Present Fishing Profits')+
        ylab('Final Spawning Biomass')+
        scale_size_continuous(range = c(5, 15)))


dev.off()

pdf(file=paste(FigureFolder,'TimeTrend.pdf',sep='/'),width=8,height=6)

CatchTrend=(ggplot(subset(MSE_ByYear,Year<= (max(Year)-1)),aes(Year,Catch))+
               geom_line(aes(color=ManagementPlan),size=1,alpha=0.6)+
               xlab('Year')+
               ylab('Catch')
             +theme(legend.position="none"))

CostsTrend=(ggplot(subset(MSE_ByYear,Year<= (max(Year)-1)),aes(Year,ManagementCosts))+
              geom_line(aes(color=ManagementPlan),size=1,alpha=0.6)+
              xlab('Year')+
              ylab('Management Costs')
            +theme(legend.position="none"))


SSBTrend=(ggplot(subset(MSE_ByYear,Year<= (max(Year)-1)),aes(Year,SSB))+
            geom_line(aes(color=ManagementPlan),size=1,alpha=0.6)+
            xlab('Year')+
            ylab('SSB')
            +theme(legend.position="none"))



DepletionTrend=(ggplot(subset(MSE_ByYear,Year<= (max(Year)-1)),aes(Year,Depletion))+
            geom_line(aes(color=ManagementPlan),size=1,alpha=0.6)+
            xlab('Year')+
            ylab('SSB / Virgin SSB') )

print(grid.arrange(CatchTrend,CostsTrend,SSBTrend,DepletionTrend,ncol=2))
dev.off()


pdf(file=paste(FigureFolder,'Profit SSB Phase Space.pdf',sep='/'),width=8,height=6)

print(ggplot(subset(MSE_ByYear,Year<=22),aes(SSB,Profits))+
        geom_line(aes(color=Year),size=1,alpha=1)+
        xlab('SSB')+
        ylab('Profits')+facet_wrap(~ManagementPlan))
dev.off()

  save.image(file=paste(ResultFolder,'/FinalResults.rdata',sep=''))


