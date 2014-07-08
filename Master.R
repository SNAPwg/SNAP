require(animation)
require(caTools)
setwd("C:/Users/Cody/Desktop/SpatialOM/")
source("lenwei.R")
source("VisualizeMovement.R")
source("movArray.R")
source("InitialPop.R")
source("Recruitment.R")
#rm(list=ls())

Graphs<-F
GraphsFish<-T
PrintLifeHistory<-T

#read.csv("LifeHistory.csv")


#==life history ==============================
kmax	 	<-20		# max age
kmat	 	<-2.63	# age at maturity
kmin	 	<-1		# minimium size
M   	 	<-0.34	# natural mortality
Linf	 	<-18.3	# VonBert Linf
K		<-.24		# VonBert K
t0  	 	<-.44		# VonBert t0
wtA		<-0.0046	# weight A
wtB		<-2.63	# weight B
mat50		<-2.5		# age at 50% maturity
mat95		<-4		# age at 95% maturity 

 MatAge<-1/(1+exp(-log(19)*((seq(1,kmax)-mat50)/(mat95-mat50))))
 lenAtAge<-rep(0,kmax)
 wgtAtAge<-rep(0,kmax)
 for(k in 1:kmax)
  {
   lenAtAge[k]<-Linf*(1-exp(-K*k-t0))
   wgtAtAge[k]<-wtA*lenAtAge[k]^wtB
  }

#==space and movement=========================
SpaceR	<-10		# Rows in the grid space
SpaceC	<-10		# cols in the grid spcae
RecDist	<-1		# How is recruitment distributed over space? (1=equally,2=determined by adults,3=determined by habitat.csv)
sdx		<-1.5		# sd of movement on x axis (in units of rows/columns)
sdy		<-1.5		# sd of movement on y axis
movP50	<-0		# probability moving by length (logistic; 50%; set both to 0 if all move)
movP95	<-0		# probability of moving at length (95%)

MoveProb<-1/(1+exp(-log(19)*((seq(1,kmax)-movP50)/(movP95-movP50))))

#==recruitment================================
detRec	<-1		# deterministic recruitment?
R0		<-1000	# Virgin recruitment
steepness	<-0.85	# steepness
det		<-1		# deterministic
sigmaR	<-0.5		# variability around recruitment curve
initDepl	<-1		# depletion at the beginning of 'fishery management' (calculate the fraction of R0 to put into the pop before allow movement to equilibrate)

#== simulation duration======================
burn		<-25		# burn in for movement equilibration (use option "Graphs" in "Initpop" to see that population is thoroughly mixed)
simTime  	<-205		# time steps for projection
yearMark	<-12		# number of time steps in a year

#==fishery characteristics===================
SizeLimit	<-0
Sel50		<-1		# selectivity by length (9 by len)
Sel95		<-2		# selectividad by len (13 by len)
q		<-1		# fishery catchability coefficient (what fraction of the exploitable biomass in a patch is catchable)
PortLc	<-c(1,1)	# location of the port of entry (x,y)
season	<-c(1,2,3)	# months in which fishing is allowed
#season	<-seq(0,11)	# no season
FishSel<-q/(1+exp(-log(19)*((seq(1,kmax)-Sel50)/(Sel95-Sel50))))

#==econ pars=================================
price		<-4		# ex-vessel price per unit harvest (kg)
costTrv	<-5		# cost to travel one patch
costFish	<-50		# cost per unit effort fishing
discRate	<-0.05	# discount rate
TimeHor	<-20		# time horizon for evaluation
Fishers	<-10		# number of fishers
maxCapac	<-50		# max capacity of a single fisherman in a timestep (kg) (can also be a 'trip limit'/permit implementation)

#==management costs==========================
MPAsunk	<-100		# start up cost of enforcing an MPA
MPAcost	<-10		# cost of maintaining a unit of MPA per unit time
SizeSunk	<-20		# start up cost of enforcing a size limit
SizeCost <-2		# cost of enforcing a size limit per unit of time per port
SeasonSunk	<-50		# start up cost of enforcing a season
SeasonCost	<-1		# cost of enforcing a season per unit of time

#==management auxilliary benefits======================
MPAtourism	<-30		# average revenue generated per unit MPA per unit time for tourism in an MPA


#==spatial matrix (0=open access, 1=TURF, 2=NTZ)
#==read in a csv with the appropriate letters denoted for the different uses
SpaceUse<-matrix(0,ncol=SpaceC,nrow=SpaceR)

#==read in csv that give habitat quality, determine movement by habitat quality
#==read csv that shows MPAs, if any
#NoTakeZone<-read.csv("notakezone40.csv",header=F)
#NoTakeZone<-read.csv("notakezone80.csv",header=F)
NoTakeZone<-read.csv("notakezoneNULL.csv",header=F)

#==add in dispersal .csv too
habitat<-read.csv("habitat.csv",header=F)
#is.numeric(habitat)
#habitat<-read.csv("habitatNULL.csv",header=F)


if(PrintLifeHistory==T)
{
 #==print out maturity, selectivity, growth, etc.
 dev.new()
 par(mfrow=c(2,1),mar=c(0.1,4,.1,4),oma=c(4,0,0,0))
 plot(MatAge,type="l",xaxt='n',las=2,ylab='')
 mtext(side=2,"Probability",line=3)
 lines(FishSel,lty=2,col=2)
 lines(MoveProb,lty=3,col=3)
 legend("bottomright",lty=c(1,2,3),col=c(1,2,3),legend=
 c("Maturity","Selectivity","Movement Probability"),bty='n')
 plot(lenAtAge,type="l",las=1,ylab='',ylim=c(0,max(lenAtAge)))
 par(new=T)
 plot(wgtAtAge,type="l",lty=2,col=2,yaxt='n',ylab='',ylim=c(0,max(wgtAtAge)))
 legend("bottomright",lty=c(1,2),col=c(1,2),legend=
 c("Length","Weight"),bty='n')
 axis(side=4,las=2)
 mtext(side=4,"Weight (kg)",line=2)
 mtext(side=2,"Length (cm)",line=2)
 mtext(side=1,"Age (yr)",line=2)
 
 dev.new()
 filled.contour(matrix(as.numeric(unlist(habitat)),ncol=SpaceC),y=seq(1,SpaceR),x=seq(1,SpaceC)) 
 mtext(side=3,"Habitat quality")

 dev.new()
 filled.contour(matrix(as.numeric(unlist(NoTakeZone)),ncol=SpaceC),y=seq(1,SpaceR),x=seq(1,SpaceC)) 
 mtext(side=3,"No take zones (0)")
 #==identify the MPAs and mark them on the graph

}

#===========================================================
#== Begin simulations====================================
#===========================================================

#==population dynamics array tracking number at length (or age) in a patch by year
SpaceNumAtLenT<-array(dim=c((simTime+burn),SpaceR,SpaceC,kmax))
CatchByFisher<-matrix(ncol=simTime-burn,nrow=Fishers)
BenefitAtTimeStep<-rep(0,simTime-burn)
CostByFisher<-matrix(ncol=simTime-burn,nrow=Fishers)
SpawningBiomass<-rep(0,simTime-burn)
CostOfManagement<-rep(0,simTime-burn)

#==costs of management
if(sum(1-NoTakeZone)>1 & timeStep==burn)
 CostOfManagement[timeStep]<- CostOfManagement[timeStep] + MPAsunk
if(SizeLimit>0)
 CostOfManagement[timeStep]<- CostOfManagement[timeStep] + SizeSunk
if(length(season) < yearMark ) 
 CostOfManagement[timeStep]<- CostOfManagement[timeStep] + SeasonSunk 

 CostOfManagement[timeStep]<- CostOfManagement[timeStep] + sum(1-NoTakeZone)*MPAcost
 CostOfManagement[timeStep]<- CostOfManagement[timeStep] + SizeCost
 CostOfManagement[timeStep]<- CostOfManagement[timeStep] + (yearMark - length(season))* SeasonCost

#==index for spatial area==========================
 coords<-NULL
 for(x in 1:SpaceR)
  {
   temp<-cbind(rep(x,SpaceC),seq(1,SpaceC))
   coords<-rbind(coords,temp)
  }

#==make movement array==============================
#=should there be two movement arrays?
#=one for adults and one for recruitment?
Movement<-movArray(SpaceR,SpaceC,sdx,sdy)

#==initialize population and allow movement to equilibrate
SpaceNumAtLenT[1:burn,,,]<-InitialPop(R0,M,kmax,SpaceR,SpaceC,MoveProb,coords,Movement,Graphs=F,burn)

#==calculate virgin spawning biomass (1843.584)================
tempSpBio<-rep(0,kmax)
for(i in 1:kmax)
 tempSpBio[i]<-sum(SpaceNumAtLenT[burn,,,i]*MatAge[i])
VirSpBio<-sum(tempSpBio)

#==begin projected harvest=================== 
#==fishers decide where to fish============
#==initial costs
 Distance<-sqrt((row(SpaceUse)-PortLc[1])^2 + (col(SpaceUse)-PortLc[2])^2)
 CostPatch<-Distance*costTrv+costFish
#==relate the cost of fishing to the number of fish in a patch??

    if(Graphs==T | GraphsFish==T)
     {
      setwd("C:/Users/Cody/Desktop/spatial gif")
      dev.new(height=7,width=7)
	ani.record(reset=TRUE)
      ani.options(interval=.01)	
     }
#for(timeStep in burn:(burn+10) )
for(timeStep in burn:simTime) 
{ 
 tempBenefit<-array(dim=c(SpaceR,SpaceC,kmax))
 for(x in 1:length(FishSel))
  tempBenefit[,,x]<-SpaceNumAtLenT[timeStep,,,x]*FishSel[x]*wgtAtAge[x]*price

 #==in dollars (through selected kg of biomass), this is what a fisher would catch by going to a given patch
  BenefitPatch<-apply(tempBenefit,c(1,2),sum)
  BenefitPatch<-BenefitPatch*NoTakeZone
  NetBenefitPatch<-BenefitPatch-CostPatch
  if(timeStep==burn)OrigMaxNetBen<-max(NetBenefitPatch)
  tempSNALT<-SpaceNumAtLenT[timeStep,,,]

#==all fishers select their patch and fish
 if(any(timeStep%%12==season))
 {for(f in 1:Fishers)
  {     
  #==find closest, highest value, fishable patch=================
  #==insert option for choosing randomly in a bivariate normal centered on the port of entry
  #==insert method to close patches here=========================
  #==allow turfs here as well====================================
  #==open access, fishers increase until profit == 0 
  maxNetBene<-which(NetBenefitPatch==max(NetBenefitPatch),arr.ind=T)
  chosenPatch<-maxNetBene[which(Distance[maxNetBene]==min(Distance[maxNetBene])),]

   if(NetBenefitPatch[chosenPatch[1],chosenPatch[2]]>0)
   {
    #===FISHHHHH===#=><">====><"">====><""">===============================
    potentialCatch<-sum(SpaceNumAtLenT[timeStep,chosenPatch[1],chosenPatch[2],]*FishSel*wgtAtAge)
    if(potentialCatch<=maxCapac)
     {
       tempSNALT[chosenPatch[1],chosenPatch[2],]<-SpaceNumAtLenT[timeStep,chosenPatch[1],chosenPatch[2],] - 
      SpaceNumAtLenT[timeStep,chosenPatch[1],chosenPatch[2],]*FishSel
      CatchByFisher[f,timeStep-burn+1]<-potentialCatch
      CostByFisher[f,timeStep-burn+1]<-CostPatch[chosenPatch[1],chosenPatch[2]]
      print(timeStep)
     }
   #==if they can catch more than capactity, they'll probably change selectivity and throw small ones back....think about that
    if(potentialCatch>maxCapac)
    {
     #===find harvest rate that would put the fisher at max capacity===
     maxHarv<-1
     minHarv<-.0000001
     for(x in 1:25)
      {
       useHarv<-(maxHarv+minHarv)/2
       tempCat<-sum(SpaceNumAtLenT[timeStep,chosenPatch[1],chosenPatch[2],]*FishSel*useHarv*wgtAtAge)
       if(tempCat<maxCapac)
        minHarv<-useHarv
       if(tempCat>maxCapac)
        maxHarv<-useHarv
       }
    #==apply harvest===============================================
      tempSNALT[chosenPatch[1],chosenPatch[2],]<-SpaceNumAtLenT[timeStep,chosenPatch[1],chosenPatch[2],] - 
       SpaceNumAtLenT[timeStep,chosenPatch[1],chosenPatch[2],]*FishSel*useHarv
      CatchByFisher[f,timeStep-burn+1]<-sum(SpaceNumAtLenT[timeStep,chosenPatch[1],chosenPatch[2],]*FishSel*useHarv*wgtAtAge)
      CostByFisher[f,timeStep-burn+1]<-CostPatch[chosenPatch[1],chosenPatch[2]]
     }

    #==update benefitPatch for next fisher=======================================
    tempBenefit<-array(dim=c(SpaceR,SpaceC,kmax))
    for(x in 1:length(FishSel))
     tempBenefit[,,x]<-tempSNALT[,,x]*FishSel[x]*wgtAtAge[x]*price
    BenefitPatch<-apply(tempBenefit,c(1,2),sum)
    BenefitPatch<-BenefitPatch*NoTakeZone
    NetBenefitPatch<-BenefitPatch-CostPatch
    if(Graphs==T)
     {
      filled.contour(NetBenefitPatch,zlim=c(-OrigMaxNetBen,OrigMaxNetBen),y=seq(1,SpaceC),x=seq(1,SpaceR))
      mtext(side=3,paste("Timestep = ",timeStep))
      ani.record()
     }
    }
   }	#end fishers
  }
BenefitAtTimeStep[(timeStep-burn+1)]<-sum(BenefitPatch)
#==natural mortality (depends on timestep size)
  tempSNALT<-tempSNALT*exp(-(M/yearMark))

#==recruitment (happens once a year)
#==calculate recruits based on total spawnign biomass, allocate them according to either the spawning biomass in a patch or evenly.
#==current SpBio
tempSpBio<-rep(0,kmax)
for(i in 1:kmax)
 tempSpBio[i]<-sum(tempSNALT[,,i]*MatAge[i])
SpawningBiomass[timeStep]<-sum(tempSpBio)

   if(timeStep%%yearMark == 0)
    {
    #==move to next year class
     tempSNALT[,,2:(kmax-1)]<-tempSNALT[,,1:(kmax-2)]
     tempSNALT[,,kmax]<-tempSNALT[,,kmax]+tempSNALT[,,kmax-1]
    #==recruitment
     Recruits<-Recruitment(SpawningBiomass[timeStep],steepness,R0,VirSpBio,detRec,sigmaR)
     if(RecDist==1) # equally spread recruitment 
       tempSNALT[,,1]<-Recruits/(SpaceR*SpaceC)
     if(RecDist==2)
      {
       #==find proportion of total SpBio in each area
	 tempSpBioMat<-matrix(0,nrow=SpaceR,ncol=SpaceC)
	 for(i in 1:kmax)
	  tempSpBioMat<-tempSpBioMat+tempSNALT[,,i]*MatAge[i]
       RecDistMat<-tempSpBioMat/sum(tempSpBioMat)   
	 tempSNALT[,,1]<-Recruits*RecDistMat
	 }
     if(RecDist==3)
      {
       ##===INSERT habitat.csv DISTRIBUTION==
      } 
     }

#==movement depends on the probability of moving at age
temp1<-tempSNALT
for(age in 1:kmax)
  {
   tempAge<-array(dim=c(SpaceR,SpaceC,SpaceR*SpaceC))
   Movers<-MoveProb[age]*tempSNALT[,,age]
   NonMovers<-(1-MoveProb[age])*tempSNALT[,,age]

   for(h in 1:nrow(coords))
    tempAge[,,h]<-Movers[coords[h,1],coords[h,2]]*Movement[,,h]

   Moved<-apply(tempAge,c(1,2),sum)
   tempSNALT[,,age]<-Moved+NonMovers
   }

if(GraphsFish==T)
{

filled.contour(temp1[,,2],main="postfishery, premovement: age 2",zlim=c(0,9))
      mtext(side=3,paste("Timestep = ",timeStep))
      ani.record()
filled.contour(tempSNALT[,,2],main="postfishery, postmovement: age 2",zlim=c(0,9))
      mtext(side=3,paste("Timestep = ",timeStep))
      ani.record()
}
 #==update popdym
 SpaceNumAtLenT[timeStep+1,,,]<-tempSNALT

} # end timestep



#===========================================================================
#== Calculate the cost of management for a given scenario
#===========================================================================
# need size of MPA, whether a size limit is imposed, length of season,
# 

totCatch<-apply(CatchByFisher,2,sum,na.rm=T)
totCost<-apply(CostByFisher,2,sum,na.rm=T)
dev.new()
par(mfrow=c(4,1),mar=c(.1,4,.1,.1))
plot(totCatch,type="b",xaxt='n',las=2)
plot(totCost,lty=2,type="b",xaxt='n',las=2)
plot(SpawningBiomass[burn:simTime],type="b",xaxt='n',las=2,ylab="SpawningBio",ylim=c(0,max(SpawningBiomass[burn:simTime],na.rm=T)))
#plot(BenefitAtTimeStep,type='b',xaxt='n',las=2,ylim=c(0,max(BenefitAtTimeStep)))
#plot(CostofManagement,type='b',las=2,ylim=c(0,max(CostofManagement)))
#plot(Profit,type='b',las=2,ylim=c(0,max(Profit)))

#ani.options(interval=.15)	
#ani.replay()
#saveHTML(ani.replay(), img.name = "record_plot_old",outdir = getwd(),interval=0.05)

	# size limit
	# 
	# identify available patches to fish in (
	# fish closest to port until density makes it unprofitable
      
	# how will that interplay with movement...
      
#==fishign mortality


