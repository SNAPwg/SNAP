require(animation)
require(caTools)
setwd("C:/SNAP")
source("lenwei.R")
source("VisualizeMovement.R")
source("movArray.R")
source("InitialPop.R")
source("Recruitment.R")
#rm(list=ls())

Graphs<-F
GraphsFish<-F
PrintLifeHistory<-F

Life<-read.csv("LifeHistory.csv")
SimCTL<-read.csv("GrandSimCtl.csv")
Fleets<-read.csv("Fleets.csv")

#==life history ==============================
kmax	 	<-Life[grep('kmax',Life[,ncol(Life)]),seq(1,ncol(Life)-1)]              # maximum age
kmat	 	<-Life[grep('kmat',Life[,ncol(Life)]),seq(1,ncol(Life)-1)]           	  # age at maturity
kmin	 	<-Life[grep('kmin',Life[,ncol(Life)]),seq(1,ncol(Life)-1)]	            # minimium size
M   	 	<-Life[grep('natural mortality',Life[,ncol(Life)]),seq(1,ncol(Life)-1)]	# natural mortality
Linf	 	<-Life[grep('Linf',Life[,ncol(Life)]),seq(1,ncol(Life)-1)]            	# VonBert Linf
K		    <-Life[grep('VonBert K',Life[,ncol(Life)]),seq(1,ncol(Life)-1)]       	# VonBert K
t0  	 	<-Life[grep('t0',Life[,ncol(Life)]),seq(1,ncol(Life)-1)]         		  	# VonBert t0
wtA	  	<-Life[grep('wtA',Life[,ncol(Life)]),seq(1,ncol(Life)-1)]	              #  weight A
wtB	  	<-Life[grep('wtB',Life[,ncol(Life)]),seq(1,ncol(Life)-1)]	              # weight B
mat50		<-Life[grep('mat50',Life[,ncol(Life)]),seq(1,ncol(Life)-1)]             # age at 50% maturity
mat95		<-Life[grep('mat95',Life[,ncol(Life)]),seq(1,ncol(Life)-1)]             # age at 95% maturity 
lenSD   <-Life[grep('lenSD',Life[,ncol(Life)]),seq(1,ncol(Life)-1)]             # standard deviation of length at age

#==recruitment================================
detRec  <-Life[grep('detRec',Life[,ncol(Life)]),seq(1,ncol(Life)-1)]        	    # deterministic recruitment?
R0	  	<-Life[grep('R0',Life[,ncol(Life)]),seq(1,ncol(Life)-1)]             		# Virgin recruitment
steepness<-Life[grep('steepness',Life[,ncol(Life)]),seq(1,ncol(Life)-1)]      		# steepness
sigmaR	<-Life[grep('sigmaR',Life[,ncol(Life)]),seq(1,ncol(Life)-1)]            	# variability around recruitment curve
RecDist <-Life[grep('RecDist',Life[,ncol(Life)]),seq(1,ncol(Life)-1)]            # How is recruitment distributed over space? (1=equally,2=determined by adults,3=determined by habitat.csv)
sdx		  <-Life[grep('sdx',Life[,ncol(Life)]),seq(1,ncol(Life)-1)]              	# sd of movement on x axis (in units of rows/columns)
sdy	  	<-Life[grep('sdy',Life[,ncol(Life)]),seq(1,ncol(Life)-1)]              	# sd of movement on y axis
movP50	<-Life[grep('movP50',Life[,ncol(Life)]),seq(1,ncol(Life)-1)]             # probability moving by length (logistic; 50%; set both to 0 if all move)
movP95	<-Life[grep('movP95',Life[,ncol(Life)]),seq(1,ncol(Life)-1)]             # probability of moving at length (95%)

MoveProb<-1/(1+exp(-log(19)*((seq(1,kmax)-movP50)/(movP95-movP50))))
MatAge<-1/(1+exp(-log(19)*((seq(1,kmax)-mat50)/(mat95-mat50))))

lenAtAge<-rep(0,kmax)
wgtAtAge<-rep(0,kmax)
 for(k in 1:kmax)
  {
   lenAtAge[k]<-Linf*(1-exp(-K*k-t0))
   wgtAtAge[k]<-wtA*lenAtAge[k]^wtB
  }

#==cross simulation ctl=========================
SpaceR	<-SimCTL[grep('SpaceR',SimCTL[,2]),1]	    	# Rows in the grid space
SpaceC	<-SimCTL[grep('SpaceC',SimCTL[,2]),1]  			# cols in the grid spcae
burn		<-SimCTL[grep('burn',SimCTL[,2]),1]    		  # burn in for movement equilibration (use option "Graphs" in "Initpop" to see that population is thoroughly mixed)
simTime <-SimCTL[grep('simTime',SimCTL[,2]),1]     	# time steps for projection
yearMark<-SimCTL[grep('yearMark',SimCTL[,2]),1]	  	# number of time steps in a year

#==management costs==========================
MPAsunk     <-SimCTL[grep('MPAsunk',SimCTL[,2]),1]  	# start up cost of enforcing an MPA
MPAcost	    <-SimCTL[grep('MPAcost',SimCTL[,2]),1]		# cost of maintaining a unit of MPA per unit time
SizeSunk	  <-SimCTL[grep('SizeSunk',SimCTL[,2]),1]		# start up cost of enforcing a size limit
SizeCost    <-SimCTL[grep('SizeCost',SimCTL[,2]),1]		# cost of enforcing a size limit per unit of time per port
SeasonSunk	<-SimCTL[grep('SeasonSunk',SimCTL[,2]),1]		# start up cost of enforcing a season
SeasonCost	<-SimCTL[grep('SeasonCost',SimCTL[,2]),1]		# cost of enforcing a season per unit of time

#==management auxilliary benefits======================
MPAtourism	<-SimCTL[grep('MPAtourism',SimCTL[,2]),1]		# average revenue generated per unit MPA per unit time for tourism in an MPA

#==fishery characteristics===================
SizeLimit	<-Fleets[grep('SizeLimit',Fleets[,ncol(Fleets)]),seq(1,ncol(Fleets)-1)] 
Sel50		  <-Fleets[grep('Sel50',Fleets[,ncol(Fleets)]),seq(1,ncol(Fleets)-1)] 		# selectivity by length (9 by len)
Sel95		  <-Fleets[grep('Sel95',Fleets[,ncol(Fleets)]),seq(1,ncol(Fleets)-1)]		# selectividad by len (13 by len)
q		      <-Fleets[grep('catchability',Fleets[,ncol(Fleets)]),seq(1,ncol(Fleets)-1)]  		# fishery catchability coefficient (what fraction of the exploitable biomass in a patch is catchable)
PortX 	  <-Fleets[grep('PortX',Fleets[,ncol(Fleets)]),seq(1,ncol(Fleets)-1)] 	# location of the port of entry (x,y)
PortY   	<-Fleets[grep('PortY',Fleets[,ncol(Fleets)]),seq(1,ncol(Fleets)-1)]   # location of the port of entry (x,y)
seasonN   <-Fleets[grep('season',Fleets[,ncol(Fleets)]),seq(1,ncol(Fleets)-1)]  	# months in which fishing is allowed

season<-seq(1,seasonN)
PortLc<-c(PortX,PortY)
FishSel<-q/(1+exp(-log(19)*((seq(1,kmax)-Sel50)/(Sel95-Sel50))))

#==econ pars=================================
price		  <-Fleets[grep('price',Fleets[,ncol(Fleets)]),seq(1,ncol(Fleets)-1)]		# ex-vessel price per unit harvest (kg)
costTrv	  <-Fleets[grep('costTrv',Fleets[,ncol(Fleets)]),seq(1,ncol(Fleets)-1)]		# cost to travel one patch
costFish	<-Fleets[grep('costFish',Fleets[,ncol(Fleets)]),seq(1,ncol(Fleets)-1)]		# cost per unit effort fishing
discRate	<-Fleets[grep('discRate',Fleets[,ncol(Fleets)]),seq(1,ncol(Fleets)-1)]	# discount rate
TimeHor	  <-Fleets[grep('TimeHor',Fleets[,ncol(Fleets)]),seq(1,ncol(Fleets)-1)]		# time horizon for evaluation
Fishers	  <-Fleets[grep('Fishers',Fleets[,ncol(Fleets)]),seq(1,ncol(Fleets)-1)]		# number of fishers
maxCapac	<-Fleets[grep('maxCapac',Fleets[,ncol(Fleets)]),seq(1,ncol(Fleets)-1)]		# max capacity of a single fisherman in a timestep (kg) (can also be a 'trip limit'/permit implementation)

#==spatial matrix (0=open access, 1=TURF, 2=NTZ)
#==read in a csv with the appropriate letters denoted for the different uses
SpaceUse<-matrix(0,ncol=SpaceC,nrow=SpaceR)

#==read in csv that give habitat quality, determine movement by habitat quality
#==read csv that shows MPAs, if any
#NoTakeZone<-read.csv("notakezone40.csv",header=F)
#NoTakeZone<-read.csv("notakezone80.csv",header=F)
NoTakeZone<-read.csv("notakezoneNULL.csv",header=F)

#==add in dispersal .csv too
habitat<-read.csv("habitatNULL.csv",header=F)
#is.numeric(habitat)
#habitat<-read.csv("habitatNULL.csv",header=F)


if(PrintLifeHistory==T)
{
 #==print out maturity, selectivity, growth, etc.
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
 
 filled.contour(matrix(as.numeric(unlist(habitat)),ncol=SpaceC),y=seq(1,SpaceR),x=seq(1,SpaceC)) 
 mtext(side=3,"Habitat quality")

 filled.contour(matrix(as.numeric(unlist(NoTakeZone)),ncol=SpaceC),y=seq(1,SpaceR),x=seq(1,SpaceC)) 
 mtext(side=3,"No take zones (0)")
 #==identify the MPAs and mark them on the graph

}

#===========================================================
#== Begin simulations====================================
#===========================================================

#==population dynamics array tracking number at length (or age) in a patch by year
SpaceNumAtAgeT<-array(dim=c((simTime),SpaceR,SpaceC,kmax))
CatchByFisher<-matrix(ncol=simTime-burn,nrow=Fishers)
ProfitByFisher<-matrix(ncol=simTime-burn,nrow=Fishers)
CostByFisher<-matrix(ncol=simTime-burn,nrow=Fishers)
SpawningBiomass<-rep(0,simTime-burn)
CostOfManagement<-rep(0,simTime-burn)

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
SpaceNumAtAgeT[1:burn,,,]<-InitialPop(R0,M,kmax,SpaceR,SpaceC,MoveProb,coords,Movement,Graphs=F,burn)
SpaceNumAtAgeT[burn,,,]

#==calculate virgin spawning biomass (1843.584)================
tempSpBio<-rep(0,kmax)
for(i in 1:kmax)
 tempSpBio[i]<-sum(SpaceNumAtAgeT[burn,,,i]*MatAge[i])
VirSpBio<-sum(tempSpBio)

#==begin projected harvest=================== 
#==fishers decide where to fish============
#==initial costs
 Distance<-sqrt((row(SpaceUse)-PortX)^2 + (col(SpaceUse)-PortY)^2)
 CostPatch<-Distance*costTrv+costFish
#==relate the cost of fishing to the number of fish in a patch??

    if(Graphs==T | GraphsFish==T)
     {
      #setwd("C:/Users/Cody/Desktop/spatial gif")
      #dev.new(height=7,width=7)
	    ani.record(reset=TRUE)
      ani.options(interval=.01)	
     }
for(timeStep in burn:simTime) 
{ 
 tempBenefit<-array(dim=c(SpaceR,SpaceC,kmax))
 for(x in 1:length(FishSel))
  tempBenefit[,,x]<-SpaceNumAtAgeT[timeStep,,,x]*FishSel[x]*wgtAtAge[x]*price

 #==in dollars (through selected kg of biomass), this is what a fisher would catch by going to a given patch
  BenefitPatch<-apply(tempBenefit,c(1,2),sum)
  BenefitPatch<-BenefitPatch*NoTakeZone
  NetBenefitPatch<-BenefitPatch-CostPatch
  if(timeStep==burn)OrigMaxNetBen<-max(NetBenefitPatch)
  tempSNALT<-SpaceNumAtAgeT[timeStep,,,]

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
    potentialCatch<-sum(SpaceNumAtAgeT[timeStep,chosenPatch[1],chosenPatch[2],]*FishSel*wgtAtAge)
    if(potentialCatch<=maxCapac)
     {
       tempSNALT[chosenPatch[1],chosenPatch[2],]<-SpaceNumAtAgeT[timeStep,chosenPatch[1],chosenPatch[2],] - 
      SpaceNumAtAgeT[timeStep,chosenPatch[1],chosenPatch[2],]*FishSel
      CatchByFisher[f,timeStep-burn+1]<-potentialCatch
      CostByFisher[f,timeStep-burn+1]<-CostPatch[chosenPatch[1],chosenPatch[2]]
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
       tempCat<-sum(SpaceNumAtAgeT[timeStep,chosenPatch[1],chosenPatch[2],]*FishSel*useHarv*wgtAtAge)
       if(tempCat<maxCapac)
        minHarv<-useHarv
       if(tempCat>maxCapac)
        maxHarv<-useHarv
       }
    #==apply harvest===============================================
      tempSNALT[chosenPatch[1],chosenPatch[2],]<-SpaceNumAtAgeT[timeStep,chosenPatch[1],chosenPatch[2],] - 
       SpaceNumAtAgeT[timeStep,chosenPatch[1],chosenPatch[2],]*FishSel*useHarv
      CatchByFisher[f,timeStep-burn+1]<-sum(SpaceNumAtAgeT[timeStep,chosenPatch[1],chosenPatch[2],]*FishSel*useHarv*wgtAtAge)
      CostByFisher[f,timeStep-burn+1]<-CostPatch[chosenPatch[1],chosenPatch[2]]
      ProfitByFisher[f,timeStep-burn+1]<-CatchByFisher[f,timeStep-burn+1]*price - CostByFisher[f,timeStep-burn+1]
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
 SpaceNumAtAgeT[timeStep+1,,,]<-tempSNALT

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
} # end timestep



#===========================================================================
#== Calculate the cost of management for a given scenario
#===========================================================================
# need size of MPA, whether a size limit is imposed, length of season,
# 

totCatch<-apply(CatchByFisher,2,sum,na.rm=T)
totCost<-apply(CostByFisher,2,sum,na.rm=T)
totProfit<-apply(ProfitByFisher,2,sum,na.rm=T)

par(mfrow=c(4,1),mar=c(.1,4,.1,.1))
plot(totCatch,type="b",xaxt='n',las=2)
plot(totCost,lty=2,type="b",xaxt='n',las=2)
plot(totProfit,lty=2,type="b",xaxt='n',las=2)
lines(CostOfManagement,lty=2,type="b",col=2)
plot(SpawningBiomass[burn:simTime],type="b",xaxt='n',las=2,ylab="SpawningBio",ylim=c(0,max(SpawningBiomass[burn:simTime],na.rm=T)))

#plot(CostofManagement,type='b',las=2,ylim=c(0,max(CostofManagement)))
#plot(Profit,type='b',las=2,ylim=c(0,max(Profit)))

#ani.options(interval=.15)	
#ani.replay()
#saveHTML(ani.replay(), img.name = "record_plot_oldf",outdir = getwd(),interval=0.05)

	# size limit
	# 
	# identify available patches to fish in (
	# fish closest to port until density makes it unprofitable
      
	# how will that interplay with movement...
      
#==fishign mortality


