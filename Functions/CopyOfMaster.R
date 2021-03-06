#==this function runs the main population and economic dynamics for simulations given several inputs (usually in the form of csvs)

Master<-function(Life,SimCTL,Fleets,season,Samp,NoTakeZoneNULL,NoTakeZoneImp,habitat,Graphs=F,GraphsFish=F,PrintLifeHistory=F)
{
FleetN<-ncol(Fleets)-1

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
   lenAtAge[k]<-Linf*(1-exp(-K*(k-t0)))
   wgtAtAge[k]<-wtA*lenAtAge[k]^wtB
  }

#==cross simulation ctl=========================
SpaceR	  <-SimCTL[grep('SpaceR',SimCTL[,2]),1]	    	# Rows in the grid space
SpaceC	  <-SimCTL[grep('SpaceC',SimCTL[,2]),1]  			# cols in the grid spcae
burn	  	<-SimCTL[grep('burn',SimCTL[,2]),1]    		  # burn in for movement equilibration (use option "Graphs" in "Initpop" to see that population is thoroughly mixed)
simTime   <-SimCTL[grep('simTime',SimCTL[,2]),1]     	# time steps for projection
yearMark  <-SimCTL[grep('yearMark',SimCTL[,2]),1]	  	# number of time steps in a year
initManage<-SimCTL[grep('initManage',SimCTL[,2]),1]   # year in which to initiate management


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
SizeLimit	<-as.numeric(Fleets[grep('SizeLimit',Fleets[,ncol(Fleets)]),seq(1,ncol(Fleets)-1)] )
Sel50		  <-as.numeric(Fleets[grep('Sel50',Fleets[,ncol(Fleets)]),seq(1,ncol(Fleets)-1)] )		# selectivity by length (9 by len)
Sel95		  <-as.numeric(Fleets[grep('Sel95',Fleets[,ncol(Fleets)]),seq(1,ncol(Fleets)-1)]	)	# selectividad by len (13 by len)
q		      <-as.numeric(Fleets[grep('catchability',Fleets[,ncol(Fleets)]),seq(1,ncol(Fleets)-1)])  		# fishery catchability coefficient (what fraction of the exploitable biomass in a patch is catchable)
PortX 	  <-as.numeric(Fleets[grep('PortX',Fleets[,ncol(Fleets)]),seq(1,ncol(Fleets)-1)]) 	# location of the port of entry (x,y)
PortY   	<-as.numeric(Fleets[grep('PortY',Fleets[,ncol(Fleets)]),seq(1,ncol(Fleets)-1)] )  # location of the port of entry (x,y)
costSteep <-as.numeric(Fleets[grep('costSteep',Fleets[,ncol(Fleets)]),seq(1,ncol(Fleets)-1)] )  # determines the influence of catch rate on fisher behavior

#==Sampling Params ==============================
DataParams <- NULL
DataParams$histCatchFD    <-Samp[grep('histCatchFD',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]              # Collect historical Catch Data?
DataParams$histEffortFD  	<-Samp[grep('histEffortFD',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]           	  # Collect historical CPUE Data?
DataParams$histStartYr	 	<-Samp[grep('histStartYr',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]	            # time step in which historical data was first collected
DataParams$Aggregate   	 	<-Samp[grep('Aggregate',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]	# Aggregate data over space? 1 for yes, 2 for no.
DataParams$SampStartYr	 	<-Samp[grep('SampStartYr',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]            	# time step in which current sampling program begins
DataParams$SampFreq		    <-Samp[grep('SampFreq',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]       	# Sampling Frequency
DataParams$sigHistCatch  	<-Samp[grep('sigHistCatch',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]         		  	# Standard deviation of observation error on total catches. If no error, sigHistCatch=0.
DataParams$catchBias	  	<-Samp[grep('catchBias',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]	              # Bias in catch reporting. Scalar between (-1, 1). 0 equals no bias.
DataParams$Sel50FI       	<-Samp[grep('Sel50FI',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]	      # Survey selectivity at age
DataParams$Sel95FI	     	<-Samp[grep('Sel95FI',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]            	# Survey selectivity at age
DataParams$Sel50VB		    <-Samp[grep('Sel50VB',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]       	#  Selectivity for von Bertlanffy surveys (usually uniform above a specific age)
DataParams$Sel95VB    	 	<-Samp[grep('Sel95VB',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]         	# Selectivity for von Bertlanffy surveys (usually uniform above a specific age)
DataParams$SurveyF	    	<-Samp[grep('SurveyF',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]	          # Fishing rate associated with survey effort
DataParams$nAgeFD         <-Samp[grep('nAgeFD',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]   # Number of fishery dependent age samples collected from each patch
DataParams$nAgeFI         <-Samp[grep('nAgeFI',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]              # Number of fishery independent age samples collected from each patch
DataParams$nSizeFD      	<-Samp[grep('nSizeFD',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]           	  # Number of fishery dependent size samples collected from each patch
DataParams$nSizeFI      	<-Samp[grep('nSizeFI',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]	            # Number of fishery independent size samples collected from each patch
DataParams$numVB      	 	<-Samp[grep('numVB',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]	# Number of VB age-size samples collected from each patch
DataParams$VBbins	      	<-Samp[grep('VBbins',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]            	#  Size bins associated with von Bertalanffy survey
DataParams$lengthObsSD_FI <-Samp[grep('lengthObsSD_FI',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]              # Observation Error on fishery independent lengths
DataParams$lengthObsSD_FD <-Samp[grep('lengthObsSD_FD',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]           	  # Observation Error on fishery dependent lengths
DataParams$AgeObsSD_FI	 	<-Samp[grep('AgeObsSD_FI',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]	            # Observation Error on fishery independent ages
DataParams$AgeObsSD_FD   	<-Samp[grep('AgeObsSD_FD',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]	# Observation Error on fishery dependent ages
DataParams$FISurvType   	<-Samp[grep('FISurvType',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]            	# Spatial pattern for fishery independent surveys (0 = input a .csv file specifying custom survey distribution, 1 = sample all patches, 2=sample no patches, 3 = sample random 20% of patches, 4= sample near and far site, 5=sample MPAs)
DataParams$FDSurvType		  <-Samp[grep('FDSurvType',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]       	# Spatial pattern for fishery dependent surveys (0 = input a .csv file specifying custom survey distribution, 1 = sample all fished patches, 2=sample no patches, 3 = sample random 20% of patches, 4= sample near and far site)
DataParams$Catch_FD       <-Samp[grep('Catch_FD',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]	            # Collect fishery dependent catch data? 1 = yes, 2 = no.
DataParams$Catch_FI   	 	<-Samp[grep('Catch_FI',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]	# Collect fishery independent catch data? 1 = yes, 2 = no.
DataParams$Effort_FI	  	<-Samp[grep('Effort_FI',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]            	# Collect fishery independent abundance index data? 1 = yes, 2 = no.
DataParams$Effort_FD      <-Samp[grep('Effort_FD',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]              # Collect fishery dependent CPUE data? 1 = yes, 2 = no.
DataParams$Ages_FD   	    <-Samp[grep('Ages_FD',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]           	  # Collect fishery dependent age composition data? 1 = yes, 2 = no.
DataParams$Ages_FI	    	<-Samp[grep('Ages_FI',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]	            # Collect fishery independent age composition data? 1 = yes, 2 = no.
DataParams$Sizes_FD   	 	<-Samp[grep('Sizes_FD',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]	# Collect fishery dependent size composition data? 1 = yes, 2 = no.
DataParams$Sizes_FI     	<-Samp[grep('Sizes_FI',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]            	# Collect fishery independent size composition data? 1 = yes, 2 = no.
DataParams$VBSurvey		    <-Samp[grep('VBSurvey',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]       	# Collect age and size data to estimate VB params? 1 = yes, 2 = no.
DataParams$Survey_q       <-Samp[grep('Survey_q',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]            	# Survey catchability. Might be different depending on mode of survey (such as dive surveys)
DataParams$sigSurvey		    <-Samp[grep('sigSurvey',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]       	# Standard deviation of observation error on total catches. If no error, sigHistCatch=0.

sampleTimeSteps <- seq(DataParams$SampStartYr,simTime-burn+1,by=DataParams$SampFreq)  #Time steps in which sampling occurs

## Set up data storage. Size is different when data is Aggregated.

if(DataParams$Aggregate == 0){
  FDCatchData <- array(dim=c(simTime-burn+1, SpaceR, SpaceC, FleetN))
  FDCPUEData <- array(dim=c(simTime-burn+1, SpaceR, SpaceC, FleetN))
  FDAgeData <- array(NA,dim=c(simTime-burn+1, DataParams$nAgeFD,SpaceR, SpaceC, FleetN))
  FDSizeData <- array(NA,dim=c(simTime-burn+1, DataParams$nSizeFD,SpaceR, SpaceC, FleetN))
  FICatchData <- array(dim=c(simTime-burn+1, SpaceR, SpaceC))
  FICPUEData <- array(dim=c(simTime-burn+1, SpaceR, SpaceC))
} else {
  FDCatchData <- array(dim=c(simTime-burn+1, FleetN))
  FDCPUEData <- array(dim=c(simTime-burn+1, FleetN))
  FDAgeData <- array(NA,dim=c(simTime-burn+1, DataParams$nAgeFD*SpaceR*SpaceC, FleetN))
  FDSizeData <- array(NA,dim=c(simTime-burn+1, DataParams$nSizeFD*SpaceR*SpaceC, FleetN))
  FICatchData <- array(dim=c(simTime-burn+1))
  FICPUEData <- array(dim=c(simTime-burn+1))
}

PortLc<-c(PortX,PortY) #Does this get used anywhere?
FishSel<-matrix(ncol=FleetN,nrow=kmax)
for(y in 1:FleetN)
 FishSel[,y]<-q[y]/(1+exp(-log(19)*((seq(1,kmax)-Sel50[y])/(Sel95[y]-Sel50[y]))))

#==econ pars=================================
price		  <-as.numeric(Fleets[grep('price',Fleets[,ncol(Fleets)]),seq(1,ncol(Fleets)-1)])		# ex-vessel price per unit harvest (kg)
costTrv	  <-as.numeric(Fleets[grep('costTrv',Fleets[,ncol(Fleets)]),seq(1,ncol(Fleets)-1)])		# cost to travel one patch
costFish	<-as.numeric(Fleets[grep('costFish',Fleets[,ncol(Fleets)]),seq(1,ncol(Fleets)-1)])		# cost per unit effort fishing
discRate	<-as.numeric(Fleets[grep('discRate',Fleets[,ncol(Fleets)]),seq(1,ncol(Fleets)-1)])	# discount rate
TimeHor	  <-as.numeric(Fleets[grep('TimeHor',Fleets[,ncol(Fleets)]),seq(1,ncol(Fleets)-1)])		# time horizon for evaluation
Fishers	  <-as.numeric(Fleets[grep('Fishers',Fleets[,ncol(Fleets)]),seq(1,ncol(Fleets)-1)])		# number of fishers
maxCapac	<-as.numeric(Fleets[grep('maxCapac',Fleets[,ncol(Fleets)]),seq(1,ncol(Fleets)-1)])		# max capacity of a single fisherman in a timestep (kg) (can also be a 'trip limit'/permit implementation)

#==dummy spatial matrix
SpaceUse<-matrix(0,ncol=SpaceC,nrow=SpaceR)

if(PrintLifeHistory==T)
{
 #==print out maturity, selectivity, growth, etc.
quartz()
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
 

quartz()
 filled.contour(matrix(as.numeric(unlist(habitat)),ncol=SpaceC),x=seq(1,SpaceR),y=seq(1,SpaceC)) 
 mtext(side=3,"Habitat quality")

quartz()
 filled.contour(matrix(as.numeric(unlist(NoTakeZone)),ncol=SpaceC),x=seq(1,SpaceR),y=seq(1,SpaceC)) 
 mtext(side=3,"No take zones (0)")
}

#===========================================================
#== Begin simulations====================================
#===========================================================

#==population dynamics array tracking number at age in a patch by year
SpaceNumAtAgeT<-array(dim=c((simTime),SpaceR,SpaceC,kmax))
SpaceCatAgeByFisher<-array(dim=c(simTime-burn+1,SpaceR,SpaceC,kmax,max(Fishers),FleetN)) #in numbers
SpaceEffort<-array(0,dim=c(simTime-burn+1,SpaceR,SpaceC,FleetN))

CatchByFisher<-array(dim=c(max(Fishers),(simTime-burn+1),FleetN))
ProfitByFisher<-array(dim=c(max(Fishers),(simTime-burn+1),FleetN))
CostByFisher<-array(dim=c(max(Fishers),simTime-burn+1,FleetN))
SpawningBiomass<-rep(0,simTime-burn)
NumbersExplt<-rep(0,simTime-burn)
CostOfManagement<-rep(0,simTime-burn)
InsideMPAspbio<-rep(0,simTime-burn)
OutsideMPAspbio<-rep(0,simTime-burn)

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

#==calculate virgin spawning biomass (1843.584)================
tempSpBio<-rep(0,kmax)
for(i in 1:kmax)
 tempSpBio[i]<-sum(SpaceNumAtAgeT[burn,,,i]*MatAge[i]*wgtAtAge[i])
VirSpBio<-sum(tempSpBio)


#==relate the cost of fishing to the number of fish in a patch?? (e.g. density is low == must fish harder)

    if(Graphs==T | GraphsFish==T)
     {
      #setwd("C:/Users/Cody/Desktop/spatial gif")
      #dev.new(height=7,width=7)
# 	    ani.record(reset=TRUE)
#       ani.options(interval=.01)	
     }

for(timeStep in burn:simTime) 
    { 
  NoTakeZone<-NoTakeZoneNULL
  if(timeStep>=initManage)
    NoTakeZone<-NoTakeZoneImp
  
  if(sum(dim(NoTakeZone)- c(SpaceR,SpaceC))!=0)
  {
    print("no take zone is of the wrong dimensions")
    break
  }
  
 for(flt in 1:FleetN)
 {
   #==fishers decide where to fish============
   #==initial costs
   Distance<-sqrt((row(SpaceUse)-PortX[flt])^2 + (col(SpaceUse)-PortY[flt])^2)
 
 tempBenefit<-array(dim=c(SpaceR,SpaceC,kmax))
 for(x in 1:nrow(FishSel))
  tempBenefit[,,x]<-SpaceNumAtAgeT[timeStep,,,x]*FishSel[x,flt]*wgtAtAge[x]
 BenefitPatch<-apply(tempBenefit,c(1,2),sum) 

 if(timeStep==burn)OrigMaxBen<-max(BenefitPatch)
 
 #==costs depend both on the distance from port and the density of fish
  #tempCostFish<-costFish[flt]*(OrigMaxBen/BenefitPatch)^costSteep
  tempCostFish<-(1+(costSteep-(costSteep)*(BenefitPatch/OrigMaxBen)))*costFish[flt]
  CostPatch<-Distance*costTrv[flt]+tempCostFish
 #==cap potential benefit at what can be removed after calculating cost of fishing
  BenefitPatch[BenefitPatch>maxCapac[flt]]<-maxCapac[flt]
 #==translate to dollars
  BenefitPatch<-BenefitPatch*price[flt]
 #==in dollars (through selected kg of biomass), this is what a fisher would catch by going to a given patch
  BenefitPatch<-BenefitPatch*NoTakeZone
  NetBenefitPatch<-BenefitPatch-CostPatch
  if(timeStep==burn)OrigMaxNetBen<-max(NetBenefitPatch)
 
 tempSNALT<-SpaceNumAtAgeT[timeStep,,,]
 if(timeStep==burn)OrigTempSNALT<-max(tempSNALT[,,3])
#==all fishers select their patch and fish if the season is 'in'

 if(any(timeStep%%12==(which(!is.na(season[,flt+1]))-1)))
 {
  for(f in 1:Fishers[flt])
  {     
  #==find closest, highest value, fishable patch=================
  #==insert option for choosing randomly in a bivariate normal centered on the port of entry
  #==insert method to close patches here=========================
  #==open access, fishers increase until profit == 0 ??
  maxNetBene<-which(NetBenefitPatch==max(NetBenefitPatch),arr.ind=T)
  chosenPatch<-maxNetBene[which(Distance[maxNetBene]==min(Distance[maxNetBene])),]

   if(NetBenefitPatch[chosenPatch[1],chosenPatch[2]]>0)
   {
    #===FISHHHHH===#=><">====><"">====><""">===============================
    potentialCatch<-sum(SpaceNumAtAgeT[timeStep,chosenPatch[1],chosenPatch[2],]*FishSel[,flt]*wgtAtAge)
    if(potentialCatch<=maxCapac[flt])
     {
       tempSNALT[chosenPatch[1],chosenPatch[2],]<-SpaceNumAtAgeT[timeStep,chosenPatch[1],chosenPatch[2],] - 
      SpaceNumAtAgeT[timeStep,chosenPatch[1],chosenPatch[2],]*FishSel[,flt]
      #==spatial catch at age by patch by fisher by fleet
      for(p in 1:kmax)
       SpaceCatAgeByFisher[timeStep-burn+1,chosenPatch[1],chosenPatch[2],p,f,flt]<-SpaceNumAtAgeT[timeStep,chosenPatch[1],chosenPatch[2],p]*FishSel[p,flt]
      CatchByFisher[f,timeStep-burn+1,flt]<-potentialCatch
      CostByFisher[f,timeStep-burn+1,flt]<-CostPatch[chosenPatch[1],chosenPatch[2]]
      ProfitByFisher[f,timeStep-burn+1,flt]<-CatchByFisher[f,timeStep-burn+1,flt]*price[flt] - CostByFisher[f,timeStep-burn+1,flt]
     }
   #==if they can catch more than capactity, they'll probably change selectivity and throw small ones back....think about that
    if(potentialCatch>maxCapac[flt])
    {
     #===find harvest rate that would put the fisher at max capacity===
     maxHarv<-1
     minHarv<-.000000001
     for(x in 1:25)
      {
       useHarv<-(maxHarv+minHarv)/2
       tempCat<-sum(SpaceNumAtAgeT[timeStep,chosenPatch[1],chosenPatch[2],]*FishSel[,flt]*useHarv*wgtAtAge)
       if(tempCat<maxCapac[flt])
        minHarv<-useHarv
       if(tempCat>maxCapac[flt])
        maxHarv<-useHarv
       }
    #==apply harvest===============================================
      tempSNALT[chosenPatch[1],chosenPatch[2],]<-SpaceNumAtAgeT[timeStep,chosenPatch[1],chosenPatch[2],] - 
       SpaceNumAtAgeT[timeStep,chosenPatch[1],chosenPatch[2],]*FishSel[,flt]*useHarv
      CatchByFisher[f,timeStep-burn+1,flt]<-sum(SpaceNumAtAgeT[timeStep,chosenPatch[1],chosenPatch[2],]*FishSel[,flt]*useHarv*wgtAtAge)
      for(p in 1:kmax)
       SpaceCatAgeByFisher[timeStep-burn+1,chosenPatch[1],chosenPatch[2],p,f,flt]<-SpaceNumAtAgeT[timeStep,chosenPatch[1],chosenPatch[2],p]*FishSel[p,flt]*useHarv

      CostByFisher[f,timeStep-burn+1,flt]<-CostPatch[chosenPatch[1],chosenPatch[2]]
      ProfitByFisher[f,timeStep-burn+1,flt]<-CatchByFisher[f,timeStep-burn+1,flt]*price[flt] - CostByFisher[f,timeStep-burn+1,flt]

     }

    #==tracks number of fishermen fishing in a patch at time
    SpaceEffort[timeStep-burn+1,chosenPatch[1],chosenPatch[2],flt]<-SpaceEffort[timeStep-burn+1,chosenPatch[1],chosenPatch[2],flt]+1

    #==update benefitPatch for next fisher=======================================
   tempBenefit<-array(dim=c(SpaceR,SpaceC,kmax))
   for(x in 1:nrow(FishSel))
     tempBenefit[,,x]<-tempSNALT[,,x]*FishSel[x,flt]*wgtAtAge[x]
   BenefitPatch<-apply(tempBenefit,c(1,2),sum) 

#    tempCostFish<-costFish[flt]*(OrigMaxBen/BenefitPatch)^costSteep
   tempCostFish<-(1+(costSteep-(costSteep)*(BenefitPatch/OrigMaxBen)))*costFish[flt]
   CostPatch<-Distance*costTrv[flt]+tempCostFish
   BenefitPatch[BenefitPatch>maxCapac[flt]]<-maxCapac[flt]   
   BenefitPatch<-BenefitPatch*price[flt]
   BenefitPatch<-BenefitPatch*NoTakeZone   
   NetBenefitPatch<-BenefitPatch-CostPatch
    if(Graphs==T)
     {
      filled.contour(as.matrix(NetBenefitPatch),y=seq(1,SpaceC),x=seq(1,SpaceR))
       #,zlim=c(-OrigMaxNetBen,OrigMaxNetBen)
      mtext(side=3,paste("Timestep = ",timeStep))
     # ani.record()
     }

    }
   }	#end fishers
  
  #Sample at the end of fishing if timeStep equals a sampling time step.
#   if(timeStep %in% (sampleTimeSteps+burn)){
#     Data <- CollectData(timeStep, simTime, burn, FleetN, DataParams, NoTakeZone, SpaceNumAtAgeT,
#                         SpaceCatAgeByFisher, SpaceEffort, wgtAtAge, lenAtAge, lenSD, PortX, PortY,kmax,Linf)
#     
#     # Assign data to storage
#     if(DataParams$Aggregate == 0){
#       if(timeStep == DataParams$SampStartYr+burn){
#         FDCatchData[DataParams$histStartYr:(timeStep-burn),,,] <- Data$histCatchDat
#         FDCPUEData[DataParams$histStartYr:(timeStep-burn),,,] <- Data$histCPUEDat
#         } else {
#         FDCatchData[(timeStep-burn),,,] <- Data$FDCatch
#         FDCPUEData[(timeStep-burn),,,] <- Data$FDCPUE
#         FICatchData[(timeStep-burn),,] <- Data$FICatch
#         FICPUEData[(timeStep-burn),,] <- Data$FICPUE
#         }
#       FDAgeData[(timeStep-burn),,,,] <- Data$AgeFD
#       FDSizeData[(timeStep-burn),,,,] <- Data$SizeFD
#     } else {
#       if(timeStep == DataParams$SampStartYr+burn){
#         FDCatchData[DataParams$histStartYr:(timeStep-burn),] <- Data$histCatchDat
#         FDCPUEData[DataParams$histStartYr:(timeStep-burn),] <- Data$histCPUEDat
#       } else {
#         FDCatchData[(timeStep-burn),] <- Data$FDCatch
#         FDCPUEData[(timeStep-burn),] <- Data$FDCPUE
#         FICatchData[(timeStep-burn)] <- Data$FICatch
#         FICPUEData[(timeStep-burn)] <- Data$FICPUE
#       }
#       FDAgeData[(timeStep-burn),,] <- Data$AgeFD
#       FDSizeData[(timeStep-burn),,] <- Data$SizeFD
#     }

#     }# if sample time steps
  
  ######## Should add if statement here with a management loop####
  ####################################################
  } #end if for season
 } #end fleet 

#==natural mortality (depends on timestep size)
  tempSNALT<-tempSNALT*exp(-(M/yearMark))

#==recruitment (happens once a year)
#==calculate recruits based on total spawnign biomass, allocate them according to either the spawning biomass in a patch or evenly.
#==current SpBio
tempNum<-rep(0,kmax)
for(i in 1:kmax)
  tempNum[i]<-sum(tempSNALT[,,i]*FishSel[i,flt])
NumbersExplt[timeStep]<-sum(tempNum)

tempSpBio<-rep(0,kmax)
for(i in 1:kmax)
 tempSpBio[i]<-sum(tempSNALT[,,i]*MatAge[i]*wgtAtAge[i])
SpawningBiomass[timeStep]<-sum(tempSpBio)

tempSpBio<-rep(0,kmax)
for(i in 1:kmax)
 tempSpBio[i]<-sum(-1*(NoTakeZoneImp-1)*tempSNALT[,,i]*MatAge[i]*wgtAtAge[i])
InsideMPAspbio[timeStep]<-sum(tempSpBio)

tempSpBio<-rep(0,kmax)
for(i in 1:kmax)
 tempSpBio[i]<-sum(NoTakeZoneImp*tempSNALT[,,i]*MatAge[i]*wgtAtAge[i])
OutsideMPAspbio[timeStep]<-sum(tempSpBio)

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
       RecDistMat<-habitat/sum(habitat)   
       tempSNALT[,,1]<-(Recruits*as.matrix(RecDistMat))
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

filled.contour(temp1[,,3],main="postfishery, premovement: age 2",zlim=c(0,2*OrigTempSNALT))
      mtext(side=3,paste("Timestep = ",timeStep))
#       ani.record()
filled.contour(tempSNALT[,,3],main="postfishery, postmovement: age 2",zlim=c(0,2*OrigTempSNALT))
      mtext(side=3,paste("Timestep = ",timeStep))
#       ani.record()
}
 #==update popdym
 if(timeStep<simTime)
  SpaceNumAtAgeT[timeStep+1,,,]<-tempSNALT


#===========================================================================
#== Calculate the cost of management for a given scenario
#===========================================================================

#==costs of management

if(sum(1-NoTakeZone)>1 & timeStep==initManage)
  CostOfManagement[timeStep]<- CostOfManagement[timeStep] + MPAsunk

if((SizeLimit[flt]>0)& timeStep==initManage)
  CostOfManagement[timeStep]<- CostOfManagement[timeStep] + SizeSunk

#==if there are any months in which fishing is not allowed, add sunk costs
FleetSeason<-season[,2]
FleetSeason[is.na(FleetSeason)]<-0
if(FleetN>1)
 FleetSeason<-apply(!is.na(season[,2:ncol(season)]),1,sum)
if(any(FleetSeason<FleetN,na.rm=T) & timeStep==initManage) 
  CostOfManagement[timeStep]<- CostOfManagement[timeStep] + SeasonSunk 

CostOfManagement[timeStep]<- CostOfManagement[timeStep] + sum(1-NoTakeZone)*MPAcost

if((SizeLimit[flt])>0)
 CostOfManagement[timeStep]<- CostOfManagement[timeStep] + SizeCost

#==if all fisheries are open, there is no cost of enforcing a season
if(FleetSeason[timeStep%%12+1]<FleetN)
 CostOfManagement[timeStep]<- CostOfManagement[timeStep] + SeasonCost

} # end timestep

list(CatchByFisher=CatchByFisher,CostByFisher=CostByFisher,ProfitByFisher=ProfitByFisher,CostOfManagement=CostOfManagement,
     SpawningBiomass=SpawningBiomass,SpaceEffort=SpaceEffort,InsideMPAspbio=InsideMPAspbio,OutsideMPAspbio=OutsideMPAspbio,
     ExploitableNumbers=NumbersExplt)
}


# outs<-Master(Life,SimCTL,Fleets,season,Samp,NoTakeZone,habitat,Graphs=F,GraphsFish=F,PrintLifeHistory=F)
# 
# totCatch<-apply(CatchByFisher,2,sum,na.rm=T)
# totCost<-apply(CostByFisher,2,sum,na.rm=T)
# totProfit<-apply(ProfitByFisher,2,sum,na.rm=T)
# 
# par(mfrow=c(4,1),mar=c(.1,4,.1,.1))
# plot(totCatch,type="b",xaxt='n',las=2)
# plot(totCost,lty=2,type="b",xaxt='n',las=2)
# plot(totProfit,lty=2,type="b",xaxt='n',las=2)
# lines(CostOfManagement,lty=2,type="b",col=2)
# plot(SpawningBiomass[burn:simTime],type="b",xaxt='n',las=2,ylab="SpawningBio",ylim=c(0,max(SpawningBiomass[burn:simTime],na.rm=T)))
# 
# dim(CatchByFisher)
# dim(CostByFisher)
# dim(ProfitByFisher)
# length(CostOfManagement)
# length(SpawningBiomass)


