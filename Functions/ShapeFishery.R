ShapeFishery<- function(Life,Fleets,SimCTL,season,NoTakeZone,Samp)
{
  #   
  
  
  FleetN<-ncol(Fleets)-1
  
  #==life history ==============================
  kmax   	<-Life[grep('kmax',Life[,ncol(Life)]),seq(1,ncol(Life)-1)]              # maximum age
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
  
  if (mat95<=mat50)
  {
    mat95<- 1.1*mat50
  }
  
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
  
  SpaceR    <-SimCTL[grep('SpaceR',SimCTL[,2]),1]	    	# Rows in the grid space
  SpaceC	  <-SimCTL[grep('SpaceC',SimCTL[,2]),1]  			# cols in the grid spcae
  burn	  	<-SimCTL[grep('burn',SimCTL[,2]),1]    		  # burn in for movement equilibration (use option "Graphs" in "Initpop" to see that population is thoroughly mixed)
  simTime   <-SimCTL[grep('simTime',SimCTL[,2]),1]     	# time steps for projection
  yearMark  <-SimCTL[grep('yearMark',SimCTL[,2]),1]	  	# number of time steps in a year
  initManage<-SimCTL[grep('initManage',SimCTL[,2]),1]   # year in which to initiate management
  
  
  #==management costs==========================
  MPAsunk     <-SimCTL[grep('MPAsunk',SimCTL[,2]),1]    # start up cost of enforcing an MPA
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
  NoTakeZone<- NoTakeZone
  season<- season
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
  
  AgeLimit<-  round((log(1-SizeLimit/Linf)/-K)+t0)
  
  FishSel<-matrix(ncol=FleetN,nrow=kmax)
  for(y in 1:FleetN)
  {
    FishSel[,y]<-q[y]/(1+exp(-log(19)*((seq(1,kmax)-Sel50[y])/(Sel95[y]-Sel50[y])))) #Can we change this to be size based?
  }
 
  if (AgeLimit>0)
  {
    FishSel[1:AgeLimit,]<- 0
  }
  #==econ pars=================================
  price		  <-as.numeric(Fleets[grep('price',Fleets[,ncol(Fleets)]),seq(1,ncol(Fleets)-1)])		# ex-vessel price per unit harvest (kg)
  costTrv	  <-as.numeric(Fleets[grep('costTrv',Fleets[,ncol(Fleets)]),seq(1,ncol(Fleets)-1)])		# cost to travel one patch
  costFish	<-as.numeric(Fleets[grep('costFish',Fleets[,ncol(Fleets)]),seq(1,ncol(Fleets)-1)])		# cost per unit effort fishing
  discRate	<-as.numeric(Fleets[grep('discRate',Fleets[,ncol(Fleets)]),seq(1,ncol(Fleets)-1)])	# discount rate
  TimeHor	  <-as.numeric(Fleets[grep('TimeHor',Fleets[,ncol(Fleets)]),seq(1,ncol(Fleets)-1)])		# time horizon for evaluation
  Fishers	  <-as.numeric(Fleets[grep('Fishers',Fleets[,ncol(Fleets)]),seq(1,ncol(Fleets)-1)])		# number of fishers
  maxCapac	<-as.numeric(Fleets[grep('maxCapac',Fleets[,ncol(Fleets)]),seq(1,ncol(Fleets)-1)])		# max capacity of a single fisherman in a timestep (kg) (can also be a 'trip limit'/permit implementation)
  
  rm(Life,Fleets,SimCTL,Samp)
  
  Things<- ls()
  
  Fishery<- list()
  
  for (t in 1:length(Things))
  {
    eval(parse(text=paste('Fishery$',Things[t],'<- ',Things[t],sep='')))
  }
  
  
  return(Fishery)
} #Close Function
