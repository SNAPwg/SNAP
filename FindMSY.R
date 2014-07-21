

Life<-read.csv("LifeHistory.csv")                 # life history characteristics
SimCTL<-read.csv("GrandSimCtl.csv",header=F)               # simulation controls
Fleets<-read.csv("Fleets.csv",header=F)                   # fleet characteristics  
season<-read.csv("seasonNULL.csv",header=F)           # fishing seasons by fleet
Samp <- read.csv("SamplingParams.csv")            # sampling controls for management
NoTakeZone<-read.csv("notakezoneNULL.csv",header=F)   # marine protected areas (0=MPA, 1=open access, 2=TURF?) (make sure that there is no MPA for MSY calcs)
NoTakeZoneNULL<-read.csv("notakezoneNULL.csv",header=F)   # marine protected areas 
habitat<-read.csv("habitatNULL.csv",header=F)         # habitat quality (recruitment suitability)

#==set costs to zero
Fleets[grep('costTrv',Fleets[,ncol(Fleets)]),seq(1,ncol(Fleets)-1)]<-0
Fleets[grep('costFish',Fleets[,ncol(Fleets)]),seq(1,ncol(Fleets)-1)]<-0
simTimet <-SimCTL[grep('simTime',SimCTL[,2]),1]       # time steps for projection
yearMarkt<-SimCTL[grep('yearMark',SimCTL[,2]),1]	  	# number of time steps in a year
burnt	<-SimCTL[grep('burn',SimCTL[,2]),1]  
#==find the Capacity that produces MSY
#==specify max capacity in function argument (10 is the max capacity at msy)
capacityRange<-seq(7,11,.5)
SustainableYield<-rep(0,length(capacityRange))

for(x in 1:length(capacityRange))
{
  Fleets[grep('maxCapac',Fleets[,ncol(Fleets)]),seq(1,ncol(Fleets)-1)]<-capacityRange[x]
  OpenAccess<-Master(Life,SimCTL,Fleets,season,Samp,NoTakeZone,habitat,Graphs=F,GraphsFish=F,PrintLifeHistory=F)
  SustainableYield[x]<-sum(apply(OpenAccess$CatchByFisher,2,sum,na.rm=T)[(simTimet-burnt-yearMarkt):(simTimet-burnt)])
  print(SustainableYield[x])
}
  
  
#==save catch
#