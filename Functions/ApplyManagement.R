ApplyManagement<- function(Strat,Management,Fishery)
{
  
  if (Strat$SizeLimit==1)
  { 
    Fishery$SizeLimit<-  Management$SizeLimit 
  }
  
  if (Strat$NTZ==1)
  {
    #     NoTakeZoneImp<- Management$NTZ
    Fishery$NoTakeZone<- Management$NTZ
  }
  
  
  if (Strat$Season==1)
  {
    Fishery$season<- Management$Season
  }
  
  if (Strat$VesselBuyback==1)
  {
    Fishery$Fishers<- Management$Effort
    
  }
  
  if (Strat$Gear==1)
  {
    Fishery$Sel50<- Fleets[grep('Sel50',Fleets[,ncol(Fleets)]),seq(1,ncol(Fleets)-1)] *Management$Gear
    
    #     Fishery$Sel50<- Fishery$Sel50*Management$Gear
    
  }
  
  if (Strat$Capacity==1)
  {
    Fishery$Price<- Fleets[grep('price',Fleets[,ncol(Fleets)]),seq(1,ncol(Fleets)-1)] * Management$Capacity[1]
    
    Fishery$maxCapac<- Fleets[grep('maxCapac',Fleets[,ncol(Fleets)]),seq(1,ncol(Fleets)-1)] *Management$Capacity[2]
    #     Fishery$Sel50<- Fishery$Sel50*Management$Gear
    
  }
  
  
  if (Strat$Tax==1)
  {
    Fishery$costFish<- Fleets[grep('costFish',Fleets[,ncol(Fleets)]),seq(1,ncol(Fleets)-1)] * Management$Tax
    
  }
  
  Fishery$FleetN<-ncol(Fleets)-1
  
  
  #   Fishery$sampleTimeSteps <- seq(Fishery$DataParams$SampStartYr,simTime-burn+1,by=DataParams$SampFreq)  #Time steps in which sampling occurs
  
  
  AgeLimit<-  round((log(1-Fishery$SizeLimit/Fishery$Linf)/-Fishery$K)+Fishery$t0)  
  AgeLimit<-0
  if (Fishery$Sel95<=Fishery$Sel50)
  {
    Fishery$Sel95<- 1.01*Fishery$Sel50
  }
  
  
  Fishery$FishSel<-matrix(ncol=Fishery$FleetN,nrow=Fishery$kmax)
  for(y in 1:Fishery$FleetN)
  {
    Fishery$FishSel[,y]<-Fishery$q[y]/(1+exp(-log(19)*((seq(1,Fishery$kmax)-Fishery$Sel50[y])/(Fishery$Sel95[y]-Fishery$Sel50[y]))))
  }
  
  if (AgeLimit>0)
  {
    Fishery$FishSel[1:AgeLimit,]<- 0
  }
  #   show(Fishery$FishSel)
  return(Fishery)
}
