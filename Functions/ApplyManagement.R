ApplyManagement<- function(Strat,Management)
{
  
    
  if (Strat$SizeLimit==1)
  {
    Fleets[grep('SizeLimit',Fleets[,ncol(Fleets)]),seq(1,ncol(Fleets)-1)]<<- Management$SizeLimit
  }
  
  if (Strat$NTZ==1)
  {
    NoTakeZoneImp<<- Management$NTZ
  }

  
  if (Strat$Season==1)
  {
    season<<- Management$Season
  }

  if (Strat$Effort==1)
  {
    Fleets[grep('Fishers',Fleets[,ncol(Fleets)]),seq(1,ncol(Fleets)-1)]<<- Management$Effort
    
  }

  if (Strat$Gear==1)
  {
    Fleets[grep('q',Fleets[,ncol(Fleets)]),seq(1,ncol(Fleets)-1)]<<- Management$Gear
    
  }

  if (Strat$Tax==1)
  {
    Fleets[grep('costFish',Fleets[,ncol(Fleets)]),seq(1,ncol(Fleets)-1)]<<- Fleets[grep('costFish',Fleets[,ncol(Fleets)]),seq(1,ncol(Fleets)-1)] * Management$Tax
    
  }
  
  
  
}
