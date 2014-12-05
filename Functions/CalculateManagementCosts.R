CalculateManagementCosts<- function(Fishery,timeStep,initManage,ManageStrats,Management)
{
  
#   ManageStrats<- ManageStrats[1,]
  
  WhereCosts<- grepl('Cost',SimCTL[,2]) 
  
  if (timeStep==initManage)
  {  
    WhereCosts<- grepl('Sunk',SimCTL[,2])
  }
  
  ManagementUsed<- colnames(ManageStrats)[ManageStrats==1]
  
  CostStrings<- strsplit(as.character(SimCTL[WhereCosts,2])," ") #Grab the second string in each variable containing the term "error"
  
  CostStrings<- ldply(CostStrings, function(Err) Err[2])
  
  SimCTLNames<- strsplit(as.character(SimCTL[,2])," " )
  
  SimCTLNames<- ldply(SimCTLNames, function(Err) Err[2])
  
  CostsToApply<- t(CostStrings)[ t(CostStrings) %in% ManagementUsed]
  
  Costs<- sum(SimCTL[(t(SimCTLNames) %in% CostsToApply) & WhereCosts==T,1])
  
  return(Costs)
  
}