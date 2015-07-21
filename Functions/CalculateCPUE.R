CalculateCPUE<- function(CPUE,Years,Weights,Form)
{
  
  CPUE$DistanceFromBorder[CPUE$DistanceFromBorder==-999]<- NA
  
  CPUE$DistanceFromBorder[is.na(CPUE$DistanceFromBorder)]<- mean(CPUE$DistanceFromBorder,na.rm=T)
  
  
  CPUEForm<- colnames(CPUE)==Form
  
  LagCPUE<- as.data.frame(matrix(NA,nrow= length(Years),ncol=4))
  
  colnames(LagCPUE)<- c('Year','MPACPUE','FishedCPUE','CPUERatio')
  
  WeightedCPUE<- as.data.frame(matrix(NA,nrow=1,ncol=4))
  
  colnames(WeightedCPUE)<- c('Year','MPACPUE','FishedCPUE','CPUERatio')
  for (y in 1:length(Years))
  {
    
    YearlyCPUE<- CPUE[CPUE$Year==Years[y],]
    
    Reserve<- YearlyCPUE$MPA==1
    
    MPACPUE<- sum(YearlyCPUE$DistanceFromBorder[Reserve]*(YearlyCPUE[Reserve,CPUEForm]/YearlyCPUE$AnglerHours[Reserve]))/sum(YearlyCPUE$DistanceFromBorder[Reserve])
    
    FishedCPUE<- sum(YearlyCPUE$DistanceFromBorder[Reserve==F]*(YearlyCPUE[Reserve==F,CPUEForm]/YearlyCPUE$AnglerHours[Reserve==F]),na.rm=T)/sum(YearlyCPUE$DistanceFromBorder[Reserve==F],na.rm=T)
    
    LagCPUE[y,]<- data.frame(Years[y],MPACPUE,FishedCPUE,FishedCPUE/MPACPUE)
    
  }
  WeightedCPUE[1,]<- data.frame(Years[length(Years)],sum(Weights*LagCPUE$MPACPUE)/sum(Weights),sum(Weights*LagCPUE$FishedCPUE)/sum(Weights),sum(Weights*LagCPUE$CPUERatio)/sum(Weights))
  
  return(WeightedCPUE)
  
}