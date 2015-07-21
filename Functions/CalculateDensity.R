CalculateDensity<- function(Densities,Years,Weights,Form)
{
  
  Densities$DistanceFromBorder[Densities$DistanceFromBorder==-999]<- NA
  
  Densities$DistanceFromBorder[is.na(Densities$DistanceFromBorder)]<- mean(Densities$DistanceFromBorder,na.rm=T)
  
  
  DensityForm<- colnames(Densities)==Form
  
  LagDensity<- as.data.frame(matrix(NA,nrow= length(Years),ncol=4))
  
  colnames(LagDensity)<- c('Year','MPADensity','FishedDensity','DensityRatio')
  
  WeightedDensity<- as.data.frame(matrix(NA,nrow=1,ncol=4))
  
  colnames(WeightedDensity)<- c('Year','MPADensity','FishedDensity','DensityRatio')	
  for (y in 1:length(Years))
  {
    
    YearlyDensity<- Densities[Densities$Year==Years[y],]
    
    Reserve<- YearlyDensity$MPA==1
    
    MPADensity<- sum(YearlyDensity$DistanceFromBorder[Reserve]*(YearlyDensity[Reserve,DensityForm]/YearlyDensity$SampleArea[Reserve]))/sum(YearlyDensity$DistanceFromBorder[Reserve])
    
    FishedDensity<- sum(YearlyDensity$DistanceFromBorder[Reserve==F]*(YearlyDensity[Reserve==F,DensityForm]/YearlyDensity$SampleArea[Reserve==F]),na.rm=T)/sum(YearlyDensity$DistanceFromBorder[Reserve==F],na.rm=T)
    
    #     MPADensity<- sum(YearlyDensity$DistanceFromBorder[Reserve]*YearlyDensity[Reserve,DensityForm])/sum(YearlyDensity$DistanceFromBorder[Reserve]*YearlyDensity$SampleArea[Reserve])
    #     
    #     FishedDensity<- sum(YearlyDensity$DistanceFromBorder[Reserve==F]*YearlyDensity[Reserve==F,DensityForm],na.rm=T)/sum(YearlyDensity$DistanceFromBorder[Reserve==F]*YearlyDensity$SampleArea[Reserve==F],na.rm=T)
    #     
    LagDensity[y,]<- data.frame(Years[y],MPADensity,FishedDensity,FishedDensity/MPADensity)
    
  }
  WeightedDensity[1,]<- data.frame(Years[length(Years)],sum(Weights*LagDensity$MPADensity)/sum(Weights),sum(Weights*LagDensity$FishedDensity)/sum(Weights),sum(Weights*LagDensity$DensityRatio)/sum(Weights))
  
  return(WeightedDensity)
  
}