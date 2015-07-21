CheckLengthSampleSize<- function(Data)
{
#     Data<- LengthData
  Data<- subset(Data,is.na(Length)==F)
  
  SampleSize<- ddply(Data,c('Year'),summarize,SampleSize=length(Length))
  
  SufficientSamples<- SampleSize$SampleSize>MinSampleSize
  
  if (sum(SufficientSamples)>0)
  {
    AvailableYears<- SampleSize$Year[SufficientSamples]
    
    Data<- Data[Data$Year %in% AvailableYears,]
  }
  
  return(list(YearsWithEnoughData=sum(SufficientSamples),ParedData=Data))
}