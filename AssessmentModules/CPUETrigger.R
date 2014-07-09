
CPUETrigger<- function(CPUEDat,TimeWindow,TimeWeight,CPUERef)
{
  # CPUE Trigger Analysis ---------------------------------------------------
  # This analysis calculates mean CPUE data across appropriate time steps, weighted as desired, against a reference CPUE level if available
  
  Output<- as.data.frame(matrix(NA,nrow=length(yr),ncol=9))
  
  colnames(Output)<- c('Year','Method','SampleSize','Value','LowerCI','UpperCI','SD','Metric','Flag')
  


  
  
  
  Output<- data.frame(yr,'DBSRA', lands.df$catch.mt,OFLStats$OFL, OFLStats$LowerQuantile, OFLStats$UpperQuantile, OFLStats$SD,'OFL',Flag)
  
  colnames(Output)<- c('Year','Method','SampleSize','Value','LowerCI','UpperCI','SD','Metric','Flag')
  
  Output$Method<- 'DBSRA'
  
  Output$Metric<- 'OFL'
  
  return(list(Output=Output,Details=list(AllResults=FlatResults,BioSeries=BioStats,BvBmsySeries=bStats,DCAC= DCAC.summary)))
}
