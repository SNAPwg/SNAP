DanHist<- function(Data,Breaks)
{
  # Data<- TempLengthDat$Age[TempLengthDat$MPA==1]
  
  # Breaks<- BinBreaks
  
  BreakStore<- as.data.frame(matrix(NA,nrow=length(Breaks),ncol=3))
  colnames(BreakStore)<- c('Age','Frequency','LogFrequency')
  for (b in 1:(length(Breaks)-1))
  {
    BreakStore[b,1:2]<- c(Breaks[b],sum(Data>= Breaks[b] & Data< Breaks[b+1],na.rm=T))
    BreakStore[b,3]<- log(BreakStore[b,2])
    
  }
  
  BreakStore$LogFrequency[is.infinite(BreakStore$LogFrequency)]<- NA
  return(BreakStore[1:(length(Breaks)-1),])
}