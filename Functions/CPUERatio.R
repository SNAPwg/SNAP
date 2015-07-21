CPUERatio<- function(CPUEDat,LagLength,Weight,Form,Iterations,BootStrap)
{
  ##################
  ###### DensityRatio ######
  ##################
  #Source:Loosely based on methods described in McGilliard et al. 2011 and Babcock & MacCall 2011
  #Summary: Estimate N/K from CPUE inside and outside of marine reserves
  
  # CPUEDat: The raw CPUE data
  # LagLength: The number of years of lagged data to use
  # Weight: The weight assigned to historic data
  #   #   
  #       CPUEDat<- DensityData
  #       
  #       LagLength=1
  #       
  #       Weight=0.2
  #       
  #       Form='Biomass'
  #       
  #       Iterations=10
  #       
  #       BootStrap=1
  #     
  ############################
  ### Process CPUE Data ###  
  ############################
  
  ### Figure out timeline of data you want ###
  
  Years<- sort(unique(CPUEDat$Year))
  
  lag<- LagLength #lagged years
  
  weight<- Weight
  
  Output<- as.data.frame(matrix(NA,nrow=length(Years),ncol=9))
  
  colnames(Output)<- c('Year','Method','SampleSize','Value','LowerCI','UpperCI','SD','Metric','Flag')

  MCOutput<- as.data.frame(matrix(NA,nrow=length(Years)*Iterations,ncol=6))
  
  colnames(MCOutput)<- c('Iteration','Year','Method','Value','Metric','Flag')
  
  MCDetails<- as.data.frame(matrix(NA,nrow=length(Years)*Iterations,ncol=4))
  
  colnames(MCDetails)<- c('Iteration','Year','FishedCPUE','MPACPUE')
  
  Flag<- 'None'
  
  c<- 0
  
  BaseFish<- Fish
  
  SampleSize<- ddply(CPUEDat,c('Year'),summarize,SampleSize=sum(Count,na.rm=T))
  for (i in 1:Iterations)
  {
    
    if (i>1)
    {
      Fish<- BaseFish
      
      Fish<- ApplyLifeHistoryError()
    }
    
    
    for (y in 1:length(Years))
    {
      
      c<- c+1
      
      Flag<- 'None'
      
      dYr<- y
      
      tempLaggedYears<- seq(from=dYr-lag,to=dYr,by=1)
      
      Lags<- t(as.matrix(tempLaggedYears[tempLaggedYears>0])) #Index of CPUEs that you want to grab for weighted calculation
      
      LaggedYears<- Years[Lags]
      
      weights<- weight^t(apply(Lags,1,rev)) #weight assigned to each year
      
      TempCPUEDat<- CPUEDat
      
      
      if (i>1 & BootStrap==1) #Resample length data inside and outside of MPAs
      {
        
        TempCPUEStorage<- as.data.frame(matrix(NA,nrow=0,ncol=dim(TempCPUEDat)[2]))
        
        colnames(TempCPUEStorage)<- colnames(TempCPUEDat)
        
        cc<- 0
        
        for (yy in 1:length(LaggedYears))
        {
          
          
          MPADat<- TempCPUEDat[TempCPUEDat$MPA==1 & TempCPUEDat$Year==LaggedYears[yy],]
          
          FishedDat<- TempCPUEDat[TempCPUEDat$MPA==0 & TempCPUEDat$Year==LaggedYears[yy],]
          
          NumPoints<- 1:dim(MPADat)[1]	
          
          BootSample<- sample(NumPoints,length(NumPoints),replace=T)
          
          TempMPADat<- MPADat[BootSample,]
          
          NumPoints<- 1:dim(FishedDat)[1]	
          
          BootSample<- sample(NumPoints,length(NumPoints),replace=T)
          
          TempFishedDat<- FishedDat[BootSample,]
          
          TempDat<- rbind(TempFishedDat,TempMPADat)
          
          Size<- dim(TempDat)[1]
          
          TempCPUEStorage[cc+1:Size,]<- TempDat
          
          cc<- cc+Size
          
        }
        
        TempCPUEDat<- TempCPUEStorage
      }
      
      #       ddply(TempCPUEDat,c('Year','MPA'),summarize,Count=sum(Count))
      WeightedCPUE<- CalculateCPUE(TempCPUEDat,LaggedYears,weights,Form)
      
      inside<- WeightedCPUE$MPACPUE
      
      outside<- WeightedCPUE$FishedCPUE
      
      WeightedRatio<- outside/inside
      
      WeightedIn<- inside
      
      WeightedOut<- outside
      
      
      MCOutput$Iteration[c]<- i
      
      MCOutput$Year[c]<- Years[y]
      
      MCOutput$Method[c]<- 'CPUERatio'
      
      MCOutput$Value[c]<- WeightedRatio
      
      MCOutput$Metric[c]<- 'Reserve/Fished'
      
      MCOutput$Flag[c]<- Form
      
      MCDetails$Iteration[c]<- i
      
      MCDetails$Year[c]<- Years[y]
      
      MCDetails$FishedCPUE[c]<- WeightedOut
      
      MCDetails$MPACPUE[c]<- WeightedIn
      
      
    } #Close year loop		
    
    
  } #Close iteration loop
  
  MCOutput$Value[is.nan(MCOutput$Value)]<- NA
  
  MCOutput$Flag[is.na(MCOutput$Value)]<- 'No CPUE Ratio Possible'
  
  MCOutput$Value[is.na(MCOutput$Value)]<- -999
  
  TrueIteration<- MCOutput$Iteration==1
  
  TrueOutput<- MCOutput[TrueIteration,]
  
  TrueDetails<- MCDetails[TrueIteration,]
  
  MCOutput<- MCOutput[TrueIteration==F,]
  
  Output$Year<- Years
  
  Output$Method<- 'CPUERatio'
  
  Output$Value<- TrueOutput$Value
  
  Output$LowerCI<- NA
  
  Output$UpperCI<- NA
  
  Output$SD<- NA
  
  Output$Metric<- 'Biomass CPUE Ratio'
  
  Output$Flag<-TrueOutput$Flag 
  
  Output$SampleSize<- SampleSize$SampleSize
  
  if (Iterations>1)
  {
    
    for (y in 1:length(Years))
    {
      Where<- MCOutput$Year==Years[y]
      
      Temp<- MCOutput[Where,]
      
      TempValue<- sort(as.numeric(Temp$Value))
      
      # pdf(file=paste(FigureFolder,Years[y],' CPUE Ratio Histogram.pdf',sep=''))
      # hist(TempValue,xlab='N/K',main=NA)
      # dev.off()
      
      Bottom<- ceiling(.025*length(TempValue))
      
      Top<- ceiling(.975*length(TempValue))
      
      MeanMetric<- mean(as.numeric(Temp$Value),na.rm=T)
      
      LowerCI<- TempValue[Bottom]
      
      UpperCI<- TempValue[Top]
      
      SD<- sd(TempValue[Bottom:Top],na.rm=T)
      Output[y,]<- I(data.frame(Years[y],'CPUERatio',SampleSize$SampleSize[y],MeanMetric,LowerCI,UpperCI,SD,'Biomass CPUE Ratio',TrueOutput$Flag[y],stringsAsFactors=F))
      
    }
    
  }
  
  MCDetails$MPACPUE[MCDetails$MPACPUE==0]<- NA

  MCDetails$CPUERatio<- MCDetails$FishedCPUE/MCDetails$MPACPUE
  
  pdf(file=paste(FigureFolder,' CPUE Ratio NvK Boxplots.pdf',sep=''))
  boxplot((MCDetails$FishedCPUE/MCDetails$MPACPUE)~MCDetails$Year,frame=F,xlab='Year',ylab='Fished CPUE/Unfished CPUE',notch=T,outline=F)
  dev.off()
  
  pdf(file=paste(FigureFolder,' CPUE Ratio Fished CPUE Boxplots.pdf',sep=''))
  boxplot((MCDetails$FishedCPUE)~MCDetails$Year,frame=F,xlab='Year',ylab='CPUE',notch=T,outline=F)
  dev.off()
  
  pdf(file=paste(FigureFolder,' CPUE Ratio MPA CPUE Boxplots.pdf',sep=''))
  boxplot((MCDetails$MPACPUE)~MCDetails$Year,frame=F,xlab='Year',ylab='CPUE',notch=T,outline=F)
  dev.off()
  
  Fish<- BaseFish
  
  return(list(Output=Output,Details=MCDetails))
  
}