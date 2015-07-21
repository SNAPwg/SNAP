DensityRatio<- function(DenDat,LagLength,Weight,Form,Iterations,BootStrap)
{
  ##################
  ###### DensityRatio ######
  ##################
  #Source:Loosely based on methods described in McGilliard et al. 2011 and Babcock & MacCall 2011
  #Summary: Estimate N/K from density inside and outside of marine reserves
  
  # DenDat: The raw density data
  # LagLength: The number of years of lagged data to use
  # Weight: The weight assigned to historic data
  #   #   
  #       DenDat<- DensityData
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
  ### Process Density Data ###	
  ############################
  
  ### Figure out timeline of data you want ###
  
  Years<- sort(unique(DenDat$Year))
  
  lag<- LagLength #lagged years
  
  weight<- Weight
  
  Output<- as.data.frame(matrix(NA,nrow=length(Years),ncol=9))
  
  colnames(Output)<- c('Year','Method','SampleSize','Value','LowerCI','UpperCI','SD','Metric','Flag')
  
  #   Output$Method<- as.character(Output$Method)
  #   
  #   Output$Metric<- as.character(Output$Metric)
  #   
  #   Output$Flag<- as.character(Output$Flag)
  #   
  MCOutput<- as.data.frame(matrix(NA,nrow=length(Years)*Iterations,ncol=6))
  
  colnames(MCOutput)<- c('Iteration','Year','Method','Value','Metric','Flag')
  
  MCDetails<- as.data.frame(matrix(NA,nrow=length(Years)*Iterations,ncol=4))
  
  colnames(MCDetails)<- c('Iteration','Year','FishedDensity','MPADensity')
  
  Flag<- 'None'
  
  c<- 0
  
  BaseFish<- Fish
  
  SampleSize<- ddply(DenDat,c('Year'),summarize,SampleSize=sum(Count,na.rm=T))
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
      
      TempDenDat<- DenDat
      
      
      if (i>1 & BootStrap==1) #Resample length data inside and outside of MPAs
      {
        
        TempDenStorage<- as.data.frame(matrix(NA,nrow=0,ncol=dim(TempDenDat)[2]))
        
        colnames(TempDenStorage)<- colnames(TempDenDat)
        
        cc<- 0
        
        for (yy in 1:length(LaggedYears))
        {
          
          
          MPADat<- TempDenDat[TempDenDat$MPA==1 & TempDenDat$Year==LaggedYears[yy],]
          
          FishedDat<- TempDenDat[TempDenDat$MPA==0 & TempDenDat$Year==LaggedYears[yy],]
          
          NumPoints<- 1:dim(MPADat)[1]	
          
          BootSample<- sample(NumPoints,length(NumPoints),replace=T)
          
          TempMPADat<- MPADat[BootSample,]
          
          NumPoints<- 1:dim(FishedDat)[1]	
          
          BootSample<- sample(NumPoints,length(NumPoints),replace=T)
          
          TempFishedDat<- FishedDat[BootSample,]
          
          TempDat<- rbind(TempFishedDat,TempMPADat)
          
          Size<- dim(TempDat)[1]
          
          TempDenStorage[cc+1:Size,]<- TempDat
          
          cc<- cc+Size
          
        }
        
        TempDenDat<- TempDenStorage
      }
      
      #       ddply(TempDenDat,c('Year','MPA'),summarize,Count=sum(Count))
      
      WeightedDensity<- CalculateDensity(TempDenDat,LaggedYears,weights,Form)
      
      inside<- WeightedDensity$MPADensity
      
      outside<- WeightedDensity$FishedDensity
      
      WeightedRatio<- outside/inside
      
      WeightedIn<- inside
      
      WeightedOut<- outside
      
      
      MCOutput$Iteration[c]<- i
      
      MCOutput$Year[c]<- Years[y]
      
      MCOutput$Method[c]<- 'DensityRatio'
      
      MCOutput$Value[c]<- WeightedRatio
      
      MCOutput$Metric[c]<- 'Reserve/Fished'
      
      MCOutput$Flag[c]<- Form
      
      MCDetails$Iteration[c]<- i
      
      MCDetails$Year[c]<- Years[y]
      
      MCDetails$FishedDensity[c]<- WeightedOut
      
      MCDetails$MPADensity[c]<- WeightedIn
      
      
    } #Close year loop		
    
    
  } #Close iteration loop
  
  MCOutput$Value[is.nan(MCOutput$Value)]<- NA
  
  MCOutput$Flag[is.na(MCOutput$Value)]<- 'No Density Ratio Possible'
  
  MCOutput$Value[is.na(MCOutput$Value)]<- -999
  
  TrueIteration<- MCOutput$Iteration==1
  
  TrueOutput<- MCOutput[TrueIteration,]
  
  TrueDetails<- MCDetails[TrueIteration,]
  
  MCOutput<- MCOutput[TrueIteration==F,]
  
  Output$Year<- Years
  
  Output$Method<- 'DensityRatio'
  
  Output$Value<- TrueOutput$Value
  
  Output$LowerCI<- NA
  
  Output$UpperCI<- NA
  
  Output$SD<- NA
  
  Output$Metric<- 'Biomass Density Ratio'
  
  Output$Flag<-TrueOutput$Flag 
  
  Output$SampleSize<- SampleSize$SampleSize
  
  if (Iterations>1)
  {
    
    for (y in 1:length(Years))
    {
      Where<- MCOutput$Year==Years[y]
      
      Temp<- MCOutput[Where,]
      
      TempValue<- sort(as.numeric(Temp$Value))
      
      # pdf(file=paste(FigureFolder,Years[y],' Density Ratio Histogram.pdf',sep=''))
      # hist(TempValue,xlab='N/K',main=NA)
      # dev.off()
      
      Bottom<- ceiling(.025*length(TempValue))
      
      Top<- ceiling(.975*length(TempValue))
      
      MeanMetric<- mean(as.numeric(Temp$Value),na.rm=T)
      
      LowerCI<- TempValue[Bottom]
      
      UpperCI<- TempValue[Top]
      
      SD<- sd(TempValue[Bottom:Top],na.rm=T)
      
      Output[y,]<- I(data.frame(Years[y],'DensityRatio',SampleSize$SampleSize[y],MeanMetric,LowerCI,UpperCI,SD,'Biomass Density Ratio',TrueOutput$Flag[y],stringsAsFactors=F))
      
    }
    
  }
  
  MCDetails$MPADensity[MCDetails$MPADensity==0]<- NA
  
  MCDetails$DensityRatio<- MCDetails$FishedDensity/MCDetails$MPADensity
  
  pdf(file=paste(FigureFolder,' Density Ratio NvK Boxplots.pdf',sep=''))
  boxplot((MCDetails$FishedDensity/MCDetails$MPADensity)~MCDetails$Year,frame=F,xlab='Year',ylab='Fished Density/Unfished Density',notch=T,outline=F)
  dev.off()
  
  pdf(file=paste(FigureFolder,' Density Ratio Fished Density Boxplots.pdf',sep=''))
  boxplot((MCDetails$FishedDensity)~MCDetails$Year,frame=F,xlab='Year',ylab='Density',notch=T,outline=F)
  dev.off()
  
  pdf(file=paste(FigureFolder,' Density Ratio MPA Density Boxplots.pdf',sep=''))
  boxplot((MCDetails$MPADensity)~MCDetails$Year,frame=F,xlab='Year',ylab='Density',notch=T,outline=F)
  dev.off()
  
  Fish<- BaseFish
  
  return(list(Output=Output,Details=MCDetails))
  
}