LBAR<-function(LengthDat,LagLength,Weight,IncludeMPA,ReserveYr,OutsideBoundYr,Iterations,BootStrap,LifeError,Lc=NA)
{
  ##################
  ###### LBAR ######
  ##################
  #Source: Based off of Ault et al. 1998
  #Summary: Calculate F based on mean length of catch, assumption of natural mortality
  
  ################
  #### Inputs ####
  ################
  
  # LengthData: The matrix of length observations
  # LagLength: The number of years of historic length data to lag when calculating Lc. A lag of 1 calculates an 
  #           independent Lc for each year. A lag that is the same number of years as there is data calculates a 
  #           a single Lc for all years. An intermediate lag functions as a moving average.
  # Weight: The weight assigned to historic length data. 1 weights all years equally.
  # IncludeMPA: 1= Include length data from MPA, 0=Exclude length data from MPA. If no data from MPA is included,
  #           M is calculated using the mean of 3 life history based methods.            
  # ReserveYr: Either year of reserve implementation or NA
  # OutsideBoundYr: Either the year mgmt changed outside reserve (for bounding), or use NA for unbounded Lbar method
  # Iterations: Number of iterations to run the model. 1 makes the model run with just the default data
  # BootStrap: 1 means data are bootstrapped, 0 means not
  # Lc: You can manually set an Lc (for example, a min size limit), or set to NA and it will be calculated for you.It
  #           can be calculated as a separate Lc for each year, a single Lc for all years, or a moving avg.
  
  
  ################
  #### Organize Data ####
  ################
  
  # LengthDat<- LengthData
  
  # LagLength<- 1
  
  # Weight<- 0.2
  
  # IncludeMPA<- 0
  
  # Iterations<- 1000
  
  # BootStrap<- 1
  
  # LifeError<- 1
  
  # ReserveYr <- 2011
  
  # OutsideBoundYr <- NA
  
  # Lc <- NA
  
  ManualLc<- Lc
  
  Years<- sort(unique(LengthDat$Year))
  
  SampleSize<- NA
  
  Output<- as.data.frame(matrix(NA,nrow=length(Years),ncol=9))
  
  colnames(Output)<- c('Year','Method','SampleSize','Value','LowerCI','UpperCI','SD','Metric','Flag')
  
  MCOutput<- as.data.frame(matrix(NA,nrow=length(Years),ncol=7))
  
  colnames(MCOutput)<- c('Iteration','Year','Method','SampleSize','Value','Metric','Flag')
  
  MCDetails<- as.data.frame(matrix(NA,nrow=length(Years)*Iterations,ncol=11))
  
  colnames(MCDetails)<- c('TotalMortality','FishingMortality','NaturalMortality','Lbar_Inside',
                          'Lbar_Outside','Llam_Inside','Llam_Outside','Lc','Iteration','Year','BeddFmsy')  
  
  
  LengthDat<- LengthDat[is.na(LengthDat$Length)==F,]
  
  if (IncludeMPA==0) #Remove MPA data if needed
  {
    LengthDat<- LengthDat[LengthDat$MPA==0,]
  }
  
  c<- 0
  
  BaseFish<- Fish #Set Fish to base life history parameters
  
  for (i in 1:Iterations) #Loop over monte carlo runs
  {
    
    if (dim(LengthDat)[1]>0) #If any data are available
    {
      
      if (i>1 & LifeError==1) #Apply life history uncertainty 
      {
        
        Fish<- BaseFish
        
        Fish<- ApplyLifeHistoryError()
        
      }
      
      for (y in 1:length(Years)) #Loop over years
      {
        
        c<- c+1
        
        Flag<- 'None'
        
        Lengths<- LengthDat$Length[LengthDat$Year==Years[y]]
        SampleSize[y]<- length(Lengths)
        
        dYr<- y #current year
        
        weight<- Weight
        
        TotYrs <- length(Years)
        
        if (dYr == 1){
          tempLaggedYears<- seq(from=1,to=LagLength,by=1) # Pull out lagged years
        } else if (dYr == length(Years)){
          tempLaggedYears<- seq(from=dYr-(LagLength-1),to=dYr,by=1) # Pull out lagged years
        } else {
          start <- floor(LagLength/2)
          end <- ceiling(LagLength/2)
          if((dYr-start)<1) {    #Calculate correction if this is below 1, add it to end.
            CorrNum1 <- 1 - (dYr-start)
            start <- start - CorrNum1
            end <- end + CorrNum1
          }
          if((dYr+(end-1))>length(Years)){
            CorrNum2 <- dYr+(end-1)-length(Years)
            start <- start + CorrNum1
            end <- end - CorrNum1
          } 
          tempLaggedYears<- seq(from=dYr-start,to=dYr+(end-1),by=1)
        }
        
        Lag<-t(as.matrix(tempLaggedYears[tempLaggedYears>0])) 
        
        weights<- as.numeric(weight^t(apply(Lag,1,rev))) #Weight assigned to each year
        
        LaggedYears<- Years[Lag]
        
        TempData<- LengthDat[LengthDat$Year %in% LaggedYears,] #Pull out all the length data in the lagged years
        
        if (i>1 & BootStrap==1) #Bootstrap data if desired
        {
          
          NumPoints<- 1:dim(TempData)[1]  
          
          BootSample<- sample(NumPoints,length(NumPoints),replace=T)
          
          TempData<- TempData[BootSample,]	
        }
        
        sizestore<- rep(NA,length(LaggedYears))
        
        ##Loop over lagged years to calculate an Lc that will work for all years
        if(is.na(ManualLc)){
          for (d in 1:length(LaggedYears)) #Loop over lagged years
          {
            
            YearlySizes<- sort(TempData$Length[TempData$Year==LaggedYears[d]])
            
            #Assigns Lc to be the mode
            x<-table(YearlySizes) #create histogram of lengths
            maxPos <- which.max(x) #find the mode of the catch distribution
            sizestore[d] <- as.numeric(names(x)[min(maxPos)]) #Store the size corresponding to the first mode
            
          }
          #Calculate weighted Lc past lagged years
          # print(sizestore)
          
          Lc<- round(sum(sizestore*weights)/sum(weights),10) 
        }
        
        # print(Lc)
        
        #Calculate Llam -- used for bounding.
        
        if(IncludeMPA == 0){
          LlamIn <- NA
        } else if (Years[y]-ReserveYr < 2){
          #Need at least 2 years since reserve went in to calculate M
          LlamIn <- NA
        } else {
          #Calculate the Expected size of fish protected since reserve.
          AgeLc <- AgeAtLength(Lc,Fish,Fish$AgeSD)
          AgeBound <- AgeLc+(Years[y]-ReserveYr)
          LlamIn = LengthAtAge(AgeBound,Fish, Fish$LengthError)     
        }
        
        if(is.na(OutsideBoundYr)){
          LlamOut <- Fish$Linf
        } else {
          AgeLc <- AgeAtLength(Lc,Fish, Fish$AgeSD)
          AgeBound <- AgeLc+(Years[y]-OutsideBoundYr)
          LlamOut = LengthAtAge(AgeBound,Fish, Fish$LengthError)
        }
        
        #Calculate mean weights 
        if (is.na(LlamIn)){
          LbarIn <- NA
        } else {
          MPAdat <- subset(TempData, MPA == 1 & Year == Years[y] & Length >= Lc & Length < LlamIn)
          LbarIn <- mean(MPAdat$Length,na.rm=TRUE)
        }
        
        if (is.na(OutsideBoundYr)){      #Should outside length data be bounded?        
          Outdat <- subset(TempData, MPA == 0 & Year == Years[y] & Length >= Lc)
          # Outdat <- subset(TempData, TempData[TempData$MPA == 0 & TempData$Year == Years[y] & TempData$Length >= Lc,])
          
          LbarOut <- mean(Outdat$Length,na.rm=TRUE)
        } else {
          Outdat <- subset(TempData, MPA == 0 & Year == Years[y] & Length >= Lc & Length > LlamOut)
          LbarOut <- mean(Outdat$Length,na.rm=TRUE)
        }
        # print(c("In",Lc,LbarIn,LlamIn,"Out",Lc, LbarOut,LlamOut))  
        
        #################
        #### Analyze Data ####
        ################
        
        Bounded.BH <- function(z,Linf,K,Lc,Llam,Lbar) #Bounded Beverton-Holt function (see Ault. et al. 1998)
        {((z*(Lc-Lbar)+K*(Linf-Lbar))/(z*(Llam-Lbar)+K*(Linf-Lbar)))-(((Linf-Llam)/(Linf-Lc))^(z/K))}
        
        #allresults=as.data.frame(matrix(NA,nrow=1,ncol=6)) 
        
        #colnames(allresults)=c('M', 'F','Z','Lbar','FvFmsy','adjust')
        EmpM <- FALSE
        
        if (!is.na(LlamIn)){
          #Test whether the M function has a root  
          test<- NULL #Test whether the Z function has a root
          s<- seq(0.0000001,10,.01)
          for (t in 1:length(s))
          {
            
            test[t]<- Bounded.BH(s[t],Fish$Linf,Fish$vbk,Lc,LlamIn,LbarIn)    
          }
          
          if ((test[1]*test[length(s)])<0){ #multiply first and last to check for opposite signs.
            EmpM <- TRUE
          }
        }
        
        if (EmpM == FALSE | IncludeMPA == 0 | length(LengthDat$MPA == 1) == 0 | Years[y]-ReserveYr < 2)
        {
          #If no MPA data, or no root, or told not to use MPA, or not enough time has passed since 
          #MPA implementation then use sum of LHI methods
          M<- mean(c((Fish$MvK*Fish$vbk),(4.22/Fish$MaxAge),(exp(1.71-1.084*log(Fish$MaxAge)))))
          
        } else {
          # Calculate Bounded to estimate M
          BOUND_M <-(uniroot(function(z) Bounded.BH(z,Fish$Linf,Fish$vbk,Lc,LlamIn,LbarIn), c(0.0000001,10)))$root  #M   
          M <- mean(c(BOUND_M,(1.2*Fish$vbk),(4.22/Fish$MaxAge),(exp(1.71-1.084*log(Fish$MaxAge)))))
        }
        
        
        
        test<- NULL #Test whether the Z function has a root
        s<- seq(0.0000001,10,.01)
        for (t in 1:length(s))
        {
          
          test[t]<- Bounded.BH(s[t],Fish$Linf,Fish$vbk,Lc, LlamOut,LbarOut)
          
        }
        
        if ((test[1]*test[length(s)])>=0) #If the function has no root (aka no possibly natural mortality)
        {
          Flag<- 'Uniroot failed, function has no zero root'
          
          MCOutput[c,]<- data.frame(i,Years[y],'LBAR',SampleSize[y],-999,'FvM',Flag)
          
        }
        #If function has a root calcualte Z/M using bounded BH function
        
        if ((test[1]*test[length(s)])<0) 
        {
          
          
          BOUND_Z <-(uniroot(function(z) Bounded.BH(z,Fish$Linf,Fish$vbk,Lc,LlamOut,LbarOut), c(0.0000001,10)))$root  #M	 
          
          #Test Ault's numbert
          # BOUND_Z <-(uniroot(function(z) Bounded.BH(z,939,0.13,400,798,493), c(0.0000001,10)))$root  #M	 
          
          
          if (M > BOUND_Z) #If the result doesn't seem possible, calculate M from life history invariants
          {
            
            Flag<- 'Real M greater than Z, using Fish$M estimate'
            M<- Fish$M
            
          }
          
          
          BOUND_F <-BOUND_Z-M #Fishing mortality
          if (BOUND_F<0)
          {
            Flag<- 'LBAR Failed, Estimated Z too low'
          }
          
          FvFmsy <- BOUND_F/M #Calculate F/Fmsy
          
          ################
          #### Store Results ####
          ################
          
          MCOutput[c,]<- c(i,Years[y],'LBAR',SampleSize[y],FvFmsy,'FvM',Flag)
          
          
          MCDetails$TotalMortality[c]<- BOUND_Z
          
          MCDetails$FishingMortality[c]<- BOUND_F
          
          MCDetails$NaturalMortality[c]<- M
          
          MCDetails$Lbar_Inside[c]<- LbarIn
          
          MCDetails$Lbar_Outside[c]<- LbarOut
          
          MCDetails$Llam_Inside[c]<- LlamIn
          
          MCDetails$Llam_Outside[c]<- LlamOut
          
          MCDetails$Lc[c]<- Lc
          
          MCDetails$Iteration[c]<- i
          
          MCDetails$Year[c]<- Years[y]
          
          BeddFmsy<- (0.6*Fish$vbk)/(0.67-Lc/Fish$Linf)
          
          MCDetails$BeddFmsy[c]<- BeddFmsy
          
          
        } #Close if there's a root
        
        
      } #Close Year Loop	
      
    }
    if (dim(LengthDat)[1]==0) #If there are no data
      
    {
      Output<- c(Years,'LBAR',NaN,'FvM',NaN,NaN,NaN,'No Usable Data')
    }
    
  }
  
  TrueIteration<- MCOutput$Iteration==1
  
  TrueOutput<- MCOutput[TrueIteration,]
  
  # MCOutput<- MCOutput[TrueIteration==F,]
  
  Output$Year<- Years
  
  Output$Method<- 'LBAR'
  
  Output$Value<- TrueOutput$Value
  
  Output$LowerCI<- NaN
  
  Output$UpperCI<- NaN
  
  Output$SD<- NaN
  
  Output$Metric<- 'FvM'
  
  Output$Flag<- TrueOutput$Flag
  
  Output$SampleSize<- SampleSize
  
  if (Iterations>1) # If there are multiple iterations, process Monte Carlo output
  {
    
    for (y in 1:length(Years))
    {
      Where<- MCOutput$Year==Years[y]
      
      Temp<- MCOutput[Where,]
      
      TempValue<- sort(as.numeric(Temp$Value))
      
      Bottom<- ceiling(.025*length(TempValue))
      
      Top<- ceiling(.975*length(TempValue))
      
      MeanMetric<- mean(as.numeric(Temp$Value),na.rm=T)
      
      LowerCI<- TempValue[Bottom]
      
      UpperCI<- TempValue[Top]
      
      SD<- sd(TempValue[Bottom:Top],na.rm=T)
      
      Output$Year<- Years
      
      Output$Method[y]<- 'LBAR'
      
      Output$Value[y]<- TrueOutput$Value[y]
      
      Output$LowerCI[y]<- LowerCI
      
      Output$UpperCI[y]<- UpperCI
      
      Output$SD[y]<- SD
      
      Output$Metric[y]<- 'FvM'
      
    }
  }
  
  ################
  #### Produce Figures ####
  ################
  
  NumNans<-(is.na(MCDetails)[,1])
  if (sum(NumNans)<length(NumNans))
  {
    
    pdf(file=paste(FigureFolder,' LBAR FvM Boxplots.pdf',sep=''))
    boxplot((MCDetails$FishingMortality/M)~MCDetails$Year,frame=F,xlab='Year',ylab='F/M',notch=T,outline=F,width=SampleSize)
    dev.off()
    
    pdf(file=paste(FigureFolder,' LBAR Fishing Mortality Boxplots.pdf',sep=''))
    boxplot((MCDetails$FishingMortality)~MCDetails$Year,frame=F,xlab='Year',ylab='F',notch=T,outline=F,width=SampleSize)
    dev.off()
    
    pdf(file=paste(FigureFolder,' LBAR Total Mortality Boxplots.pdf',sep=''))
    boxplot((MCDetails$TotalMortality)~MCDetails$Year,frame=F,xlab='Year',ylab='Z',notch=T,outline=F, width=SampleSize)
    dev.off()
    
    pdf(file=paste(FigureFolder,' LBAR Mean Length Boxplots.pdf',sep=''))
    if(IncludeMPA == 1 & sum(is.na(MCDetails$Lbar_Inside))<length(MCDetails$Lbar_Inside)){
      par(mfrow=c(1,2))
      boxplot((MCDetails$Lbar_Inside)~MCDetails$Year,frame=F,xlab='Year',ylab='Mean Length',notch=T,main = "Mean Lengths Inside",outline=F, width=SampleSize)
      boxplot((MCDetails$Lbar_Outside)~MCDetails$Year,frame=F,xlab='Year',ylab='Mean Length',notch=T,main = "Mean Lengths Outside",outline=F, width=SampleSize)
    } else {
      boxplot((MCDetails$Lbar_Outside)~MCDetails$Year,frame=F,xlab='Year',ylab='Mean Length',notch=T,main = "Mean Lengths Outside",outline=F, width=SampleSize)
    }
    dev.off()
  }
  
  Fish<- BaseFish
  
  return(list(Output=Output,Details=MCDetails))	
}
