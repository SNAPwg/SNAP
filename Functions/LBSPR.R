LBSPR<-function(LengthDat,EstimateM,Iterations,BootStrap,LifeError,LengthBins,ReserveYear,SL50Min,SL50Max,
                DeltaMin,DeltaMax,IncludeReserve)
{
  
  
  
  ######################
  ###### LBSPR #########
  ######################
  #Source: Based on the length based SPR methods developed by Prince, Valencia, Adrian
  #Summary: Estimate SPR by estimating selectivity, F using observed and predicted length frequency data
  
  ######################
  ###### Inputs #########
  ######################
  
  #   LengthDat: LengthData
  #   EstimateM: 1 estiamtes M using catch curve 
  #   Iterations: Number of iterations to run (1 runs only on the default data)
  #   BootStrap: 1 bootstaps the length data in the monte carlo run
  
  
  ######################
  ###### Process Data #########
  ######################
  LengthDat<- LengthDat[is.na(LengthDat$Length)==F,]
  
  # LengthDat$Length[LengthDat$Length>Fish$Linf]<- Fish$Linf * 0.98
  
  # LengthDat<- LengthDat[LengthDat$Length<Fish$Linf,]
  
  SampleSize<- NA
  
  Years<- sort(unique(LengthDat$Year))
  
  Output<- as.data.frame(matrix(NA,nrow=length(Years),ncol=9))
  
  colnames(Output)<- c('Year','Method','SampleSize','Value','LowerCI','UpperCI','SD','Metric','Flag')
  
  MCOutput<- as.data.frame(matrix(NA,nrow=length(Years)*Iterations,ncol=9))
  
  colnames(MCOutput)<- c('Iteration','Year','Method','Value','LowerCI','UpperCI','SD','Metric','Flag')
  
  Flag<- 'None'
  
  Details<- as.data.frame(matrix(NA,nrow=length(Years),ncol=4))
  
  colnames(Details)<- c('Year','FvM','SelL50','SelL95')
  
  MCDetails<- as.data.frame(matrix(NA,nrow=length(Years)*Iterations,ncol=7))
  
  colnames(MCDetails)<- c('Iteration','Year','FvM','F','SelL50','SelL95','SPR')
  
  ###################################
  # Source LBSPR function #    
  #################################### 
  
  WD<- getwd()
  
  LBSPRDir <- paste(WD,"LBSPR_ADMB_Code",sep="/")
  
  
  #   setwd(LBSPRDir)
  
  #      compile_admb("lbspr",verbose=FALSE)
  
  #   setwd(WD)  #changes directory back to main folder.
  
  ##Set Length bins
  # LengthBins <- 5 #1cm
  
  BaseFish<- Fish
  
  c<- 0
  
  CohortDeviates<- as.data.frame(matrix(NA,nrow=0,ncol=2))
  
  AgeDeviates<- as.data.frame(matrix(NA,nrow=0,ncol=2))
  
  #   library(dplyr)
  
  for (i in 1:Iterations)
  {
    
    if (i>1 & LifeError==1) #Apply life history error
    {
      Fish<- BaseFish
      Fish<- ApplyLifeHistoryError()
    }
    for (y in 1:length(Years)) #loop over years
    {
      #Read in the Size data and Species parameter file
      c<- c+1
      
      WhereFished<- LengthDat$Year==Years[y] & LengthDat$MPA==0
      
      WhereAll<- LengthDat$Year==Years[y] 
      
      AnyMPA<- sum(LengthDat$MPA[WhereAll])>0
      
      ## Convert lenth data into appropriate format-- vector of raw size data
      if (IncludeReserve==F)
      {
        CatchatLength <- LengthDat$Length[WhereFished]
      }
      if (IncludeReserve==T)
      {
        CatchatLength <- LengthDat$Length[WhereAll]
      }
      SampleSize[y]<- length(CatchatLength)
      
      if (i>1 & BootStrap==1) #Bootstrap length data 
      {
        
        NumPoints<- 1:length(CatchatLength)
        
        BootSample<- sample(NumPoints,length(NumPoints),replace=T)
        
        CatchatLength<- CatchatLength[BootSample]	
      }
      
      
      EstimatedM<-'None'
      
      if (EstimateM==1 & AnyMPA==T) #Estimate M using catch curve if you want
      {
        AllCatchAtLength<- LengthDat[WhereAll,]
        
        EstimatedM<- CatchCurve(AllCatchAtLength,'AgeBased',1,ReserveYear,NA,0,10,1,1,1)$Details$NaturalMortality
      }
      
      AssessPars<- BuildLBSPRPars(Fish=Fish,Mpow=0,NGTG=100,MaxSD=3,FecB=3,
                                  SL50Min=SL50Min,SL50Max=SL50Max,DeltaMin=DeltaMin,DeltaMax=DeltaMax)
      
      LenHist<- hist(CatchatLength,breaks=seq(min(CatchatLength,na.rm=T),
                                              max(CatchatLength,na.rm=T)+LengthBins,by=LengthBins),plot=F)
      
      LenFreq <- LenHist$counts/sum(LenHist$counts,na.rm=T)
      LenMids <- LenHist$mids 
      ADMBDir <- paste(WD,"/LBSPR_ADMB_Code",sep='')
      
      runMod <- RunLBSPRAssess(AssessPars, LenFreq, LenMids, ADMBDir, ExName="lbspr", MaxCount=10, ADMBRead=NULL)
      
      
      if (runMod$ModelFailed==T)
      {
        pdf(paste(FigureFolder,'LBSPR Fit.pdf',sep=''))
        plot(LenMids,LenFreq)
        text(x=mean(LenMids),y=mean(LenFreq),'LBSPR Model Failed',cex=3)      
        dev.off()
        
        MCOutput[c,]<- data.frame(i,Years[y],'LBSPR',NA,NA,NA,
                                  NA,'SPR',Flag,stringsAsFactors=F)
        
        MCDetails[c,]<- data.frame(i,Years[y],NA,
                                   NA,
                                   NA,
                                   NA,
                                   NA,stringsAsFactors=F)
      }
      if (runMod$ModelFailed==F)
      {
        
        FitDiagnostic<- data.frame(runMod$Bins,runMod$Obs,runMod$Pred,runMod$Unfished)
        colnames(FitDiagnostic)<- c('Length','Observed','Predicted','Unfished')
        FitDiagnostic<- gather(FitDiagnostic,Model,Proportion,2:4)
        
        Selectivity<- runMod$Estimates %>% subset(Par=='SL50' | Par=='SL95')
        pdf(paste(FigureFolder,'LBSPR Fit.pdf',sep=''))
        print(ggplot(FitDiagnostic,aes(Length,Proportion,color=Model))+geom_point(size=3,alpha=0.9)
              + geom_vline(xintercept=Selectivity$Est))
        dev.off()
        ResidualAnalysis<- AssessLBSPRResiduals(runMod,Fish,Years[y])
        
        CohortDeviates<- rbind(CohortDeviates,ResidualAnalysis$CohortDeviates)
        
        AgeDeviates<- rbind(AgeDeviates,ResidualAnalysis$AgeDeviates)
        
        MCOutput[c,]<- data.frame(i,Years[y],'LBSPR',runMod$Estimates$Est[runMod$Estimates$Par=='SPR'],NA,NA,
                                  runMod$Estimates$SD[runMod$Estimates$Par=='SPR'],'SPR',Flag,stringsAsFactors=F)
        
        MCDetails[c,]<- data.frame(i,Years[y],runMod$Estimates$Est[runMod$Estimates$Par=='FMpar'],
                                   runMod$Estimates$Est[runMod$Estimates$Par=='FMpar']*Fish$M,
                                   runMod$Estimates$Est[runMod$Estimates$Par=='SL50'],
                                   runMod$Estimates$Est[runMod$Estimates$Par=='SL95'],
                                   runMod$Estimates$Est[runMod$Estimates$Par=='SPR'],stringsAsFactors=F)
      }
      # Sel50, Sel95, F/M, and SPR stored in Estimates.
    } #Close year loop	
  } #Close iteration loop
  #   detach("package:dplyr", unload=TRUE)
  ########################################
  ###### Process Monte Carlo Data #########
  #########################################
  MCOutput<- subset(MCOutput,is.na(Iteration)==F)
  
  SampleSize<- ddply(LengthDat,c('Year'),summarize,Samples=length(Length))
  
  TrueIteration<- MCOutput$Iteration==1 & is.na(MCOutput$Iteration)==F
  
  TrueOutput<- MCOutput[TrueIteration,]
  
  TrueDetails<- MCDetails[TrueIteration,]
  
  # MCOutput<- MCOutput[TrueIteration==F,]
  Output$Year<- Years
  
  Output$Method<- 'LBSPR'
  
  Output$Value<- TrueOutput$Value
  
  Output$LowerCI<- NA
  
  Output$UpperCI<- NA
  
  Output$SD<- NA
  
  Output$Metric<- 'SPR'
  
  Output$Flag<-TrueOutput$Flag 
  
  Output$SampleSize<- SampleSize$Samples
  
  MCOutput<- join(MCOutput,SampleSize,by='Year')  
  
  FOutput<- ddply(MCDetails,c('Year'),summarize,Method='LBSPR',SampleSize=NA,Value=mean(FvM,na.rm=T),LowerCI=NA,UpperCI=NA,
                  SD=NA,Metric='FvM',Flag=NA)
  
  if (Iterations>1)
  {
    
    
    MCOutput$value<- MCOutput$Value
    
    Output<- ddply(MCOutput,c('Year'),summarize,Method='LBSPR',SampleSize=unique(Samples),Value=mean(value,na.rm=T),
                   LowerCI=quantile(value,0.025,na.rm=T),UpperCI=quantile(value,0.975,na.rm=T),SD=sd(value,na.rm=T),
                   Metric='SPR',Flag='None')
    
    FOutput<- ddply(MCDetails,c('Year'),summarize,Method='LBSPR',SampleSize=NA,Value=mean(FvM,na.rm=T),
                    LowerCI=quantile(FvM,0.025,na.rm=T),UpperCI=quantile(FvM,0.975,na.rm=T),
                    SD=sd(FvM,na.rm=T),Metric='FvM',Flag=NA)
    
    Output$Flag<-TrueOutput$Flag 
    
    MCOutput$value<- NULL
    
  }
  
  FOutput$SampleSize<- Output$SampleSize
  
  
  ######################
  ###### Make Plots #########
  ######################
  
  if (any(!is.na(MCOutput$Value)))
  {
    pdf(file=paste(FigureFolder,'Age Residuals Boxplots.pdf',sep=''))
    print(ggplot(data=AgeDeviates,aes(factor(Age),Residuals))+geom_boxplot(varwidth=F)
          +xlab('Age')+ylab('Residuals')+geom_hline(yintercept=0))
    dev.off()
    
    pdf(file=paste(FigureFolder,'Cohort Residuals Boxplots.pdf',sep=''))
    print(ggplot(data=CohortDeviates,aes(factor(Cohort),Residuals))+geom_boxplot(varwidth=F)
          +xlab('Cohort')+ylab('Residuals')+geom_hline(yintercept=0)+theme(axis.text.x=element_text(angle=45)))
    
    dev.off()
    
    if (length(unique(MCOutput$Year))>1)
    {
      
      pdf(file=paste(FigureFolder,' LBSPR SPR Boxplots.pdf',sep=''))
      print(ggplot(data=MCOutput,aes(factor(Year),Value,fill=Samples))+geom_boxplot() +
            scale_fill_gradient(low = 'red', high = 'green')
            + xlab('Year')+ylab('SPR'))
      dev.off()
      
      MCDetails<- join(MCDetails,SampleSize,by='Year')
      pdf(file=paste(FigureFolder,' LBSPR FvM Boxplots.pdf',sep=''))
      print(ggplot(data=MCDetails,aes(factor(Year),FvM,fill=Samples))+geom_boxplot(varwidth=F) + 
              scale_fill_gradient(low = 'red', high = 'green')
            +xlab('Year')+ylab('F/M'))  
      dev.off()

            pdf(file=paste(FigureFolder,' LBSPR F Boxplots.pdf',sep=''))
      print(ggplot(data=MCDetails,aes(factor(Year),F,fill=Samples))+geom_boxplot(varwidth=F) + 
              scale_fill_gradient(low = 'red', high = 'green') +
            xlab('Year')+ylab('F')+geom_hline(yintercept=Fish$M,linetype='longdash'))  
      dev.off()
    }
    if (length(unique(MCOutput$Year))==1)
    {
      pdf(file=paste(FigureFolder,' LBSPR SPR Boxplots.pdf',sep=''))
      print(ggplot(data=MCOutput,aes(factor(Year),Value))+geom_boxplot()
            +xlab('Year')+ylab('SPR'))
      dev.off()
      
      MCDetails<- join(MCDetails,SampleSize,by='Year')
      pdf(file=paste(FigureFolder,' LBSPR FvM Boxplots.pdf',sep=''))
      print(ggplot(data=MCDetails,aes(factor(Year),FvM))+geom_boxplot(varwidth=F)
            +xlab('Year')+ylab('F/M'))  
      dev.off()
      
      pdf(file=paste(FigureFolder,' LBSPR F Boxplots.pdf',sep=''))
      print(ggplot(data=MCDetails,aes(factor(Year),F))+geom_boxplot(varwidth=F)
            +xlab('Year')+ylab('F')+geom_hline(yintercept=Fish$M,linetype='longdash'))  
      dev.off()
    }
    
  }
  Fish<- BaseFish
  Output<- rbind(Output,FOutput)
  return(list(Output=Output,Details=MCDetails))
}