
# Assess ------------------------------------------------------------------
# The assessment module of the 'Mam package.
# This function takes life history and data and performs a suite of
# assessments

Assess <- function(Fishery, Data, Assessments, BatchFolder, FigureFolder, MakeFigs = T,
                   NumIterations = 1, MinSampleSize = 150, DefaultSD = 0.05,
                   ReserveYear = NA,
                   Theme = theme_update())
{

  library(plyr, quietly = T)
  library(gridExtra, quietly = T)
  library(ggplot2, quietly = T)
  library(ggthemes, quietly = T)
  library(dplyr, quietly = T)
  library(ggmap, quietly = T)
  library(R2admb, quietly = T)
  library(tidyr, quietly = T)

  sapply(list.files(pattern="[.]R$", path="Functions", full.names=TRUE), source)

  # High Level Assessment Options -------------------------------------------
  #   Assessment <- 'Test Generic Module'

  #   dir.create(Assessment)

  #   NumIterations <- 10

  #   SPRRef <- 0.4

  #   CPUERef <- 0.4

  #   RunAssessments <- TRUE

  #   NumberOfSpecies <- 5

  #   ReserveYear <- 2007

  #   Assessments <- c('CatchCurve','DBSRA')

  MonteCarloNames  <- c('Site','Species','Assessment','Itertation','Metric','Value')
  #   DefaultSD <- 0.05

  #   MinSampleSize <- 150

  MPAColor <- "#1b9e77"

  FishedColor <- "#d95f02"

  Font <- 'Helvetica'

  FontColor <- 'Black'

  PlotFontSize <- 11

  # Load Data ---------------------------------------------------------------

  Files <- names(Data)

  if ( any(grepl('_LengthData',Files)) )
  {
    LengthData <- Data[[which(grepl('_LengthData',Files))]]
  }
  if ( any(grepl('_CatchData',Files)) )
  {
    CatchData <-  Data[[which(grepl('_CatchData',Files))]]
  }
  if ( any(grepl('_EffortData',Files)) )
  {
    EffortData <- Data[[which(grepl('_EffortData',Files))]]
  }
  if ( any(grepl('_CPUEData',Files)) )
  {
    CPUEData <- Data[[which(grepl('_CPUEData',Files))]]
  }
  if ( any(grepl('_DensityData',Files)) )
  {
    DensityData <- Data[[which(grepl('_DensityData',Files))]]
  }
  if ( any(grepl('_MapData',Files)) )
  {
    Locations <- Data[[which(grepl('_MapData',Files))]]
  }

  #   LifeHistory<-read.csv(paste(Assessment,'/',Files[grepl('_LifeHistory',Files)],sep=''), stringsAsFactors = F)

  #   Fishery$Species <- 'P. Gracilis'

  Fish <- with(Fishery,
               list(CommName = CommName, Taxa = Taxa,
                    HasLifeHistory = 1, MvK = MvK, Linf = Linf, vbk = K, t0 = t0,
                    MaxAge = kmax, AgeMat50 = kmat, LengthMatRatio = NA, Mat50 = NA,
                    WeightA = wtA, WeightB = wtB , GrowthType = PrinceType, VBErrorSlope = 0.1,
                    stringsAsFactors = F))

  Fishes <- unique(Fish$CommName)

  LHI<- read.csv('data/Prince et al 2014 LHI table.csv',stringsAsFactors=F)

  #   LifeData<- colnames(LifeHistory)[5:dim(LifeHistory)[2]]

  #   Fish <- CreateFish()
  #GFD now stands for generic fishery data since I don't feel like changing it

  Sites<- c('All')

  AssessmentResults<- list()

  MonteResults<- list()

  Counter<- 0

  # Run Assessments ---------------------------------------------------------


  #   if (RunAssessments==T)
  #   {

  for (f in seq_len(length(Fishes)))
    #     for (f in 9)
  {

    show(Fishes[f])

    Species<- Fishes[f]

    #       iGFD<- GFD[WhereSite & GFD$CommName==Fishes[f],]

    #       AssessmentName <- paste(Assessment,Sites[s],Fishes[f],sep='_')
    #
    #       Directory<- paste(Assessment, "/", Sites[s],'/',Fishes[f],'/',sep='')
    #
    #       if (file.exists(Directory)==F)
    #       {
    #         dir.create(paste(Assessment, "/", Sites[s],sep=''))
    #         dir.create( paste(Assessment, "/", Sites[s],'/',Fishes[f],'/',sep=''))
    #       }

    # Prepare Fish Object -------------------------------------------------------------

    #     SpeciesLifeHistory<- Fish[Fish$CommName == Species[f],colnames(Fish) %in% LifeData]

    sLHI <- filter(LHI, Taxa %in% Fish$Taxa )

    if (is.na(Fish$MvK)==F ){

      ClosestMvK<- which((Fish$MvK-sLHI$MeanMvK)^2==min((Fish$MvK-sLHI$MeanMvK)^2,na.rm=T))[1]

      Fish$LengthMatRatio<- sLHI$MeanLMATvMLinf[ClosestMvK]

      Fish$MinMvK<- sLHI$MinMvK[ClosestMvK]

      Fish$MaxMvK<- sLHI$MaxMvK[ClosestMvK]

      Fish$MinLengthMatRatio<- sLHI$MinLMATvMLinf[ClosestMvK]

      Fish$MaxLengthMatRatio<- sLHI$MaxLMATvMLinf[ClosestMvK]
    } else {

      ClosestMvK <- sLHI$GrowthType == Fish$GrowthType

      Fish$MvK <- (sLHI$MeanMvK)[ClosestMvK]

      Fish$LengthMatRatio<- sLHI$MeanLMATvMLinf[ClosestMvK]

      Fish$MinMvK<- sLHI$MinMvK[ClosestMvK]

      Fish$MaxMvK<- sLHI$MaxMvK[ClosestMvK]

      Fish$MinLengthMatRatio<- sLHI$MinLMATvMLinf[ClosestMvK]

      Fish$MaxLengthMatRatio<- sLHI$MaxLMATvMLinf[ClosestMvK]
    }

    Fish$LengthError<- DefaultSD

    Fish$M<- Fish$vbk*Fish$MvK

    Fish$MaxAge<- ceiling(-log(0.005)/Fish$M)
    if (is.na(Fish$Mat50) | Fish$Mat50 > (Fish$Linf* Fish$LengthMatRatio) ) #Use Prince et al. 2014 LHI
    {
      Fish$Mat50<- Fish$Linf*Fish$LengthMatRatio

      Fish$Mat95<- as.numeric(1.1*Fish$Mat50)
    }

    if (is.na(Fish$AgeMat50)==F)
    {
      Fish$Mat50<- LengthAtAge(Fish$AgeMat50,Fish,0)

      Fish$Mat95<- as.numeric(1.1*Fish$Mat50)
    }

    #       FigureFolder<- paste(Directory,'Figures/',sep='')
    #
    #       ResultFolder<- paste(Directory,'Results/',sep='')
    #
    #       if (file.exists(FigureFolder)==F)
    #       {
    #         dir.create(FigureFolder,recursive=T)
    #
    #         dir.create(ResultFolder,recursive=T)
    #       }

    #       write.csv(file=paste(ResultFolder,Species,' Life History.csv',sep=''),as.data.frame(Fish))


    LifePlot <- PlotLifeHistory(Fish)

    Theme<- theme(legend.position='top',plot.background=element_rect(color=NA),
                  rect=element_rect(fill='transparent',color=NA)
                  ,text=element_text(size=12,family=Font,color=FontColor),
                  axis.text.x=element_text(color=FontColor),
                  axis.text.y=element_text(color=FontColor),legend.key.size=unit(1,'cm'))

    Site <- 'All'

#     if (exists('LengthData')) {PlotLengthData(LengthData,FigureFolder,Fish,Species,Sites,Theme)}
#
#     if (exists('DensityData')) { PlotDensityData(DensityData,FigureFolder,Fish,Species,Sites,Theme)}
#
#     if (exists('CPUEData')) {PlotCPUEData(CPUEData,FigureFolder,Fish,Species,Sites,Theme)}
#
#     if (exists('CatchData')) {PlotCatchData(CatchData,FigureFolder)}

    #       MapCCFRP(ReformData)


    Fish$LHITol<- 0.99


    for (a in 1:length(Assessments)) #Loop over possible assessments, store in Assessment results. Many assessments have more detailed outputs than can also be accessed
    {

      Counter<- Counter+1
      if (Assessments[a]=='LBAR') #Run LBAR assessment
      {

        SampleCheck<- CheckLengthSampleSize(LengthData)

        if (SampleCheck$YearsWithEnoughData>0)
        {


          Temp<- LBAR(SampleCheck$ParedData,LagLength=1,Weight=0.2,IncludeMPA=0,ReserveYr=ReserveYear,OutsideBoundYr=NA,Iterations=NumIterations,
                      BootStrap=1,LifeError=1,Lc=NA)$Output

          StoreAssess<- data.frame(Species,Sites[s],Assessments[a],Temp,stringsAsFactors=F) %>%
            rename(Site=Sites.s.,Assessment=Assessments.a.)

          AssessmentResults[[Counter]]<-StoreAssess
        }
      }

      if (Assessments[a]=='CatchCurve') #Run Catch Curve analysis
      {

        SampleCheck<- CheckLengthSampleSize(LengthData)

        if (SampleCheck$YearsWithEnoughData>0)
        {


          Temp<- CatchCurve(SampleCheck$ParedData,CatchCurveWeight=0,WeightedRegression=1,
                            ReserveYr=ReserveYear,OutsideBoundYr=NA,ManualM=0,GroupMPA=1,Iterations=NumIterations,BootStrap=1,LifeError=1,HistInterval=1)

          MonteCarlo<- Temp$Details

          Temp<- Temp$Output

          StoreAssess<- data.frame(Species,Sites[s],Assessments[a],Temp,stringsAsFactors=F) %>%
            rename(Site=Sites.s.,Assessment=Assessments.a.)

          StoreMonte<- data.frame(Species,Sites[s],Assessments[a],MonteCarlo[,c('Iteration','Year','FvM')],stringsAsFactors=F) %>%
            rename(Site=Sites.s.,Assessment=Assessments.a.,Value=FvM) %>%
            mutate(Metric='F/M (CC)')

          MonteResults[[Counter]]<- StoreMonte

          AssessmentResults[[Counter]]<-StoreAssess


        }
      }


      if (Assessments[a]=='DensityRatio') #Run density ratio analysis
      {

        Temp<- DensityRatio(DensityData,LagLength=1,Weight=1,Form='Biomass',Iterations=NumIterations,BootStrap=1)

        MonteCarlo<- Temp$Details

        Temp<- Temp$Output

        StoreMonte<- data.frame(Species,Sites[s],Assessments[a],MonteCarlo[,c('Iteration','Year','DensityRatio')],stringsAsFactors=F) %>%
          rename(Site=Sites.s.,Assessment=Assessments.a.,Value=DensityRatio) %>%
          mutate(Metric='DensityRatio')

        MonteResults[[Counter]]<- StoreMonte

        StoreAssess<- data.frame(Species,Sites[s],Assessments[a],Temp,stringsAsFactors=F) %>%
          rename(Site=Sites.s.,Assessment=Assessments.a.)

        AssessmentResults[[Counter]]<-StoreAssess
      }

      if (Assessments[a]=='CPUERatio') #Run density ratio analysis
      {

        Temp<- CPUERatio(CPUEData,LagLength=1,Weight=1,Form='Biomass',Iterations=NumIterations,BootStrap=1)

        MonteCarlo<- Temp$Details

        Temp<- Temp$Output

        StoreMonte<- data.frame(Species,Sites[s],Assessments[a],MonteCarlo[,c('Iteration','Year','CPUERatio')],stringsAsFactors=F) %>%
          rename(Site=Sites.s.,Assessment=Assessments.a.,Value=CPUERatio) %>%
          mutate(Metric='CPUERatio')

        MonteResults[[Counter]]<- StoreMonte


        StoreAssess<- data.frame(Species,Sites[s],Assessments[a],Temp,stringsAsFactors=F) %>%
          rename(Site=Sites.s.,Assessment=Assessments.a.)

        AssessmentResults[[Counter]]<-StoreAssess
      }

      if (Assessments[a]=='CatchMSY')
      {
        Temp2<- CatchMSY(CatchData,1000,0.05,0,0,1,0,0,1,NA,c(0.75,0.99),NA,NA,c(0.25,0.65))

        Temp<- Temp2$Output

        StoreAssess<- data.frame(Species,Sites[s],Assessments[a],Temp,stringsAsFactors=F) %>%
          rename(Site=Sites.s.,Assessment=Assessments.a.)

        AssessmentResults[[Counter]]<-StoreAssess

      }

      if (Assessments[a]=='LBSPR') #Run LBSPR Assessment
      {

        SampleCheck<- CheckLengthSampleSize(LengthData)

        if (SampleCheck$YearsWithEnoughData>0)
        {

          LengthQuantile<- quantile(SampleCheck$ParedData$Length,na.rm=T)


          #           Temp2<- LBSPR(SampleCheck$ParedData,EstimateM=0,Iterations=1,BootStrap=1,
          #                         LifeError=1,LengthBins=1,ReserveYear=ReserveYear,SL50Min=LengthQuantile[1],
          #                         SL50Max=LengthQuantile[2],DeltaMin=NA,DeltaMax=NA,IncludeReserve=TRUE)

          Temp2<- LBSPR(SampleCheck$ParedData,EstimateM=0,Iterations=NumIterations,BootStrap=1,
                        LifeError=1,LengthBins=1,ReserveYear=ReserveYear,SL50Min=LengthQuantile[1],
                        SL50Max=LengthQuantile[2],DeltaMin=0.01,DeltaMax=.5*Fish$Linf,IncludeReserve=TRUE)

          MonteCarlo<- Temp2$Details

          StoreMonte<- data.frame(Species,Sites[s],Assessments[a],MonteCarlo[,c('Iteration','Year','FvM','SPR')],stringsAsFactors=F) %>%
            gather('Metric','Value',FvM,SPR,convert=T) %>%
            rename(Site=Sites.s.,Assessment=Assessments.a.) %>%
            select(Species,Site,Assessment,Iteration,Year,Value,Metric)

          StoreMonte$Metric[StoreMonte$Metric=='FvM']<- 'F/M (LBSR)'

          MonteResults[[Counter]]<- StoreMonte

          Temp<- Temp2$Output

          StoreAssess<- data.frame(Species,Sites[s],Assessments[a],Temp,stringsAsFactors=F) %>%
            rename(Site=Sites.s.,Assessment=Assessments.a.)

          AssessmentResults[[Counter]]<-StoreAssess
        }
      }


      if (Assessments[a]=='DBSRA') #Run DBSRA Assessment
      {
        DCAC.start.yr <- CatchData$Year[1] #start of the catch period
        DCAC.end.yr<- CatchData$Year[length(CatchData$Year)] #end of the catch period
        delta.yr<- CatchData$Year[length(CatchData$Year)] #Year that current depletion is fit to
        DBSRA.OFL.yr<- CatchData$Year[length(CatchData$Year)] #Year to calculate DBSRA OFL outputs
        FMSYtoMratio <- 0.8 #ratio of Fmsy to M
        SD.FMSYtoMratio<- 0.05
        Delta<- 0.7
        SD.Delta<- 0.1
        DeltaLowerBound<- 0.5
        DeltaUpperBound<- 0.9
        BMSYtoB0ratio <- 0.3
        SD.BMSYtoB0ratio<- 0.1
        BMSYtoB0LowerBound<- 0.2
        BMSYtoB0UpperBound<- 0.5
        CatchInterp<-1
        NIter<- 500


        Temp2<- DBSRA(CatchData, DCAC.start.yr, DCAC.end.yr, delta.yr, DBSRA.OFL.yr, FMSYtoMratio, SD.FMSYtoMratio, Delta, SD.Delta, DeltaLowerBound, DeltaUpperBound, BMSYtoB0ratio, SD.BMSYtoB0ratio, BMSYtoB0LowerBound, BMSYtoB0UpperBound, CatchInterp, NIter)

        Temp<- Temp2$Output

        StoreAssess<- data.frame(Species,Sites[s],Assessments[a],Temp,stringsAsFactors=F) %>%
          rename(Site=Sites.s.,Assessment=Assessments.a.)

        AssessmentResults[[Counter]]<-StoreAssess
      }


      #         show(paste('Finished ',Assessments[a],'-',round(100*a/length(Assessments),2),'% Done',sep=''))

    }

    CurrentResults<- ldply(AssessmentResults) %>% subset(Species==Fishes[f] & Site==Sites[s])

    #       save.image(file=paste(ResultFolder,AssessmentName,'_Settings.RData',sep='')) #Save settings used to produce current results
    #       write.csv(file=paste(ResultFolder,AssessmentName,'_Results.csv',sep=''),CurrentResults) #Save current results

  } #Close species  (f) loop

  # Process Assessments -----------------------------------------------------

  FlatAssessments<- ldply(AssessmentResults)

  Ramp<- 'RdYlGn'

  return(list(Results = Results, LengthStructure = LengthStructure, Diagnostics = Diagnostics))
}