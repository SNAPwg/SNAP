
# Assess ------------------------------------------------------------------
# The assessment module of the 'Mam package.
# This function takes life history and data and performs a suite of
# assessments

Assess <- function(Life,Data,NumIterations = 1, ReserveYear = NA, Assessments = 'All', Theme = theme_update())
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

  Assessment <- 'Test Generic Module'

  dir.create(Assessment)

  NumIterations <- 10

  SPRRef <- 0.4

  CPUERef <- 0.4

  RunAssessments <- TRUE

  NumberOfSpecies <- 5

  ReserveYear <- 2007

  Assessments <- c('CatchCurve','DBSRA')

  MonteCarloNames  <- c('Site','Species','Assessment','Itertation','Metric','Value')

  DefaultSD <- 0.05

  MinSampleSize <- 150

  MPAColor <- "#1b9e77"

  FishedColor <- "#d95f02"

  Font <- 'Helvetica'

  FontColor <- 'Black'

  PlotFontSize <- 11

  # Load Data ---------------------------------------------------------------

  Files <- list.files(Assessment)

  if ( any(grepl('_LengthData',Files)) )
  {
    LengthData <- read.csv(paste(Assessment,'/',Files[grepl('_LengthData',Files)],sep=''), stringsAsFactors = F)
  }
  if ( any(grepl('_CatchData',Files)) )
  {
    CatchData <- read.csv(paste(Assessment,'/',Files[grepl('_CatchData',Files)],sep=''), stringsAsFactors = F)
  }
  if ( any(grepl('_EffortData',Files)) )
  {
    CatchData <- read.csv(paste(Assessment,'/',Files[grepl('_EffortData',Files)],sep=''), stringsAsFactors = F)
  }
  if ( any(grepl('_CPUEData',Files)) )
  {
    CPUEData <- read.csv(paste(Assessment,'/',Files[grepl('_CPUEData',Files)],sep=''), stringsAsFactors = F)
  }
  if ( any(grepl('_DensityData',Files)) )
  {
    DensityData <- read.csv(paste(Assessment,'/',Files[grepl('_DensityData',Files)],sep=''), stringsAsFactors = F)
  }
  if ( any(grepl('_MapData',Files)) )
  {
    Locations <- read.csv(paste(Assessment,'/',Files[grepl('_MapData',Files)],sep=''), stringsAsFactors = F)
  }

  LifeHistory<-read.csv(paste(Assessment,'/',Files[grepl('_LifeHistory',Files)],sep=''), stringsAsFactors = F)

  Fishes <- unique(LifeHistory$CommName)

  LHI<- read.csv('data/Prince et al 2014 LHI table.csv',stringsAsFactors=F)

  LifeData<- colnames(LifeHistory)[5:dim(LifeHistory)[2]]

  source('data/Default_Controlfile.R')
  #GFD now stands for generic fishery data since I don't feel like changing it


  Fish$AgeMat50<- NA

  Fish$AgeMatSource<- NA


  Sites<- c('All')

  AssessmentResults<- list()

  MonteResults<- list()

  Counter<- 0


  # Run Assessments ---------------------------------------------------------


  if (RunAssessments==T)
  {

    for (f in seq_len(length(Fishes)))
      #     for (f in 9)
    {

      show(Fishes[f])

      Species<- Fishes[f]

      #       iGFD<- GFD[WhereSite & GFD$CommName==Fishes[f],]

      AssessmentName <- paste(Assessment,Sites[s],Fishes[f],sep='_')

      Directory<- paste(Assessment, "/", Sites[s],'/',Fishes[f],'/',sep='')

      if (file.exists(Directory)==F)
      {
        dir.create(paste(Assessment, "/", Sites[s],sep=''))
        dir.create( paste(Assessment, "/", Sites[s],'/',Fishes[f],'/',sep=''))
      }

      # Prepare Fish Object -------------------------------------------------------------

      SpeciesLifeHistory<- LifeHistory[LifeHistory$CommName == Species[s],colnames(LifeHistory) %in% LifeData]

      sTaxa<- LifeHistory$Taxa[LifeHistory$CommName == Species[s]]

      HasLifeHistory<- SpeciesLifeHistory[which(is.na(SpeciesLifeHistory)==F)]

      HasLife<- colnames(HasLifeHistory)

      Fish$CommName <- Species[s]

      for (l in 1:length(HasLife))
      {
        WhereLife<- which(names(Fish)==HasLife[l])

        Fish[[WhereLife]]<- as.numeric(HasLifeHistory[l])
      }

      Fish$M<- Fish$vbk*Fish$MvK

      Fish$MaxAge<- ceiling(-log(0.005)/Fish$M)

      sLHI <- filter(LHI, Taxa %in% sTaxa )

      ClosestMvK<- which((Fish$MvK-LHI$MeanMvK)^2==min((Fish$MvK-LHI$MeanMvK)^2,na.rm=T))[1]

      Fish$LengthMatRatio<- LHI$MeanLMATvMLinf[ClosestMvK]

      Fish$MinMvK<- LHI$MinMvK[ClosestMvK]

      Fish$MaxMvK<- LHI$MaxMvK[ClosestMvK]

      Fish$MinLengthMatRatio<- LHI$MinLMATvMLinf[ClosestMvK]

      Fish$MaxLengthMatRatio<- LHI$MaxLMATvMLinf[ClosestMvK]

      if (is.na(SpeciesLifeHistory$Mat50) | SpeciesLifeHistory$Mat50>(SpeciesLifeHistory$Linf* SpeciesLifeHistory$LengthMatRatio) ) #Use Prince et al. 2014 LHI
      {
        Fish$Mat50<- Fish$Linf*Fish$LengthMatRatio

        Fish$Mat95<- as.numeric(1.1*Fish$Mat50)
      }

      if (is.na(SpeciesLifeHistory$AgeMat50)==F)
      {
        Fish$Mat50<- LengthAtAge(SpeciesLifeHistory$AgeMat50,Fish,0)

        Fish$Mat95<- as.numeric(1.1*Fish$Mat50)
      }

      Fish$Mat95<- as.numeric(1.1*Fish$Mat50)

      #       ReformData<- FormatCCFRPData(iGFD)

      #       LengthData<- ReformData$LengthData
      #
      #       DensityData<- ReformData$DensityData
      #
      #       CPUEData<- ReformData$CPUEData
      #
      #       write.csv(file=paste(Directory,AssessmentName,'_LengthData.csv',sep=''),LengthData)
      #
      #       write.csv(file=paste(Directory,AssessmentName,'_DensityData.csv',sep=''),DensityData)
      #
      #       write.csv(file=paste(Directory,AssessmentName,'_CPUEData.csv',sep=''),DensityData)
      #
      FigureFolder<- paste(Directory,'Figures/',sep='')

      ResultFolder<- paste(Directory,'Results/',sep='')

      if (file.exists(FigureFolder)==F)
      {
        dir.create(FigureFolder,recursive=T)

        dir.create(ResultFolder,recursive=T)
      }

      #       write.csv(file=paste(ResultFolder,Species,' Life History.csv',sep=''),as.data.frame(Fish))

      PlotLifeHistory()

      Theme<- theme(legend.position='top',plot.background=element_rect(color=NA),
                    rect=element_rect(fill='transparent',color=NA)
                    ,text=element_text(size=12,family=Font,color=FontColor),
                    axis.text.x=element_text(color=FontColor),
                    axis.text.y=element_text(color=FontColor),legend.key.size=unit(1,'cm'))

      if (exists('LengthData')) {PlotLengthData(LengthData,FigureFolder,Fish,Species,Sites[s],Theme)}


      if (exists('DensityData')) { PlotDensityData(DensityData,FigureFolder,Fish,Species,Sites[s],Theme)}

      if (exists('CPUEData')) {PlotCPUEData(CPUEData,FigureFolder,Fish,Species,Sites[s],Theme)}

      if (exists('CatchData')) {PlotCatchData(CatchData,FigureFolder)}

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


        show(paste('Finished ',Assessments[a],'-',round(100*a/length(Assessments),2),'% Done',sep=''))

      }

      CurrentResults<- ldply(AssessmentResults) %>% subset(Species==Fishes[f] & Site==Sites[s])

      if (any(CurrentResults$Assessment=='CatchCurve') & any(CurrentResults$Assessment=='CPUERatio')
          & any(CurrentResults$Assessment=='LBSPR') )
      {

        Theme<- theme(legend.position='top',plot.background=element_rect(color=NA),
                      rect=element_rect(fill='transparent',color=NA)
                      ,text=element_text(size=12,family=Font,color=FontColor),
                      axis.text.x=element_text(color=FontColor),
                      axis.text.y=element_text(color=FontColor))

        SummaryPanel(CurrentResults,LengthData,CPUEData,Species,Sites[s],YearsToSmooth=3,Theme)

      }
      save.image(file=paste(ResultFolder,AssessmentName,'_Settings.RData',sep='')) #Save settings used to produce current results
      write.csv(file=paste(ResultFolder,AssessmentName,'_Results.csv',sep=''),CurrentResults) #Save current results

    } #Close species  (f) loop

    save(file=paste(Assessment,'/Assessment Results.Rdata',sep=''),AssessmentResults,MonteResults)

    save(file=paste(Assessment,'/Raw Data.Rdata',sep=''),GFD)
  } #Close Run Assessments
  if (RunAssessments==F)
  {
    try(load(file=paste(Assessment,'/Assessment Results.Rdata',sep='')))
    try(load(file=paste(Assessment,'/Raw Data.Rdata',sep='')))
  }

  # Process Assessments -----------------------------------------------------

  FlatAssessments<- ldply(AssessmentResults)

  Ramp<- 'RdYlGn'

  LineKeynoteTheme<- theme(plot.background=element_rect(color=NA),rect=element_rect(fill='transparent',color=NA)
                           ,text=element_text(size=16,family=Font,color=FontColor))


  PanelTheme<- theme(text=element_text(size=16,family=Font,color=FontColor))

  KeynoteTheme<- theme(plot.background=element_rect(color=NA),rect=element_rect(fill='transparent',color=NA)
                       ,text=element_text(size=16,family=Font,color=FontColor),
                       axis.text.x=element_text(angle=45,hjust=0.9,vjust=0.9,color=FontColor),
                       axis.text.y=element_text(color=FontColor),legend.position='none')


  AggPlot<- function(PlotDat,PlotName,Theme,Order)
  {

    PlotDat$Species<- as.factor(PlotDat$Species)

    PlotDat$Species<- reorder((PlotDat$Species),PlotDat$Value)

    if (Order=='Descending')
    {
      PlotDat$Species<- reorder((PlotDat$Species),1/(PlotDat$Value+abs(min(PlotDat$Value,na.rm=T))))
    }



    print(ggplot(data=PlotDat,
                 aes(x=(Species),y=Value,fill=Species))+scale_fill_brewer(palette=Ramp)
          +ylab(PlotName)+xlab(NULL)
          +geom_bar(stat='identity',color='black')+
            Theme+geom_hline(yintercept=0,size=2))

  }

  OrderDat<- function(PlotDat,OrdVar,OrdBy,Order)
  {
    PlotDat[,OrdVar]<- as.factor(PlotDat[,OrdVar])

    PlotDat[,OrdVar]<- reorder(PlotDat[,OrdVar],PlotDat[,OrdBy])

    if (Order=='Descending')
    {
      PlotDat$Species<- reorder(PlotDat[,OrdVar],1/(PlotDat[,OrdBy]+abs(min(PlotDat[,OrdBy],na.rm=T))))
    }

    return(PlotDat)
  }

  pdf(paste(Assessment,'/CC FvM.pdf',sep=''))

  PlotDat<- subset(FlatAssessments,Year==2014 & Metric=='FvM' & Method=='CatchCurve' & Site=='All')

  CCFvMPlot<- AggPlot(PlotDat,'FvM',KeynoteTheme,'Descending')$plot

  print(CCFvMPlot)
  dev.off()

  pdf(paste(Assessment,'/LBSPR FvM.pdf',sep=''))

  PlotDat<- subset(FlatAssessments,Year==2014 & Metric=='FvM' & Method=='LBSPR' & Site=='All')

  LBSPRFvMPlot<- AggPlot(PlotDat,'FvM',KeynoteTheme,'Descending')$plot

  print(LBSPRFvMPlot)

  dev.off()

  pdf(paste(Assessment,'/LBSPR SPR.pdf',sep=''))

  PlotDat<- subset(FlatAssessments,Year==2014 & Metric=='SPR' & Site=='All')

  SPRPlot<- AggPlot(PlotDat,'SPR',KeynoteTheme,'Ascending')$plot

  print(SPRPlot)

  dev.off()

  pdf(paste(Assessment,'/CPUE Ratio.pdf',sep=''))

  PlotDat<- subset(FlatAssessments,Year==2014 & Metric=='Biomass CPUE Ratio' & Site=='All')

  CPUEPlot<- AggPlot(PlotDat,'CPUE Ratio',KeynoteTheme,'Ascending')$plot

  print(CPUEPlot)

  dev.off()

  AllAssess<- subset(FlatAssessments,Site=='All')  %>% mutate(SA=paste(Species,Method,Metric,sep='-'))

  SAS<- unique(AllAssess$SA)

  AllAssess$FSlope<- NA

  AllAssess$Value[!is.finite(AllAssess$Value)]<- NA

  AllAssess$LowerCI[!is.finite(AllAssess$LowerCI)]<- NA

  AllAssess$UpperCI[!is.finite(AllAssess$UpperCI)]<- NA

  for (i in 1:length(SAS))
  {

    FSlope<- lm(Value~Year,data=subset(AllAssess,SA==SAS[i]))$coefficient[2]

    AllAssess$Slope[AllAssess$SA==SAS[i]]<- FSlope
  }

  write.csv(file=paste(Assessment,'/Complete Assessment Results.csv',sep=''),AllAssess)

  Ordered<- OrderDat(subset(AllAssess,Year==2014 & Metric=='FvM' & Method=='CatchCurve'),'Species','Slope','Descending')

  pdf(paste(Assessment,'/CC FvM Trend.pdf',sep=''))

  CCFvMTrendPlot<- (ggplot(data=,Ordered,aes(x=factor(Species),y=Slope,fill=Species))+scale_fill_brewer(palette=Ramp)+ylab('Trend in F/M')
                    +geom_bar(stat='identity',color='black')+KeynoteTheme+xlab(NULL)+geom_hline(yintercept=0,size=2))

  print(CCFvMTrendPlot)

  dev.off()

  Ordered<- OrderDat(subset(AllAssess,Year==2014 & Metric=='SPR'),'Species','Slope','Increasing')

  pdf(paste(Assessment,'/SPR Trend.pdf',sep=''))


  SPRTrend<- (ggplot(data=,Ordered,aes(x=factor(Species),y=Slope,fill=Species))+scale_fill_brewer(palette=Ramp)+ylab('Trend in SPR')
              +geom_bar(stat='identity',color='black')+KeynoteTheme+xlab(NULL)+geom_hline(yintercept=0,size=2))

  print(SPRTrend)

  dev.off()

  Ordered<- OrderDat(subset(AllAssess,Year==2014 & Metric=='Biomass CPUE Ratio'),'Species','Slope','Increasing')

  pdf(paste(Assessment,'/CPUE Ratio Trend.pdf',sep=''))

  CPUETrend<- (ggplot(data=,Ordered,aes(x=factor(Species),y=Slope,fill=Species))+scale_fill_brewer(palette=Ramp)+ylab('Trend in CR')
               +geom_bar(stat='identity',color='black')+KeynoteTheme+xlab(NULL)+geom_hline(yintercept=0,size=2))
  print(CPUETrend)

  dev.off()

  # More Summary Plots ------------------------------------------------------


  AllMonte<- ldply(MonteResults) %>%
    filter(is.na(Value)==F & Value<2000)

  AllMonte$RefPoint[grepl('F/M',AllMonte$Metric)]<- AllMonte$Value[grepl('F/M',AllMonte$Metric)]

  AllMonte$RefPoint[AllMonte$Metric=='SPR']<- (AllMonte$Value/SPRRef )[AllMonte$Metric=='SPR']

  AllMonte$RefPoint[AllMonte$Metric=='CPUERatio']<- (AllMonte$Value/CPUERef )[AllMonte$Metric=='CPUERatio']

  AllMonte$GoodThing[grepl('F/M',AllMonte$Metric) & AllMonte$RefPoint>1] <- 'Bad Thing'

  AllMonte$GoodThing[grepl('F/M',AllMonte$Metric) & AllMonte$RefPoint<=1] <- 'Good Thing'

  AllMonte$GoodThing[AllMonte$Metric == 'SPR' & AllMonte$RefPoint<=1] <- 'Bad Thing'

  AllMonte$GoodThing[AllMonte$Metric == 'SPR' & AllMonte$RefPoint>1] <- 'Good Thing'

  AllMonte$GoodThing[AllMonte$Metric == 'CPUERatio' & AllMonte$RefPoint<=1] <- 'Bad Thing'

  AllMonte$GoodThing[AllMonte$Metric == 'CPUERatio' & AllMonte$RefPoint>1] <- 'Good Thing'

  AllMonte$RelRefPoint[grepl('F/M',AllMonte$Metric)]<- 100*(1-AllMonte$Value[grepl('F/M',AllMonte$Metric)])

  AllMonte$RelRefPoint[AllMonte$Metric=='SPR']<- 100*((AllMonte$Value/SPRRef )-1)[AllMonte$Metric=='SPR']

  AllMonte$RelRefPoint[AllMonte$Metric=='CPUERatio']<- 100*((AllMonte$Value/CPUERef)-1)[AllMonte$Metric=='CPUERatio']

  AllMonte$RelRefPoint<- pmax(-200,pmin(AllMonte$RelRefPoint,200))


  AllMonte<- AllMonte %>%
    group_by(Species) %>%
    mutate(RanAllAssess=length(unique(Metric))==4) %>%
    ungroup() %>%
    group_by(Species,Metric,Year) %>%
    mutate(MeanRelRef=mean(RelRefPoint,na.rm=T))

  PaperBarTheme<- theme(text=element_text(size=PlotFontSize,family = Font,color = FontColor),
                        axis.text.x=element_text(angle=45,hjust=0.9,vjust=0.9))
  # PaperTheme<- theme(text=element_text(size=PlotFontSize,family = Font,color = FontColor),panel.background = element_blank())

  PaperTheme<- theme(text=element_text(size=PlotFontSize,family = Font,color = FontColor))


  AllMonte<- AllMonte %>%
    group_by(Species,Metric) %>%
    mutate(Lag1RefPoint=lag(RelRefPoint),Lag2RefPoint=lag(RelRefPoint,2),Lag1Value=lag(Value,1),Lag2Value=lag(Value,2))

  AllMonte$SmoothRelRefPoint=apply(AllMonte[,c('RelRefPoint','Lag1RefPoint','Lag2RefPoint')],1,weighted.mean,c(1,1,1))

  AllMonte$SmoothValue=apply(AllMonte[,c('Value','Lag1Value','Lag2Value')],1,weighted.mean,c(1,1,1))


  pdf(file=paste(Assessment,'Percent Deviation Plot.pdf',sep='/'))
  PercDevPlot<- (ggplot(data=subset(AllMonte,Year==max(Year) & RanAllAssess==T),aes(factor(Metric),SmoothRelRefPoint,fill=MeanRelRef))
                 +geom_boxplot()+facet_wrap(~Species)+PaperBarTheme+geom_hline(yintercept=1)+
                   ylab('% Deviation From Target')
                 +theme(legend.position='none',axis.title.x=element_blank())+scale_fill_gradient2(low='red',mid='white',high='green'))
  (PercDevPlot)
  dev.off()

  pdf(file=paste(Assessment,'Reference Point Plot.pdf',sep='/'))
  ReferencePlot<- (ggplot(data=subset(AllMonte,Year==max(Year) & RanAllAssess==T),aes(factor(Metric),RefPoint,fill=GoodThing))
                   +geom_violin()+facet_wrap(~Species,scales='free_y')+PaperBarTheme+geom_hline(yintercept=1)+
                     xlab("")+ylab('Relative to Reference Point')
                   +scale_fill_manual(name='',values=c('red','blue')))
  (ReferencePlot)
  dev.off()


  AllMonte<- AllMonte %>%
    group_by(Species,Assessment) %>%
    mutate(NumAssessYears=length(unique(Year)))

  pdf(file=paste(Assessment,'Monte Carlo Summary Plots.pdf',sep='/'))

  FMCCSumPlot<- (ggplot(data=subset(AllMonte,Metric=="F/M (CC)" & NumAssessYears>=1),aes(Year,Value))
                 +geom_point(aes(Year,Value),alpha=0.1)+ylab('F/M (CC)')+
                   geom_smooth(method='loess',size=2)+facet_wrap(~Species,scales='free_y')
                 +geom_hline(yintercept=1,linetype='longdash')+PaperTheme+theme(axis.title.x=element_blank()))

  FMCCSumPlot

  FMLBSPRSumPlot<- (ggplot(data=subset(AllMonte,Metric=="F/M (LBSR)" & NumAssessYears>=1),aes(Year,Value))
                    +geom_point(aes(Year,Value),alpha=0.1) +ylab('F/M (LBSPR)') +
                      geom_smooth(method='loess',size=2) +facet_wrap(~Species,scales='free_y')
                    +geom_hline(yintercept=1,linetype='longdash')
                    +PaperTheme+theme(axis.title.x=element_blank()))

  FMLBSPRSumPlot

  SPRSumPlot<- (ggplot(data=subset(AllMonte,Metric=="SPR" & NumAssessYears>=1),aes(Year,Value))
                +geom_point(aes(Year,Value),alpha=0.1)+ylab('SPR')+
                  geom_smooth(method='loess',size=2)+facet_wrap(~Species,scales='free_y')+geom_hline(yintercept=0.4,linetype='longdash')
                +PaperTheme+theme(axis.title.x=element_blank()))

  SPRSumPlot

  CPUESumPlot<- (ggplot(data=subset(AllMonte,Metric=="CPUERatio" & NumAssessYears>=1),aes(Year,Value))
                 +geom_point(aes(Year,Value),alpha=0.1)+ylab('CPUE Ratio')+
                   geom_smooth(method='loess',size=2)+facet_wrap(~Species,scales='free_y')+geom_hline(yintercept=0.4,linetype='longdash')
                 +PaperTheme+theme(axis.title.x=element_blank()))

  CPUESumPlot

  dev.off()

  pdf(file=paste(Assessment,'2014 Status Plots.pdf',sep='/'))
  RecentStatusPlot<- (ggplot(data=subset(AllMonte,Year==2014),aes(Species,SmoothValue,fill=Species))+geom_violin()+facet_wrap(~Metric,scales='free_y')+
                        theme(legend.position='none',axis.title.x=element_blank())
                      +PaperBarTheme+ylab('Value'))
  RecentStatusPlot
  dev.off()

  pdf(file=paste(Assessment,'2014 Status BoxPlots.pdf',sep='/'))
  RecentStatusBoxPlot<- (ggplot(data=subset(AllMonte,Year==2014),aes(Species,SmoothValue,fill=Species))+geom_boxplot()+facet_wrap(~Metric,scales='free_y')+
                           theme(legend.position='none',axis.title.x=element_blank())
                         +PaperBarTheme+ylab('Value'))
  RecentStatusBoxPlot
  dev.off()

  # ggplot(data=subset(AllMonte,Metric=="F/M (CC)" & Species=='Black Rockfish'),aes((Year),Value))+geom_smooth()

  save(file=paste(Assessment,'/AllMonte.rdata',sep=''),AllMonte)

  save(RecentStatusPlot,RecentStatusBoxPlot,FMCCSumPlot,FMLBSPRSumPlot,SPRSumPlot,CPUESumPlot,PercDevPlot,ReferencePlot,CPUETrend,SPRTrend,CCFvMTrendPlot,CPUEPlot,SPRPlot,CCFvMPlot,LBSPRFvMPlot,file=paste(Assessment,'/MonteCarlo Plots.rdata',sep=''))





  return(list(Results = Results, LengthStructure = LengthStructure, Diagnostics = Diagnostics))
}