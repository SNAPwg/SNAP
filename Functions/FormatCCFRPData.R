FormatCCFRPData<- function(Data)
{
  # Format Length Data ------------------------------------------------------
  
  #   Data<- iGFD
  LengthDataNames<- c('Fish','Year','Month','Site','Length','LengthType','Sex','Special','MPA','FisheryDependent','MeanLongitude','MeanLatitude')
  
  LengthData<- as.data.frame(matrix(NA,nrow=dim(Data)[1],ncol=length(LengthDataNames)))
  
  colnames(LengthData)<- LengthDataNames
  
  LengthData$Fish<- Data$CommName 
  
  LengthData$Year<- Data$Year
  
  LengthData$Month<- Data$Month
  
  LengthData$SiteType<- Data$SiteId
  
  LengthData$Length<- Data$length_cm
  
  LengthData$Length[LengthData$Length==0]<- NA
  
  LengthData$LengthType<- 'cm'
  
  LengthData$Sex<- 'Unknown'
  
  LengthData$Special<- paste('Gear is ',Data$Sample_Type,sep='')
  
  LengthData$MPA<- Data$MPA_or_REF
  
  LengthData$MPA[LengthData$MPA=='REF']<- 0
  
  LengthData$MPA[LengthData$MPA=='MPA']<- 1
  
  LengthData$MPA<- as.numeric(LengthData$MPA)
  
  LengthData$FisheryDependent<- 1
  
  LengthData$MeanLongitude<- Data$MeanLon
  
  LengthData$MeanLatitude<- Data$MeanLat
  
  # Format Density Data -----------------------------------------------------
  Data$Weight<- Fish$WeightA* Data$length_cm ^ Fish$WeightB
  
  #   DensityData<- ddply(Data,c('Year','Month','sample_Idcellday'),plyr::summarize,Site='All',Count=sum(length_cm>0 | is.na(length_cm),na.rm=T),Biomass=sum(Weight,na.rm=T)
  #                       ,SampleArea= mean(Sample_Area,na.rm=T),AreaUnits=unique(Area_units),DistanceFromBorder=mean(Meters.to.MPA.border,na.rm=T)
  #                       ,SampleType=unique(Sample_Type),MPA=unique(MPA_or_REF),
  #                       DistanceProtected=mean(Meters.to.MPA.border,na.rm=T),MeanLongitude=mean(MeanLon,na.rm=T),
  #                       MeanLatitude=mean(MeanLat,na.rm=T),SiteType=unique(SiteId),Fish=unique(CommName))
  #   
  DensityData<- Data %>%
    group_by(Year,Month,sample_Idcellday) %>%
    summarize(Count=sum(length_cm>0 | is.na(length_cm),na.rm=T),Biomass=sum(Weight,na.rm=T)
              ,SampleArea= mean(Sample_Area,na.rm=T),AreaUnits=unique(Area_units),DistanceFromBorder=mean(Meters.to.MPA.border,na.rm=T)
              ,SampleType=unique(Sample_Type),MPA=unique(MPA_or_REF),
              DistanceProtected=mean(Meters.to.MPA.border,na.rm=T),MeanLongitude=mean(MeanLon,na.rm=T),
              MeanLatitude=mean(MeanLat,na.rm=T),SiteType=unique(SiteId),Fish=unique(CommName))
  
  DensityData$MPA[DensityData$MPA=='REF']<- 0
  
  DensityData$MPA[DensityData$MPA=='MPA']<- 1
  
  DensityData$MPA<- as.numeric(DensityData$MPA)
  
  DensityData$DistanceProtected[DensityData$MPA==0]<- (DensityData$DistanceProtected*-1)[DensityData$MPA==0]
  
  
  # Format CPUE Data -----------------------------------------------------
  
  CPUEData<- ddply(Data,c('Year','Month','sample_Idcellday'),plyr::summarize,Site='All',Count=sum(length_cm>0 | is.na(length_cm),na.rm=T),Biomass=sum(Weight,na.rm=T)
                   ,AnglerHours= sum(Angler_hours,na.rm=T),DistanceFromBorder=mean(Meters.to.MPA.border,na.rm=T)
                   ,SampleType=unique(Sample_Type),MPA=unique(MPA_or_REF),
                   DistanceProtected=mean(Meters.to.MPA.border,na.rm=T),MeanLongitude=mean(MeanLon,na.rm=T),
                   MeanLatitude=mean(MeanLat,na.rm=T),SiteType=unique(SiteId),Fish=unique(CommName))
  
  CPUEData<- Data %>%
    group_by(Year,Month,sample_Idcellday) %>%
    summarize(Site='All',Count=sum(length_cm>0 | is.na(length_cm),na.rm=T),Biomass=sum(Weight,na.rm=T)
              ,AnglerHours= sum(Angler_hours,na.rm=T),DistanceFromBorder=mean(Meters.to.MPA.border,na.rm=T)
              ,SampleType=unique(Sample_Type),MPA=unique(MPA_or_REF),
              DistanceProtected=mean(Meters.to.MPA.border,na.rm=T),MeanLongitude=mean(MeanLon,na.rm=T),
              MeanLatitude=mean(MeanLat,na.rm=T),SiteType=unique(SiteId),Fish=unique(CommName))
  
  CPUEData$MPA[CPUEData$MPA=='REF']<- 0
  
  CPUEData$MPA[CPUEData$MPA=='MPA']<- 1
  
  CPUEData$MPA<- as.numeric(CPUEData$MPA)
  
  CPUEData$DistanceProtected[CPUEData$MPA==0]<- CPUEData$DistanceProtected[CPUEData$MPA==0]*-1
  
  return(list(LengthData=LengthData,DensityData=DensityData,CPUEData=CPUEData))
  
}




