PlotDensityData<- function(DensityDat,FigureFolder,Fish,Species,Site,Theme)
{
#   DensityDat<- DensityData
  DensitySummary<- ddply(DensityDat,c('Year','MPA'),plyr::summarize,NumberDensity=mean(Count/SampleArea,na.rm=T),BiomassDensity=mean(Biomass/SampleArea,na.rm=T))
  
  DensitySummary$SiteType[DensitySummary$MPA==0]<- 'Fished'
  
  DensitySummary$SiteType[DensitySummary$MPA==1]<- 'MPA'

  DensityDat$SiteType[DensitySummary$MPA==0]<- 'Fished'
  
  DensityDat$SiteType[DensitySummary$MPA==1]<- 'MPA'
  
  pdf(file=paste(FigureFolder,Species,'-',Site,' Density Data.pdf',sep=''),height=5,width=5)
  print(ggplot(data=DensitySummary,aes(Year,NumberDensity,color=SiteType))+geom_line(size=3)+
          scale_color_manual(name='',values=c(FishedColor,MPAColor))+
    ylab(expression(paste('Number/',km^{2},sep='')))+ggtitle(paste(Species,Site,sep='-'))+Theme)
  
  print(ggplot(data=DensitySummary,aes(Year,BiomassDensity,color=SiteType))+geom_line(size=3)+
          scale_color_manual(name='',values=c(FishedColor,MPAColor))+
     ylab(expression(paste('kg/',km^{2},sep='')))+ggtitle(paste(Species,Site,sep='-'))+Theme)
  
  print((ggplot(data=DensityDat,aes(Year,Biomass/SampleArea,color=SiteType)))+geom_point(size=6,alpha=0.6)
  +ylab(expression(paste('kg/',km^{2}, ' per trip',sep='')))+ggtitle(paste(Species,Site,sep='-'))
  +scale_color_manual(name='',values=c(FishedColor,MPAColor))+Theme)
  
  dev.off()
  write.csv(file=paste(ResultFolder,AssessmentName,' Density Data Summary.csv',sep=''),DensitySummary)
  
}