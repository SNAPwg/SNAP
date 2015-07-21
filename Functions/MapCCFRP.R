MapCCFRP<- function(Data)
{
  
  ## Map length data ##
  
  LengthData<- Data$LengthData %>% subset(is.na(SiteType)==F)
  
  MapDat<- ddply(LengthData,c('SiteType','Year','Fish'),summarize,MeanLength=mean(Length,na.rm=T),Lon=mean(MeanLongitude),
                 Lat=mean(MeanLatitude),Reserve=unique(MPA))
  
  MapDat$MPA[MapDat$Reserve==0]<- 'Reference'
  
  MapDat$MPA[MapDat$Reserve==1]<- 'Reserve'
  
  CentralCoast<- get_map(location=c(lon=mean(MapDat$Lon),lat=mean(MapDat$Lat)),zoom=8,maptype='hybrid')
  
  Years<- unique(MapDat$Year)
  
  BaseDir<- getwd()
  
#   setwd(FigureFolder)
#   
#   saveGIF(
# {
#   for (y in Years)
#   {
#     print(ggmap(CentralCoast, fullpage = TRUE)+
#             geom_point(aes(x=Lon,y=Lat,color=factor(MPA),size=MeanLength),alpha=0.8,
#                        data=subset(MapDat,Year==y))+ggtitle(y)
#           +scale_size_continuous(range=c(2,10),name='Mean Length (cm)')+
#             scale_color_manual(name='',values=c(FishedColor,MPAColor)))
#   }
#   
# },movie.name='MPA Average Size Movie.gif',imgdir=FigureFolder)
# 
# setwd(BaseDir)


pdf(file=paste(FigureFolder,'Study Sites.pdf',sep=''))
print(ggmap(CentralCoast, fullpage = TRUE)+
        geom_point(aes(x=Lon,y=Lat,color=factor(MPA),shape=factor(MPA)),size=7,alpha=0.9,data=subset(MapDat,Year==2014))+
        scale_color_manual(name='',values=c(FishedColor,MPAColor)))
dev.off()



pdf(file=paste(FigureFolder,'MPA Average Size.pdf',sep=''))
print(ggmap(CentralCoast, fullpage = TRUE)+
        geom_point(aes(x=Lon,y=Lat,color=factor(MPA),size=MeanLength),alpha=0.7,data=(MapDat))
      +ggtitle(unique(MapDat$Fish))+scale_size_continuous(range=c(2,10),name='Mean Length (cm)')+
        scale_color_manual(name='',values=c(FishedColor,MPAColor))+facet_wrap(~Year,as.table=F))

dev.off()

pdf(file=paste(FigureFolder,'Final MPA Average Size.pdf',sep=''))
print(ggmap(CentralCoast, fullpage = TRUE)+
        geom_point(aes(x=Lon,y=Lat,color=factor(MPA),size=MeanLength),alpha=0.7,data=subset(MapDat,Year==max(Year)))
      +ggtitle(unique(MapDat$Fish))+scale_size_continuous(range=c(2,10),name='Mean Length (cm)')+
        scale_color_manual(name='',values=c(FishedColor,MPAColor)))

dev.off()




DensityData<- Data$DensityData %>% subset(is.na(SiteType)==F)

MapDat<- ddply(DensityData,c('SiteType','Year','Fish'),summarize,
               Density=sum(Biomass,na.rm=T)/sum(SampleArea,na.rm=T),
               WeightedDensity=sum(DistanceFromBorder*Biomass/SampleArea,na.rm=T)/sum(DistanceFromBorder,na.rm=T),
               Lon=mean(MeanLongitude),
               Lat=mean(MeanLatitude),Reserve=unique(MPA))

MapDat$MPA[MapDat$Reserve==0]<- 'Reference'

MapDat$MPA[MapDat$Reserve==1]<- 'Reserve'

Years<- unique(MapDat$Year)


# setwd(FigureFolder)
# 
# saveGIF(
# {
#   for (y in Years)
#   {
#     print(ggmap(CentralCoast, fullpage = TRUE)+
#             geom_point(aes(x=Lon,y=Lat,color=factor(MPA),size=Density),alpha=0.7,
#                        data=subset(MapDat,Year==y))+ggtitle(y)
#           +scale_size_continuous(range=c(2,10),name='Density (kg/km2)')+
#             scale_color_manual(name='',values=c(FishedColor,MPAColor)))
#   }
#   
# },movie.name='MPA Density Movie.gif',imgdir=FigureFolder)

# setwd(BaseDir)

pdf(file=paste(FigureFolder,'MPA Density.pdf',sep=''))
print(ggmap(CentralCoast, fullpage = TRUE)+
        geom_point(aes(x=Lon,y=Lat,color=factor(MPA),size=Density),alpha=0.7,data=(MapDat))
      +ggtitle(unique(MapDat$Fish))+scale_size_continuous(range=c(2,10),name='Density (kg/km2)')+
        scale_color_manual(name='',values=c(FishedColor,MPAColor))+facet_wrap(~Year,as.table=F))

dev.off()


pdf(file=paste(FigureFolder,'Distance Weighted MPA Density.pdf',sep=''))
print(ggmap(CentralCoast, fullpage = TRUE)+
        geom_point(aes(x=Lon,y=Lat,color=factor(MPA),size=WeightedDensity),alpha=0.7,data=(MapDat))
      +ggtitle(unique(MapDat$Fish))+scale_size_continuous(range=c(2,10),name='Density (kg/km2)')+
        scale_color_manual(name='',values=c(FishedColor,MPAColor))+facet_wrap(~Year,as.table=F))

dev.off()

pdf(file=paste(FigureFolder,'Final Distance Weighted MPA Density.pdf',sep=''))
print(ggmap(CentralCoast, fullpage = TRUE)+
        geom_point(aes(x=Lon,y=Lat,color=factor(MPA),size=WeightedDensity),alpha=0.7,data=subset(MapDat,Year==max(Year)))
      +ggtitle(unique(MapDat$Fish))+scale_size_continuous(range=c(2,10),name='Density (kg/km2)')+
        scale_color_manual(name='',values=c(FishedColor,MPAColor)))

dev.off()

}