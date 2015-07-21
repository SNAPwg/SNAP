PlotCPUEData<- function(CPUEDat,FigureFolder,Fish,Species,Site,Theme)
{
  
  #   DensityDat<- DensityData
#   CPUESummary<- ddply(CPUEDat,c('Year','MPA'),summarize,NumberCPUE=mean(Count/AnglerHours,na.rm=T),BiomassCPUE=mean(Biomass/AnglerHours,na.rm=T))
  CPUESummary<- CPUEDat %>%
    group_by(Year,MPA) %>%
      summarize(NumberCPUE=mean(Count/AnglerHours,na.rm=T),BiomassCPUE=mean(Biomass/AnglerHours,na.rm=T))
    
#     summarize(NumberCPUE=sum(Count,na.rm=T)/sum(AnglerHours,na.rm=T),BiomassCPUE=sum(Biomass,na.rm=T)/sum(AnglerHours,na.rm=T))
  
#     ddply(CPUEDat,c('Year','MPA'),summarize,NumberCPUE=mean(Count/AnglerHours,na.rm=T),BiomassCPUE=mean(Biomass/AnglerHours,na.rm=T))
  
  CPUESummary$SiteType[CPUESummary$MPA==0]<- 'Fished'
  
  CPUESummary$SiteType[CPUESummary$MPA==1]<- 'MPA'
  
  CPUEDat$SiteType[CPUESummary$MPA==0]<- 'Fished'
  
  CPUEDat$SiteType[CPUESummary$MPA==1]<- 'MPA'
  
  pdf(file=paste(FigureFolder,Species,'-',Site,' CPUE Data.pdf',sep=''),height=5,width=5)
  print(ggplot(data=CPUESummary,aes(Year,NumberCPUE,color=SiteType))+geom_line(size=3)+
          scale_color_manual(name='',values=c(FishedColor,MPAColor))+
          ylab(paste('Number/Angler Hour',sep=''))+ggtitle(paste(Species,Site,sep='-'))+Theme)
  
  BioCPUE<- (ggplot(data=CPUESummary,aes(Year,BiomassCPUE,color=SiteType))+geom_line(size=3)+
               scale_color_manual(name='',values=c(FishedColor,MPAColor))+
               ylab(paste('kg/Angler Hour',sep=''))+ggtitle(paste(Species,Site,sep='-'))+Theme)
  print(BioCPUE)
  
  DotCPUE<- ((ggplot(data=CPUEDat,aes(Year,Biomass/AnglerHours,color=SiteType)))+geom_point(size=6,alpha=0.6)
             +ylab(paste('kg/Angler Hour per trip',sep=''))+ggtitle(paste(Species,Site,sep='-'))
             +scale_color_manual(name='',values=c(FishedColor,MPAColor))+Theme)
  
  BoxCPUE<- ((ggplot(data=CPUEDat,aes(factor(Year),Biomass/AnglerHours,fill=SiteType)))+geom_boxplot()
             +ylab(paste('kg/Angler Hour per trip',sep=''))+ggtitle(paste(Species,Site,sep='-'))
             +scale_fill_manual(name='',values=c(FishedColor,MPAColor))+Theme)
  
  print(DotCPUE)
  
  dev.off()
  write.csv(file=paste(ResultFolder,AssessmentName,' CPUE Data Summary.csv',sep=''),CPUESummary)
  
  save(BioCPUE,DotCPUE,BoxCPUE,file=paste(FigureFolder,Species,'-',Site,'CPUE Data.Rdata',sep=''))
  
  return(BioCPUE)
}