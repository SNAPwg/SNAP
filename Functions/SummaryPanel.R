SummaryPanel<- function(AssessData,LengthDat,CPUEDat,Species,Site,YearsToSmooth,Theme)
{
  
  
  LengthDat$MPA<- as.factor(LengthData$MPA)
  
  MaxYear<- max(LengthDat$Year,na.rm=T)
  
  LengthDat$Year<- as.factor(LengthData$Year)
  
  levels(LengthDat$MPA)<- c('Fished','MPA')
  
  MeanLength<- ddply(LengthDat,c('Year','MPA'),summarize,MeanLength=mean(Length,na.rm=T))
  
  
  CPUESummary<- ddply(CPUEDat,c('Year','MPA'),summarize,NumberCPUE=mean(Count/AnglerHours,na.rm=T),BiomassCPUE=mean(Biomass/AnglerHours,na.rm=T))
  
  CPUESummary$SiteType[CPUESummary$MPA==0]<- 'Fished'
  
  CPUESummary$SiteType[CPUESummary$MPA==1]<- 'MPA'
  
  CPUEDat$SiteType[CPUESummary$MPA==0]<- 'Fished'
  
  CPUEDat$SiteType[CPUESummary$MPA==1]<- 'MPA'
  
  BioCPUE<- (ggplot(data=CPUESummary,aes(Year,BiomassCPUE,color=SiteType))+geom_smooth(size=2,se=F)+
               scale_color_manual(name='',values=c(FishedColor,MPAColor))+
               ylab(paste('kg/Angler Hour',sep=''))+Theme)
  
  
#   pdf(file=paste(FigureFolder,'Assessment Summary.pdf',sep=''),width=8,height=6)
  
  #   theme(legend.position='top')
  LengthPlot<- (ggplot(data=subset(LengthDat,Year==MaxYear),aes(Length,fill=MPA))+
                  geom_density(alpha=0.7,aes(y=..count..))+xlab('Length (cm)')
                +geom_vline(xintercept=Fish$Mat50,linetype='longdash',size=2)
                +scale_fill_manual(name='',values=c(FishedColor,MPAColor))+Theme+ylab('Count'))
  
  limits <- aes(ymax = UpperCI, ymin=LowerCI)
  
  AssessData$LowerCI[is.na(AssessData$LowerCI)]<- AssessData$Value[is.na(AssessData$LowerCI)]
  
  AssessData$UpperCI[is.na(AssessData$UpperCI)]<- AssessData$Value[is.na(AssessData$UpperCI)]
  #   FPlot<- (ggplot(data=subset(AssessData,Method=='CatchCurve'),aes(x=Year,y=Value))+
  #     geom_smooth(se=F,size=2)+geom_errorbar(limits)+ylab('F/Fmsy')+geom_hline(yintercept=1))
  FPlot<- (ggplot(data=subset(AssessData,(Method=='CatchCurve' | Method=='LBSPR') & Metric=='FvM'),aes(x=Year,y=Value,color=Method))+
             geom_smooth(se=F,aes(ymin=LowerCI,ymax=UpperCI,fill=Method),stat='identity',size=2,alpha=0.4,size=2)+ylab('F/M')+
             geom_hline(yintercept=1,size=2,linetype='longdash')
           +scale_color_manual(name='',values=c("#1f78b4","#33a02c"))+scale_fill_manual(name='',values=c("#1f78b4","#33a02c"))+Theme)
  
  
  CPUEPlot<- (ggplot(data=subset(AssessData,Method=='CPUERatio'),aes(x=Year,y=Value))+
                geom_smooth(se=F,stat='identity',aes(ymin=LowerCI,ymax=UpperCI),fill='#253494',size=2,color='#253494')+ylab('CPUE Ratio')+geom_hline(yintercept=1,size=2,linetype='longdash')
              +Theme)
  
  SPRPlot<- (ggplot(data=subset(AssessData,Method=='LBSPR' & Metric=='SPR'),aes(x=Year,y=Value))+
               geom_smooth(se=F,stat='identity',aes(ymin=LowerCI,ymax=UpperCI),size=2,fill='#b30000',color='#b30000')+ylab('SPR')+geom_hline(yintercept=0.4,size=2,linetype='longdash')
             +Theme)
  
  #   grid.arrange(LengthPlot,FPlot,BioCPUE,SPRPlot
  #                ,nrow=2,ncol=2,main=paste(Species,'-',Sites[s],sep=''))
  SummaryGrid=arrangeGrob(LengthPlot,FPlot,BioCPUE,SPRPlot
               ,nrow=2,ncol=2)
  
  ggsave(SummaryGrid,file=paste(FigureFolder,'Assessment Summary.pdf',sep=''),width=8,height=6)
  
#   dev.off()
  
  
#   pdf(file=paste(FigureFolder,'Assessment Summary Smooth.pdf',sep=''),width=8,height=6)
  
  #   theme(legend.position='top')
  LengthPlot<- (ggplot(data=subset(LengthDat,Year==MaxYear),aes(Length,fill=MPA))+
                  geom_density(alpha=0.7,aes(y=..count..))+xlab('Length (cm)')
                +geom_vline(xintercept=Fish$Mat50,linetype='longdash',size=2)
                +scale_fill_manual(name='',values=c(FishedColor,MPAColor))+Theme+ylab('Count'))
  
  limits <- aes(ymax = UpperCI, ymin=LowerCI)
  
  AssessData$LowerCI[is.na(AssessData$LowerCI)]<- AssessData$Value[is.na(AssessData$LowerCI)]
  
  AssessData$UpperCI[is.na(AssessData$UpperCI)]<- AssessData$Value[is.na(AssessData$UpperCI)]
  #   FPlot<- (ggplot(data=subset(AssessData,Method=='CatchCurve'),aes(x=Year,y=Value))+
  #     geom_smooth(se=F,size=2)+geom_errorbar(limits)+ylab('F/Fmsy')+geom_hline(yintercept=1))
  FPlot<- (ggplot(data=subset(AssessData,(Method=='CatchCurve' | Method=='LBSPR') & Metric=='FvM'),aes(x=Year,y=Value,color=Method))+
             geom_smooth(se=F,size=2,size=2)+ylab('F/M')+
             geom_hline(yintercept=1,size=2,linetype='longdash')
           +scale_color_manual(name='',values=c("#1f78b4","#33a02c"))+Theme)
  
  
  CPUEPlot<- (ggplot(data=subset(AssessData,Method=='CPUERatio'),aes(x=Year,y=Value))+
                geom_smooth(se=F,size=2,color='#253494')+ylab('CPUE Ratio')+geom_hline(yintercept=1,size=2,linetype='longdash')
              +Theme)
  
  SPRPlot<- (ggplot(data=subset(AssessData,Method=='LBSPR' & Metric=='SPR'),aes(x=Year,y=Value))+
               geom_smooth(se=F,size=2,color='#b30000')+ylab('SPR')+geom_hline(yintercept=0.4,size=2,linetype='longdash')
             +Theme)
  
  #   grid.arrange(LengthPlot,FPlot,BioCPUE,SPRPlot
  #                ,nrow=2,ncol=2,main=paste(Species,'-',Sites[s],sep=''))
  SmoothSummary<- arrangeGrob(LengthPlot,FPlot,BioCPUE,SPRPlot
               ,nrow=2,ncol=2)
  
  ggsave(SmoothSummary,file=paste(FigureFolder,'Assessment Summary Smooth.pdf',sep=''),width=8,height=6)
  
#   dev.off()
  
#   pdf(file=paste(FigureFolder,'Assessment Summary 2.pdf',sep=''),width=8,height=6)
  
  #   theme(legend.position='top')
  LengthPlot<- (ggplot(data=subset(LengthDat,Year==MaxYear),aes(Length,fill=MPA))+
                  geom_density(alpha=0.7,aes(y=..count..))+xlab('Length (cm)')
                +geom_vline(xintercept=Fish$Mat50,linetype='longdash',size=2)
                +scale_fill_manual(name='',values=c(FishedColor,MPAColor))+Theme+ylab('Count'))
  
  limits <- aes(ymax = UpperCI, ymin=LowerCI)
  
  AssessData$LowerCI[is.na(AssessData$LowerCI)]<- AssessData$Value[is.na(AssessData$LowerCI)]
  AssessData$UpperCI[is.na(AssessData$UpperCI)]<- AssessData$Value[is.na(AssessData$UpperCI)]
  #   FPlot<- (ggplot(data=subset(AssessData,Method=='CatchCurve'),aes(x=Year,y=Value))+
  #     geom_smooth(se=F,size=2)+geom_errorbar(limits)+ylab('F/Fmsy')+geom_hline(yintercept=1))
  FPlot<- (ggplot(data=subset(AssessData,(Method=='CatchCurve') & Metric=='FvM'),aes(x=Year,y=Value))+
             geom_smooth(se=F,aes(ymin=LowerCI,ymax=UpperCI),stat='identity',size=2,alpha=0.4,size=2,
                         color="#336600",fill="#336600")+ylab('F/M')+
             geom_hline(yintercept=1,size=2,linetype='longdash')+Theme)
  
  CPUEPlot<- (ggplot(data=subset(AssessData,Method=='CPUERatio'),aes(x=Year,y=Value))+
                geom_smooth(se=F,stat='identity',aes(ymin=LowerCI,ymax=UpperCI),fill='#253494',size=2,color='#253494')+ylab('CPUE Ratio')+geom_hline(yintercept=1,size=2,linetype='longdash')
              +Theme)
  
  SPRPlot<- (ggplot(data=subset(AssessData,Method=='LBSPR' & Metric=='SPR'),aes(x=Year,y=Value))+
               geom_smooth(se=F,stat='identity',aes(ymin=LowerCI,ymax=UpperCI),size=2,fill='#b30000',color='#b30000')+ylab('SPR')+geom_hline(yintercept=0.4,size=2,linetype='longdash')
             +Theme)
  
  #   grid.arrange(LengthPlot,FPlot,BioCPUE,SPRPlot
  #                ,nrow=2,ncol=2,main=paste(Species,'-',Sites[s],sep=''))
  SummaryGrid2=arrangeGrob(LengthPlot,FPlot,BioCPUE,SPRPlot
                       ,nrow=2,ncol=2)
  
  ggsave(SummaryGrid2,file=paste(FigureFolder,'Assessment Summary 2.pdf',sep=''),width=8,height=6)

  save(SummaryGrid,SummaryGrid2,SmoothSummary,file=paste(FigureFolder,'Assessment Summaries.Rdata',sep=''))
  
#   print(GridPlot)
#   dev.off()
}